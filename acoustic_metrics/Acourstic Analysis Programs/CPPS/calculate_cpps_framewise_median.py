#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frame-wise CPPS calculation for WAV files using Praat/Parselmouth.

This script calculates CPPS for each 2048-sample frame of a WAV file and
summarizes each file by the median of the frame-wise CPPS values.

Primary design:
    - Analyze the whole WAV file, not only a fixed vowel-midpoint segment.
    - Split the signal into 2048-sample non-overlapping frames.
    - Retain the final frame even if it is shorter than 2048 samples.
    - Apply a Hann window to each actual frame.
    - Zero-pad short final frames to 2048 samples for FFT computation.
    - Convert each frame to a Praat PowerCepstrum via Parselmouth.
    - Smooth the PowerCepstrum along the quefrency axis.
    - Calculate cepstral peak prominence for each frame.
    - Use the median of all finite frame-wise CPPS values as the file-level CPPS.
    - Flag likely unvoiced or low-periodicity frames.

Important:
    CPPS values depend on implementation details. This script uses Praat's
    PowerCepstrum "Get peak prominence" command for each fixed 2048-sample
    frame. It is not the same as Praat's whole-signal "PowerCepstrogram: Get
    CPPS" command, which uses its own analysis-window settings.

Outputs:
    - One frame-wise CSV per WAV file.
    - One summary CSV for all processed WAV files.

Dependencies:
    numpy, scipy, praat-parselmouth
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import parselmouth
from parselmouth.praat import call
from scipy.io import wavfile


EPS = 1e-20


@dataclass
class CPPSParameters:
    frame_size: int = 2048
    hop_size: int = 2048
    pitch_floor_hz: float = 60.0
    pitch_ceiling_hz: float = 1000.0
    quefrency_smoothing_s: float = 0.0005
    smoothing_iterations: int = 1
    interpolation: str = "parabolic"
    trend_q_min_s: float = 0.001
    trend_q_max_s: float = 0.05
    trend_type: str = "straight"
    fit_method: str = "robust slow"
    rms_unvoiced_threshold_dbfs: float = -60.0
    autocorr_voicing_threshold: float = 0.30


@dataclass
class FrameResult:
    frame_index: int
    start_sample: int
    end_sample: int
    frame_length_samples: int
    start_time_s: float
    end_time_s: float
    short_final_frame: bool
    rms_dbfs: float
    autocorr_voicing_score: float
    is_voiced: bool
    voicing_warning: str
    cpps_db: float
    peak_quefrency_s: float
    peak_frequency_hz: float
    calculation_warning: str


@dataclass
class FileSummary:
    input_file: str
    output_frame_csv: str
    sampling_rate_hz: int
    duration_s: float
    frame_size_samples: int
    hop_size_samples: int
    frame_duration_ms: float
    n_frames_total: int
    n_voiced_frames: int
    n_unvoiced_frames: int
    unvoiced_frame_indices: str
    median_cpps_all_frames_db: float
    mean_cpps_all_frames_db: float
    sd_cpps_all_frames_db: float
    median_cpps_voiced_frames_db: float | None
    mean_cpps_voiced_frames_db: float | None
    sd_cpps_voiced_frames_db: float | None
    warning: str


def read_wav_as_mono_float(wav_path: Path) -> tuple[int, np.ndarray]:
    """Read a WAV file and return sample rate and mono float64 signal."""
    sr, x = wavfile.read(str(wav_path))

    if x.ndim == 2:
        x = x.mean(axis=1)

    if np.issubdtype(x.dtype, np.integer):
        max_abs = max(abs(np.iinfo(x.dtype).min), np.iinfo(x.dtype).max)
        x = x.astype(np.float64) / max_abs
    else:
        x = x.astype(np.float64)

    return sr, x


def iter_wav_files(input_path: Path, recursive: bool = False) -> list[Path]:
    """Return WAV files from a single file or directory."""
    if input_path.is_file():
        if input_path.suffix.lower() != ".wav":
            raise ValueError(f"Input file is not a .wav file: {input_path}")
        return [input_path]

    if input_path.is_dir():
        pattern = "**/*.wav" if recursive else "*.wav"
        wav_files = sorted(input_path.glob(pattern))
        if not wav_files:
            raise FileNotFoundError(f"No .wav files found in: {input_path}")
        return wav_files

    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def make_output_stem(wav_path: Path, root_path: Path) -> str:
    """Create a safe output stem for a WAV file."""
    try:
        if root_path.is_dir():
            rel = wav_path.relative_to(root_path)
        else:
            rel = wav_path.name
    except ValueError:
        rel = wav_path.name

    rel_str = str(rel)
    if rel_str.lower().endswith(".wav"):
        rel_str = rel_str[:-4]

    safe = rel_str.replace("\\", "__").replace("/", "__")
    safe = safe.replace(":", "_").replace(" ", "_")
    return safe


def split_frames_include_last(
    x: np.ndarray,
    frame_size: int,
    hop_size: int,
) -> Iterable[tuple[int, np.ndarray]]:
    """Yield frames while retaining the final short frame."""
    if frame_size <= 0:
        raise ValueError("frame_size must be positive.")
    if hop_size <= 0:
        raise ValueError("hop_size must be positive.")
    if len(x) == 0:
        return

    start = 0
    while start < len(x):
        frame = x[start:start + frame_size]
        if len(frame) == 0:
            break
        yield start, frame
        start += hop_size


def rms_dbfs(frame: np.ndarray) -> float:
    """Return frame RMS in dBFS for normalized audio."""
    if len(frame) == 0:
        return float("nan")
    rms = float(np.sqrt(np.mean(frame.astype(np.float64) ** 2)))
    return 20.0 * math.log10(rms + EPS)


def autocorr_voicing_score(
    frame: np.ndarray,
    sr: int,
    pitch_floor_hz: float,
    pitch_ceiling_hz: float,
) -> float:
    """
    Estimate frame periodicity using the normalized autocorrelation peak.

    This value is used only as a warning flag for likely unvoiced or
    low-periodicity frames. CPPS is still calculated for all frames when
    possible and included in the all-frame median.
    """
    if len(frame) < 4:
        return float("nan")

    y = frame.astype(np.float64) - np.mean(frame)
    energy = float(np.dot(y, y))
    if energy <= EPS:
        return 0.0

    if len(y) > 1:
        y = y * np.hanning(len(y))

    corr = np.correlate(y, y, mode="full")[len(y) - 1:]
    if corr[0] <= EPS:
        return 0.0
    corr = corr / (corr[0] + EPS)

    min_lag = int(math.floor(sr / pitch_ceiling_hz))
    max_lag = int(math.ceil(sr / pitch_floor_hz))
    min_lag = max(1, min_lag)
    max_lag = min(len(corr) - 1, max_lag)

    if max_lag <= min_lag:
        return float("nan")

    return float(np.max(corr[min_lag:max_lag + 1]))


def prepare_frame_for_praat(frame: np.ndarray, frame_size: int) -> np.ndarray:
    """
    Apply DC removal, Hann window, and zero-padding to a frame.

    The Hann window is applied to the actual frame length. If the final frame
    is shorter than frame_size, the windowed frame is then zero-padded to
    frame_size samples.
    """
    if len(frame) == 0:
        return np.zeros(frame_size, dtype=np.float64)

    y = frame.astype(np.float64) - np.mean(frame)
    if len(y) > 1:
        y = y * np.hanning(len(y))

    if len(y) < frame_size:
        padded = np.zeros(frame_size, dtype=np.float64)
        padded[:len(y)] = y
        return padded

    return y[:frame_size]


def calculate_praat_cpps_for_frame(
    frame: np.ndarray,
    sr: int,
    params: CPPSParameters,
) -> tuple[float, float, float, str]:
    """
    Calculate CPPS for one frame using Praat PowerCepstrum commands.

    Returns:
        cpps_db, peak_quefrency_s, peak_frequency_hz, calculation_warning
    """
    warnings: list[str] = []

    if len(frame) == 0:
        return float("nan"), float("nan"), float("nan"), "empty_frame"

    if rms_dbfs(frame) < -120.0:
        return float("nan"), float("nan"), float("nan"), "near_silent_frame"

    y = prepare_frame_for_praat(frame, params.frame_size)
    if len(frame) < params.frame_size:
        warnings.append("zero_padded_short_frame")

    try:
        sound = parselmouth.Sound(y, sampling_frequency=sr)
        spectrum = call(sound, "To Spectrum", "yes")
        power_cepstrum = call(spectrum, "To PowerCepstrum")

        if params.quefrency_smoothing_s > 0 and params.smoothing_iterations > 0:
            power_cepstrum = call(
                power_cepstrum,
                "Smooth",
                params.quefrency_smoothing_s,
                params.smoothing_iterations,
            )

        cpps_db = call(
            power_cepstrum,
            "Get peak prominence",
            params.pitch_floor_hz,
            params.pitch_ceiling_hz,
            params.interpolation,
            params.trend_q_min_s,
            params.trend_q_max_s,
            params.trend_type,
            params.fit_method,
        )

        peak_q = call(
            power_cepstrum,
            "Get quefrency of peak",
            params.pitch_floor_hz,
            params.pitch_ceiling_hz,
            params.interpolation,
        )
        peak_freq = 1.0 / peak_q if peak_q and peak_q > 0 else float("nan")

        return float(cpps_db), float(peak_q), float(peak_freq), ";".join(warnings)

    except Exception as exc:
        warnings.append(f"praat_calculation_failed:{exc}")
        return float("nan"), float("nan"), float("nan"), ";".join(warnings)


def analyze_wav_file(
    wav_path: Path,
    root_path: Path,
    output_dir: Path,
    params: CPPSParameters,
) -> FileSummary:
    """Analyze one WAV file and write its frame-wise CSV."""
    sr, x = read_wav_as_mono_float(wav_path)
    duration_s = len(x) / sr if sr > 0 else float("nan")

    output_stem = make_output_stem(wav_path, root_path)
    frame_csv = output_dir / f"{output_stem}_cpps_framewise.csv"

    frame_results: list[FrameResult] = []

    for frame_index, (start_sample, frame) in enumerate(
        split_frames_include_last(x, params.frame_size, params.hop_size)
    ):
        end_sample = min(start_sample + len(frame), len(x))
        short_final = len(frame) < params.frame_size

        rms_value = rms_dbfs(frame)
        ac_score = autocorr_voicing_score(
            frame,
            sr,
            pitch_floor_hz=params.pitch_floor_hz,
            pitch_ceiling_hz=params.pitch_ceiling_hz,
        )

        voicing_warnings: list[str] = []
        if not np.isfinite(rms_value) or rms_value < params.rms_unvoiced_threshold_dbfs:
            voicing_warnings.append("low_rms")
        if not np.isfinite(ac_score) or ac_score < params.autocorr_voicing_threshold:
            voicing_warnings.append("low_periodicity")
        if short_final:
            voicing_warnings.append("short_final_frame")

        is_voiced = not any(w in voicing_warnings for w in ["low_rms", "low_periodicity"])

        cpps_db, peak_q, peak_freq, calc_warning = calculate_praat_cpps_for_frame(
            frame, sr, params
        )

        frame_results.append(
            FrameResult(
                frame_index=frame_index,
                start_sample=start_sample,
                end_sample=end_sample,
                frame_length_samples=len(frame),
                start_time_s=start_sample / sr,
                end_time_s=end_sample / sr,
                short_final_frame=short_final,
                rms_dbfs=rms_value,
                autocorr_voicing_score=ac_score,
                is_voiced=is_voiced,
                voicing_warning=";".join(voicing_warnings),
                cpps_db=cpps_db,
                peak_quefrency_s=peak_q,
                peak_frequency_hz=peak_freq,
                calculation_warning=calc_warning,
            )
        )

    if not frame_results:
        raise ValueError(f"No frames were produced for {wav_path}.")

    write_framewise_csv(frame_csv, frame_results)

    all_cpps = np.array([r.cpps_db for r in frame_results], dtype=np.float64)
    all_cpps = all_cpps[np.isfinite(all_cpps)]

    voiced_cpps = np.array(
        [r.cpps_db for r in frame_results if r.is_voiced and np.isfinite(r.cpps_db)],
        dtype=np.float64,
    )

    unvoiced_indices = [str(r.frame_index) for r in frame_results if not r.is_voiced]
    n_unvoiced = len(unvoiced_indices)
    n_voiced = len(frame_results) - n_unvoiced

    warning_parts = []
    if n_unvoiced > 0:
        warning_parts.append(f"unvoiced_or_low_periodicity_frames_detected:{n_unvoiced}")
    if len(all_cpps) < len(frame_results):
        warning_parts.append(f"nonfinite_cpps_frames:{len(frame_results) - len(all_cpps)}")

    return FileSummary(
        input_file=str(wav_path),
        output_frame_csv=str(frame_csv),
        sampling_rate_hz=sr,
        duration_s=duration_s,
        frame_size_samples=params.frame_size,
        hop_size_samples=params.hop_size,
        frame_duration_ms=1000.0 * params.frame_size / sr,
        n_frames_total=len(frame_results),
        n_voiced_frames=n_voiced,
        n_unvoiced_frames=n_unvoiced,
        unvoiced_frame_indices=";".join(unvoiced_indices),
        median_cpps_all_frames_db=safe_median(all_cpps),
        mean_cpps_all_frames_db=safe_mean(all_cpps),
        sd_cpps_all_frames_db=safe_sd(all_cpps),
        median_cpps_voiced_frames_db=safe_median_or_none(voiced_cpps),
        mean_cpps_voiced_frames_db=safe_mean_or_none(voiced_cpps),
        sd_cpps_voiced_frames_db=safe_sd_or_none(voiced_cpps),
        warning=";".join(warning_parts),
    )


def safe_median(x: np.ndarray) -> float:
    if len(x) == 0:
        return float("nan")
    return float(np.median(x))


def safe_mean(x: np.ndarray) -> float:
    if len(x) == 0:
        return float("nan")
    return float(np.mean(x))


def safe_sd(x: np.ndarray) -> float:
    if len(x) <= 1:
        return float("nan")
    return float(np.std(x, ddof=1))


def safe_median_or_none(x: np.ndarray) -> float | None:
    if len(x) == 0:
        return None
    return float(np.median(x))


def safe_mean_or_none(x: np.ndarray) -> float | None:
    if len(x) == 0:
        return None
    return float(np.mean(x))


def safe_sd_or_none(x: np.ndarray) -> float | None:
    if len(x) <= 1:
        return None
    return float(np.std(x, ddof=1))


def fmt_float(value: float | None, digits: int = 6) -> str:
    if value is None:
        return ""
    if not np.isfinite(value):
        return "nan"
    return f"{value:.{digits}f}"


def write_framewise_csv(csv_path: Path, rows: list[FrameResult]) -> None:
    """Write frame-wise CPPS results."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow([
            "frame_index",
            "start_sample",
            "end_sample",
            "frame_length_samples",
            "start_time_s",
            "end_time_s",
            "short_final_frame",
            "rms_dbfs",
            "autocorr_voicing_score",
            "is_voiced",
            "voicing_warning",
            "cpps_db",
            "peak_quefrency_s",
            "peak_frequency_hz",
            "calculation_warning",
        ])
        for r in rows:
            writer.writerow([
                r.frame_index,
                r.start_sample,
                r.end_sample,
                r.frame_length_samples,
                f"{r.start_time_s:.9f}",
                f"{r.end_time_s:.9f}",
                int(r.short_final_frame),
                fmt_float(r.rms_dbfs),
                fmt_float(r.autocorr_voicing_score),
                int(r.is_voiced),
                r.voicing_warning,
                fmt_float(r.cpps_db),
                fmt_float(r.peak_quefrency_s, 9),
                fmt_float(r.peak_frequency_hz),
                r.calculation_warning,
            ])


def write_summary_csv(csv_path: Path, summaries: list[FileSummary]) -> None:
    """Write summary CSV for all analyzed WAV files."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow([
            "input_file",
            "output_frame_csv",
            "sampling_rate_hz",
            "duration_s",
            "frame_size_samples",
            "hop_size_samples",
            "frame_duration_ms",
            "n_frames_total",
            "n_voiced_frames",
            "n_unvoiced_frames",
            "unvoiced_frame_indices",
            "median_cpps_all_frames_db",
            "mean_cpps_all_frames_db",
            "sd_cpps_all_frames_db",
            "median_cpps_voiced_frames_db",
            "mean_cpps_voiced_frames_db",
            "sd_cpps_voiced_frames_db",
            "warning",
        ])
        for s in summaries:
            writer.writerow([
                s.input_file,
                s.output_frame_csv,
                s.sampling_rate_hz,
                f"{s.duration_s:.9f}",
                s.frame_size_samples,
                s.hop_size_samples,
                f"{s.frame_duration_ms:.6f}",
                s.n_frames_total,
                s.n_voiced_frames,
                s.n_unvoiced_frames,
                s.unvoiced_frame_indices,
                fmt_float(s.median_cpps_all_frames_db),
                fmt_float(s.mean_cpps_all_frames_db),
                fmt_float(s.sd_cpps_all_frames_db),
                fmt_float(s.median_cpps_voiced_frames_db),
                fmt_float(s.mean_cpps_voiced_frames_db),
                fmt_float(s.sd_cpps_voiced_frames_db),
                s.warning,
            ])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate frame-wise CPPS and summarize each WAV file by median CPPS."
    )
    parser.add_argument(
        "input_path",
        nargs="?",
        default="input.wav",
        help="Input WAV file or directory containing WAV files. Default: input.wav",
    )
    parser.add_argument(
        "--output-dir",
        default="cpps_results",
        help="Directory for output CSV files. Default: cpps_results",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="When input_path is a directory, search recursively for WAV files.",
    )
    parser.add_argument(
        "--frame-size",
        type=int,
        default=2048,
        help="Frame size in samples. Default: 2048",
    )
    parser.add_argument(
        "--hop-size",
        type=int,
        default=2048,
        help="Hop size in samples. Default: 2048, i.e., non-overlapping frames",
    )
    parser.add_argument(
        "--pitch-floor",
        type=float,
        default=60.0,
        help="Pitch floor in Hz for cepstral peak search and voicing check. Default: 60",
    )
    parser.add_argument(
        "--pitch-ceiling",
        type=float,
        default=1000.0,
        help="Pitch ceiling in Hz for cepstral peak search and voicing check. Default: 1000",
    )
    parser.add_argument(
        "--quefrency-smoothing",
        type=float,
        default=0.0005,
        help="Praat PowerCepstrum smoothing window in quefrency seconds. Default: 0.0005",
    )
    parser.add_argument(
        "--smoothing-iterations",
        type=int,
        default=1,
        help="Number of Praat PowerCepstrum smoothing iterations. Default: 1",
    )
    parser.add_argument(
        "--interpolation",
        choices=["none", "parabolic", "cubic", "sinc70", "sinc700"],
        default="parabolic",
        help="Praat peak interpolation method. Default: parabolic",
    )
    parser.add_argument(
        "--trend-q-min",
        type=float,
        default=0.001,
        help="Lower quefrency limit for trend-line fitting in seconds. Default: 0.001",
    )
    parser.add_argument(
        "--trend-q-max",
        type=float,
        default=0.05,
        help="Upper quefrency limit for trend-line fitting in seconds. Default: 0.05",
    )
    parser.add_argument(
        "--trend-type",
        choices=["straight", "exponential"],
        default="straight",
        help="Praat trend type. Default: straight",
    )
    parser.add_argument(
        "--fit-method",
        choices=["least squares", "robust", "robust slow"],
        default="robust slow",
        help="Praat trend-line fit method. Default: robust slow",
    )
    parser.add_argument(
        "--rms-unvoiced-threshold-dbfs",
        type=float,
        default=-60.0,
        help="Frames below this RMS level are flagged as low_rms. Default: -60 dBFS",
    )
    parser.add_argument(
        "--autocorr-voicing-threshold",
        type=float,
        default=0.30,
        help="Frames below this autocorrelation peak are flagged as low_periodicity. Default: 0.30",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_path)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    params = CPPSParameters(
        frame_size=args.frame_size,
        hop_size=args.hop_size,
        pitch_floor_hz=args.pitch_floor,
        pitch_ceiling_hz=args.pitch_ceiling,
        quefrency_smoothing_s=args.quefrency_smoothing,
        smoothing_iterations=args.smoothing_iterations,
        interpolation=args.interpolation,
        trend_q_min_s=args.trend_q_min,
        trend_q_max_s=args.trend_q_max,
        trend_type=args.trend_type,
        fit_method=args.fit_method,
        rms_unvoiced_threshold_dbfs=args.rms_unvoiced_threshold_dbfs,
        autocorr_voicing_threshold=args.autocorr_voicing_threshold,
    )

    wav_files = iter_wav_files(input_path, recursive=args.recursive)
    summaries: list[FileSummary] = []

    print("Frame-wise CPPS calculation")
    print("---------------------------")
    print(f"Input path : {input_path}")
    print(f"Output dir : {output_dir}")
    print(f"WAV files  : {len(wav_files)}")
    print(f"Frame size : {params.frame_size} samples")
    print(f"Hop size   : {params.hop_size} samples")
    print()

    for wav_path in wav_files:
        summary = analyze_wav_file(wav_path, input_path, output_dir, params)
        summaries.append(summary)

        print(f"File: {wav_path}")
        print(f"  sampling rate             : {summary.sampling_rate_hz} Hz")
        print(f"  frame duration            : {summary.frame_duration_ms:.6f} ms")
        print(f"  median CPPS, all frames   : {fmt_float(summary.median_cpps_all_frames_db)} dB")
        print(f"  voiced frames             : {summary.n_voiced_frames}/{summary.n_frames_total}")
        if summary.n_unvoiced_frames > 0:
            print(f"  WARNING: {summary.n_unvoiced_frames} unvoiced/low-periodicity frame(s) detected.")
            print(f"           Frame indices: {summary.unvoiced_frame_indices}")
        print(f"  frame-wise CSV            : {summary.output_frame_csv}")
        print()

    summary_csv = output_dir / "cpps_summary.csv"
    write_summary_csv(summary_csv, summaries)

    print(f"Summary CSV: {summary_csv}")
    print("Primary value for each file: median_cpps_all_frames_db")


if __name__ == "__main__":
    main()
