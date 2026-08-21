#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Frame-wise SFR calculation based on the Excel definition.

This script extends the original one-frame SFR calculation to the whole WAV file.
It preserves the Excel-derived spectral calculation for each 2048-sample frame:

    SFR (%) = 100 * sum(amplitude in 2.4-4.0 kHz)
                    / sum(amplitude in 0-4.0 kHz)

Main features
-------------
- Processes a single WAV file or all WAV files in a folder.
- Uses non-overlapping 2048-sample frames.
- Includes the final frame in the frame-wise CSV even when it is shorter than
  2048 samples; the final short frame is zero-padded for FFT calculation.
- Detects silent or near-silent frames by RMS level.
- Excludes silent/near-silent frames and short final frames from the file-level
  median SFR calculation.
- Writes frame-wise CSV files and a summary CSV.

Requirements
------------
Python 3.10+ and NumPy.
The input WAV files are expected to be mono, uncompressed, 16-bit PCM, 44.1 kHz,
matching the original Excel definition.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import wave
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

N_FFT = 2048
EXPECTED_SR = 44100
LOW_HZ = 2400.0
HIGH_HZ = 4000.0
FULL_SCALE_INT16 = 32768.0


@dataclass
class FrameSFRResult:
    frame_index: int
    start_sample: int
    end_sample: int
    start_time_s: float
    end_time_s: float
    frame_length_samples: int
    is_short_frame: bool
    rms_dbfs: float
    max_abs_sample: float
    is_silent_frame: bool
    numerator_sum_2p4_4k: float | None
    denominator_sum_0_4k: float | None
    sfr_percent: float | None
    include_in_median: bool
    exclusion_reason: str
    calculation_warning: str


@dataclass
class FileSummary:
    file: str
    status: str
    error_message: str
    sample_rate_hz: int | None
    total_samples: int | None
    duration_s: float | None
    frame_size_samples: int
    frame_duration_ms: float | None
    n_frames_total: int
    n_full_frames: int
    n_short_frames: int
    short_frame_indices: str
    n_silent_frames: int
    silent_frame_indices: str
    n_frames_included_in_median: int
    median_sfr_percent: float | None
    mean_sfr_percent: float | None
    sd_sfr_percent: float | None
    min_sfr_percent: float | None
    max_sfr_percent: float | None
    framewise_csv: str
    warning: str


def read_mono_pcm16_wav(path: Path) -> tuple[np.ndarray, int]:
    """Read a mono, 16-bit PCM WAV without resampling or normalization."""
    with wave.open(str(path), "rb") as wf:
        channels = wf.getnchannels()
        sample_width = wf.getsampwidth()
        sample_rate = wf.getframerate()
        n_samples = wf.getnframes()
        compression = wf.getcomptype()
        raw = wf.readframes(n_samples)

    if channels != 1:
        raise ValueError(f"Mono WAV required; found {channels} channels.")
    if sample_width != 2:
        raise ValueError(f"16-bit PCM required; found {sample_width * 8}-bit samples.")
    if compression != "NONE":
        raise ValueError(f"Uncompressed PCM WAV required; found compression type {compression!r}.")
    if sample_rate != EXPECTED_SR:
        raise ValueError(
            f"The Excel definition assumes {EXPECTED_SR} Hz; found {sample_rate} Hz. "
            "No resampling is performed."
        )

    samples = np.frombuffer(raw, dtype="<i2").astype(np.float64)
    return samples, sample_rate


def periodic_hann_window(n_fft: int = N_FFT) -> np.ndarray:
    """Return the exact periodic Hann window used in the Excel definition."""
    n = np.arange(n_fft, dtype=np.float64)
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * n / n_fft))


def rms_dbfs(frame: np.ndarray) -> float:
    """Return RMS level in dBFS, using 32768 as full-scale for int16 PCM."""
    if frame.size == 0:
        return -math.inf
    rms = float(np.sqrt(np.mean(frame.astype(np.float64) ** 2)))
    if rms <= 0.0:
        return -math.inf
    return 20.0 * math.log10(rms / FULL_SCALE_INT16)


def calculate_sfr_for_frame(
    frame: np.ndarray,
    sample_rate: int,
    n_fft: int = N_FFT,
    window: np.ndarray | None = None,
) -> tuple[float | None, float | None, float | None, str]:
    """Calculate SFR for one frame using the Excel-derived formula.

    If the frame is shorter than n_fft, it is zero-padded for FFT calculation.
    The caller should decide whether the short frame is included in file-level
    statistics. This script excludes short frames from the median by default.

    Returns
    -------
    sfr_percent, numerator_sum, denominator_sum, calculation_warning
    """
    calculation_warning = ""

    if frame.size == 0:
        return None, None, None, "empty_frame"

    if frame.size < n_fft:
        padded = np.zeros(n_fft, dtype=np.float64)
        padded[: frame.size] = frame
        frame_for_fft = padded
        calculation_warning = "zero_padded_short_frame"
    elif frame.size == n_fft:
        frame_for_fft = frame.astype(np.float64, copy=False)
    else:
        # This should not occur when frames are generated by this script.
        frame_for_fft = frame[:n_fft].astype(np.float64, copy=False)
        calculation_warning = "frame_truncated_to_n_fft"

    if window is None:
        window = periodic_hann_window(n_fft)

    spectrum = np.fft.fft(frame_for_fft * window, n=n_fft)

    # Reproduce the Excel columns used in the original script:
    # H = 10*LOG10((Re^2+Im^2)/2048) - 20*LOG10(2^16)
    # J = 10^(H/20)
    magnitude_squared = spectrum.real**2 + spectrum.imag**2
    with np.errstate(divide="ignore"):
        spectrum_db = (
            10.0 * np.log10(magnitude_squared / n_fft)
            - 20.0 * np.log10(2.0**16)
        )
    amplitude = np.power(10.0, spectrum_db / 20.0)

    frequencies = np.arange(n_fft, dtype=np.float64) * sample_rate / n_fft
    denominator_mask = (frequencies >= 0.0) & (frequencies <= HIGH_HZ)
    numerator_mask = (frequencies >= LOW_HZ) & (frequencies <= HIGH_HZ)

    denominator = float(np.sum(amplitude[denominator_mask]))
    numerator = float(np.sum(amplitude[numerator_mask]))

    if denominator <= 0.0 or not np.isfinite(denominator):
        return None, numerator, denominator, _join_warnings(calculation_warning, "zero_or_invalid_denominator")

    sfr_percent = 100.0 * numerator / denominator
    if not np.isfinite(sfr_percent):
        return None, numerator, denominator, _join_warnings(calculation_warning, "invalid_sfr")

    return float(sfr_percent), numerator, denominator, calculation_warning


def _join_warnings(*parts: str) -> str:
    return "; ".join([p for p in parts if p])


def iter_frames_include_last(samples: np.ndarray, frame_size: int = N_FFT) -> Iterable[tuple[int, int, np.ndarray]]:
    """Yield non-overlapping frames and retain the final short frame."""
    if samples.size == 0:
        return
    frame_index = 0
    for start in range(0, samples.size, frame_size):
        end = min(start + frame_size, samples.size)
        yield frame_index, start, samples[start:end]
        frame_index += 1


def process_wav_file(
    wav_path: Path,
    output_dir: Path,
    silence_threshold_dbfs: float,
    frame_size: int = N_FFT,
) -> FileSummary:
    """Process one WAV file and write its frame-wise CSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    framewise_csv = output_dir / f"{wav_path.stem}_sfr_framewise.csv"

    try:
        samples, sample_rate = read_mono_pcm16_wav(wav_path)
        duration_s = samples.size / sample_rate if sample_rate else None
        frame_duration_ms = frame_size / sample_rate * 1000.0
        window = periodic_hann_window(frame_size)

        frame_results: list[FrameSFRResult] = []

        for frame_index, start_sample, frame in iter_frames_include_last(samples, frame_size):
            frame_len = int(frame.size)
            end_sample = start_sample + frame_len
            is_short = frame_len < frame_size
            rms = rms_dbfs(frame)
            max_abs = float(np.max(np.abs(frame))) if frame_len > 0 else 0.0
            is_silent = bool(rms <= silence_threshold_dbfs or max_abs == 0.0)

            sfr, numerator, denominator, calc_warning = calculate_sfr_for_frame(
                frame=frame,
                sample_rate=sample_rate,
                n_fft=frame_size,
                window=window,
            )

            exclusion_reasons: list[str] = []
            if is_short:
                exclusion_reasons.append("short_final_frame_less_than_2048_samples")
            if is_silent:
                exclusion_reasons.append(f"silent_or_near_silent_rms_le_{silence_threshold_dbfs:g}_dBFS")
            if sfr is None:
                exclusion_reasons.append("sfr_not_calculated")

            include_in_median = len(exclusion_reasons) == 0

            frame_results.append(
                FrameSFRResult(
                    frame_index=frame_index,
                    start_sample=start_sample,
                    end_sample=end_sample,
                    start_time_s=start_sample / sample_rate,
                    end_time_s=end_sample / sample_rate,
                    frame_length_samples=frame_len,
                    is_short_frame=is_short,
                    rms_dbfs=rms,
                    max_abs_sample=max_abs,
                    is_silent_frame=is_silent,
                    numerator_sum_2p4_4k=numerator,
                    denominator_sum_0_4k=denominator,
                    sfr_percent=sfr,
                    include_in_median=include_in_median,
                    exclusion_reason="; ".join(exclusion_reasons),
                    calculation_warning=calc_warning,
                )
            )

        write_framewise_csv(framewise_csv, wav_path, frame_results)

        included_values = np.array(
            [r.sfr_percent for r in frame_results if r.include_in_median and r.sfr_percent is not None],
            dtype=np.float64,
        )

        short_indices = [str(r.frame_index) for r in frame_results if r.is_short_frame]
        silent_indices = [str(r.frame_index) for r in frame_results if r.is_silent_frame]

        warnings: list[str] = []
        if short_indices:
            warnings.append("short_final_frame_excluded_from_median")
        if silent_indices:
            warnings.append("silent_or_near_silent_frames_excluded_from_median")
        if included_values.size == 0:
            warnings.append("no_valid_frames_for_median")

        return FileSummary(
            file=str(wav_path),
            status="ok",
            error_message="",
            sample_rate_hz=sample_rate,
            total_samples=int(samples.size),
            duration_s=duration_s,
            frame_size_samples=frame_size,
            frame_duration_ms=frame_duration_ms,
            n_frames_total=len(frame_results),
            n_full_frames=sum(1 for r in frame_results if not r.is_short_frame),
            n_short_frames=sum(1 for r in frame_results if r.is_short_frame),
            short_frame_indices=";".join(short_indices),
            n_silent_frames=sum(1 for r in frame_results if r.is_silent_frame),
            silent_frame_indices=";".join(silent_indices),
            n_frames_included_in_median=int(included_values.size),
            median_sfr_percent=float(np.median(included_values)) if included_values.size else None,
            mean_sfr_percent=float(np.mean(included_values)) if included_values.size else None,
            sd_sfr_percent=float(np.std(included_values, ddof=1)) if included_values.size > 1 else (0.0 if included_values.size == 1 else None),
            min_sfr_percent=float(np.min(included_values)) if included_values.size else None,
            max_sfr_percent=float(np.max(included_values)) if included_values.size else None,
            framewise_csv=str(framewise_csv),
            warning="; ".join(warnings),
        )

    except Exception as exc:  # keep folder batch processing robust
        return FileSummary(
            file=str(wav_path),
            status="error",
            error_message=str(exc),
            sample_rate_hz=None,
            total_samples=None,
            duration_s=None,
            frame_size_samples=frame_size,
            frame_duration_ms=None,
            n_frames_total=0,
            n_full_frames=0,
            n_short_frames=0,
            short_frame_indices="",
            n_silent_frames=0,
            silent_frame_indices="",
            n_frames_included_in_median=0,
            median_sfr_percent=None,
            mean_sfr_percent=None,
            sd_sfr_percent=None,
            min_sfr_percent=None,
            max_sfr_percent=None,
            framewise_csv="",
            warning="processing_failed",
        )


def format_float(value: float | None, digits: int = 12) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        if math.isinf(value) and value < 0:
            return "-inf"
        if math.isinf(value) and value > 0:
            return "inf"
        return "nan"
    return f"{value:.{digits}f}"


def write_framewise_csv(csv_path: Path, source_wav: Path, frame_results: list[FrameSFRResult]) -> None:
    """Write one frame-wise CSV file."""
    fieldnames = [
        "source_file",
        "frame_index",
        "start_sample",
        "end_sample",
        "start_time_s",
        "end_time_s",
        "frame_length_samples",
        "is_short_frame",
        "rms_dbfs",
        "max_abs_sample",
        "is_silent_frame",
        "numerator_sum_2p4_4k",
        "denominator_sum_0_4k",
        "sfr_percent",
        "include_in_median",
        "exclusion_reason",
        "calculation_warning",
    ]
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in frame_results:
            writer.writerow(
                {
                    "source_file": source_wav.name,
                    "frame_index": r.frame_index,
                    "start_sample": r.start_sample,
                    "end_sample": r.end_sample,
                    "start_time_s": format_float(r.start_time_s, 9),
                    "end_time_s": format_float(r.end_time_s, 9),
                    "frame_length_samples": r.frame_length_samples,
                    "is_short_frame": int(r.is_short_frame),
                    "rms_dbfs": format_float(r.rms_dbfs, 6),
                    "max_abs_sample": format_float(r.max_abs_sample, 6),
                    "is_silent_frame": int(r.is_silent_frame),
                    "numerator_sum_2p4_4k": format_float(r.numerator_sum_2p4_4k, 15),
                    "denominator_sum_0_4k": format_float(r.denominator_sum_0_4k, 15),
                    "sfr_percent": format_float(r.sfr_percent, 15),
                    "include_in_median": int(r.include_in_median),
                    "exclusion_reason": r.exclusion_reason,
                    "calculation_warning": r.calculation_warning,
                }
            )


def write_summary_csv(csv_path: Path, summaries: list[FileSummary]) -> None:
    """Write the folder-level summary CSV."""
    fieldnames = [
        "file",
        "status",
        "error_message",
        "sample_rate_hz",
        "total_samples",
        "duration_s",
        "frame_size_samples",
        "frame_duration_ms",
        "n_frames_total",
        "n_full_frames",
        "n_short_frames",
        "short_frame_indices",
        "n_silent_frames",
        "silent_frame_indices",
        "n_frames_included_in_median",
        "median_sfr_percent",
        "mean_sfr_percent",
        "sd_sfr_percent",
        "min_sfr_percent",
        "max_sfr_percent",
        "framewise_csv",
        "warning",
    ]
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in summaries:
            writer.writerow(
                {
                    "file": s.file,
                    "status": s.status,
                    "error_message": s.error_message,
                    "sample_rate_hz": "" if s.sample_rate_hz is None else s.sample_rate_hz,
                    "total_samples": "" if s.total_samples is None else s.total_samples,
                    "duration_s": format_float(s.duration_s, 9),
                    "frame_size_samples": s.frame_size_samples,
                    "frame_duration_ms": format_float(s.frame_duration_ms, 6),
                    "n_frames_total": s.n_frames_total,
                    "n_full_frames": s.n_full_frames,
                    "n_short_frames": s.n_short_frames,
                    "short_frame_indices": s.short_frame_indices,
                    "n_silent_frames": s.n_silent_frames,
                    "silent_frame_indices": s.silent_frame_indices,
                    "n_frames_included_in_median": s.n_frames_included_in_median,
                    "median_sfr_percent": format_float(s.median_sfr_percent, 15),
                    "mean_sfr_percent": format_float(s.mean_sfr_percent, 15),
                    "sd_sfr_percent": format_float(s.sd_sfr_percent, 15),
                    "min_sfr_percent": format_float(s.min_sfr_percent, 15),
                    "max_sfr_percent": format_float(s.max_sfr_percent, 15),
                    "framewise_csv": s.framewise_csv,
                    "warning": s.warning,
                }
            )


def collect_wav_files(input_path: Path, recursive: bool = False) -> list[Path]:
    """Collect WAV files from a file or folder path."""
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        pattern = "**/*.wav" if recursive else "*.wav"
        return sorted(input_path.glob(pattern))
    raise FileNotFoundError(f"Input path not found: {input_path}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate frame-wise Excel-defined SFR for one WAV file or a folder of WAV files. "
            "The file-level SFR is the median of valid frame-wise SFR values."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input WAV file or folder containing WAV files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("sfr_results"),
        help="Output directory. Default: sfr_results",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="When input is a folder, process WAV files recursively.",
    )
    parser.add_argument(
        "--silence-threshold-dbfs",
        type=float,
        default=-60.0,
        help=(
            "Frames with RMS at or below this level are marked as silent/near-silent "
            "and excluded from the median. Default: -60 dBFS."
        ),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    wav_files = collect_wav_files(args.input, recursive=args.recursive)

    if not wav_files:
        print(f"No WAV files found in {args.input}", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[FileSummary] = []

    for wav_path in wav_files:
        print(f"Processing: {wav_path}")
        summary = process_wav_file(
            wav_path=wav_path,
            output_dir=args.output_dir,
            silence_threshold_dbfs=args.silence_threshold_dbfs,
            frame_size=N_FFT,
        )
        summaries.append(summary)
        if summary.status == "ok":
            median_text = format_float(summary.median_sfr_percent, 6)
            print(f"  median SFR: {median_text} %")
            if summary.warning:
                print(f"  warning: {summary.warning}")
        else:
            print(f"  ERROR: {summary.error_message}", file=sys.stderr)

    summary_csv = args.output_dir / "sfr_summary.csv"
    write_summary_csv(summary_csv, summaries)
    print(f"\nSummary CSV: {summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
