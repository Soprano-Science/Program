#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calculate_spectral_centroid_framewise_median_iqr.py

Frame-wise Spectral Centroid calculation for WAV files.

Default method:
- Whole WAV file is analyzed.
- Frame size: 2048 samples.
- Hop size: 2048 samples, i.e., non-overlapping frames.
- FFT size: 2048 points.
- Final short frame is retained and zero-padded to 2048 samples.
- Hann window is applied to each actual frame before zero-padding.
- Spectral centroid is calculated for each frame.
- File-level representative value: median of finite frame-wise centroids.
- Variability: IQR = 75th percentile - 25th percentile.
- Likely unvoiced or low-periodicity frames are flagged but not excluded from
  the primary all-frame median/IQR.

Dependencies:
    pip install numpy scipy
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
from scipy.io import wavfile

EPS = 1e-20


def read_wav_as_mono_float(path: Path) -> tuple[int, np.ndarray]:
    """Read WAV and return sampling rate and mono float64 signal."""
    sr, x = wavfile.read(str(path))
    if x.ndim == 2:
        x = x.mean(axis=1)
    if np.issubdtype(x.dtype, np.integer):
        info = np.iinfo(x.dtype)
        scale = max(abs(info.min), abs(info.max))
        x = x.astype(np.float64) / scale
    else:
        x = x.astype(np.float64)
    if len(x) == 0:
        raise ValueError(f"Empty WAV file: {path}")
    return sr, x


def iter_wav_files(input_path: Path, recursive: bool) -> list[Path]:
    """Return one WAV file or all WAV files in a folder."""
    if input_path.is_file():
        if input_path.suffix.lower() != ".wav":
            raise ValueError(f"Input file is not a WAV file: {input_path}")
        return [input_path]
    if input_path.is_dir():
        pattern = "**/*.wav" if recursive else "*.wav"
        files = sorted(input_path.glob(pattern))
        if not files:
            raise FileNotFoundError(f"No WAV files found in: {input_path}")
        return files
    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def safe_stem(path: Path, root: Path) -> str:
    """Create safe output stem."""
    if root.is_dir():
        try:
            stem = str(path.relative_to(root).with_suffix(""))
        except ValueError:
            stem = path.stem
    else:
        stem = path.stem
    for old, new in [("\\", "__"), ("/", "__"), (":", "_"), (" ", "_")]:
        stem = stem.replace(old, new)
    return stem


def frame_generator(x: np.ndarray, frame_size: int, hop_size: int) -> Iterable[tuple[int, np.ndarray]]:
    """Yield frames and include the final short frame."""
    if frame_size <= 0 or hop_size <= 0:
        raise ValueError("frame_size and hop_size must be positive.")
    start = 0
    while start < len(x):
        frame = x[start : start + frame_size]
        if len(frame) == 0:
            break
        yield start, frame
        if start + frame_size >= len(x):
            break
        start += hop_size


def rms_dbfs(frame: np.ndarray) -> float:
    """Return RMS in dBFS for a frame scaled to approximately [-1, 1]."""
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
    """Return normalized autocorrelation peak used only as a voicing diagnostic."""
    if len(frame) < 4:
        return float("nan")
    y = frame.astype(np.float64) - float(np.mean(frame))
    if len(y) > 1:
        y = y * np.hanning(len(y))
    energy = float(np.dot(y, y))
    if energy <= EPS:
        return 0.0
    corr = np.correlate(y, y, mode="full")[len(y) - 1 :]
    if len(corr) < 2 or corr[0] <= EPS:
        return 0.0
    corr = corr / (corr[0] + EPS)
    min_lag = max(1, int(math.floor(sr / pitch_ceiling_hz)))
    max_lag = min(len(corr) - 1, int(math.ceil(sr / pitch_floor_hz)))
    if max_lag <= min_lag:
        return float("nan")
    return float(np.max(corr[min_lag : max_lag + 1]))


def spectral_centroid_one_frame(
    frame: np.ndarray,
    sr: int,
    frame_size: int,
    fft_size: int,
    weighting: str,
    min_frequency_hz: float,
    max_frequency_hz: Optional[float],
) -> tuple[float, str]:
    """Calculate one frame's spectral centroid in Hz."""
    warnings: list[str] = []
    if len(frame) == 0:
        return float("nan"), "empty_frame"
    if len(frame) < frame_size:
        warnings.append("zero_padded_short_frame")

    y = frame.astype(np.float64) - float(np.mean(frame))
    if len(y) > 1:
        y = y * np.hanning(len(y))

    if len(y) < fft_size:
        y_fft = np.zeros(fft_size, dtype=np.float64)
        y_fft[: len(y)] = y
    else:
        y_fft = y[:fft_size]
        if len(y) > fft_size:
            warnings.append("frame_truncated_to_fft_size")

    spectrum = np.fft.rfft(y_fft, n=fft_size)
    magnitude = np.abs(spectrum)
    if weighting == "magnitude":
        weights = magnitude
    elif weighting == "power":
        weights = magnitude ** 2
    else:
        raise ValueError("weighting must be either 'magnitude' or 'power'.")

    freqs = np.fft.rfftfreq(fft_size, d=1.0 / sr)
    nyquist = sr / 2.0
    max_freq = nyquist if max_frequency_hz is None else max_frequency_hz
    if min_frequency_hz < 0:
        raise ValueError("min_frequency_hz must be >= 0.")
    if max_freq > nyquist:
        raise ValueError(f"max_frequency_hz={max_freq} exceeds Nyquist frequency {nyquist}.")
    if max_freq <= min_frequency_hz:
        raise ValueError("max_frequency_hz must be greater than min_frequency_hz.")

    mask = (freqs >= min_frequency_hz) & (freqs <= max_freq)
    if not np.any(mask):
        return float("nan"), "; ".join(warnings + ["no_frequency_bins_in_selected_range"])

    w = weights[mask]
    f = freqs[mask]
    wsum = float(np.sum(w))
    if wsum <= EPS:
        return float("nan"), "; ".join(warnings + ["zero_or_near_zero_spectral_energy"])

    centroid = float(np.sum(f * w) / (wsum + EPS))
    return centroid, "; ".join(warnings)


def finite_array(values: Iterable[float]) -> np.ndarray:
    arr = np.array(list(values), dtype=np.float64)
    return arr[np.isfinite(arr)]


def stats(values: np.ndarray) -> dict[str, float]:
    """Return median, mean, sd, q25, q75, iqr, min, max."""
    if values.size == 0:
        return {k: float("nan") for k in ["median", "mean", "sd", "q25", "q75", "iqr", "min", "max"]}
    q25 = float(np.percentile(values, 25))
    q75 = float(np.percentile(values, 75))
    return {
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
        "q25": q25,
        "q75": q75,
        "iqr": q75 - q25,
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def fmt(value) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return f"{value:.9f}"
    return str(value)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k, "")) for k in fieldnames})


def process_file(args, wav_path: Path, root_path: Path) -> dict:
    sr, x = read_wav_as_mono_float(wav_path)
    nyquist = sr / 2.0
    max_frequency_used = nyquist if args.max_frequency is None else args.max_frequency

    output_stem = safe_stem(wav_path, root_path)
    framewise_csv = args.output_dir / f"{output_stem}_spectral_centroid_framewise.csv"

    frame_rows: list[dict] = []
    for idx, (start, frame) in enumerate(frame_generator(x, args.frame_size, args.hop_size)):
        end = start + len(frame)
        short_final = len(frame) < args.frame_size
        rdb = rms_dbfs(frame)
        voice_score = autocorr_voicing_score(frame, sr, args.pitch_floor_hz, args.pitch_ceiling_hz)
        voice_warnings: list[str] = []
        if short_final:
            voice_warnings.append("short_final_frame")
        if math.isnan(rdb) or rdb < args.rms_unvoiced_threshold_dbfs:
            voice_warnings.append("low_rms")
        if math.isnan(voice_score) or voice_score < args.autocorr_voicing_threshold:
            voice_warnings.append("low_periodicity")
        is_voiced = len([w for w in voice_warnings if w in {"low_rms", "low_periodicity"}]) == 0

        centroid, calc_warning = spectral_centroid_one_frame(
            frame=frame,
            sr=sr,
            frame_size=args.frame_size,
            fft_size=args.fft_size,
            weighting=args.weighting,
            min_frequency_hz=args.min_frequency,
            max_frequency_hz=args.max_frequency,
        )
        if short_final and "short_final_frame" not in calc_warning:
            calc_warning = "; ".join([w for w in [calc_warning, "short_final_frame"] if w])

        frame_rows.append({
            "frame_index": idx,
            "start_sample": start,
            "end_sample": end,
            "frame_length_samples": len(frame),
            "start_time_s": start / sr,
            "end_time_s": end / sr,
            "short_final_frame": short_final,
            "rms_dbfs": rdb,
            "autocorr_voicing_score": voice_score,
            "is_likely_voiced": is_voiced,
            "voicing_warning": "; ".join(voice_warnings),
            "spectral_centroid_hz": centroid,
            "spectral_centroid_khz": centroid / 1000.0 if math.isfinite(centroid) else float("nan"),
            "calculation_warning": calc_warning,
        })

    frame_fields = [
        "frame_index", "start_sample", "end_sample", "frame_length_samples",
        "start_time_s", "end_time_s", "short_final_frame", "rms_dbfs",
        "autocorr_voicing_score", "is_likely_voiced", "voicing_warning",
        "spectral_centroid_hz", "spectral_centroid_khz", "calculation_warning",
    ]
    write_csv(framewise_csv, frame_rows, frame_fields)

    all_vals = finite_array(row["spectral_centroid_hz"] for row in frame_rows)
    voiced_vals = finite_array(row["spectral_centroid_hz"] for row in frame_rows if row["is_likely_voiced"])
    st_all = stats(all_vals)
    st_voiced = stats(voiced_vals)
    unvoiced_indices = [str(row["frame_index"]) for row in frame_rows if not row["is_likely_voiced"]]
    invalid_count = len(frame_rows) - int(all_vals.size)

    warnings: list[str] = []
    if unvoiced_indices:
        warnings.append(f"likely_unvoiced_or_low_periodicity_frames_detected: {','.join(unvoiced_indices)}")
    if invalid_count:
        warnings.append(f"undefined_centroid_frames_excluded_from_statistics: {invalid_count}")
    if any(row["frame_length_samples"] < args.frame_size for row in frame_rows):
        warnings.append("short_final_frame_retained_and_zero_padded")

    return {
        "input_file": str(wav_path),
        "framewise_csv": str(framewise_csv),
        "sampling_rate_hz": sr,
        "duration_s": len(x) / sr,
        "frame_size_samples": args.frame_size,
        "hop_size_samples": args.hop_size,
        "fft_size_points": args.fft_size,
        "frame_duration_ms": args.frame_size / sr * 1000.0,
        "window": "hann",
        "weighting": args.weighting,
        "min_frequency_hz": args.min_frequency,
        "max_frequency_hz": max_frequency_used,
        "centroid_definition": "sum(frequency_bin_hz * spectral_weight) / sum(spectral_weight)",
        "n_frames_total": len(frame_rows),
        "n_valid_centroid_frames": int(all_vals.size),
        "n_invalid_centroid_frames": invalid_count,
        "n_likely_voiced_frames": len(frame_rows) - len(unvoiced_indices),
        "n_likely_unvoiced_frames": len(unvoiced_indices),
        "likely_unvoiced_frame_indices": ",".join(unvoiced_indices),
        "median_spectral_centroid_all_frames_hz": st_all["median"],
        "q25_spectral_centroid_all_frames_hz": st_all["q25"],
        "q75_spectral_centroid_all_frames_hz": st_all["q75"],
        "iqr_spectral_centroid_all_frames_hz": st_all["iqr"],
        "mean_spectral_centroid_all_frames_hz": st_all["mean"],
        "sd_spectral_centroid_all_frames_hz": st_all["sd"],
        "min_spectral_centroid_all_frames_hz": st_all["min"],
        "max_spectral_centroid_all_frames_hz": st_all["max"],
        "median_spectral_centroid_likely_voiced_frames_hz": st_voiced["median"],
        "iqr_spectral_centroid_likely_voiced_frames_hz": st_voiced["iqr"],
        "warning": "; ".join(warnings),
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Calculate frame-wise Spectral Centroid and summarize by median and IQR."
    )
    p.add_argument("input", type=Path, help="Input WAV file or folder containing WAV files.")
    p.add_argument("--output-dir", type=Path, default=Path("spectral_centroid_results"), help="Output directory.")
    p.add_argument("--recursive", action="store_true", help="Process WAV files recursively when input is a folder.")
    p.add_argument("--frame-size", type=int, default=2048, help="Frame size in samples. Default: 2048")
    p.add_argument("--hop-size", type=int, default=2048, help="Hop size in samples. Default: 2048")
    p.add_argument("--fft-size", type=int, default=2048, help="FFT size in points. Default: 2048")
    p.add_argument("--weighting", choices=["magnitude", "power"], default="magnitude", help="Spectral weighting. Default: magnitude")
    p.add_argument("--min-frequency", type=float, default=0.0, help="Minimum frequency in Hz. Default: 0")
    p.add_argument("--max-frequency", type=float, default=None, help="Maximum frequency in Hz. Default: Nyquist")
    p.add_argument("--rms-unvoiced-threshold-dbfs", type=float, default=-60.0, help="Low-RMS warning threshold. Default: -60 dBFS")
    p.add_argument("--autocorr-voicing-threshold", type=float, default=0.30, help="Low-periodicity warning threshold. Default: 0.30")
    p.add_argument("--pitch-floor-hz", type=float, default=60.0, help="Pitch floor for voicing diagnostic. Default: 60 Hz")
    p.add_argument("--pitch-ceiling-hz", type=float, default=1000.0, help="Pitch ceiling for voicing diagnostic. Default: 1000 Hz")
    return p


def main() -> None:
    args = build_parser().parse_args()
    if args.fft_size < args.frame_size:
        raise ValueError("fft_size must be greater than or equal to frame_size.")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    wavs = iter_wav_files(args.input, args.recursive)
    summaries: list[dict] = []
    for wav in wavs:
        print(f"Processing: {wav}")
        s = process_file(args, wav, args.input)
        summaries.append(s)
        print(
            f"  median={s['median_spectral_centroid_all_frames_hz']:.6f} Hz, "
            f"IQR={s['iqr_spectral_centroid_all_frames_hz']:.6f} Hz, "
            f"likely_unvoiced_frames={s['n_likely_unvoiced_frames']}"
        )

    summary_csv = args.output_dir / "spectral_centroid_summary.csv"
    summary_fields = [
        "input_file", "framewise_csv", "sampling_rate_hz", "duration_s",
        "frame_size_samples", "hop_size_samples", "fft_size_points", "frame_duration_ms",
        "window", "weighting", "min_frequency_hz", "max_frequency_hz", "centroid_definition",
        "n_frames_total", "n_valid_centroid_frames", "n_invalid_centroid_frames",
        "n_likely_voiced_frames", "n_likely_unvoiced_frames", "likely_unvoiced_frame_indices",
        "median_spectral_centroid_all_frames_hz", "q25_spectral_centroid_all_frames_hz",
        "q75_spectral_centroid_all_frames_hz", "iqr_spectral_centroid_all_frames_hz",
        "mean_spectral_centroid_all_frames_hz", "sd_spectral_centroid_all_frames_hz",
        "min_spectral_centroid_all_frames_hz", "max_spectral_centroid_all_frames_hz",
        "median_spectral_centroid_likely_voiced_frames_hz",
        "iqr_spectral_centroid_likely_voiced_frames_hz", "warning",
    ]
    write_csv(summary_csv, summaries, summary_fields)
    print(f"Summary CSV: {summary_csv}")
    print(f"Processed WAV files: {len(summaries)}")


if __name__ == "__main__":
    main()
