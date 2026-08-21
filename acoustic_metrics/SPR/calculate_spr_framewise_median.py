#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frame-wise Singing Power Ratio (SPR) calculation.

This script calculates SPR for one WAV file or for all WAV files in a folder.
It is designed for reproducible acoustic analysis of singing-voice recordings.

Default method:
    SPR is calculated frame by frame as the dB ratio between
    the greatest spectral peak in the 2-4 kHz band and the greatest
    spectral peak in the 0-2 kHz band:

        SPR_dB = 10 * log10(P_peak_2_4kHz / P_peak_0_2kHz)

    where P denotes spectral power. This is equivalent to the dB difference
    between peak amplitudes when power is magnitude squared.

Default analysis settings:
    - frame_size = 2048 samples
    - hop_size   = 2048 samples (non-overlapping frames)
    - n_fft      = 2048 points
    - window     = Hann window
    - final frame is retained even if shorter than 2048 samples
    - final short frame is zero-padded to 2048 samples for FFT computation
    - representative file-level SPR = median of frame-wise SPR values

A band-power variant is also available:
    --method bandpower

This variant calculates the dB ratio between summed spectral power in
2-4 kHz and summed spectral power in 0-2 kHz for each frame.

Example:
    python calculate_spr_framewise_median.py input.wav
    python calculate_spr_framewise_median.py ./audio --recursive
    python calculate_spr_framewise_median.py input.wav --method bandpower
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.io import wavfile


@dataclass(frozen=True)
class SPRConfig:
    """Analysis settings for frame-wise SPR."""

    frame_size: int = 2048
    hop_size: int = 2048
    n_fft: int = 2048
    low_band: tuple[float, float] = (0.0, 2000.0)
    high_band: tuple[float, float] = (2000.0, 4000.0)
    method: str = "peak"  # "peak" or "bandpower"
    min_rms: float = 0.0  # default: include all frames
    eps: float = 1e-20


def read_wav_as_mono_float(path: str | Path) -> tuple[int, np.ndarray]:
    """
    Read a WAV file and return sampling rate and mono float64 signal.

    Integer PCM is scaled to approximately [-1, 1].
    Stereo or multi-channel audio is converted to mono by averaging channels.
    """
    sr, x = wavfile.read(str(path))

    if x.ndim > 1:
        x = x.mean(axis=1)

    if np.issubdtype(x.dtype, np.integer):
        x = x.astype(np.float64) / np.iinfo(x.dtype).max
    else:
        x = x.astype(np.float64)

    return sr, x


def iter_frames_include_last(
    x: np.ndarray,
    frame_size: int,
    hop_size: int,
) -> Iterable[tuple[int, np.ndarray]]:
    """
    Yield frames and their start sample positions.

    The final frame is retained even if it is shorter than frame_size.
    """
    if frame_size <= 0:
        raise ValueError("frame_size must be positive.")
    if hop_size <= 0:
        raise ValueError("hop_size must be positive.")

    for start in range(0, len(x), hop_size):
        frame = x[start : start + frame_size]
        if len(frame) == 0:
            continue
        yield start, frame
        if start + frame_size >= len(x):
            break


def prepare_frame_for_fft(frame: np.ndarray, n_fft: int) -> tuple[np.ndarray, float]:
    """
    Remove frame-level DC offset, apply a Hann window, and zero-pad to n_fft.

    For the final short frame, the Hann window is generated using the actual
    frame length and the windowed frame is then zero-padded to n_fft.

    Returns:
        padded_frame, frame_rms_before_windowing
    """
    if len(frame) == 0:
        raise ValueError("Empty frame was given.")
    if len(frame) > n_fft:
        raise ValueError(
            f"Frame length ({len(frame)}) is longer than n_fft ({n_fft}). "
            "Use n_fft >= frame_size."
        )

    frame = frame.astype(np.float64, copy=False)
    frame = frame - np.mean(frame)
    rms = float(np.sqrt(np.mean(frame**2)))

    if len(frame) > 1:
        windowed = frame * np.hanning(len(frame))
    else:
        windowed = frame.copy()

    padded = np.zeros(n_fft, dtype=np.float64)
    padded[: len(windowed)] = windowed

    return padded, rms


def spectral_power(frame_for_fft: np.ndarray, sr: int, n_fft: int) -> tuple[np.ndarray, np.ndarray]:
    """Return frequency bins and power spectrum for a prepared frame."""
    spectrum = np.fft.rfft(frame_for_fft, n=n_fft)
    power = np.abs(spectrum) ** 2
    freqs = np.fft.rfftfreq(n_fft, d=1.0 / sr)
    return freqs, power


def calculate_spr_for_frame(
    frame: np.ndarray,
    sr: int,
    config: SPRConfig,
) -> dict[str, float | str]:
    """
    Calculate SPR for one frame.

    Default peak method:
        10 * log10(max_power_2_4kHz / max_power_0_2kHz)

    Bandpower method:
        10 * log10(sum_power_2_4kHz / sum_power_0_2kHz)
    """
    frame_for_fft, rms = prepare_frame_for_fft(frame, config.n_fft)

    if rms < config.min_rms:
        return {
            "spr_db": np.nan,
            "low_value": np.nan,
            "high_value": np.nan,
            "low_peak_freq_hz": np.nan,
            "high_peak_freq_hz": np.nan,
            "rms": rms,
            "status": "excluded_low_rms",
        }

    freqs, power = spectral_power(frame_for_fft, sr, config.n_fft)

    low_min, low_max = config.low_band
    high_min, high_max = config.high_band

    # DC offset is removed above; bin 0 is excluded to avoid DC dominance.
    low_mask = (freqs > low_min) & (freqs < low_max)
    high_mask = (freqs >= high_min) & (freqs <= high_max)

    if not np.any(low_mask):
        raise ValueError("No FFT bins found in the low band. Check sampling rate and n_fft.")
    if not np.any(high_mask):
        raise ValueError("No FFT bins found in the high band. Check sampling rate and n_fft.")

    low_power = power[low_mask]
    high_power = power[high_mask]
    low_freqs = freqs[low_mask]
    high_freqs = freqs[high_mask]

    if config.method == "peak":
        low_idx = int(np.argmax(low_power))
        high_idx = int(np.argmax(high_power))
        low_value = float(low_power[low_idx])
        high_value = float(high_power[high_idx])
        low_peak_freq_hz = float(low_freqs[low_idx])
        high_peak_freq_hz = float(high_freqs[high_idx])

    elif config.method == "bandpower":
        # Frequency-bin width. Multiplication by df cancels in the ratio,
        # but it keeps the values interpretable as approximate band power.
        df = float(freqs[1] - freqs[0])
        low_value = float(np.sum(low_power) * df)
        high_value = float(np.sum(high_power) * df)

        low_idx = int(np.argmax(low_power))
        high_idx = int(np.argmax(high_power))
        low_peak_freq_hz = float(low_freqs[low_idx])
        high_peak_freq_hz = float(high_freqs[high_idx])

    else:
        raise ValueError("method must be either 'peak' or 'bandpower'.")

    spr_db = 10.0 * np.log10((high_value + config.eps) / (low_value + config.eps))

    return {
        "spr_db": float(spr_db),
        "low_value": low_value,
        "high_value": high_value,
        "low_peak_freq_hz": low_peak_freq_hz,
        "high_peak_freq_hz": high_peak_freq_hz,
        "rms": rms,
        "status": "included",
    }


def analyze_wav(
    wav_path: str | Path,
    output_dir: str | Path,
    config: SPRConfig,
) -> dict[str, str | int | float]:
    """Analyze one WAV file and write a frame-wise CSV."""
    wav_path = Path(wav_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sr, x = read_wav_as_mono_float(wav_path)

    if sr < 8000:
        raise ValueError(
            f"Sampling rate is {sr} Hz. Analysis up to 4 kHz requires sr >= 8000 Hz."
        )
    if config.n_fft < config.frame_size:
        raise ValueError("n_fft must be greater than or equal to frame_size.")

    frame_rows: list[dict[str, str | int | float]] = []
    spr_values: list[float] = []

    for frame_index, (start_sample, frame) in enumerate(
        iter_frames_include_last(x, config.frame_size, config.hop_size)
    ):
        result = calculate_spr_for_frame(frame, sr, config)
        spr_db = result["spr_db"]

        if isinstance(spr_db, float) and not np.isnan(spr_db):
            spr_values.append(spr_db)

        end_sample = start_sample + len(frame)
        row = {
            "file": wav_path.name,
            "frame_index": frame_index,
            "start_sample": start_sample,
            "end_sample": end_sample,
            "frame_length_samples": len(frame),
            "start_time_s": start_sample / sr,
            "end_time_s": end_sample / sr,
            "rms": result["rms"],
            "spr_db": result["spr_db"],
            "low_band_value": result["low_value"],
            "high_band_value": result["high_value"],
            "low_peak_freq_hz": result["low_peak_freq_hz"],
            "high_peak_freq_hz": result["high_peak_freq_hz"],
            "status": result["status"],
        }
        frame_rows.append(row)

    if not spr_values:
        median_spr = np.nan
        mean_spr = np.nan
        sd_spr = np.nan
        min_spr = np.nan
        max_spr = np.nan
    else:
        spr_array = np.asarray(spr_values, dtype=np.float64)
        median_spr = float(np.median(spr_array))
        mean_spr = float(np.mean(spr_array))
        sd_spr = float(np.std(spr_array, ddof=1)) if len(spr_array) > 1 else 0.0
        min_spr = float(np.min(spr_array))
        max_spr = float(np.max(spr_array))

    frame_csv = output_dir / f"{wav_path.stem}_spr_framewise.csv"
    fieldnames = [
        "file",
        "frame_index",
        "start_sample",
        "end_sample",
        "frame_length_samples",
        "start_time_s",
        "end_time_s",
        "rms",
        "spr_db",
        "low_band_value",
        "high_band_value",
        "low_peak_freq_hz",
        "high_peak_freq_hz",
        "status",
    ]
    with frame_csv.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in frame_rows:
            writer.writerow(row)

    return {
        "file": wav_path.name,
        "path": str(wav_path),
        "sampling_rate_hz": sr,
        "duration_s": len(x) / sr,
        "frame_size_samples": config.frame_size,
        "hop_size_samples": config.hop_size,
        "n_fft": config.n_fft,
        "frame_duration_ms": 1000.0 * config.frame_size / sr,
        "method": config.method,
        "low_band_hz": f"{config.low_band[0]}-{config.low_band[1]}",
        "high_band_hz": f"{config.high_band[0]}-{config.high_band[1]}",
        "n_frames_total": len(frame_rows),
        "n_frames_included": len(spr_values),
        "median_spr_db": median_spr,
        "mean_spr_db": mean_spr,
        "sd_spr_db": sd_spr,
        "min_spr_db": min_spr,
        "max_spr_db": max_spr,
        "framewise_csv": str(frame_csv),
    }


def find_wav_files(input_path: str | Path, recursive: bool = False) -> list[Path]:
    """Return WAV files from a file or directory path."""
    input_path = Path(input_path)

    if input_path.is_file():
        if input_path.suffix.lower() != ".wav":
            raise ValueError(f"Input file is not a WAV file: {input_path}")
        return [input_path]

    if input_path.is_dir():
        pattern = "**/*.wav" if recursive else "*.wav"
        return sorted(input_path.glob(pattern))

    raise FileNotFoundError(f"Input path was not found: {input_path}")


def write_summary_csv(summary_rows: list[dict[str, str | int | float]], output_dir: str | Path) -> Path:
    """Write file-level summary CSV."""
    output_dir = Path(output_dir)
    summary_csv = output_dir / "spr_summary.csv"

    fieldnames = [
        "file",
        "path",
        "sampling_rate_hz",
        "duration_s",
        "frame_size_samples",
        "hop_size_samples",
        "n_fft",
        "frame_duration_ms",
        "method",
        "low_band_hz",
        "high_band_hz",
        "n_frames_total",
        "n_frames_included",
        "median_spr_db",
        "mean_spr_db",
        "sd_spr_db",
        "min_spr_db",
        "max_spr_db",
        "framewise_csv",
    ]

    with summary_csv.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    return summary_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate frame-wise Singing Power Ratio (SPR) and use the median "
            "as the representative file-level SPR value."
        )
    )
    parser.add_argument(
        "input",
        help="Input WAV file or folder containing WAV files.",
    )
    parser.add_argument(
        "--output-dir",
        default="spr_results",
        help="Output directory. Default: spr_results",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Process WAV files recursively when input is a folder.",
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
        help="Hop size in samples. Default: 2048, non-overlapping frames",
    )
    parser.add_argument(
        "--n-fft",
        type=int,
        default=2048,
        help="FFT size. Default: 2048",
    )
    parser.add_argument(
        "--method",
        choices=["peak", "bandpower"],
        default="peak",
        help=(
            "SPR calculation method. 'peak' uses the ratio of maximum spectral "
            "peaks in 2-4 kHz and 0-2 kHz. 'bandpower' uses summed spectral "
            "power in those bands. Default: peak"
        ),
    )
    parser.add_argument(
        "--min-rms",
        type=float,
        default=0.0,
        help=(
            "Optional RMS threshold for excluding very low-level frames. "
            "Default: 0.0, meaning all frames are included."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    config = SPRConfig(
        frame_size=args.frame_size,
        hop_size=args.hop_size,
        n_fft=args.n_fft,
        method=args.method,
        min_rms=args.min_rms,
    )

    wav_files = find_wav_files(args.input, recursive=args.recursive)
    if not wav_files:
        raise FileNotFoundError("No WAV files were found.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for wav_path in wav_files:
        print(f"Analyzing: {wav_path}")
        row = analyze_wav(wav_path, output_dir, config)
        summary_rows.append(row)
        print(f"  median SPR: {row['median_spr_db']:.6f} dB")

    summary_csv = write_summary_csv(summary_rows, output_dir)
    print("\nDone.")
    print(f"Summary CSV: {summary_csv}")
    print(f"Frame-wise CSV files are in: {output_dir}")


if __name__ == "__main__":
    main()
