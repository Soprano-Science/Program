#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calculate_alpha_ratio_framewise_median.py

Frame-wise alpha ratio calculation for WAV files.

This script calculates alpha ratio for each frame of a WAV file and uses
the median of the frame-wise values as the representative alpha ratio for
that audio file.

Default method
--------------
- Input frame length: 2048 samples
- Hop size: 2048 samples, i.e., non-overlapping frames
- The final frame is retained even when shorter than 2048 samples.
- For FFT calculation, the final short frame is zero-padded to 2048 samples.
- Window: Hann window applied to each actual frame before zero-padding
- Low-frequency band:  50-1000 Hz
- High-frequency band: 1000-5000 Hz
- Alpha ratio definition:

      alpha_ratio_db = 10 * log10(E_50_1000 / E_1000_5000)

  where E_50_1000 and E_1000_5000 are summed spectral power values
  in the corresponding frequency bands.

Interpretation of this sign convention
--------------------------------------
With the default definition above, a smaller alpha ratio value means that
energy in the 1000-5000 Hz band is relatively stronger than energy in the
50-1000 Hz band. If your previous implementation used the inverse ratio,
use the column named "alpha_ratio_db_high_over_low_check" for checking the
sign convention.

Dependencies
------------
    pip install numpy scipy

Basic usage
-----------
    python calculate_alpha_ratio_framewise_median.py input.wav

This creates:
    alpha_ratio_results/input_alpha_ratio_framewise.csv
    alpha_ratio_results/alpha_ratio_summary.csv

Folder usage
------------
    python calculate_alpha_ratio_framewise_median.py ./wav_files --recursive

Author note
-----------
This implementation is intended for research workflows in which an acoustic
feature is calculated frame by frame and summarized by the median. It is not
an LTAS implementation. In a manuscript or README, describe it as a
"frame-wise alpha ratio summarized by the median".
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
from scipy.io import wavfile


@dataclass
class FrameResult:
    """Frame-wise result for one frame."""

    frame_index: int
    start_sample: int
    end_sample: int
    frame_length_samples: int
    start_time_s: float
    end_time_s: float
    rms: float
    energy_50_1000: float
    energy_1000_5000: float
    alpha_ratio_db_low_over_high: float
    alpha_ratio_db_high_over_low_check: float


@dataclass
class FileSummary:
    """Summary result for one WAV file."""

    input_file: str
    sampling_rate_hz: int
    duration_s: float
    frame_size_samples: int
    hop_size_samples: int
    n_frames: int
    n_short_frames: int
    low_band_hz: str
    high_band_hz: str
    alpha_definition: str
    median_alpha_ratio_db: float
    mean_alpha_ratio_db: float
    sd_alpha_ratio_db: float
    min_alpha_ratio_db: float
    max_alpha_ratio_db: float
    framewise_csv: str


def read_wav_as_mono_float(wav_path: Path) -> Tuple[int, np.ndarray]:
    """
    Read a WAV file and return sampling rate and mono float64 signal.

    Stereo files are converted to mono by averaging channels.
    Integer PCM data are scaled to approximately [-1, 1].
    """
    sr, x = wavfile.read(str(wav_path))

    if x.size == 0:
        raise ValueError(f"Empty WAV file: {wav_path}")

    # Convert multi-channel to mono by averaging channels.
    if x.ndim == 2:
        x = x.mean(axis=1)
    elif x.ndim != 1:
        raise ValueError(f"Unsupported WAV shape {x.shape}: {wav_path}")

    # Normalize integer PCM data.
    if np.issubdtype(x.dtype, np.integer):
        info = np.iinfo(x.dtype)
        x = x.astype(np.float64)

        if np.issubdtype(np.dtype(info.dtype), np.unsignedinteger):
            # Example: uint8 WAV. Center around zero.
            midpoint = (info.max + info.min) / 2.0
            scale = (info.max - info.min + 1) / 2.0
            x = (x - midpoint) / scale
        else:
            # Example: int16 WAV. Use the larger absolute bound.
            scale = max(abs(info.min), abs(info.max))
            x = x / scale
    else:
        x = x.astype(np.float64)

    return sr, x


def iter_frames_include_last(
    x: np.ndarray,
    frame_size: int,
    hop_size: int,
) -> Iterable[Tuple[int, int, np.ndarray]]:
    """
    Yield frames while retaining the final short frame.

    Returns tuples:
        start_sample, end_sample, frame
    """
    if frame_size <= 0:
        raise ValueError("frame_size must be positive.")
    if hop_size <= 0:
        raise ValueError("hop_size must be positive.")

    n = len(x)
    start = 0

    while start < n:
        end = min(start + frame_size, n)
        frame = x[start:end]
        yield start, end, frame
        start += hop_size


def calculate_alpha_ratio_for_frame(
    frame: np.ndarray,
    sr: int,
    n_fft: int,
    low_band: Tuple[float, float],
    high_band: Tuple[float, float],
    eps: float,
    remove_dc: bool,
) -> Tuple[float, float, float, float, float]:
    """
    Calculate alpha ratio for one frame.

    The actual frame may be shorter than n_fft. The frame is windowed using
    its actual length and then zero-padded to n_fft.

    Returns:
        alpha_low_over_high_db,
        alpha_high_over_low_db,
        energy_low,
        energy_high,
        rms
    """
    if len(frame) == 0:
        raise ValueError("Empty frame was given.")

    frame = frame.astype(np.float64, copy=True)

    if remove_dc:
        frame -= np.mean(frame)

    rms = float(np.sqrt(np.mean(frame ** 2)))

    # Apply a Hann window to the actual frame length.
    # For a single-sample frame, windowing is skipped.
    if len(frame) > 1:
        frame *= np.hanning(len(frame))

    # Zero-pad or trim to n_fft.
    if len(frame) < n_fft:
        padded = np.zeros(n_fft, dtype=np.float64)
        padded[: len(frame)] = frame
    else:
        padded = frame[:n_fft]

    spectrum = np.fft.rfft(padded, n=n_fft)
    power = np.abs(spectrum) ** 2

    freqs = np.fft.rfftfreq(n_fft, d=1.0 / sr)
    df = freqs[1] - freqs[0]

    low_min, low_max = low_band
    high_min, high_max = high_band

    # Put 1000 Hz into the high band, not both bands.
    low_mask = (freqs >= low_min) & (freqs < low_max)
    high_mask = (freqs >= high_min) & (freqs <= high_max)

    if not np.any(low_mask):
        raise ValueError(
            f"No FFT bins found in low band {low_band}. "
            f"Check sampling rate={sr} and n_fft={n_fft}."
        )
    if not np.any(high_mask):
        raise ValueError(
            f"No FFT bins found in high band {high_band}. "
            f"Check sampling rate={sr} and n_fft={n_fft}."
        )

    energy_low = float(power[low_mask].sum() * df)
    energy_high = float(power[high_mask].sum() * df)

    alpha_low_over_high = 10.0 * math.log10((energy_low + eps) / (energy_high + eps))
    alpha_high_over_low = 10.0 * math.log10((energy_high + eps) / (energy_low + eps))

    return alpha_low_over_high, alpha_high_over_low, energy_low, energy_high, rms


def calculate_file(
    wav_path: Path,
    output_dir: Path,
    frame_size: int,
    hop_size: int,
    low_band: Tuple[float, float],
    high_band: Tuple[float, float],
    eps: float,
    remove_dc: bool,
    save_framewise_csv: bool,
) -> FileSummary:
    """Calculate frame-wise alpha ratios for one WAV file."""
    sr, x = read_wav_as_mono_float(wav_path)

    if sr < 2 * high_band[1]:
        raise ValueError(
            f"Sampling rate is {sr} Hz, but the high band extends to "
            f"{high_band[1]} Hz. A sampling rate of at least "
            f"{2 * high_band[1]} Hz is required. File: {wav_path}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    frame_results: List[FrameResult] = []

    for frame_index, (start_sample, end_sample, frame) in enumerate(
        iter_frames_include_last(x, frame_size=frame_size, hop_size=hop_size)
    ):
        alpha_low_high, alpha_high_low, energy_low, energy_high, rms = calculate_alpha_ratio_for_frame(
            frame=frame,
            sr=sr,
            n_fft=frame_size,
            low_band=low_band,
            high_band=high_band,
            eps=eps,
            remove_dc=remove_dc,
        )

        frame_results.append(
            FrameResult(
                frame_index=frame_index,
                start_sample=start_sample,
                end_sample=end_sample,
                frame_length_samples=len(frame),
                start_time_s=start_sample / sr,
                end_time_s=end_sample / sr,
                rms=rms,
                energy_50_1000=energy_low,
                energy_1000_5000=energy_high,
                alpha_ratio_db_low_over_high=alpha_low_high,
                alpha_ratio_db_high_over_low_check=alpha_high_low,
            )
        )

    if not frame_results:
        raise ValueError(f"No frames were generated for file: {wav_path}")

    alpha_values = np.array(
        [r.alpha_ratio_db_low_over_high for r in frame_results],
        dtype=np.float64,
    )

    framewise_csv_path = output_dir / f"{wav_path.stem}_alpha_ratio_framewise.csv"

    if save_framewise_csv:
        write_framewise_csv(framewise_csv_path, frame_results)
    else:
        framewise_csv_path = Path("")

    n_short_frames = sum(1 for r in frame_results if r.frame_length_samples < frame_size)

    sd = float(np.std(alpha_values, ddof=1)) if len(alpha_values) > 1 else 0.0

    return FileSummary(
        input_file=str(wav_path),
        sampling_rate_hz=sr,
        duration_s=len(x) / sr,
        frame_size_samples=frame_size,
        hop_size_samples=hop_size,
        n_frames=len(frame_results),
        n_short_frames=n_short_frames,
        low_band_hz=f"{low_band[0]:g}-{low_band[1]:g}",
        high_band_hz=f"{high_band[0]:g}-{high_band[1]:g}",
        alpha_definition="10log10(E_50_1000/E_1000_5000)",
        median_alpha_ratio_db=float(np.median(alpha_values)),
        mean_alpha_ratio_db=float(np.mean(alpha_values)),
        sd_alpha_ratio_db=sd,
        min_alpha_ratio_db=float(np.min(alpha_values)),
        max_alpha_ratio_db=float(np.max(alpha_values)),
        framewise_csv=str(framewise_csv_path),
    )


def write_framewise_csv(path: Path, frame_results: List[FrameResult]) -> None:
    """Write frame-wise results to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8-sig") as f:
        fieldnames = list(asdict(frame_results[0]).keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in frame_results:
            d = asdict(row)
            for key, value in d.items():
                if isinstance(value, float):
                    d[key] = f"{value:.10g}"
            writer.writerow(d)


def write_summary_csv(path: Path, summaries: List[FileSummary]) -> None:
    """Write file-level summary results to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8-sig") as f:
        fieldnames = list(asdict(summaries[0]).keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for summary in summaries:
            d = asdict(summary)
            for key, value in d.items():
                if isinstance(value, float):
                    d[key] = f"{value:.10g}"
            writer.writerow(d)


def collect_wav_files(input_path: Path, recursive: bool) -> List[Path]:
    """Collect WAV files from a file or directory path."""
    if input_path.is_file():
        if input_path.suffix.lower() != ".wav":
            raise ValueError(f"Input file is not a .wav file: {input_path}")
        return [input_path]

    if input_path.is_dir():
        pattern = "**/*.wav" if recursive else "*.wav"
        wav_files = sorted(input_path.glob(pattern))
        if not wav_files:
            raise FileNotFoundError(f"No .wav files found in directory: {input_path}")
        return wav_files

    raise FileNotFoundError(f"Input path was not found: {input_path}")


def build_arg_parser() -> argparse.ArgumentParser:
    """Build command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate frame-wise alpha ratio for WAV files and summarize "
            "each file by the median of frame-wise values."
        )
    )

    parser.add_argument(
        "input_path",
        nargs="?",
        default="input.wav",
        help="Input WAV file or directory. Default: input.wav",
    )
    parser.add_argument(
        "--output_dir",
        default="alpha_ratio_results",
        help="Output directory. Default: alpha_ratio_results",
    )
    parser.add_argument(
        "--summary_csv",
        default=None,
        help=(
            "Summary CSV path. Default: <output_dir>/alpha_ratio_summary.csv"
        ),
    )
    parser.add_argument(
        "--frame_size",
        type=int,
        default=2048,
        help="Frame size in samples. Default: 2048",
    )
    parser.add_argument(
        "--hop_size",
        type=int,
        default=2048,
        help=(
            "Hop size in samples. Default: 2048, i.e., non-overlapping frames. "
            "Use 1024 for 50%% overlap if needed."
        ),
    )
    parser.add_argument(
        "--low_band",
        nargs=2,
        type=float,
        default=(50.0, 1000.0),
        metavar=("LOW_MIN", "LOW_MAX"),
        help="Low band in Hz. Default: 50 1000",
    )
    parser.add_argument(
        "--high_band",
        nargs=2,
        type=float,
        default=(1000.0, 5000.0),
        metavar=("HIGH_MIN", "HIGH_MAX"),
        help="High band in Hz. Default: 1000 5000",
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=1e-20,
        help="Small constant to avoid division by zero. Default: 1e-20",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="If input_path is a directory, search WAV files recursively.",
    )
    parser.add_argument(
        "--keep_dc",
        action="store_true",
        help="Do not remove the DC offset from each frame. Default: DC is removed.",
    )
    parser.add_argument(
        "--no_framewise_csv",
        action="store_true",
        help="Do not save per-frame CSV files. Summary CSV is still saved.",
    )

    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    input_path = Path(args.input_path)
    output_dir = Path(args.output_dir)
    summary_csv = Path(args.summary_csv) if args.summary_csv else output_dir / "alpha_ratio_summary.csv"

    low_band = (float(args.low_band[0]), float(args.low_band[1]))
    high_band = (float(args.high_band[0]), float(args.high_band[1]))

    if low_band[0] < 0 or low_band[1] <= low_band[0]:
        raise ValueError(f"Invalid low_band: {low_band}")
    if high_band[0] < 0 or high_band[1] <= high_band[0]:
        raise ValueError(f"Invalid high_band: {high_band}")

    wav_files = collect_wav_files(input_path, recursive=args.recursive)

    summaries: List[FileSummary] = []

    for wav_path in wav_files:
        try:
            summary = calculate_file(
                wav_path=wav_path,
                output_dir=output_dir,
                frame_size=args.frame_size,
                hop_size=args.hop_size,
                low_band=low_band,
                high_band=high_band,
                eps=args.eps,
                remove_dc=not args.keep_dc,
                save_framewise_csv=not args.no_framewise_csv,
            )
            summaries.append(summary)

            print(
                f"OK: {wav_path} | "
                f"median alpha ratio = {summary.median_alpha_ratio_db:.6f} dB | "
                f"frames = {summary.n_frames}"
            )
        except Exception as exc:
            print(f"ERROR: {wav_path} | {exc}")

    if not summaries:
        raise RuntimeError("No WAV files were successfully processed.")

    write_summary_csv(summary_csv, summaries)
    print()
    print(f"Summary CSV written to: {summary_csv}")


if __name__ == "__main__":
    main()
