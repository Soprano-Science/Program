#!/usr/bin/env python3
"""Y_Q_value: frame-wise LPC-based Q value after Yamashita's Excel example.

This script generalizes the one-frame Q-value workflow illustrated in
``1_Cessa_mono_Q.xlsx`` to an entire WAV file or to all WAV files in a folder.

Definition for one valid 2048-sample frame
------------------------------------------
1. A 2048-sample frame is multiplied by the periodic Hann window
   ``0.5 * (1 - cos(2*pi*n/2048))``.
2. LPC coefficients are estimated using the autocorrelation / Yule-Walker method.
3. The LPC spectral envelope is evaluated on a normal 2048-point one-sided
   frequency grid. The output frequencies are bin-centered values:
   ``(k + 0.5) * sample_rate / 2048`` for ``k = 0 ... 1023``.
4. Candidate LPC-envelope peaks are searched in 2400--4000 Hz.
5. For each candidate peak, the left and right -3 dB crossing frequencies
   are estimated by linear interpolation.
6. Q is defined as

       Q = fa / (fb - fc)

   where ``fa`` is the peak frequency, ``fc`` is the lower -3 dB crossing,
   and ``fb`` is the upper -3 dB crossing.

The program writes frame-wise CSV files and a summary CSV. The file-level
representative Q value is the median of valid frame-wise Q values.

Short frames (<2048 samples) and silent/near-silent frames are recorded in the
CSV but excluded from Q calculation and from the median.
"""

from __future__ import annotations

import argparse
import csv
import math
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import scipy
from scipy.linalg import solve_toeplitz
from scipy.signal import find_peaks, freqz
import soundfile as sf


FRAME_SIZE = 2048
HOP_SIZE = 2048
EXPECTED_SR = 44_100
DEFAULT_LPC_ORDER = 12
PEAK_LOW_HZ = 2_400.0
PEAK_HIGH_HZ = 4_000.0
DEFAULT_SILENCE_THRESHOLD_DBFS = -60.0
EPS = np.finfo(np.float64).tiny


@dataclass(frozen=True)
class QCandidate:
    fa_hz: float
    fc_hz: float
    fb_hz: float
    bandwidth_hz: float
    q_value: float
    peak_db: float
    target_db: float
    peak_index: int


def read_wav_mono_float(path: Path) -> tuple[np.ndarray, int]:
    """Read a WAV file and return a mono float64 signal in approximately [-1, 1]."""
    audio, sr = sf.read(str(path), always_2d=False)
    audio = np.asarray(audio)

    if audio.ndim == 2:
        audio = audio.mean(axis=1)
    elif audio.ndim != 1:
        raise ValueError(f"Unsupported audio shape {audio.shape} in {path}")

    audio = audio.astype(np.float64, copy=False)
    return audio, int(sr)


def periodic_hann(n: int) -> np.ndarray:
    """Return Excel-style periodic Hann window: 0.5*(1-cos(2*pi*n/N))."""
    idx = np.arange(n, dtype=np.float64)
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * idx / n))


def rms_dbfs(frame: np.ndarray) -> float:
    """Return RMS in dBFS for a floating-point signal normalized near [-1, 1]."""
    if frame.size == 0:
        return -math.inf
    rms = float(np.sqrt(np.mean(np.square(frame, dtype=np.float64))))
    if rms <= 0.0 or not math.isfinite(rms):
        return -math.inf
    return 20.0 * math.log10(rms)


def lpc_autocorrelation(frame: np.ndarray, order: int) -> np.ndarray:
    """Return LPC denominator coefficients by autocorrelation/Yule-Walker.

    The returned polynomial is A(z) = 1 + a1 z^-1 + ... + ap z^-p.
    No pre-emphasis is applied.
    """
    frame = np.asarray(frame, dtype=np.float64)
    if frame.ndim != 1:
        raise ValueError("frame must be one-dimensional")
    if frame.size <= order:
        raise ValueError("frame is too short for the requested LPC order")

    x = frame - np.mean(frame)
    corr_full = np.correlate(x, x, mode="full")
    mid = frame.size - 1
    r = corr_full[mid : mid + order + 1] / frame.size

    if not np.isfinite(r[0]) or r[0] <= np.finfo(np.float64).eps:
        raise ValueError("silent or numerically degenerate frame")

    a_tail = solve_toeplitz((r[:-1], r[:-1]), -r[1:])
    return np.concatenate(([1.0], a_tail))


def interpolate_crossing(f1: float, y1: float, f2: float, y2: float, target: float) -> float:
    """Linearly interpolate the frequency where a segment crosses target dB."""
    if y2 == y1:
        return float((f1 + f2) / 2.0)
    return float(f1 + (target - y1) * (f2 - f1) / (y2 - y1))


def lpc_envelope_db(
    a: np.ndarray,
    sample_rate: int,
    frame_size: int = FRAME_SIZE,
) -> tuple[np.ndarray, np.ndarray]:
    """Return bin-centered one-sided frequency grid and LPC envelope in dB.

    The grid matches the frequency layout shown in Yamashita's Excel example:
    10.766602, 32.299805, ... Hz at 44.1 kHz, corresponding to
    (k + 0.5) * sample_rate / 2048.
    """
    freqs = (np.arange(frame_size // 2, dtype=np.float64) + 0.5) * sample_rate / frame_size
    _, response = freqz(b=[1.0], a=a, worN=freqs, fs=sample_rate)
    magnitude_db = 20.0 * np.log10(np.maximum(np.abs(response), 1e-15))
    return freqs, magnitude_db


def q_candidates_from_lpc_envelope(
    freqs: np.ndarray,
    magnitude_db: np.ndarray,
    peak_low_hz: float = PEAK_LOW_HZ,
    peak_high_hz: float = PEAK_HIGH_HZ,
) -> list[QCandidate]:
    """Find valid Q candidates in the specified LPC peak-search range."""
    band = (freqs >= peak_low_hz) & (freqs <= peak_high_hz)
    band_indices = np.flatnonzero(band)
    if band_indices.size < 3:
        return []

    local_peaks, _ = find_peaks(magnitude_db[band_indices])
    peak_indices = band_indices[local_peaks]
    results: list[QCandidate] = []

    for peak_idx in peak_indices:
        peak_db = float(magnitude_db[peak_idx])
        target_db = peak_db - 3.0

        left_idx = int(peak_idx)
        while left_idx > 0 and magnitude_db[left_idx] > target_db:
            left_idx -= 1
        if left_idx == 0 and magnitude_db[left_idx] > target_db:
            continue

        right_idx = int(peak_idx)
        last_idx = len(magnitude_db) - 1
        while right_idx < last_idx and magnitude_db[right_idx] > target_db:
            right_idx += 1
        if right_idx == last_idx and magnitude_db[right_idx] > target_db:
            continue

        fc_hz = interpolate_crossing(
            float(freqs[left_idx]),
            float(magnitude_db[left_idx]),
            float(freqs[left_idx + 1]),
            float(magnitude_db[left_idx + 1]),
            target_db,
        )
        fb_hz = interpolate_crossing(
            float(freqs[right_idx - 1]),
            float(magnitude_db[right_idx - 1]),
            float(freqs[right_idx]),
            float(magnitude_db[right_idx]),
            target_db,
        )

        fa_hz = float(freqs[peak_idx])
        bandwidth_hz = float(fb_hz - fc_hz)
        if bandwidth_hz <= 0.0 or not math.isfinite(bandwidth_hz):
            continue
        q_value = float(fa_hz / bandwidth_hz)
        if math.isfinite(q_value):
            results.append(
                QCandidate(
                    fa_hz=fa_hz,
                    fc_hz=float(fc_hz),
                    fb_hz=float(fb_hz),
                    bandwidth_hz=bandwidth_hz,
                    q_value=q_value,
                    peak_db=peak_db,
                    target_db=float(target_db),
                    peak_index=int(peak_idx),
                )
            )
    return results


def select_candidate(candidates: list[QCandidate], policy: str) -> QCandidate:
    """Select one candidate when multiple valid peaks are present."""
    if not candidates:
        raise ValueError("no candidates")
    if policy == "largest-q":
        return max(candidates, key=lambda item: item.q_value)
    if policy == "highest-peak":
        return max(candidates, key=lambda item: item.peak_db)
    raise ValueError(f"Unknown peak selection policy: {policy}")


def iter_wav_files(input_path: Path, recursive: bool = False) -> list[Path]:
    """Return WAV files from a single path or folder."""
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        pattern = "**/*.wav" if recursive else "*.wav"
        return sorted(input_path.glob(pattern))
    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def safe_stem(path: Path) -> str:
    """Return a filename-safe stem for output files."""
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in path.stem)


def analyze_wav(
    wav_path: Path,
    output_dir: Path,
    *,
    frame_size: int = FRAME_SIZE,
    hop_size: int = HOP_SIZE,
    expected_sr: int | None = EXPECTED_SR,
    lpc_order: int = DEFAULT_LPC_ORDER,
    peak_low_hz: float = PEAK_LOW_HZ,
    peak_high_hz: float = PEAK_HIGH_HZ,
    silence_threshold_dbfs: float = DEFAULT_SILENCE_THRESHOLD_DBFS,
    peak_selection: str = "largest-q",
) -> dict[str, str | int | float]:
    """Analyze one WAV file and write its frame-wise CSV."""
    audio, sr = read_wav_mono_float(wav_path)
    if expected_sr is not None and sr != expected_sr:
        raise ValueError(f"{wav_path.name}: expected {expected_sr} Hz, found {sr} Hz")

    output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = output_dir / f"{safe_stem(wav_path)}_Y_Q_value_framewise.csv"

    window = periodic_hann(frame_size)
    rows: list[dict[str, str | int | float]] = []
    valid_q_values: list[float] = []

    total_samples = int(audio.size)
    frame_index = 0
    for start in range(0, total_samples, hop_size):
        frame = audio[start : start + frame_size]
        frame_length = int(frame.size)
        end = start + frame_length
        is_short = frame_length < frame_size
        frame_rms_dbfs = rms_dbfs(frame)
        is_silent = (not math.isfinite(frame_rms_dbfs)) or frame_rms_dbfs <= silence_threshold_dbfs

        base_row: dict[str, str | int | float] = {
            "file_name": wav_path.name,
            "frame_index": frame_index,
            "start_sample": int(start),
            "end_sample_exclusive": int(end),
            "frame_length_samples": frame_length,
            "start_time_s": start / sr,
            "center_time_s": (start + frame_length / 2.0) / sr,
            "end_time_s": end / sr,
            "rms_dbfs": frame_rms_dbfs if math.isfinite(frame_rms_dbfs) else "-inf",
            "is_short_frame": str(is_short),
            "is_silent_frame": str(is_silent),
            "include_in_median": "False",
            "exclusion_reason": "",
            "fa_hz": "",
            "fc_hz": "",
            "fb_hz": "",
            "bandwidth_hz": "",
            "q_value": "",
            "peak_db": "",
            "minus_3db_target_db": "",
            "n_valid_candidates": 0,
            "selected_peak_policy": peak_selection,
            "status": "",
        }

        exclusion_reasons = []
        if is_short:
            exclusion_reasons.append(f"short_frame_less_than_{frame_size}_samples")
        if is_silent:
            exclusion_reasons.append(f"silent_or_near_silent_RMS_le_{silence_threshold_dbfs:g}_dBFS")

        if exclusion_reasons:
            base_row["exclusion_reason"] = "; ".join(exclusion_reasons)
            base_row["status"] = "excluded_before_Q_calculation"
            rows.append(base_row)
            frame_index += 1
            continue

        try:
            frame_windowed = frame * window
            a = lpc_autocorrelation(frame_windowed, lpc_order)
            freqs, envelope_db = lpc_envelope_db(a, sr, frame_size=frame_size)
            candidates = q_candidates_from_lpc_envelope(
                freqs,
                envelope_db,
                peak_low_hz=peak_low_hz,
                peak_high_hz=peak_high_hz,
            )
            base_row["n_valid_candidates"] = len(candidates)
            if not candidates:
                base_row["exclusion_reason"] = "no_valid_2p4_4p0_kHz_LPC_peak_or_3dB_crossings"
                base_row["status"] = "no_valid_Q"
                rows.append(base_row)
                frame_index += 1
                continue
            selected = select_candidate(candidates, peak_selection)
        except Exception as exc:  # Keep batch processing robust and transparent.
            base_row["exclusion_reason"] = "LPC_or_Q_calculation_failed"
            base_row["status"] = f"error: {type(exc).__name__}: {exc}"
            rows.append(base_row)
            frame_index += 1
            continue

        base_row.update(
            {
                "include_in_median": "True",
                "fa_hz": selected.fa_hz,
                "fc_hz": selected.fc_hz,
                "fb_hz": selected.fb_hz,
                "bandwidth_hz": selected.bandwidth_hz,
                "q_value": selected.q_value,
                "peak_db": selected.peak_db,
                "minus_3db_target_db": selected.target_db,
                "exclusion_reason": "",
                "status": "valid",
            }
        )
        valid_q_values.append(selected.q_value)
        rows.append(base_row)
        frame_index += 1

    fieldnames = [
        "file_name",
        "frame_index",
        "start_sample",
        "end_sample_exclusive",
        "frame_length_samples",
        "start_time_s",
        "center_time_s",
        "end_time_s",
        "rms_dbfs",
        "is_short_frame",
        "is_silent_frame",
        "include_in_median",
        "exclusion_reason",
        "fa_hz",
        "fc_hz",
        "fb_hz",
        "bandwidth_hz",
        "q_value",
        "peak_db",
        "minus_3db_target_db",
        "n_valid_candidates",
        "selected_peak_policy",
        "status",
    ]
    with out_csv.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    q_arr = np.asarray(valid_q_values, dtype=np.float64)
    if q_arr.size:
        median_q = float(np.median(q_arr))
        mean_q = float(np.mean(q_arr))
        sd_q = float(np.std(q_arr, ddof=1)) if q_arr.size > 1 else 0.0
        min_q = float(np.min(q_arr))
        max_q = float(np.max(q_arr))
        q1 = float(np.percentile(q_arr, 25))
        q3 = float(np.percentile(q_arr, 75))
        iqr = float(q3 - q1)
    else:
        median_q = mean_q = sd_q = min_q = max_q = q1 = q3 = iqr = math.nan

    n_short = sum(row["is_short_frame"] == "True" for row in rows)
    n_silent = sum(row["is_silent_frame"] == "True" for row in rows)
    n_no_valid_q = sum(row["status"] == "no_valid_Q" for row in rows)

    return {
        "file_name": wav_path.name,
        "sample_rate_hz": sr,
        "duration_s": total_samples / sr,
        "total_samples": total_samples,
        "frame_size_samples": frame_size,
        "hop_size_samples": hop_size,
        "lpc_order": lpc_order,
        "peak_search_low_hz": peak_low_hz,
        "peak_search_high_hz": peak_high_hz,
        "silence_threshold_dbfs": silence_threshold_dbfs,
        "total_frames_recorded": len(rows),
        "n_short_frames_excluded": int(n_short),
        "n_silent_frames_excluded": int(n_silent),
        "n_no_valid_q_frames": int(n_no_valid_q),
        "n_valid_q_frames": int(q_arr.size),
        "median_Y_Q_value": median_q,
        "mean_Y_Q_value": mean_q,
        "sd_Y_Q_value": sd_q,
        "min_Y_Q_value": min_q,
        "max_Y_Q_value": max_q,
        "q1_Y_Q_value": q1,
        "q3_Y_Q_value": q3,
        "iqr_Y_Q_value": iqr,
        "framewise_csv": str(out_csv),
        "status": "ok" if q_arr.size else "no_valid_q_values",
    }


def write_summary(rows: list[dict[str, str | int | float]], output_csv: Path) -> None:
    """Write summary rows to CSV."""
    if not rows:
        return
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with output_csv.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate Yamashita-style frame-wise LPC-based Q values from WAV files."
    )
    parser.add_argument("input", type=Path, help="Input WAV file or folder containing WAV files.")
    parser.add_argument("--output-dir", type=Path, default=Path("Y_Q_value_results"))
    parser.add_argument("--recursive", action="store_true", help="Recursively process WAV files in folders.")
    parser.add_argument("--frame-size", type=int, default=FRAME_SIZE)
    parser.add_argument("--hop-size", type=int, default=HOP_SIZE)
    parser.add_argument("--lpc-order", type=int, default=DEFAULT_LPC_ORDER)
    parser.add_argument("--peak-low-hz", type=float, default=PEAK_LOW_HZ)
    parser.add_argument("--peak-high-hz", type=float, default=PEAK_HIGH_HZ)
    parser.add_argument("--silence-threshold-dbfs", type=float, default=DEFAULT_SILENCE_THRESHOLD_DBFS)
    parser.add_argument(
        "--peak-selection",
        choices=("largest-q", "highest-peak"),
        default="largest-q",
        help="How to select a candidate when multiple valid 2.4-4.0 kHz peaks exist.",
    )
    parser.add_argument(
        "--allow-non-44100",
        action="store_true",
        help="Allow sampling rates other than 44.1 kHz. No resampling is performed.",
    )
    parser.add_argument("--version", action="version", version=f"Python {platform.python_version()}, NumPy {np.__version__}, SciPy {scipy.__version__}, SoundFile {sf.__version__}")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        wav_files = iter_wav_files(args.input, recursive=args.recursive)
        if not wav_files:
            raise ValueError(f"No WAV files found in {args.input}")

        rows = []
        for wav_path in wav_files:
            print(f"Processing: {wav_path}")
            row = analyze_wav(
                wav_path,
                args.output_dir,
                frame_size=args.frame_size,
                hop_size=args.hop_size,
                expected_sr=None if args.allow_non_44100 else EXPECTED_SR,
                lpc_order=args.lpc_order,
                peak_low_hz=args.peak_low_hz,
                peak_high_hz=args.peak_high_hz,
                silence_threshold_dbfs=args.silence_threshold_dbfs,
                peak_selection=args.peak_selection,
            )
            rows.append(row)

        summary_csv = args.output_dir / "Y_Q_value_summary.csv"
        write_summary(rows, summary_csv)
        print(f"Wrote summary: {summary_csv}")
        for row in rows:
            print(f"{row['file_name']}: median_Y_Q_value={row['median_Y_Q_value']}, valid_frames={row['n_valid_q_frames']}")
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
