#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frame-wise H1-H2 analysis for WAV files.

This script calculates H1-H2 for an entire WAV file using non-overlapping
2048-sample frames by default. The final frame is retained even when it is
shorter than 2048 samples; it is zero-padded to the FFT size for spectral
analysis.

Primary output:
    median_h1_minus_h2_db

Optional/diagnostic output:
    median_h1_star_minus_h2_star_db

Important methodological note:
    H1-H2 is defined here as the dB amplitude difference between the first
    harmonic (H1) and the second harmonic (H2):

        H1-H2 = H1_dB - H2_dB

    The uncorrected H1-H2 median is the default representative value.
    An approximate F1-corrected value H1*-H2* can also be requested with
    --with-correction, but this correction is only as reliable as the formant
    estimate. For high-pitched singing, especially high vowels or frames where
    F0 is close to F1, the corrected value can be unstable.

Dependencies:
    numpy, scipy

Examples:
    python calculate_h1_h2_framewise_median.py input.wav
    python calculate_h1_h2_framewise_median.py ./audio --recursive
    python calculate_h1_h2_framewise_median.py input.wav --min-f0 150 --max-f0 1100
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.io import wavfile
from scipy.linalg import solve_toeplitz, toeplitz
from scipy.signal import find_peaks


EPS = 1e-20


@dataclass
class FrameResult:
    frame_index: int
    start_sample: int
    end_sample: int
    frame_length_samples: int
    start_time_s: float
    end_time_s: float
    rms_dbfs: float
    valid_h1h2: bool
    invalid_reason: str
    f0_hz: float
    autocorr_clarity: float
    h1_frequency_hz: float
    h2_frequency_hz: float
    h1_db: float
    h2_db: float
    h1_minus_h2_db: float
    f1_hz: float
    b1_hz: float
    h1_star_db: float
    h2_star_db: float
    h1_star_minus_h2_star_db: float
    formant_warning: str


@dataclass
class FileSummary:
    input_file: str
    sampling_rate_hz: int
    duration_s: float
    frame_size_samples: int
    hop_size_samples: int
    fft_size_samples: int
    n_frames_total: int
    n_valid_h1h2_frames: int
    n_valid_corrected_frames: int
    n_corrected_warning_frames: int
    n_h2_above_f1_warning_frames: int
    median_h1_minus_h2_db: float
    mean_h1_minus_h2_db: float
    sd_h1_minus_h2_db: float
    median_h1_star_minus_h2_star_db: float
    mean_h1_star_minus_h2_star_db: float
    sd_h1_star_minus_h2_star_db: float
    framewise_csv: str


def read_wav_as_mono_float(path: str | Path) -> tuple[int, np.ndarray]:
    """Read a WAV file and convert it to mono float64."""
    sr, x = wavfile.read(str(path))

    if x.ndim == 2:
        x = x.mean(axis=1)

    if np.issubdtype(x.dtype, np.integer):
        x = x.astype(np.float64) / max(abs(np.iinfo(x.dtype).min), np.iinfo(x.dtype).max)
    else:
        x = x.astype(np.float64)

    # Replace non-finite values, if any.
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return sr, x


def iter_frames_include_last(
    x: np.ndarray,
    frame_size: int,
    hop_size: int,
) -> Iterable[tuple[int, int, np.ndarray]]:
    """
    Yield frames with their start/end sample positions.

    The final frame is retained even when shorter than frame_size.
    """
    if len(x) == 0:
        return

    start = 0
    while start < len(x):
        end = min(start + frame_size, len(x))
        yield start, end, x[start:end]
        start += hop_size


def rms_dbfs(frame: np.ndarray) -> float:
    """Return RMS in dBFS for a frame normalized to approximately [-1, 1]."""
    if len(frame) == 0:
        return float("nan")
    rms = float(np.sqrt(np.mean(np.square(frame))))
    return 20.0 * math.log10(rms + EPS)


def parabolic_interpolation(y_minus: float, y0: float, y_plus: float) -> float:
    """
    Return the fractional peak offset for three adjacent samples.

    The result is in the range roughly [-0.5, 0.5] when the middle sample is a
    local maximum. If the parabola is ill-conditioned, return 0.
    """
    denom = y_minus - 2.0 * y0 + y_plus
    if abs(denom) < EPS:
        return 0.0
    return 0.5 * (y_minus - y_plus) / denom


def estimate_f0_autocorr(
    frame: np.ndarray,
    sr: int,
    min_f0: float,
    max_f0: float,
    min_clarity: float,
    min_rms_dbfs: float,
) -> tuple[float, float, str]:
    """
    Estimate F0 from one frame using normalized autocorrelation.

    Returns:
        f0_hz, clarity, invalid_reason

    If invalid_reason is not empty, f0_hz and clarity are NaN.
    """
    if len(frame) < 3:
        return float("nan"), float("nan"), "frame_too_short"

    frame_rms_db = rms_dbfs(frame)
    if frame_rms_db < min_rms_dbfs:
        return float("nan"), float("nan"), "below_min_rms"

    x = frame.astype(np.float64)
    x = x - np.mean(x)

    if len(x) > 1:
        x = x * np.hanning(len(x))

    if np.max(np.abs(x)) < EPS:
        return float("nan"), float("nan"), "near_zero_signal"

    min_lag = int(math.floor(sr / max_f0))
    max_lag = int(math.ceil(sr / min_f0))

    min_lag = max(1, min_lag)
    max_lag = min(max_lag, len(x) - 2)

    # Need enough data to cover at least the maximum lag.
    if max_lag <= min_lag:
        return float("nan"), float("nan"), "frame_too_short_for_f0_range"

    ac = np.correlate(x, x, mode="full")[len(x) - 1:]
    if ac[0] <= EPS:
        return float("nan"), float("nan"), "zero_autocorrelation_energy"

    ac_norm = ac / (ac[0] + EPS)
    search = ac_norm[min_lag:max_lag + 1]
    if len(search) == 0:
        return float("nan"), float("nan"), "empty_f0_search_range"

    # Use local peaks and prefer the earliest strong peak. This reduces the risk
    # of octave errors in which a later autocorrelation peak is selected.
    peaks, _ = find_peaks(search)
    if len(peaks) > 0:
        peak_heights = search[peaks]
        max_peak = float(np.max(peak_heights))
        strong_threshold = max(float(min_clarity), 0.85 * max_peak)
        strong_peaks = peaks[peak_heights >= strong_threshold]
        if len(strong_peaks) > 0:
            best_offset = int(strong_peaks[0])
        else:
            best_offset = int(peaks[int(np.argmax(peak_heights))])
    else:
        best_offset = int(np.argmax(search))

    best_lag = min_lag + best_offset
    clarity = float(ac_norm[best_lag])

    if clarity < min_clarity:
        return float("nan"), clarity, "low_autocorrelation_clarity"

    # Fractional lag refinement.
    refined_lag = float(best_lag)
    if 1 <= best_lag < len(ac_norm) - 1:
        delta = parabolic_interpolation(
            float(ac_norm[best_lag - 1]),
            float(ac_norm[best_lag]),
            float(ac_norm[best_lag + 1]),
        )
        if math.isfinite(delta) and abs(delta) <= 1.0:
            refined_lag += delta

    if refined_lag <= 0:
        return float("nan"), clarity, "invalid_refined_lag"

    f0 = sr / refined_lag
    if not (min_f0 <= f0 <= max_f0):
        return float("nan"), clarity, "f0_out_of_range"

    return float(f0), clarity, ""


def prepare_fft_frame(
    frame: np.ndarray,
    fft_size: int,
) -> np.ndarray:
    """Remove DC, apply Hann window, and zero-pad/truncate to fft_size."""
    x = frame.astype(np.float64)
    x = x - np.mean(x)

    if len(x) > 1:
        x = x * np.hanning(len(x))

    y = np.zeros(fft_size, dtype=np.float64)
    n = min(len(x), fft_size)
    y[:n] = x[:n]
    return y


def amplitude_peak_near(
    magnitude: np.ndarray,
    freqs: np.ndarray,
    target_hz: float,
    half_width_hz: float,
) -> tuple[float, float]:
    """
    Find a local spectral amplitude peak near target_hz.

    Returns:
        peak_frequency_hz, peak_amplitude_db
    """
    if target_hz <= 0 or target_hz >= freqs[-1]:
        return float("nan"), float("nan")

    mask = (freqs >= target_hz - half_width_hz) & (freqs <= target_hz + half_width_hz)
    idxs = np.flatnonzero(mask)
    if len(idxs) == 0:
        return float("nan"), float("nan")

    local_mags = magnitude[idxs]
    local_best = int(np.argmax(local_mags))
    best_idx = int(idxs[local_best])

    # Parabolic interpolation in log magnitude for frequency and amplitude.
    mag_db = 20.0 * np.log10(magnitude + EPS)
    peak_idx_float = float(best_idx)
    peak_db = float(mag_db[best_idx])

    if 1 <= best_idx < len(mag_db) - 1:
        delta = parabolic_interpolation(
            float(mag_db[best_idx - 1]),
            float(mag_db[best_idx]),
            float(mag_db[best_idx + 1]),
        )
        if math.isfinite(delta) and abs(delta) <= 1.0:
            peak_idx_float += delta
            # Value of fitted parabola at the vertex.
            y_minus = float(mag_db[best_idx - 1])
            y0 = float(mag_db[best_idx])
            y_plus = float(mag_db[best_idx + 1])
            peak_db = y0 - 0.25 * (y_minus - y_plus) * delta

    df = freqs[1] - freqs[0]
    peak_freq = peak_idx_float * df
    return float(peak_freq), float(peak_db)


def lpc_coefficients_autocorr(x: np.ndarray, order: int) -> np.ndarray | None:
    """Estimate LPC coefficients using the autocorrelation method."""
    if len(x) <= order + 1:
        return None

    r_full = np.correlate(x, x, mode="full")
    r = r_full[len(x) - 1:len(x) + order]
    if len(r) < order + 1 or r[0] <= EPS:
        return None

    try:
        # Solve Toeplitz system: R a = -r[1:]
        a_rest = solve_toeplitz((r[:order], r[:order]), -r[1:order + 1])
    except Exception:
        try:
            R = toeplitz(r[:order])
            a_rest = np.linalg.solve(R, -r[1:order + 1])
        except Exception:
            return None

    a = np.concatenate(([1.0], np.asarray(a_rest, dtype=np.float64)))
    if not np.all(np.isfinite(a)):
        return None
    return a


def estimate_f1_lpc(
    frame: np.ndarray,
    sr: int,
    lpc_order: int,
    preemphasis: float,
    min_f1: float,
    max_f1: float,
    max_bandwidth: float,
) -> tuple[float, float, str]:
    """
    Estimate F1 and B1 using LPC.

    Returns:
        f1_hz, b1_hz, warning_or_reason
    """
    if len(frame) < lpc_order + 2:
        return float("nan"), float("nan"), "frame_too_short_for_lpc"

    x = frame.astype(np.float64)
    x = x - np.mean(x)

    if np.max(np.abs(x)) < EPS:
        return float("nan"), float("nan"), "near_zero_signal"

    # Pre-emphasis helps LPC formant estimation.
    if len(x) > 1:
        x = np.append(x[0], x[1:] - preemphasis * x[:-1])
        x = x * np.hanning(len(x))

    a = lpc_coefficients_autocorr(x, lpc_order)
    if a is None:
        return float("nan"), float("nan"), "lpc_failed"

    roots = np.roots(a)
    roots = roots[np.imag(roots) >= 0]

    formants = []
    for root in roots:
        radius = abs(root)
        if radius <= 0 or radius >= 1.0:
            continue
        angle = math.atan2(root.imag, root.real)
        freq = angle * sr / (2.0 * math.pi)
        bandwidth = -math.log(radius) * sr / math.pi

        if (
            min_f1 <= freq <= max_f1
            and 20.0 <= bandwidth <= max_bandwidth
            and math.isfinite(freq)
            and math.isfinite(bandwidth)
        ):
            formants.append((float(freq), float(bandwidth)))

    if not formants:
        return float("nan"), float("nan"), "no_valid_f1_candidate"

    formants.sort(key=lambda pair: pair[0])
    f1, b1 = formants[0]
    return f1, b1, ""


def formant_boost_db(freq_hz: float, formant_hz: float, bandwidth_hz: float) -> float:
    """
    Approximate all-pole resonance boost in dB at freq_hz.

    This is a simple one-formant correction model. It is provided as a
    diagnostic approximation, not as an exact clone of VoiceSauce.
    """
    if (
        freq_hz <= 0
        or formant_hz <= 0
        or bandwidth_hz <= 0
        or not all(math.isfinite(v) for v in [freq_hz, formant_hz, bandwidth_hz])
    ):
        return float("nan")

    x = freq_hz / formant_hz
    b = bandwidth_hz / formant_hz
    denom = (1.0 - x * x) ** 2 + (b * x) ** 2
    return -10.0 * math.log10(denom + EPS)


def analyze_frame(
    frame: np.ndarray,
    sr: int,
    frame_index: int,
    start_sample: int,
    end_sample: int,
    frame_size: int,
    fft_size: int,
    min_f0: float,
    max_f0: float,
    min_clarity: float,
    min_rms_dbfs: float,
    harmonic_search_half_width_hz: float | None,
    estimate_corrected: bool,
    lpc_order: int,
    preemphasis: float,
    min_f1: float,
    max_f1: float,
    max_f1_bandwidth: float,
    formant_proximity_warning_hz: float,
) -> FrameResult:
    """Analyze one frame and return all frame-level values."""
    start_time_s = start_sample / sr
    end_time_s = end_sample / sr
    frame_rms_db = rms_dbfs(frame)

    f0, clarity, reason = estimate_f0_autocorr(
        frame=frame,
        sr=sr,
        min_f0=min_f0,
        max_f0=max_f0,
        min_clarity=min_clarity,
        min_rms_dbfs=min_rms_dbfs,
    )

    nan = float("nan")
    if reason:
        return FrameResult(
            frame_index=frame_index,
            start_sample=start_sample,
            end_sample=end_sample,
            frame_length_samples=len(frame),
            start_time_s=start_time_s,
            end_time_s=end_time_s,
            rms_dbfs=frame_rms_db,
            valid_h1h2=False,
            invalid_reason=reason,
            f0_hz=f0,
            autocorr_clarity=clarity,
            h1_frequency_hz=nan,
            h2_frequency_hz=nan,
            h1_db=nan,
            h2_db=nan,
            h1_minus_h2_db=nan,
            f1_hz=nan,
            b1_hz=nan,
            h1_star_db=nan,
            h2_star_db=nan,
            h1_star_minus_h2_star_db=nan,
            formant_warning="",
        )

    fft_frame = prepare_fft_frame(frame, fft_size=fft_size)
    spectrum = np.fft.rfft(fft_frame, n=fft_size)
    magnitude = np.abs(spectrum)
    freqs = np.fft.rfftfreq(fft_size, d=1.0 / sr)

    if 2.0 * f0 >= sr / 2.0:
        return FrameResult(
            frame_index=frame_index,
            start_sample=start_sample,
            end_sample=end_sample,
            frame_length_samples=len(frame),
            start_time_s=start_time_s,
            end_time_s=end_time_s,
            rms_dbfs=frame_rms_db,
            valid_h1h2=False,
            invalid_reason="h2_above_nyquist",
            f0_hz=f0,
            autocorr_clarity=clarity,
            h1_frequency_hz=nan,
            h2_frequency_hz=nan,
            h1_db=nan,
            h2_db=nan,
            h1_minus_h2_db=nan,
            f1_hz=nan,
            b1_hz=nan,
            h1_star_db=nan,
            h2_star_db=nan,
            h1_star_minus_h2_star_db=nan,
            formant_warning="",
        )

    if harmonic_search_half_width_hz is None:
        # A moderate adaptive window around the expected harmonic.
        # It is wide enough to tolerate small F0 errors but narrow enough to avoid
        # adjacent harmonics for the default singing range.
        half_width_hz = min(60.0, max(20.0, 0.15 * f0))
    else:
        half_width_hz = harmonic_search_half_width_hz

    h1_freq, h1_db = amplitude_peak_near(
        magnitude=magnitude,
        freqs=freqs,
        target_hz=f0,
        half_width_hz=half_width_hz,
    )
    h2_freq, h2_db = amplitude_peak_near(
        magnitude=magnitude,
        freqs=freqs,
        target_hz=2.0 * f0,
        half_width_hz=half_width_hz,
    )

    if not (math.isfinite(h1_db) and math.isfinite(h2_db)):
        return FrameResult(
            frame_index=frame_index,
            start_sample=start_sample,
            end_sample=end_sample,
            frame_length_samples=len(frame),
            start_time_s=start_time_s,
            end_time_s=end_time_s,
            rms_dbfs=frame_rms_db,
            valid_h1h2=False,
            invalid_reason="harmonic_peak_not_found",
            f0_hz=f0,
            autocorr_clarity=clarity,
            h1_frequency_hz=h1_freq,
            h2_frequency_hz=h2_freq,
            h1_db=h1_db,
            h2_db=h2_db,
            h1_minus_h2_db=nan,
            f1_hz=nan,
            b1_hz=nan,
            h1_star_db=nan,
            h2_star_db=nan,
            h1_star_minus_h2_star_db=nan,
            formant_warning="",
        )

    h1_h2 = h1_db - h2_db

    f1 = nan
    b1 = nan
    h1_star = nan
    h2_star = nan
    h1star_h2star = nan
    formant_warning = ""

    if estimate_corrected:
        f1, b1, f1_reason = estimate_f1_lpc(
            frame=frame,
            sr=sr,
            lpc_order=lpc_order,
            preemphasis=preemphasis,
            min_f1=min_f1,
            max_f1=max_f1,
            max_bandwidth=max_f1_bandwidth,
        )

        if f1_reason:
            formant_warning = f1_reason
        else:
            correction_warnings: list[str] = []

            # When H2 is above F1, H2 lies on the descending side of the
            # first-formant resonance. In that configuration, a small F1 or B1
            # estimation error can produce a large change in the corrected H2
            # amplitude. This warning is therefore written explicitly to CSV.
            if math.isfinite(h2_freq) and h2_freq > f1:
                correction_warnings.append(
                    "H2_above_F1_correction_sensitive_to_F1_estimate"
                )

            if abs(f1 - f0) < formant_proximity_warning_hz:
                correction_warnings.append("F1_close_to_H1_correction_may_be_unstable")
            if abs(f1 - 2.0 * f0) < formant_proximity_warning_hz:
                correction_warnings.append("F1_close_to_H2_correction_may_be_unstable")

            formant_warning = ";".join(correction_warnings)

            boost_h1 = formant_boost_db(f0, f1, b1)
            boost_h2 = formant_boost_db(2.0 * f0, f1, b1)
            if math.isfinite(boost_h1) and math.isfinite(boost_h2):
                h1_star = h1_db - boost_h1
                h2_star = h2_db - boost_h2
                h1star_h2star = h1_star - h2_star
            else:
                formant_warning = formant_warning or "correction_failed"

    return FrameResult(
        frame_index=frame_index,
        start_sample=start_sample,
        end_sample=end_sample,
        frame_length_samples=len(frame),
        start_time_s=start_time_s,
        end_time_s=end_time_s,
        rms_dbfs=frame_rms_db,
        valid_h1h2=True,
        invalid_reason="",
        f0_hz=f0,
        autocorr_clarity=clarity,
        h1_frequency_hz=h1_freq,
        h2_frequency_hz=h2_freq,
        h1_db=h1_db,
        h2_db=h2_db,
        h1_minus_h2_db=h1_h2,
        f1_hz=f1,
        b1_hz=b1,
        h1_star_db=h1_star,
        h2_star_db=h2_star,
        h1_star_minus_h2_star_db=h1star_h2star,
        formant_warning=formant_warning,
    )


def safe_output_stem(input_path: Path, root: Path | None = None) -> str:
    """Create a safe unique-ish output stem from a file path."""
    try:
        if root is not None:
            rel = input_path.resolve().relative_to(root.resolve())
            raw = str(rel.with_suffix(""))
        else:
            raw = input_path.stem
    except Exception:
        raw = input_path.stem
    return re.sub(r"[^A-Za-z0-9._-]+", "_", raw)


def float_for_csv(value: float, digits: int = 6) -> str:
    """Format floats for CSV, leaving NaN as empty."""
    if value is None or not math.isfinite(float(value)):
        return ""
    return f"{float(value):.{digits}f}"


def save_framewise_csv(path: Path, rows: list[FrameResult]) -> None:
    """Write frame-wise results to CSV."""
    fieldnames = [
        "frame_index",
        "start_sample",
        "end_sample",
        "frame_length_samples",
        "start_time_s",
        "end_time_s",
        "rms_dbfs",
        "valid_h1h2",
        "invalid_reason",
        "f0_hz",
        "autocorr_clarity",
        "h1_frequency_hz",
        "h2_frequency_hz",
        "h1_db",
        "h2_db",
        "h1_minus_h2_db",
        "f1_hz",
        "b1_hz",
        "h1_star_db",
        "h2_star_db",
        "h1_star_minus_h2_star_db",
        "formant_warning",
    ]

    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({
                "frame_index": r.frame_index,
                "start_sample": r.start_sample,
                "end_sample": r.end_sample,
                "frame_length_samples": r.frame_length_samples,
                "start_time_s": float_for_csv(r.start_time_s),
                "end_time_s": float_for_csv(r.end_time_s),
                "rms_dbfs": float_for_csv(r.rms_dbfs),
                "valid_h1h2": "TRUE" if r.valid_h1h2 else "FALSE",
                "invalid_reason": r.invalid_reason,
                "f0_hz": float_for_csv(r.f0_hz),
                "autocorr_clarity": float_for_csv(r.autocorr_clarity),
                "h1_frequency_hz": float_for_csv(r.h1_frequency_hz),
                "h2_frequency_hz": float_for_csv(r.h2_frequency_hz),
                "h1_db": float_for_csv(r.h1_db),
                "h2_db": float_for_csv(r.h2_db),
                "h1_minus_h2_db": float_for_csv(r.h1_minus_h2_db),
                "f1_hz": float_for_csv(r.f1_hz),
                "b1_hz": float_for_csv(r.b1_hz),
                "h1_star_db": float_for_csv(r.h1_star_db),
                "h2_star_db": float_for_csv(r.h2_star_db),
                "h1_star_minus_h2_star_db": float_for_csv(r.h1_star_minus_h2_star_db),
                "formant_warning": r.formant_warning,
            })


def nanmedian(values: list[float]) -> float:
    arr = np.asarray([v for v in values if math.isfinite(float(v))], dtype=np.float64)
    if arr.size == 0:
        return float("nan")
    return float(np.median(arr))


def nanmean(values: list[float]) -> float:
    arr = np.asarray([v for v in values if math.isfinite(float(v))], dtype=np.float64)
    if arr.size == 0:
        return float("nan")
    return float(np.mean(arr))


def nansd(values: list[float]) -> float:
    arr = np.asarray([v for v in values if math.isfinite(float(v))], dtype=np.float64)
    if arr.size <= 1:
        return float("nan")
    return float(np.std(arr, ddof=1))


def analyze_wav_file(
    wav_path: Path,
    output_dir: Path,
    root_for_output_name: Path | None,
    frame_size: int,
    hop_size: int,
    fft_size: int,
    min_f0: float,
    max_f0: float,
    min_clarity: float,
    min_rms_dbfs: float,
    harmonic_search_half_width_hz: float | None,
    estimate_corrected: bool,
    lpc_order: int,
    preemphasis: float,
    min_f1: float,
    max_f1: float,
    max_f1_bandwidth: float,
    formant_proximity_warning_hz: float,
) -> FileSummary:
    """Analyze one WAV file and return a summary."""
    sr, x = read_wav_as_mono_float(wav_path)

    if sr <= 0:
        raise ValueError(f"Invalid sampling rate for {wav_path}: {sr}")
    if max_f0 * 2.0 >= sr / 2.0:
        raise ValueError(
            f"max_f0={max_f0} Hz is too high for sampling rate {sr} Hz. "
            "H2 must be below Nyquist. Lower --max-f0 or use a higher sampling rate."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = safe_output_stem(wav_path, root_for_output_name)
    framewise_csv = output_dir / f"{stem}_h1_h2_framewise.csv"

    rows: list[FrameResult] = []
    for frame_index, (start, end, frame) in enumerate(
        iter_frames_include_last(x, frame_size=frame_size, hop_size=hop_size)
    ):
        rows.append(
            analyze_frame(
                frame=frame,
                sr=sr,
                frame_index=frame_index,
                start_sample=start,
                end_sample=end,
                frame_size=frame_size,
                fft_size=fft_size,
                min_f0=min_f0,
                max_f0=max_f0,
                min_clarity=min_clarity,
                min_rms_dbfs=min_rms_dbfs,
                harmonic_search_half_width_hz=harmonic_search_half_width_hz,
                estimate_corrected=estimate_corrected,
                lpc_order=lpc_order,
                preemphasis=preemphasis,
                min_f1=min_f1,
                max_f1=max_f1,
                max_f1_bandwidth=max_f1_bandwidth,
                formant_proximity_warning_hz=formant_proximity_warning_hz,
            )
        )

    save_framewise_csv(framewise_csv, rows)

    h1h2_values = [r.h1_minus_h2_db for r in rows if r.valid_h1h2 and math.isfinite(r.h1_minus_h2_db)]
    h1h2c_values = [
        r.h1_star_minus_h2_star_db
        for r in rows
        if r.valid_h1h2 and math.isfinite(r.h1_star_minus_h2_star_db)
    ]
    corrected_warning_rows = [
        r
        for r in rows
        if r.valid_h1h2
        and math.isfinite(r.h1_star_minus_h2_star_db)
        and bool(r.formant_warning)
    ]
    h2_above_f1_warning_rows = [
        r
        for r in corrected_warning_rows
        if "H2_above_F1" in r.formant_warning
    ]

    return FileSummary(
        input_file=str(wav_path),
        sampling_rate_hz=sr,
        duration_s=len(x) / sr,
        frame_size_samples=frame_size,
        hop_size_samples=hop_size,
        fft_size_samples=fft_size,
        n_frames_total=len(rows),
        n_valid_h1h2_frames=len(h1h2_values),
        n_valid_corrected_frames=len(h1h2c_values),
        n_corrected_warning_frames=len(corrected_warning_rows),
        n_h2_above_f1_warning_frames=len(h2_above_f1_warning_rows),
        median_h1_minus_h2_db=nanmedian(h1h2_values),
        mean_h1_minus_h2_db=nanmean(h1h2_values),
        sd_h1_minus_h2_db=nansd(h1h2_values),
        median_h1_star_minus_h2_star_db=nanmedian(h1h2c_values),
        mean_h1_star_minus_h2_star_db=nanmean(h1h2c_values),
        sd_h1_star_minus_h2_star_db=nansd(h1h2c_values),
        framewise_csv=str(framewise_csv),
    )


def discover_wav_files(input_path: Path, recursive: bool) -> list[Path]:
    """Return WAV files to analyze."""
    if input_path.is_file():
        if input_path.suffix.lower() != ".wav":
            raise ValueError(f"Input file is not a WAV file: {input_path}")
        return [input_path]

    if input_path.is_dir():
        pattern = "**/*.wav" if recursive else "*.wav"
        files = sorted(input_path.glob(pattern))
        if not files:
            raise FileNotFoundError(f"No WAV files found in {input_path}")
        return files

    raise FileNotFoundError(f"Input path not found: {input_path}")


def save_summary_csv(path: Path, summaries: list[FileSummary]) -> None:
    """Write file-level summary CSV."""
    fieldnames = [
        "input_file",
        "sampling_rate_hz",
        "duration_s",
        "frame_size_samples",
        "hop_size_samples",
        "fft_size_samples",
        "n_frames_total",
        "n_valid_h1h2_frames",
        "n_valid_corrected_frames",
        "n_corrected_warning_frames",
        "n_h2_above_f1_warning_frames",
        "median_h1_minus_h2_db",
        "mean_h1_minus_h2_db",
        "sd_h1_minus_h2_db",
        "median_h1_star_minus_h2_star_db",
        "mean_h1_star_minus_h2_star_db",
        "sd_h1_star_minus_h2_star_db",
        "framewise_csv",
    ]

    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in summaries:
            writer.writerow({
                "input_file": s.input_file,
                "sampling_rate_hz": s.sampling_rate_hz,
                "duration_s": float_for_csv(s.duration_s),
                "frame_size_samples": s.frame_size_samples,
                "hop_size_samples": s.hop_size_samples,
                "fft_size_samples": s.fft_size_samples,
                "n_frames_total": s.n_frames_total,
                "n_valid_h1h2_frames": s.n_valid_h1h2_frames,
                "n_valid_corrected_frames": s.n_valid_corrected_frames,
                "n_corrected_warning_frames": s.n_corrected_warning_frames,
                "n_h2_above_f1_warning_frames": s.n_h2_above_f1_warning_frames,
                "median_h1_minus_h2_db": float_for_csv(s.median_h1_minus_h2_db),
                "mean_h1_minus_h2_db": float_for_csv(s.mean_h1_minus_h2_db),
                "sd_h1_minus_h2_db": float_for_csv(s.sd_h1_minus_h2_db),
                "median_h1_star_minus_h2_star_db": float_for_csv(s.median_h1_star_minus_h2_star_db),
                "mean_h1_star_minus_h2_star_db": float_for_csv(s.mean_h1_star_minus_h2_star_db),
                "sd_h1_star_minus_h2_star_db": float_for_csv(s.sd_h1_star_minus_h2_star_db),
                "framewise_csv": s.framewise_csv,
            })


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate frame-wise H1-H2 from one WAV file or all WAV files in a folder. "
            "The file-level H1-H2 value is the median of valid frame-wise H1-H2 values."
        )
    )

    parser.add_argument(
        "input",
        nargs="?",
        default="input.wav",
        help="Input WAV file or folder. Default: input.wav",
    )
    parser.add_argument(
        "--output-dir",
        default="h1_h2_results",
        help="Output directory. Default: h1_h2_results",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="If input is a folder, analyze WAV files recursively.",
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
        "--fft-size",
        type=int,
        default=2048,
        help="FFT size in samples. Default: 2048",
    )
    parser.add_argument(
        "--min-f0",
        type=float,
        default=70.0,
        help="Minimum F0 for autocorrelation search in Hz. Default: 70",
    )
    parser.add_argument(
        "--max-f0",
        type=float,
        default=1200.0,
        help="Maximum F0 for autocorrelation search in Hz. Default: 1200",
    )
    parser.add_argument(
        "--min-clarity",
        type=float,
        default=0.30,
        help="Minimum normalized autocorrelation peak for voiced frames. Default: 0.30",
    )
    parser.add_argument(
        "--min-rms-dbfs",
        type=float,
        default=-60.0,
        help="Frames below this RMS level are marked invalid. Default: -60 dBFS",
    )
    parser.add_argument(
        "--harmonic-search-half-width-hz",
        type=float,
        default=None,
        help=(
            "Half-width around expected H1/H2 frequency for local peak search. "
            "Default: adaptive width, min(60 Hz, max(20 Hz, 0.15*F0))."
        ),
    )
    parser.set_defaults(estimate_corrected=False)
    parser.add_argument(
        "--with-correction",
        dest="estimate_corrected",
        action="store_true",
        help=(
            "Also estimate an approximate F1-corrected H1*-H2* value. "
            "This is diagnostic only and is not the primary output."
        ),
    )
    parser.add_argument(
        "--no-correction",
        dest="estimate_corrected",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--lpc-order",
        type=int,
        default=14,
        help="LPC order for optional F1 estimation. Default: 14",
    )
    parser.add_argument(
        "--preemphasis",
        type=float,
        default=0.97,
        help="Pre-emphasis coefficient for optional LPC formant estimation. Default: 0.97",
    )
    parser.add_argument(
        "--min-f1",
        type=float,
        default=150.0,
        help="Minimum acceptable F1 candidate in Hz. Default: 150",
    )
    parser.add_argument(
        "--max-f1",
        type=float,
        default=1500.0,
        help="Maximum acceptable F1 candidate in Hz. Default: 1500",
    )
    parser.add_argument(
        "--max-f1-bandwidth",
        type=float,
        default=1000.0,
        help="Maximum acceptable F1 bandwidth in Hz. Default: 1000",
    )
    parser.add_argument(
        "--formant-proximity-warning-hz",
        type=float,
        default=100.0,
        help="Warn when F1 is within this distance of H1 or H2. Default: 100 Hz",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    wav_files = discover_wav_files(input_path, recursive=args.recursive)
    root_for_output_name = input_path if input_path.is_dir() else None

    summaries: list[FileSummary] = []
    for wav in wav_files:
        try:
            summary = analyze_wav_file(
                wav_path=wav,
                output_dir=output_dir,
                root_for_output_name=root_for_output_name,
                frame_size=args.frame_size,
                hop_size=args.hop_size,
                fft_size=args.fft_size,
                min_f0=args.min_f0,
                max_f0=args.max_f0,
                min_clarity=args.min_clarity,
                min_rms_dbfs=args.min_rms_dbfs,
                harmonic_search_half_width_hz=args.harmonic_search_half_width_hz,
                estimate_corrected=args.estimate_corrected,
                lpc_order=args.lpc_order,
                preemphasis=args.preemphasis,
                min_f1=args.min_f1,
                max_f1=args.max_f1,
                max_f1_bandwidth=args.max_f1_bandwidth,
                formant_proximity_warning_hz=args.formant_proximity_warning_hz,
            )
            summaries.append(summary)
            print(
                f"OK: {wav} | median H1-H2 = "
                f"{float_for_csv(summary.median_h1_minus_h2_db)} dB | "
                f"valid frames = {summary.n_valid_h1h2_frames}/{summary.n_frames_total}"
            )
        except Exception as exc:
            print(f"ERROR: {wav} | {exc}")

    summary_csv = output_dir / "h1_h2_summary.csv"
    save_summary_csv(summary_csv, summaries)

    print("\nDone.")
    print(f"Analyzed files: {len(summaries)}")
    print(f"Summary CSV   : {summary_csv}")
    print(f"Output folder : {output_dir}")


if __name__ == "__main__":
    main()
