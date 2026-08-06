from __future__ import annotations

import csv
import platform
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import solve_toeplitz
from scipy.signal import find_peaks, freqz
import soundfile as sf

try:
    import librosa
except Exception:  # Librosa is not required for the calculation.
    librosa = None


WAV_FILE = Path("kawano-after_2022_20190408tanto_tanto_02_mono(2).wav")

# Conditions stated in the cited study
EXPECTED_SR = 44_100
STABLE_START_SEC = 0.430
STABLE_END_SEC = 0.930          # 500-ms stable /a/ interval
FRAME_MS = 30.0
HOP_MS = 10.0
LPC_ORDER = 12
PEAK_LOW_HZ = 2_000.0
PEAK_HIGH_HZ = 4_000.0
FREQZ_POINTS = 131_072          # dense grid for interpolated -3 dB crossings


@dataclass(frozen=True)
class QResult:
    center_sec: float
    fa_hz: float
    fc_hz: float
    fb_hz: float
    bandwidth_hz: float
    q_value: float


def lpc_autocorrelation(frame: np.ndarray, order: int) -> np.ndarray:
    """Return LPC denominator coefficients by the autocorrelation/Yule-Walker method.

    The returned polynomial is A(z) = 1 + a1 z^-1 + ... + ap z^-p.
    """
    frame = np.asarray(frame, dtype=np.float64)
    if frame.ndim != 1:
        raise ValueError("frame must be one-dimensional")
    if len(frame) <= order:
        raise ValueError("frame is too short for the requested LPC order")

    # Remove only the frame DC component; no pre-emphasis is applied because
    # it is not stated in the cited analysis conditions.
    x = frame - np.mean(frame)

    # Biased autocorrelation for lags 0..order.
    corr_full = np.correlate(x, x, mode="full")
    mid = len(x) - 1
    r = corr_full[mid : mid + order + 1] / len(x)

    if not np.isfinite(r[0]) or r[0] <= np.finfo(float).eps:
        raise ValueError("silent or numerically degenerate frame")

    # Yule-Walker: R a = -r[1:].
    a_tail = solve_toeplitz((r[:-1], r[:-1]), -r[1:])
    return np.concatenate(([1.0], a_tail))


def interpolate_crossing(
    f1: float,
    y1: float,
    f2: float,
    y2: float,
    target: float,
) -> float:
    """Linearly interpolate a frequency where y crosses target."""
    if y2 == y1:
        return float((f1 + f2) / 2.0)
    return float(f1 + (target - y1) * (f2 - f1) / (y2 - y1))


def q_candidates_from_lpc(
    a: np.ndarray,
    sample_rate: int,
    center_sec: float,
) -> list[QResult]:
    """Find valid Q candidates whose LPC-envelope peak lies in 2-4 kHz."""
    freq, response = freqz(
        b=[1.0],
        a=a,
        worN=FREQZ_POINTS,
        fs=sample_rate,
    )
    magnitude_db = 20.0 * np.log10(np.maximum(np.abs(response), 1e-15))

    band = (freq >= PEAK_LOW_HZ) & (freq <= PEAK_HIGH_HZ)
    band_indices = np.flatnonzero(band)
    if band_indices.size < 3:
        return []

    local_peaks, _ = find_peaks(magnitude_db[band_indices])
    peak_indices = band_indices[local_peaks]
    results: list[QResult] = []

    for peak_idx in peak_indices:
        peak_db = magnitude_db[peak_idx]
        target_db = peak_db - 3.0

        # Search outward from the peak for the first -3 dB crossing on each side.
        left_idx = peak_idx
        while left_idx > 0 and magnitude_db[left_idx] > target_db:
            left_idx -= 1
        if left_idx == 0 and magnitude_db[left_idx] > target_db:
            continue

        right_idx = peak_idx
        last_idx = len(magnitude_db) - 1
        while right_idx < last_idx and magnitude_db[right_idx] > target_db:
            right_idx += 1
        if right_idx == last_idx and magnitude_db[right_idx] > target_db:
            continue

        fc_hz = interpolate_crossing(
            freq[left_idx],
            magnitude_db[left_idx],
            freq[left_idx + 1],
            magnitude_db[left_idx + 1],
            target_db,
        )
        fb_hz = interpolate_crossing(
            freq[right_idx - 1],
            magnitude_db[right_idx - 1],
            freq[right_idx],
            magnitude_db[right_idx],
            target_db,
        )
        fa_hz = float(freq[peak_idx])
        bandwidth_hz = fb_hz - fc_hz

        if bandwidth_hz <= 0.0:
            continue

        q_value = fa_hz / bandwidth_hz
        if np.isfinite(q_value):
            results.append(
                QResult(
                    center_sec=center_sec,
                    fa_hz=fa_hz,
                    fc_hz=fc_hz,
                    fb_hz=fb_hz,
                    bandwidth_hz=bandwidth_hz,
                    q_value=float(q_value),
                )
            )

    return results


def main() -> None:
    audio, sample_rate = sf.read(WAV_FILE, always_2d=False)
    audio = np.asarray(audio, dtype=np.float64)

    if audio.ndim != 1:
        raise ValueError("A mono WAV file is required.")
    if sample_rate != EXPECTED_SR:
        raise ValueError(
            f"Expected {EXPECTED_SR} Hz, but input is {sample_rate} Hz."
        )

    start_sample = int(round(STABLE_START_SEC * sample_rate))
    end_sample = int(round(STABLE_END_SEC * sample_rate))
    stable_vowel = audio[start_sample:end_sample]

    frame_length = int(round(FRAME_MS * sample_rate / 1000.0))
    hop_length = int(round(HOP_MS * sample_rate / 1000.0))
    window = np.hamming(frame_length)

    frame_results: list[QResult] = []
    frame_rows: list[dict[str, float | str]] = []
    total_frames = 0

    for offset in range(0, len(stable_vowel) - frame_length + 1, hop_length):
        total_frames += 1
        frame = stable_vowel[offset : offset + frame_length]
        center_sample = start_sample + offset + frame_length / 2.0
        center_sec = center_sample / sample_rate

        try:
            a = lpc_autocorrelation(frame * window, LPC_ORDER)
        except (ValueError, np.linalg.LinAlgError):
            frame_rows.append({
                "frame_number": total_frames,
                "center_sec": center_sec,
                "fa_hz": "",
                "fc_hz": "",
                "fb_hz": "",
                "bandwidth_hz": "",
                "q_value": "",
                "status": "LPC calculation failed",
            })
            continue

        candidates = q_candidates_from_lpc(a, sample_rate, center_sec)
        if not candidates:
            frame_rows.append({
                "frame_number": total_frames,
                "center_sec": center_sec,
                "fa_hz": "",
                "fc_hz": "",
                "fb_hz": "",
                "bandwidth_hz": "",
                "q_value": "",
                "status": "No valid 2-4 kHz peak / -3 dB crossings",
            })
            continue

        # If more than one valid high-frequency peak exists in a frame,
        # retain the sharpest peak (largest Q), as in the user's prior Q workflow.
        selected = max(candidates, key=lambda item: item.q_value)
        frame_results.append(selected)
        frame_rows.append({
            "frame_number": total_frames,
            "center_sec": selected.center_sec,
            "fa_hz": selected.fa_hz,
            "fc_hz": selected.fc_hz,
            "fb_hz": selected.fb_hz,
            "bandwidth_hz": selected.bandwidth_hz,
            "q_value": selected.q_value,
            "status": "valid",
        })

    if not frame_results:
        raise RuntimeError("No valid 2-4 kHz LPC peak with two -3 dB crossings was found.")

    q_values = np.array([item.q_value for item in frame_results], dtype=np.float64)
    mean_q = float(np.mean(q_values))
    median_q = float(np.median(q_values))
    std_q = float(np.std(q_values, ddof=1)) if len(q_values) > 1 else 0.0

    representative = min(frame_results, key=lambda item: abs(item.q_value - median_q))


    csv_path = WAV_FILE.with_name("kawano_after_q_frames.csv")
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "frame_number",
                "center_sec",
                "fa_hz",
                "fc_hz",
                "fb_hz",
                "bandwidth_hz",
                "q_value",
                "status",
            ],
        )
        writer.writeheader()
        writer.writerows(frame_rows)

    print("[Versions]")
    print(f"Python:    {platform.python_version()}")
    print(f"NumPy:     {np.__version__}")
    print(f"SciPy:     {scipy.__version__}")
    print(f"SoundFile: {sf.__version__}")
    print(
        "Librosa:   "
        + (librosa.__version__ if librosa is not None else "not installed")
        + " (not used in the primary calculation)"
    )
    print()

    print("[Analysis conditions]")
    print(f"Stable interval: {STABLE_START_SEC:.3f}-{STABLE_END_SEC:.3f} s")
    print(f"Sampling rate:   {sample_rate} Hz")
    print(f"Frame length:    {frame_length} samples ({FRAME_MS:.1f} ms)")
    print(f"Frame shift:     {hop_length} samples ({HOP_MS:.1f} ms)")
    print(f"Window:          Hamming")
    print(f"LPC order:       {LPC_ORDER}")
    print(f"Peak search:     {PEAK_LOW_HZ:.0f}-{PEAK_HIGH_HZ:.0f} Hz")
    print()

    print("[Q results]")
    print(f"Total frames:    {total_frames}")
    print(f"Valid frames:    {len(frame_results)}")
    print(f"Mean Q:          {mean_q:.12f}")
    print(f"Median Q:        {median_q:.12f}")
    print(f"SD Q:            {std_q:.12f}")
    print(f"Min-Max Q:       {np.min(q_values):.12f} - {np.max(q_values):.12f}")
    print()

    print("[Representative valid frame: closest to the median Q]")
    print(f"Center time:     {representative.center_sec:.6f} s")
    print(f"fa (peak):       {representative.fa_hz:.9f} Hz")
    print(f"fc (left -3 dB): {representative.fc_hz:.9f} Hz")
    print(f"fb (right -3dB): {representative.fb_hz:.9f} Hz")
    print(f"Bandwidth:       {representative.bandwidth_hz:.9f} Hz")
    print(f"Q:               {representative.q_value:.12f}")
    print()
    print(f"Frame CSV:       {csv_path}")


if __name__ == "__main__":
    main()
