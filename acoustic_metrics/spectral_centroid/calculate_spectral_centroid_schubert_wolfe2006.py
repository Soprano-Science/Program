from __future__ import annotations

import argparse
import platform
from pathlib import Path

import librosa
import numpy as np
import soundfile as sf


DEFAULT_N_FFT = 1024
DEFAULT_HOP_LENGTH = 256


def power_weighted_spectral_centroid(
    audio: np.ndarray,
    sample_rate: int,
    *,
    n_fft: int = DEFAULT_N_FFT,
    hop_length: int = DEFAULT_HOP_LENGTH,
) -> dict[str, float | int]:
    """Calculate Schubert & Wolfe's power-spectrum centroid.

    The paper defines spectral centroid as

        fc = sum(f_i * P_i) / sum(P_i),

    where P_i is spectral power at frequency f_i.

    For a multi-frame recording, this function applies the same equation to
    all time-frequency bins together. This is equivalent to calculating each
    frame centroid and averaging the frame values with frame power as weights.
    """

    if audio.ndim != 1:
        raise ValueError("A mono audio signal is required.")
    if len(audio) < n_fft:
        raise ValueError(
            f"The signal has {len(audio)} samples, fewer than n_fft={n_fft}."
        )
    if n_fft <= 0 or hop_length <= 0:
        raise ValueError("n_fft and hop_length must be positive integers.")

    # center=False prevents zero-padding before and after the fragment.
    # librosa's 'hann' window is a periodic Hann window for FFT analysis.
    stft = librosa.stft(
        y=np.asarray(audio, dtype=np.float64),
        n_fft=n_fft,
        hop_length=hop_length,
        win_length=n_fft,
        window="hann",
        center=False,
    )

    # Schubert & Wolfe explicitly define the centroid from the POWER spectrum.
    power = np.abs(stft) ** 2
    frequencies = librosa.fft_frequencies(sr=sample_rate, n_fft=n_fft)

    numerator = float(np.sum(frequencies[:, np.newaxis] * power))
    denominator = float(np.sum(power))

    if not np.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("Spectral power is zero or invalid.")

    centroid_hz = numerator / denominator

    frame_power = np.sum(power, axis=0)
    frame_centroids = np.divide(
        np.sum(frequencies[:, np.newaxis] * power, axis=0),
        frame_power,
        out=np.full(frame_power.shape, np.nan, dtype=np.float64),
        where=frame_power > 0.0,
    )

    valid = np.isfinite(frame_centroids) & (frame_power > 0.0)
    simple_frame_mean_hz = float(np.mean(frame_centroids[valid]))
    power_weighted_frame_mean_hz = float(
        np.average(frame_centroids[valid], weights=frame_power[valid])
    )

    frame_count = int(power.shape[1])
    used_samples = (frame_count - 1) * hop_length + n_fft

    return {
        "centroid_hz": float(centroid_hz),
        "numerator": numerator,
        "denominator": denominator,
        "frame_count": frame_count,
        "used_samples": int(used_samples),
        "unused_samples": int(len(audio) - used_samples),
        "frame_duration_ms": 1000.0 * n_fft / sample_rate,
        "hop_duration_ms": 1000.0 * hop_length / sample_rate,
        "frequency_resolution_hz": sample_rate / n_fft,
        "simple_frame_mean_hz": simple_frame_mean_hz,
        "power_weighted_frame_mean_hz": power_weighted_frame_mean_hz,
    }


def whole_fragment_fft_centroid(
    audio: np.ndarray,
    sample_rate: int,
) -> float:
    """Sensitivity check: one FFT of the complete fragment, no window."""

    spectrum = np.fft.rfft(np.asarray(audio, dtype=np.float64))
    power = np.abs(spectrum) ** 2
    frequencies = np.fft.rfftfreq(len(audio), d=1.0 / sample_rate)

    denominator = float(np.sum(power))
    if not np.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("Spectral power is zero or invalid.")

    return float(np.sum(frequencies * power) / denominator)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate power-spectrum spectral centroid using the definition "
            "in Schubert & Wolfe (2006)."
        )
    )
    parser.add_argument("wav_path", type=Path, help="Path to a mono WAV file")
    parser.add_argument("--n-fft", type=int, default=DEFAULT_N_FFT)
    parser.add_argument("--hop-length", type=int, default=DEFAULT_HOP_LENGTH)
    args = parser.parse_args()

    audio, sample_rate = sf.read(args.wav_path, always_2d=False)
    audio = np.asarray(audio, dtype=np.float64)

    if audio.ndim != 1:
        raise ValueError("The WAV file must be mono.")

    result = power_weighted_spectral_centroid(
        audio,
        sample_rate,
        n_fft=args.n_fft,
        hop_length=args.hop_length,
    )
    direct_centroid = whole_fragment_fft_centroid(audio, sample_rate)

    print("[Software]")
    print(f"Python:    {platform.python_version()}")
    print(f"NumPy:     {np.__version__}")
    print(f"Librosa:   {librosa.__version__}")
    print(f"SoundFile: {sf.__version__}")
    print()

    print("[Audio]")
    print(f"Sampling rate: {sample_rate} Hz")
    print(f"Samples:       {len(audio)}")
    print(f"Duration:      {len(audio) / sample_rate:.12f} s")
    print()

    print("[STFT settings]")
    print(f"n_fft:                {args.n_fft}")
    print(f"hop_length:           {args.hop_length}")
    print("window:               Hann")
    print("center:               False")
    print("spectral weighting:   power = |STFT|^2")
    print(f"frames:               {result['frame_count']}")
    print(f"frame duration:       {result['frame_duration_ms']:.12f} ms")
    print(f"hop duration:         {result['hop_duration_ms']:.12f} ms")
    print(f"frequency resolution: {result['frequency_resolution_hz']:.12f} Hz")
    print(f"unused tail samples:  {result['unused_samples']}")
    print()

    print("[Schubert-Wolfe spectral centroid]")
    print(f"Numerator sum(f*P):   {result['numerator']:.15f}")
    print(f"Denominator sum(P):   {result['denominator']:.15f}")
    print(f"Spectral centroid:    {result['centroid_hz']:.12f} Hz")
    print()

    print("[Checks]")
    print(
        "Power-weighted mean of frame centroids: "
        f"{result['power_weighted_frame_mean_hz']:.12f} Hz"
    )
    print(
        "Unweighted mean of frame centroids:     "
        f"{result['simple_frame_mean_hz']:.12f} Hz"
    )
    print(
        "One FFT of complete fragment:           "
        f"{direct_centroid:.12f} Hz"
    )


if __name__ == "__main__":
    main()
