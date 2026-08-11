"""Calculate SFR exactly as defined in the SFR definition workbook.

Definition for one FFT frame:
    SFR = 100 * sum(amplitude in 2.4-4.0 kHz)
                / sum(amplitude in 0-4.0 kHz)

The workbook uses:
- sampling rate: 44,100 Hz
- FFT length: 2,048 samples
- first 2,048 samples of the WAV as the analysis frame
- periodic Hann window: 0.5 * (1 - cos(2*pi*n/2048))
- frequency-bin spacing: 44,100 / 2,048 = 21.533203125 Hz
- numerator bins: 112..185 (2411.71875..3983.642578125 Hz)
- denominator bins: 0..185 (0..3983.642578125 Hz)

No librosa is used. Dependencies: Python standard library + NumPy.
"""

from __future__ import annotations

import argparse
import wave
from pathlib import Path

import numpy as np

N_FFT = 2048
EXPECTED_SR = 44100
LOW_HZ = 2400.0
HIGH_HZ = 4000.0


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
        raise ValueError(f"Uncompressed PCM WAV required; found {compression}.")

    samples = np.frombuffer(raw, dtype="<i2").astype(np.float64)
    return samples, sample_rate


def calculate_sfr_excel_definition(
    samples: np.ndarray,
    sample_rate: int,
    n_fft: int = N_FFT,
) -> tuple[float, float, float]:
    """Return (SFR percent, numerator sum, denominator sum) for the first frame."""
    if sample_rate != EXPECTED_SR:
        raise ValueError(
            f"The Excel definition assumes {EXPECTED_SR} Hz; found {sample_rate} Hz. "
            "No resampling is performed."
        )
    if samples.size < n_fft:
        raise ValueError(f"At least {n_fft} samples are required; found {samples.size}.")

    # Excel uses the first 2,048 samples for its one-frame example.
    frame = samples[:n_fft]

    # Exact Excel formula: 0.5*(1-COS(2*PI()*n/2048)).
    # This is a periodic Hann window. np.hanning(n_fft) is not used because
    # np.hanning uses n_fft-1 in the denominator and would not match exactly.
    n = np.arange(n_fft, dtype=np.float64)
    window = 0.5 * (1.0 - np.cos(2.0 * np.pi * n / n_fft))
    spectrum = np.fft.fft(frame * window, n=n_fft)

    # Reproduce Excel columns H and J exactly:
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
    if denominator == 0.0:
        raise ZeroDivisionError("The 0-4.0 kHz denominator is zero.")

    sfr_percent = 100.0 * numerator / denominator
    return sfr_percent, numerator, denominator


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate Excel-defined one-frame SFR.")
    parser.add_argument("wav_file", type=Path, help="Path to a mono 16-bit, 44.1-kHz PCM WAV")
    args = parser.parse_args()

    samples, sample_rate = read_mono_pcm16_wav(args.wav_file)
    sfr, numerator, denominator = calculate_sfr_excel_definition(samples, sample_rate)

    print(f"File: {args.wav_file}")
    print(f"Sampling rate: {sample_rate} Hz")
    print(f"Total samples: {samples.size}")
    print(f"Analyzed samples: 0-{N_FFT - 1} ({N_FFT} samples)")
    print(f"Frequency resolution: {sample_rate / N_FFT:.12f} Hz/bin")
    print(f"Numerator sum (2.4-4.0 kHz): {numerator:.15f}")
    print(f"Denominator sum (0-4.0 kHz): {denominator:.15f}")
    print(f"SFR: {sfr:.15f}")


if __name__ == "__main__":
    main()
