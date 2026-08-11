#!/usr/bin/env python3
"""Calculate Singing Power Ratio (SPR) following Omori et al. (1996).

Definition used:
- Select a steady 4096-sample (92.9 ms at 44.1 kHz) vowel segment.
- Apply a Hann/Hanning window.
- Compute a 4096-point FFT and its power spectrum.
- Find the greatest power-spectrum bin from 0 to <2000 Hz.
- Find the greatest power-spectrum bin from 2000 to 4000 Hz (SPP).
- SPR [dB] = 10*log10(P_high / P_low).

For the supplied 'tanto' recording, the analysis window is centered at
0.680 s, in the stable /a/ portion.
"""

from __future__ import annotations

import argparse
import platform

import numpy as np
import scipy
from scipy.signal.windows import hann
import soundfile as sf

N_FFT = 4096
LOW_BAND_MAX_HZ = 2000.0
HIGH_BAND_MAX_HZ = 4000.0


def calculate_spr(
    wav_path: str,
    center_sec: float,
) -> dict[str, float | int]:
    audio, sample_rate = sf.read(wav_path, always_2d=False)

    if audio.ndim != 1:
        raise ValueError("A mono WAV file is required.")
    if sample_rate != 44100:
        raise ValueError(
            f"Omori et al. used 44.1 kHz; input is {sample_rate} Hz. "
            "This script does not resample automatically."
        )

    start_sample = int(round(center_sec * sample_rate - N_FFT / 2))
    end_sample = start_sample + N_FFT
    if start_sample < 0 or end_sample > len(audio):
        raise ValueError("The requested 4096-sample window is outside the file.")

    segment = np.asarray(audio[start_sample:end_sample], dtype=np.float64)

    # Periodic Hann window, appropriate for FFT spectral analysis.
    # Omori et al. state 'Hanning window' but do not specify the endpoint convention.
    window = hann(N_FFT, sym=False)
    spectrum = np.fft.rfft(segment * window, n=N_FFT)
    power = np.abs(spectrum) ** 2
    frequencies = np.fft.rfftfreq(N_FFT, d=1.0 / sample_rate)

    low_mask = (frequencies >= 0.0) & (frequencies < LOW_BAND_MAX_HZ)
    high_mask = (
        (frequencies >= LOW_BAND_MAX_HZ)
        & (frequencies <= HIGH_BAND_MAX_HZ)
    )

    low_indices = np.flatnonzero(low_mask)
    high_indices = np.flatnonzero(high_mask)
    low_peak_bin = int(low_indices[np.argmax(power[low_mask])])
    high_peak_bin = int(high_indices[np.argmax(power[high_mask])])

    low_peak_power = float(power[low_peak_bin])
    high_peak_power = float(power[high_peak_bin])
    spr_db = float(10.0 * np.log10(high_peak_power / low_peak_power))

    return {
        "sample_rate": sample_rate,
        "start_sample": start_sample,
        "end_sample": end_sample,
        "start_sec": start_sample / sample_rate,
        "end_sec": end_sample / sample_rate,
        "duration_ms": 1000.0 * N_FFT / sample_rate,
        "frequency_resolution_hz": sample_rate / N_FFT,
        "low_peak_frequency_hz": float(frequencies[low_peak_bin]),
        "high_peak_frequency_hz": float(frequencies[high_peak_bin]),
        "low_peak_power": low_peak_power,
        "high_peak_power": high_peak_power,
        "power_ratio": high_peak_power / low_peak_power,
        "spr_db": spr_db,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("wav_path", help="Path to the mono 44.1-kHz WAV file")
    parser.add_argument(
        "--center-sec",
        type=float,
        default=0.680,
        help="Center time of the steady 4096-sample vowel segment (default: 0.680)",
    )
    args = parser.parse_args()

    result = calculate_spr(args.wav_path, args.center_sec)

    print(f"Python:    {platform.python_version()}")
    print(f"NumPy:     {np.__version__}")
    print(f"SciPy:     {scipy.__version__}")
    print(f"SoundFile: {sf.__version__}")
    print("Librosa:   not used")
    print()
    print(f"Sampling rate:       {result['sample_rate']} Hz")
    print(f"Analysis interval:   {result['start_sec']:.12f}-{result['end_sec']:.12f} s")
    print(f"Window duration:     {result['duration_ms']:.9f} ms")
    print(f"Frequency resolution:{result['frequency_resolution_hz']:.12f} Hz/bin")
    print(f"0-2 kHz peak:        {result['low_peak_frequency_hz']:.9f} Hz")
    print(f"2-4 kHz SPP:         {result['high_peak_frequency_hz']:.9f} Hz")
    print(f"Low peak power:      {result['low_peak_power']:.15g}")
    print(f"SPP power:           {result['high_peak_power']:.15g}")
    print(f"Power ratio:         {result['power_ratio']:.15g}")
    print(f"SPR:                 {result['spr_db']:.12f} dB")


if __name__ == "__main__":
    main()
