#!/usr/bin/env python3
"""
Run the CPPS procedure described by Baker et al.
through Praat-Parselmouth.

Example:
    python calculate_cpps_baker2024.py input.wav --start 0.430 --end 0.930

Main settings:
- Sound: To PowerCepstrogram
    pitch floor = 60 Hz
    time step = 0.002 s
    maximum frequency = 5000 Hz
    pre-emphasis from = 50 Hz
- PowerCepstrogram: Get CPPS
    subtract trend before smoothing = yes
    time averaging window = 0.02 s
    quefrency averaging window = 0.0005 s
    peak search pitch range = 60-1000 Hz
    tolerance = 0.05
    interpolation = Parabolic
    trend line quefrency range = 0.001-0.05 s
    trend type = Straight
    fit method = Robust
"""

from __future__ import annotations

import argparse
import platform
from pathlib import Path

import parselmouth
from parselmouth.praat import call


def calculate_cpps(sound: parselmouth.Sound) -> tuple[float, int]:
    """Calculate CPPS using Praat PowerCepstrogram and Get CPPS."""

    power_cepstrogram = call(
        sound,
        "To PowerCepstrogram",
        60.0,       # Pitch floor (Hz)
        0.002,      # Time step (s)
        5000.0,     # Maximum frequency (Hz)
        50.0,       # Pre-emphasis from (Hz)
    )

    number_of_frames = int(
        call(power_cepstrogram, "Get number of frames")
    )

    cpps_db = float(
        call(
            power_cepstrogram,
            "Get CPPS",
            True,           # Subtract trend before smoothing
            0.02,           # Time averaging window (s)
            0.0005,         # Quefrency averaging window (s)
            60.0,           # Peak-search pitch floor (Hz)
            1000.0,         # Peak-search pitch ceiling (Hz)
            0.05,           # Tolerance
            "Parabolic",    # Peak interpolation
            0.001,          # Trend-line lower quefrency (s)
            0.05,           # Trend-line upper quefrency (s)
            "Straight",     # Trend type
            "Robust",       # Fit method
        )
    )

    return cpps_db, number_of_frames


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate CPPS using the procedure described by Baker et al."
    )
    parser.add_argument("wav_file", type=Path)
    parser.add_argument(
        "--start",
        type=float,
        default=0.430,
        help="Start time of the primary analysis interval in seconds (default: 0.430).",
    )
    parser.add_argument(
        "--end",
        type=float,
        default=0.930,
        help="End time of the primary analysis interval in seconds (default: 0.930).",
    )
    args = parser.parse_args()

    if not args.wav_file.is_file():
        raise FileNotFoundError(args.wav_file)

    sound = parselmouth.Sound(str(args.wav_file))

    if sound.n_channels != 1:
        raise ValueError("The input audio file must be mono.")

    if not (0.0 <= args.start < args.end <= sound.duration):
        raise ValueError(
            f"The analysis interval must be within 0 to {sound.duration:.6f} seconds."
        )

    stable_sound = sound.extract_part(
        from_time=args.start,
        to_time=args.end,
        preserve_times=False,
    )

    stable_cpps, stable_frames = calculate_cpps(stable_sound)
    whole_cpps, whole_frames = calculate_cpps(sound)

    print("[Software environment]")
    print(f"Python:             {platform.python_version()}")
    print(f"Praat-Parselmouth:  {parselmouth.__version__}")
    print(f"Embedded Praat:     {parselmouth.PRAAT_VERSION}")
    print(f"Praat version date: {parselmouth.PRAAT_VERSION_DATE}")
    print("Librosa:            Not used")
    print()

    print("[Audio file]")
    print(f"File:                {args.wav_file.name}")
    print(f"Sampling frequency:  {sound.sampling_frequency:.0f} Hz")
    print(f"Total duration:      {sound.duration:.12f} seconds")
    print()

    print("[Primary result: stable /a/ interval]")
    print(f"Analysis interval:   {args.start:.3f}-{args.end:.3f} seconds")
    print(f"Interval duration:   {stable_sound.duration:.6f} seconds")
    print(f"Number of frames:    {stable_frames}")
    print(f"CPPS:                 {stable_cpps:.12f} dB")
    print()

    print("[Reference result: entire audio file]")
    print(f"Number of frames:    {whole_frames}")
    print(f"CPPS:                 {whole_cpps:.12f} dB")


if __name__ == "__main__":
    main()
