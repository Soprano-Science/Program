# H1-H2

This directory contains `calculate_h1_h2_holmberg1995.py`, a Python
implementation of H1-H2 analysis based on Holmberg et al. (1995).

## Measure

```text
H1-H2 = amplitude of the first harmonic - amplitude of the second harmonic
```

The script also reports an F1-corrected value:

```text
H1*-H2*
```

## Analysis settings

- Original sampling rate: 44100 Hz
- Analysis sampling rate: 10000 Hz
- Low-pass cutoff: 4500 Hz
- Analysis frame: 51.2 ms at the vowel midpoint
- Window: Hamming
- F0 and F1 estimation: Praat through Praat-Parselmouth

## Requirements

- Python 3.10 or later
- NumPy
- SciPy
- SoundFile
- praat-parselmouth

```bash
pip install numpy scipy soundfile praat-parselmouth
```

## Usage

The example WAV filename and analysis center time are set near the beginning
of the script.

```bash
python calculate_h1_h2_holmberg1995.py
```

For a different recording, edit `WAV_FILE` and `CENTER_SEC`.

## Important note

The F1-corrected value should be interpreted cautiously when H2 is above F1,
because the correction is sensitive to the estimated F1 frequency.
