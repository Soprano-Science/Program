# Smoothed Cepstral Peak Prominence (CPPS)

This directory contains `calculate_cpps_baker2024.py`, a Praat-Parselmouth
implementation of the CPPS procedure used by Baker et al. for singing-voice
analysis.

## Main Praat settings

### To PowerCepstrogram

- Pitch floor: 60 Hz
- Time step: 0.002 s
- Maximum frequency: 5000 Hz
- Pre-emphasis from: 50 Hz

### Get CPPS

- Subtract trend before smoothing: yes
- Time averaging window: 0.02 s
- Quefrency averaging window: 0.0005 s
- Peak-search pitch range: 60-1000 Hz
- Peak interpolation: Parabolic
- Trend-line quefrency range: 0.001-0.05 s
- Trend type: Straight
- Fit method: Robust

## Requirements

- Python 3.10 or later
- praat-parselmouth

```bash
pip install praat-parselmouth
```

## Usage

```bash
python calculate_cpps_baker2024.py input.wav --start 0.430 --end 0.930
```

The selected interval should contain a stable vowel portion and exclude onset,
offset, consonants, and vowel transitions.
