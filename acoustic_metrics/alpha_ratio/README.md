# Alpha Ratio

This directory contains `calculate_alpha_ratio_patel2010.py`, a reproducible
Python implementation of Alpha Ratio from a long-term average spectrum
estimate.

## Definition used in the script

```text
Alpha Ratio (dB) = 10 * log10(E_high / E_low)
```

- `E_low`: energy from 50 to below 1000 Hz
- `E_high`: energy from 1000 to 5000 Hz

## LTAS implementation

- Stable interval: 0.430-0.930 seconds
- Welch power spectral density estimate
- Hann window
- FFT length: 2048
- Overlap: 1024 samples

## Requirements

- Python 3.10 or later
- NumPy
- SciPy
- SoundFile

```bash
pip install numpy scipy soundfile
```

## Usage

The WAV filename is set by the `WAV_FILE` constant near the beginning of the
script. Place the WAV file in the working directory or edit that constant.

```bash
python calculate_alpha_ratio_patel2010.py
```

Librosa is not used in the calculation; its installed version is reported only
for documentation.
