# Alpha Ratio

This directory contains `calculate_alpha_ratio_patel2010.py`, a reproducible Python implementation of the alpha-ratio measure using the 50–1000 Hz and 1–5 kHz frequency bands used by Patel et al. (2010). In this repository, the ratio is explicitly defined using a high-to-low energy convention.

## Definition used in the script

Alpha Ratio (dB) = 10 * log10(E_high / E_low)

- `E_low`: summed energy from 50 to below 1000 Hz
- `E_high`: summed energy from 1000 to 5000 Hz

Under this convention, a higher Alpha Ratio indicates relatively greater high-frequency energy (1–5 kHz) compared with low-frequency energy (50–1000 Hz), whereas a lower value indicates relatively less high-frequency energy.

Because the direction of the ratio is not always stated consistently across the literature, the numerator and denominator are explicitly defined here. The measure is used descriptively and is not interpreted as a direct physiological measure.

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
