# Alpha Ratio

This directory contains `calculate_alpha_ratio_patel2010.py`, a reproducible Python implementation of the alpha-ratio measure using the 50–1000 Hz and 1–5 kHz frequency bands described by Patel et al. (2010). In this repository, alpha ratio is explicitly defined as the low-frequency-to-high-frequency spectral-energy ratio.

## Definition used in the script

`Alpha Ratio (dB) = 10 * log10(E_low / E_high)`

- `E_low`: summed spectral energy from 50 Hz to below 1000 Hz
- `E_high`: summed spectral energy from 1000 Hz to 5000 Hz

Under this convention, a lower alpha-ratio value indicates relatively greater high-frequency energy in the 1–5 kHz band compared with low-frequency energy in the 50–1000 Hz band, whereas a higher value indicates relatively greater low-frequency energy.

Because the direction of the alpha-ratio convention is not always stated consistently across the literature, the numerator and denominator are explicitly defined here. The measure is used descriptively as a measure of spectral energy balance and is not interpreted as a direct measure of glottal adduction or any other physiological state.

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
