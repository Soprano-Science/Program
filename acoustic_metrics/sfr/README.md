# Singer's Formant Ratio (SFR)

This directory contains `calculate_sfr_excel_definition.py`, a Python script
that reproduces the SFR calculation defined in the SFR workbook.

## Definition

```text
SFR = 100 * sum(amplitude from 2.4 to 4.0 kHz)
            / sum(amplitude from 0 to 4.0 kHz)
```

## Analysis settings

- Sampling rate: 44100 Hz
- FFT length: 2048 samples
- Analysis frame: the first 2048 samples
- Window: periodic Hann window
- Numerator bins: 112-185
- Denominator bins: 0-185
- No resampling or normalization

## Requirements

- Python 3.10 or later
- NumPy

```bash
pip install numpy
```

## Usage

```bash
python calculate_sfr_excel_definition.py input.wav
```

The input must be a mono, uncompressed, 16-bit PCM WAV file sampled at
44100 Hz. Librosa is not used.
