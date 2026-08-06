# Spectral Centroid

This directory contains
`calculate_spectral_centroid_schubert_wolfe2006.py`, a Python implementation
of the power-weighted Spectral Centroid.

## Definition

```text
Spectral Centroid = sum(frequency * power) / sum(power)
```

The script uses spectral power, not spectral amplitude, as the weight.

## Default implementation

- STFT length: 1024 samples
- Hop length: 256 samples
- Window: Hann
- Boundary zero-padding: disabled (`center=False`)
- Spectrum: one-sided power spectrum
- Aggregation: all time-frequency power values are pooled

## Requirements

- Python 3.10 or later
- NumPy
- Librosa
- SoundFile

```bash
pip install numpy librosa soundfile
```

## Usage

```bash
python calculate_spectral_centroid_schubert_wolfe2006.py input.wav
```

Optional settings:

```bash
python calculate_spectral_centroid_schubert_wolfe2006.py input.wav \
  --n-fft 1024 \
  --hop-length 256
```
