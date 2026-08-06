# Singing Power Ratio (SPR)

This directory contains `calculate_spr_omori1996.py`, a Python implementation
of the Singing Power Ratio described by Omori et al. (1996).

## Definition

```text
SPR (dB) = 10 * log10(P_high / P_low)
```

- `P_low`: greatest power-spectrum peak from 0 to below 2000 Hz
- `P_high`: greatest power-spectrum peak from 2000 to 4000 Hz
- `P_high` is the Singing Power Peak (SPP)

## Analysis settings

- Sampling rate: 44100 Hz
- FFT length: 4096 samples
- Window: periodic Hann window
- Segment length: approximately 92.9 ms

## Requirements

- Python 3.10 or later
- NumPy
- SciPy
- SoundFile

```bash
pip install numpy scipy soundfile
```

## Usage

```bash
python calculate_spr_omori1996.py input.wav --center-sec 0.680
```

Set `--center-sec` so that the 4096-sample window lies within a steady vowel
segment. Librosa is not used.
