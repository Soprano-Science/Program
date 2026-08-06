# LPC-Based Q Value

This directory contains `calculate_q_lpc_jstage2014.py`, a Python script for
calculating a resonance Q value from an LPC spectral envelope.

## Definition

```text
Q = fa / (fb - fc)
```

- `fa`: peak frequency
- `fc`: left -3 dB crossing
- `fb`: right -3 dB crossing

## Analysis settings in the script

- Sampling rate: 44100 Hz
- Stable interval: 0.430-0.930 seconds
- Frame length: 30 ms
- Frame shift: 10 ms
- Window: Hamming
- LPC order: 12
- Peak search range: 2000-4000 Hz
- LPC method: autocorrelation and Yule-Walker

## Requirements

- Python 3.10 or later
- NumPy
- SciPy
- SoundFile
- Librosa is optional and is not used in the primary calculation

```bash
pip install numpy scipy soundfile
```

## Usage

The WAV filename is set by the `WAV_FILE` constant near the beginning of the
script. Place the WAV file in the working directory or edit that constant.

```bash
python calculate_q_lpc_jstage2014.py
```

The script also writes frame-level results to
`kawano_after_q_frames.csv`.
