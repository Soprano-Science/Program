# Frame-wise Spectral Centroid Median and IQR Calculation

This repository provides a Python script for calculating **Spectral Centroid** from WAV files using a fixed frame-based procedure.

The script was prepared for short singing-voice tokens such as `tanto` or `cessa`, where the aim is to calculate spectral centroid over the **entire WAV file** rather than only a manually selected vowel-midpoint segment.

## Overview

For each WAV file, the script:

1. Reads the whole WAV file.
2. Converts stereo audio to mono by averaging channels.
3. Splits the signal into **2048-sample non-overlapping frames**.
4. Retains the final frame even when it is shorter than 2048 samples.
5. Removes the DC component from each frame.
6. Applies a Hann window to each actual frame.
7. Zero-pads short final frames to 2048 samples for FFT computation.
8. Calculates the spectral centroid for each frame.
9. Writes all frame-wise spectral centroid values to CSV.
10. Uses the **median of all finite frame-wise spectral centroid values** as the representative spectral centroid value for that WAV file.
11. Calculates the **IQR** as the 75th percentile minus the 25th percentile.
12. Flags likely unvoiced, low-periodicity, or very low-amplitude frames in both the frame-wise CSV and the summary CSV.

At a sampling rate of 44.1 kHz, a 2048-sample frame corresponds to approximately:

```text
2048 / 44100 = 0.04644 s = 46.4 ms
```

## Primary output values

The primary file-level value is:

```text
median_spectral_centroid_all_frames_hz
```

The primary file-level variability measure is:

```text
iqr_spectral_centroid_all_frames_hz
```

where:

```text
IQR = q75_spectral_centroid_all_frames_hz - q25_spectral_centroid_all_frames_hz
```

The median and IQR are calculated from all finite frame-wise spectral centroid values. Frames flagged as likely unvoiced or low-periodicity are **not automatically excluded** from the primary values. They are retained because the intended design is to represent the entire WAV token.

For diagnostic purposes, the script also reports:

```text
median_spectral_centroid_likely_voiced_frames_hz
iqr_spectral_centroid_likely_voiced_frames_hz
```

These diagnostic values use only frames that pass the simple voicing check. They are provided for inspection, but they are not the primary values unless the study explicitly decides to exclude unvoiced frames.

## Spectral centroid definition

For each frame, the spectral centroid is calculated as the magnitude-weighted mean frequency by default:

```text
spectral_centroid = sum(frequency_bin_hz * magnitude_bin) / sum(magnitude_bin)
```

In other words, each magnitude spectrum is treated as a distribution over frequency bins, and the centroid is the frequency-weighted average.

The default spectral weighting is:

```text
weighting = magnitude
```

A power-weighted version can be selected for diagnostic comparison:

```bash
--weighting power
```

If the weighting method is changed, it must be reported because magnitude-weighted and power-weighted centroids are not numerically identical.

## Frequency range

The default frequency range is:

```text
0 Hz to Nyquist frequency
```

You can optionally restrict the frequency range with:

```bash
--min-frequency 50 --max-frequency 5000
```

If the frequency range is changed, it must be reported because spectral centroid values depend strongly on the selected analysis band.

## Important methodological note

This script uses a **fixed 2048-sample frame-wise method**. It does not use a vowel-midpoint segment, and it does not use a longer analysis frame such as 4096 samples.

Therefore, spectral centroid values from this script should be compared only with values calculated using the same script and the same parameters:

```text
Frame size:       2048 samples
Hop size:         2048 samples
FFT size:         2048 points
Window:           Hann
Weighting:        magnitude
Final frame:      included even if shorter than 2048 samples
Summary value:    median of all finite frame-wise centroid values
Variability:      IQR = 75th percentile - 25th percentile
```

## Unvoiced or low-periodicity frame warnings

Spectral centroid can be calculated for voiced vowels, unvoiced consonants, and noise-like frames. Therefore, a frame being unvoiced does not automatically make the centroid invalid.

However, in singing-voice research, unvoiced or very low-amplitude frames may have a very different spectral structure from sustained voiced singing. This script therefore flags such frames for inspection.

The default warning criteria are:

```text
low_rms:         RMS < -60 dBFS
low_periodicity: normalized autocorrelation peak < 0.30
short_final_frame: final frame shorter than 2048 samples
```

The frame-wise CSV includes:

```text
rms_dbfs
autocorr_voicing_score
is_likely_voiced
voicing_warning
calculation_warning
```

The summary CSV includes:

```text
n_likely_unvoiced_frames
likely_unvoiced_frame_indices
warning
```

These warnings are intended to make the analysis transparent. They do not automatically remove frames from the primary all-frame median or IQR.

## Installation

Install the required packages:

```bash
pip install -r requirements.txt
```

or:

```bash
pip install numpy scipy
```

## Basic usage

### Single WAV file

```bash
python calculate_spectral_centroid_framewise_median_iqr.py input.wav
```

This creates an output folder named:

```text
spectral_centroid_results/
```

and writes:

```text
spectral_centroid_results/input_spectral_centroid_framewise.csv
spectral_centroid_results/spectral_centroid_summary.csv
```

### Folder of WAV files

```bash
python calculate_spectral_centroid_framewise_median_iqr.py ./audio
```

### Folder of WAV files, including subfolders

```bash
python calculate_spectral_centroid_framewise_median_iqr.py ./audio --recursive
```

### Specify an output directory

```bash
python calculate_spectral_centroid_framewise_median_iqr.py input.wav --output-dir results_centroid
```

### Restrict the frequency range

```bash
python calculate_spectral_centroid_framewise_median_iqr.py input.wav --min-frequency 50 --max-frequency 5000
```

### Change warning thresholds

```bash
python calculate_spectral_centroid_framewise_median_iqr.py input.wav \
  --rms-unvoiced-threshold-dbfs -55 \
  --autocorr-voicing-threshold 0.35
```

## Frame-wise CSV columns

Each file produces one frame-wise CSV file with the following columns:

| Column | Meaning |
|---|---|
| `frame_index` | Frame number, starting from 0 |
| `start_sample` | Start sample of the frame |
| `end_sample` | End sample of the frame |
| `frame_length_samples` | Actual frame length; the final frame may be shorter than 2048 samples |
| `start_time_s` | Frame start time in seconds |
| `end_time_s` | Frame end time in seconds |
| `short_final_frame` | Whether the frame is shorter than 2048 samples |
| `rms_dbfs` | RMS amplitude in dBFS |
| `autocorr_voicing_score` | Normalized autocorrelation peak used as a voicing diagnostic |
| `is_likely_voiced` | `True` if the frame passes the RMS and periodicity checks |
| `voicing_warning` | Warning labels such as `low_rms`, `low_periodicity`, or `short_final_frame` |
| `spectral_centroid_hz` | Frame-wise spectral centroid in Hz |
| `spectral_centroid_khz` | Frame-wise spectral centroid in kHz |
| `calculation_warning` | Calculation warning, such as `zero_padded_short_frame` |

## Summary CSV columns

The summary CSV contains one row per WAV file. Important columns include:

| Column | Meaning |
|---|---|
| `median_spectral_centroid_all_frames_hz` | Primary representative spectral centroid value |
| `q25_spectral_centroid_all_frames_hz` | 25th percentile of frame-wise spectral centroid |
| `q75_spectral_centroid_all_frames_hz` | 75th percentile of frame-wise spectral centroid |
| `iqr_spectral_centroid_all_frames_hz` | IQR = q75 - q25 |
| `mean_spectral_centroid_all_frames_hz` | Mean of finite frame-wise centroid values |
| `sd_spectral_centroid_all_frames_hz` | Standard deviation of finite frame-wise centroid values |
| `min_spectral_centroid_all_frames_hz` | Minimum finite frame-wise spectral centroid |
| `max_spectral_centroid_all_frames_hz` | Maximum finite frame-wise spectral centroid |
| `median_spectral_centroid_likely_voiced_frames_hz` | Diagnostic median using only likely voiced frames |
| `iqr_spectral_centroid_likely_voiced_frames_hz` | Diagnostic IQR using only likely voiced frames |
| `n_likely_unvoiced_frames` | Number of frames flagged as likely unvoiced or low-periodicity |
| `likely_unvoiced_frame_indices` | Indices of flagged frames |
| `warning` | File-level warnings |
| `framewise_csv` | Path to the frame-wise CSV file |

## Suggested Methods wording

For a manuscript or supplementary material, the method can be described as follows:

> Spectral centroid was calculated frame by frame over the entire WAV file. Each audio file was divided into non-overlapping 2048-sample frames. The final frame was retained even when shorter than 2048 samples and was zero-padded to 2048 samples for FFT computation. For each frame, the DC component was removed, a Hann window was applied, and the magnitude spectrum was computed using a 2048-point FFT. Spectral centroid was defined as the magnitude-weighted mean frequency, i.e., sum(frequency_bin_hz × magnitude_bin) / sum(magnitude_bin). The median of all finite frame-wise spectral centroid values was used as the representative value for each WAV file. The interquartile range (IQR) was calculated as the 75th percentile minus the 25th percentile of the frame-wise values. Frames with low RMS or low periodicity were flagged for inspection but were not automatically excluded from the primary all-frame median and IQR.

## Interpretation

Spectral centroid is often interpreted as a correlate of spectral brightness. Higher values generally indicate that spectral energy is distributed toward higher frequencies, whereas lower values generally indicate a lower-frequency spectral balance.

For short singing-voice tokens, however, the value may be influenced by:

- vowels and consonants included in the token,
- breath noise or fricative-like components,
- room noise,
- microphone response,
- sound level differences,
- the selected frequency range,
- the inclusion of short final frames.

For this reason, the script reports both the frame-wise data and the file-level median/IQR, and it flags frames that may be unvoiced or low-periodicity.

## Reproducibility parameters to report

When reporting results, include at least the following:

```text
Sampling rate
Frame size
Hop size
FFT size
Window type
Spectral weighting: magnitude or power
Frequency range used for centroid calculation
Whether the final short frame was included
Whether zero padding was used
Whether unvoiced/low-periodicity frames were retained
Summary statistic: median
Variability statistic: IQR
Software and script version
```
