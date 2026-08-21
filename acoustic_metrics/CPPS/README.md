# Frame-wise CPPS Median Calculation

This repository provides a Python script for calculating **Cepstral Peak Prominence Smoothed (CPPS)** from WAV files using a fixed frame-based procedure.

The script was prepared for short singing-voice tokens such as `tanto` or `cessa`, where the aim is to calculate CPPS over the **entire WAV file** rather than only a manually selected vowel-midpoint segment.

## Overview

For each WAV file, the script:

1. Reads the whole WAV file.
2. Converts stereo audio to mono by averaging channels.
3. Splits the signal into **2048-sample non-overlapping frames**.
4. Retains the final frame even when it is shorter than 2048 samples.
5. Applies a Hann window to each actual frame.
6. Zero-pads short final frames to 2048 samples for FFT computation.
7. Calculates CPPS for each frame using Praat's `PowerCepstrum: Get peak prominence` command through Parselmouth.
8. Writes all frame-wise CPPS values to CSV.
9. Uses the **median of all finite frame-wise CPPS values** as the representative CPPS value for that WAV file.
10. Flags likely unvoiced or low-periodicity frames in both the frame-wise CSV and the summary CSV.

At a sampling rate of 44.1 kHz, a 2048-sample frame corresponds to approximately:

```text
2048 / 44100 = 0.04644 s = 46.4 ms
```

## Primary output value

The primary file-level CPPS value is:

```text
median_cpps_all_frames_db
```

This value is the median of all finite frame-wise CPPS values, including frames that are flagged as likely unvoiced or low-periodicity. This follows the intended design of representing the entire WAV token.

For diagnostic purposes, the script also reports:

```text
median_cpps_voiced_frames_db
```

This value uses only frames that pass the voicing check. It is provided for inspection, but it is not the primary value unless the study explicitly decides to exclude unvoiced frames.

## Important methodological note

This script uses a **fixed 2048-sample frame-wise method**. It does not reproduce Praat's whole-signal `PowerCepstrogram: Get CPPS` calculation exactly.

Praat's `Sound: To PowerCepstrogram` uses an analysis window whose effective length is determined by the pitch floor. Praat's `PowerCepstrogram: Get CPPS` then returns the average of cepstral peak prominences over the frames in the selected PowerCepstrogram. In contrast, the present script deliberately fixes each analysis frame to 2048 samples and summarizes the resulting frame-wise CPPS values by the median.

Therefore, CPPS values from this script should be compared only with values calculated using the same script and the same parameters. They should not be directly compared with published clinical cutoff values or with CPPS values from Praat, ADSV, VoiceSauce, or other software unless the implementation differences are explicitly addressed.

## CPPS calculation in each frame

For each frame, the script performs the following steps:

1. Remove the DC component.
2. Apply a Hann window.
3. Zero-pad to 2048 samples when necessary.
4. Convert the frame to a Praat `Sound` object.
5. Convert the sound to a `Spectrum`.
6. Convert the spectrum to a `PowerCepstrum`.
7. Smooth the `PowerCepstrum` in the quefrency domain.
8. Search for the cepstral peak within the specified pitch range.
9. Fit a trend line to the cepstral background.
10. Calculate CPPS as the difference between the cepstral peak and the fitted trend line at the same quefrency.

The default settings are:

```text
Frame size:                 2048 samples
Hop size:                   2048 samples
Window:                     Hann
Final short frame:          included
Zero padding:               to 2048 samples
Pitch floor:                60 Hz
Pitch ceiling:              1000 Hz
Quefrency smoothing window: 0.0005 s
Smoothing iterations:       1
Peak interpolation:         parabolic
Trend-line quefrency range: 0.001–0.05 s
Trend type:                 straight
Fit method:                 robust slow
Representative value:       median of all finite frame-wise CPPS values
```

## Unvoiced or low-periodicity frame detection

The script flags likely unvoiced or low-periodicity frames. A frame is flagged when either of the following conditions is met:

```text
RMS level < -60 dBFS
```

or

```text
normalized autocorrelation peak < 0.30
```

The default warning labels are:

```text
low_rms
low_periodicity
short_final_frame
```

A flagged frame is still analyzed when possible, and its CPPS value is still included in `median_cpps_all_frames_db`. The purpose of the flag is to make the presence of unvoiced or unstable frames transparent.

This is especially important for short sung words because consonants, releases, breathy onsets, and very low-amplitude endings may produce lower or unstable CPPS values. The summary CSV therefore reports:

```text
n_unvoiced_frames
unvoiced_frame_indices
warning
```

If unvoiced frames are detected, the console output also prints a warning.

## Installation

Install the required packages:

```bash
pip install -r requirements.txt
```

The required packages are:

```text
numpy
scipy
praat-parselmouth
```

## Usage

### Analyze a single WAV file

```bash
python calculate_cpps_framewise_median.py input.wav
```

By default, output files are written to:

```text
cpps_results/
```

### Analyze a single WAV file and specify the output directory

```bash
python calculate_cpps_framewise_median.py input.wav --output-dir cpps_results
```

### Analyze all WAV files in a folder

```bash
python calculate_cpps_framewise_median.py ./audio
```

### Analyze all WAV files in a folder recursively

```bash
python calculate_cpps_framewise_median.py ./audio --recursive
```

### Change pitch search range

For high female classical singing, the default pitch ceiling of 1000 Hz usually covers the target range. If needed, it can be changed:

```bash
python calculate_cpps_framewise_median.py input.wav --pitch-floor 80 --pitch-ceiling 1200
```

### Change unvoiced-frame thresholds

```bash
python calculate_cpps_framewise_median.py input.wav \
  --rms-unvoiced-threshold-dbfs -55 \
  --autocorr-voicing-threshold 0.35
```

## Output files

### 1. Summary CSV

The summary file is:

```text
cpps_results/cpps_summary.csv
```

Main columns:

| Column | Description |
|---|---|
| `input_file` | Input WAV file path |
| `sampling_rate_hz` | Sampling rate |
| `duration_s` | Duration of the WAV file |
| `frame_size_samples` | Frame size, default 2048 samples |
| `hop_size_samples` | Hop size, default 2048 samples |
| `frame_duration_ms` | Frame duration in milliseconds |
| `n_frames_total` | Total number of frames, including the final short frame |
| `n_voiced_frames` | Number of frames not flagged as low-rms or low-periodicity |
| `n_unvoiced_frames` | Number of flagged frames |
| `unvoiced_frame_indices` | Frame indices of flagged frames |
| `median_cpps_all_frames_db` | Primary representative CPPS value |
| `mean_cpps_all_frames_db` | Mean of all finite frame-wise CPPS values |
| `sd_cpps_all_frames_db` | Standard deviation of all finite frame-wise CPPS values |
| `median_cpps_voiced_frames_db` | Median using only frames that passed the voicing check |
| `warning` | Summary warning |

### 2. Frame-wise CSV

For each WAV file, a frame-wise CSV is created:

```text
<filename>_cpps_framewise.csv
```

Main columns:

| Column | Description |
|---|---|
| `frame_index` | Frame number |
| `start_sample` | Start sample of the frame |
| `end_sample` | End sample of the frame |
| `frame_length_samples` | Actual frame length before zero-padding |
| `start_time_s` | Frame start time in seconds |
| `end_time_s` | Frame end time in seconds |
| `short_final_frame` | 1 if the frame is shorter than 2048 samples |
| `rms_dbfs` | RMS level in dBFS |
| `autocorr_voicing_score` | Normalized autocorrelation peak |
| `is_voiced` | 1 if the frame passed the voicing check |
| `voicing_warning` | Warning labels such as `low_rms` or `low_periodicity` |
| `cpps_db` | Frame-wise CPPS value |
| `peak_quefrency_s` | Quefrency of the detected cepstral peak |
| `peak_frequency_hz` | Frequency corresponding to the detected quefrency |
| `calculation_warning` | Praat/Parselmouth or zero-padding warnings |

## Suggested Methods wording

The following text may be adapted for a manuscript or supplementary material:

> CPPS was calculated frame by frame over the entire WAV file. Each recording was divided into non-overlapping 2048-sample frames. The final frame was retained even when shorter than 2048 samples and was zero-padded to 2048 samples for FFT-based cepstral computation. A Hann window was applied to each frame. For each frame, CPPS was obtained using Praat's PowerCepstrum peak-prominence calculation through Parselmouth. The cepstral peak was searched within the specified pitch range, and peak prominence was defined as the difference between the cepstral peak and the fitted cepstral trend line at the same quefrency. The median of all finite frame-wise CPPS values was used as the representative CPPS value for each WAV file. Frames with low RMS level or low autocorrelation-based periodicity were flagged as likely unvoiced or low-periodicity frames and reported in the output CSV files.

## Interpretation notes

Higher CPPS generally indicates a stronger harmonic structure and clearer periodic organization of the voice signal. Lower CPPS may reflect breathiness, noise, weak periodicity, consonantal segments, unvoiced frames, or low-amplitude signal portions.

For singing-voice research, CPPS should not be interpreted in isolation. It is best interpreted together with other measures such as HNR, H1-H2, alpha ratio, SPR, spectral centroid, SFR/Q, F0 stability, and expert auditory-perceptual evaluation.

## References

- Praat manual: `PowerCepstrum: Get peak prominence` describes CPP as the difference between the cepstral peak and the corresponding value on the trend line below the peak.
- Praat manual: `PowerCepstrogram: Get CPPS` describes CPPS as the average of cepstral peak prominences across frames in a PowerCepstrogram.
- Praat manual: `Sound: To PowerCepstrogram` notes that the effective analysis-window length is determined by the pitch floor.
- Maryn et al. and related clinical voice literature discuss CPPS as a widely used objective voice-quality measure, but implementation details differ across software packages.
