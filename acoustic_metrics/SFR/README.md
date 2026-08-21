# Frame-wise SFR Calculation Using the Excel Definition

This repository calculates **SFR (Singing Formant Ratio)** for WAV files using the same spectral definition as the original Excel-based implementation, but applies it to **all frames** in the audio file rather than only the first frame.

The script writes frame-wise SFR values to CSV and uses the **median of valid frame-wise SFR values** as the representative SFR value for each WAV file.

## Definition

For each frame, SFR is defined as:

```text
SFR (%) = 100 × sum(amplitude in 2.4–4.0 kHz)
                / sum(amplitude in 0–4.0 kHz)
```

The implementation follows the Excel-derived calculation:

```text
Sampling rate: 44,100 Hz
Frame size: 2048 samples
FFT length: 2048 points
Window: periodic Hann window, 0.5 × (1 − cos(2πn/2048))
Frequency-bin spacing: 44,100 / 2048 = 21.533203125 Hz/bin
Numerator bins: 112–185, approximately 2411.71875–3983.64258 Hz
Denominator bins: 0–185, approximately 0–3983.64258 Hz
```

At a sampling rate of 44.1 kHz, one 2048-sample frame corresponds to approximately:

```text
2048 / 44100 = 0.0464399 s = 46.44 ms
```

## What changed from the original one-frame script?

The original script reproduced the Excel workbook for **one analysis frame**, using the first 2048 samples of the WAV file. This revised script keeps the same frame-level formula but applies it to the entire audio file.

The new procedure is:

```text
Read the whole WAV file
↓
Split the waveform into non-overlapping 2048-sample frames
↓
Calculate SFR for each frame
↓
Write all frame-wise SFR values to CSV
↓
Exclude silent/near-silent frames and short final frames from the median
↓
Use the median of the remaining valid frame-wise SFR values as the file-level SFR
```

## Input requirements

The input files should match the original Excel definition:

```text
WAV format: uncompressed PCM
Channels: mono
Bit depth: 16-bit
Sampling rate: 44,100 Hz
```

The script does not resample, normalize, or convert stereo files.

## Installation

```bash
pip install -r requirements.txt
```

The only external dependency is NumPy.

## Usage

### Single WAV file

```bash
python calculate_sfr_framewise_median.py input.wav
```

### Folder of WAV files

```bash
python calculate_sfr_framewise_median.py ./audio
```

### Folder including subfolders

```bash
python calculate_sfr_framewise_median.py ./audio --recursive
```

### Change the output folder

```bash
python calculate_sfr_framewise_median.py input.wav --output-dir sfr_results
```

### Change the silence threshold

By default, frames with RMS level at or below `-60 dBFS` are marked as silent or near-silent and excluded from the median.

```bash
python calculate_sfr_framewise_median.py input.wav --silence-threshold-dbfs -70
```

## Output files

The default output directory is:

```text
sfr_results/
```

For each WAV file, the script creates a frame-wise CSV file:

```text
<wav_name>_sfr_framewise.csv
```

It also creates a summary file:

```text
sfr_summary.csv
```

## Frame-wise CSV columns

| Column | Description |
|---|---|
| `source_file` | Input WAV filename |
| `frame_index` | Frame number, starting from 0 |
| `start_sample` | Start sample of the frame |
| `end_sample` | End sample of the frame |
| `start_time_s` | Start time in seconds |
| `end_time_s` | End time in seconds |
| `frame_length_samples` | Number of samples in the frame |
| `is_short_frame` | `1` if the frame is shorter than 2048 samples |
| `rms_dbfs` | Frame RMS level in dBFS |
| `max_abs_sample` | Maximum absolute sample value in the frame |
| `is_silent_frame` | `1` if the frame is silent or near-silent |
| `numerator_sum_2p4_4k` | Sum of spectral amplitude in 2.4–4.0 kHz |
| `denominator_sum_0_4k` | Sum of spectral amplitude in 0–4.0 kHz |
| `sfr_percent` | Frame-wise SFR value |
| `include_in_median` | `1` if the frame is included in the file-level median |
| `exclusion_reason` | Reason why the frame was excluded from the median |
| `calculation_warning` | Additional warning, such as zero-padding of a short frame |

## Summary CSV columns

The most important output column is:

```text
median_sfr_percent
```

This is the representative SFR value for the WAV file.

Other useful columns include:

| Column | Description |
|---|---|
| `n_frames_total` | Total number of frames, including the final short frame if present |
| `n_full_frames` | Number of full 2048-sample frames |
| `n_short_frames` | Number of frames shorter than 2048 samples |
| `short_frame_indices` | Indices of short frames |
| `n_silent_frames` | Number of silent or near-silent frames |
| `silent_frame_indices` | Indices of silent or near-silent frames |
| `n_frames_included_in_median` | Number of frames used for the median |
| `median_sfr_percent` | Median SFR of valid frames |
| `mean_sfr_percent` | Mean SFR of valid frames |
| `sd_sfr_percent` | Standard deviation of valid frame-wise SFR values |
| `warning` | File-level warnings |

## Treatment of silent frames and the final short frame

The script intentionally reports but excludes the following frames from the file-level median:

1. **Silent or near-silent frames**  
   Frames whose RMS level is at or below the specified silence threshold, by default `-60 dBFS`, are marked as silent or near-silent.

2. **Short final frames**  
   If the final frame contains fewer than 2048 samples, it is still written to the frame-wise CSV. For diagnostic purposes, it is zero-padded to 2048 samples before FFT calculation. However, it is excluded from the median because its spectral estimate is less stable than that of a full frame.

This means that the primary file-level SFR value is calculated as:

```text
median_sfr_percent = median(SFR values from full, non-silent frames)
```

## Recommended Methods text

The following wording can be used in a manuscript or supplementary material:

> SFR was calculated frame by frame using the Excel-derived spectral definition. Each WAV file was divided into non-overlapping 2048-sample frames. For each frame, a periodic Hann window was applied, followed by a 2048-point FFT. SFR was defined as 100 × the summed spectral amplitude in the 2.4–4.0 kHz band divided by the summed spectral amplitude in the 0–4.0 kHz band. The final frame was retained in the frame-wise CSV even when shorter than 2048 samples and was zero-padded for diagnostic FFT calculation, but it was excluded from the file-level median. Silent or near-silent frames were also excluded from the median. The representative SFR value for each WAV file was defined as the median of valid frame-wise SFR values.

## Notes on interpretation

SFR is a relative spectral-amplitude ratio. In this implementation, larger values indicate a larger proportion of amplitude in the 2.4–4.0 kHz band relative to the 0–4.0 kHz band.

Because this revised version calculates SFR across the entire WAV file, its representative value is not expected to be identical to the original one-frame result. The new value reflects the median of the time-varying SFR trajectory over the whole token.

For short words or sung syllables containing consonants, the frame-wise SFR trajectory may include rapid local changes. The median is used to reduce the influence of local outliers, silent frames, and unstable final short frames.
