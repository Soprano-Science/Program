# Frame-wise Alpha Ratio Calculation

This repository contains a Python implementation for calculating **alpha ratio** from WAV files using a **frame-wise median approach**.

The script calculates alpha ratio for each audio frame and then uses the **median of the frame-wise values** as the representative alpha ratio for each audio file.

This implementation is intended for research workflows in singing-voice analysis, especially when a short sung token such as a word, syllable, or phrase contains time-varying spectral characteristics.

---

## Script

```text
calculate_alpha_ratio_framewise_median.py
```

---

## Method Overview

The program processes each WAV file as follows:

1. Read the WAV file.
2. Convert stereo audio to mono by averaging channels.
3. Convert integer PCM audio to floating-point values.
4. Divide the whole audio signal into frames.
5. Calculate alpha ratio for each frame.
6. Use the median of all frame-wise alpha ratio values as the representative alpha ratio for the file.
7. Save both file-level summary results and optional frame-level results as CSV files.

---

## Default Parameters

| Parameter | Default value |
|---|---:|
| Frame size | 2048 samples |
| Hop size | 2048 samples |
| Overlap | None |
| Window | Hann window |
| Low-frequency band | 50-1000 Hz |
| High-frequency band | 1000-5000 Hz |
| Alpha ratio definition | `10log10(E_50_1000 / E_1000_5000)` |
| Representative value | Median of frame-wise alpha ratio values |
| DC offset removal | Enabled by default |

The final frame is retained even when it is shorter than 2048 samples. For FFT computation, the final short frame is zero-padded to the frame size.

---

## Alpha Ratio Definition

For each frame, alpha ratio is calculated as:

```text
alpha_ratio_db = 10 * log10(E_50_1000 / E_1000_5000)
```

where:

```text
E_50_1000   = summed spectral power from 50 Hz to 1000 Hz
E_1000_5000 = summed spectral power from 1000 Hz to 5000 Hz
```

With this sign convention, a **smaller alpha ratio value** indicates that the 1000-5000 Hz band is relatively stronger than the 50-1000 Hz band.

For sign checking, the program also outputs the inverse ratio:

```text
alpha_ratio_db_high_over_low_check = 10 * log10(E_1000_5000 / E_50_1000)
```

This inverse value is provided only for checking compatibility with other definitions or previous implementations.

---

## Difference from an LTAS Implementation

This script does **not** calculate alpha ratio from a single long-term average spectrum.

Instead, it calculates alpha ratio separately for each short frame and then summarizes the frame-wise values using the median.

Therefore, this method should be described as:

```text
frame-wise alpha ratio summarized by the median
```

rather than:

```text
LTAS-based alpha ratio
```

---

## Why Use the Median?

The median is robust to local outliers. In singing-voice recordings, short sections such as consonants, onsets, offsets, or low-amplitude frames may produce unusually high or low spectral ratios.

Using the median helps reduce the influence of such local fluctuations while still representing the time-varying spectral balance of the whole audio file.

This is especially useful for short sung words or phrases that include both consonants and vowels.

---

## Requirements

Python 3.9 or later is recommended.

Install the required packages with:

```bash
pip install numpy scipy
```

---

## Basic Usage

### Single WAV file

```bash
python calculate_alpha_ratio_framewise_median.py input.wav
```

This creates an output folder named:

```text
alpha_ratio_results/
```

Default output files:

```text
alpha_ratio_results/input_alpha_ratio_framewise.csv
alpha_ratio_results/alpha_ratio_summary.csv
```

---

## Folder Processing

### Process all WAV files in a folder

```bash
python calculate_alpha_ratio_framewise_median.py ./wav_files
```

### Process all WAV files recursively

```bash
python calculate_alpha_ratio_framewise_median.py ./wav_files --recursive
```

---

## Optional Arguments

| Argument | Description | Default |
|---|---|---:|
| `--output_dir` | Output directory | `alpha_ratio_results` |
| `--summary_csv` | Path for summary CSV | `<output_dir>/alpha_ratio_summary.csv` |
| `--frame_size` | Frame size in samples | `2048` |
| `--hop_size` | Hop size in samples | `2048` |
| `--low_band LOW_MIN LOW_MAX` | Low-frequency band in Hz | `50 1000` |
| `--high_band HIGH_MIN HIGH_MAX` | High-frequency band in Hz | `1000 5000` |
| `--eps` | Small constant to avoid division by zero | `1e-20` |
| `--recursive` | Search WAV files recursively | Off |
| `--keep_dc` | Do not remove DC offset | Off |
| `--no_framewise_csv` | Do not save frame-wise CSV files | Off |

---

## Examples

### Use 50% overlap

```bash
python calculate_alpha_ratio_framewise_median.py input.wav --hop_size 1024
```

### Change the output directory

```bash
python calculate_alpha_ratio_framewise_median.py input.wav --output_dir results_alpha_ratio
```

### Save only the summary CSV

```bash
python calculate_alpha_ratio_framewise_median.py input.wav --no_framewise_csv
```

### Use a custom summary CSV path

```bash
python calculate_alpha_ratio_framewise_median.py ./wav_files --recursive --summary_csv alpha_ratio_all_files.csv
```

---

## Output Files

### 1. File-level summary CSV

Default path:

```text
alpha_ratio_results/alpha_ratio_summary.csv
```

Main columns:

| Column | Description |
|---|---|
| `input_file` | Input WAV file path |
| `sampling_rate_hz` | Sampling rate |
| `duration_s` | Duration in seconds |
| `frame_size_samples` | Frame size |
| `hop_size_samples` | Hop size |
| `n_frames` | Number of analyzed frames |
| `n_short_frames` | Number of frames shorter than the frame size |
| `low_band_hz` | Low-frequency band |
| `high_band_hz` | High-frequency band |
| `alpha_definition` | Alpha ratio formula |
| `median_alpha_ratio_db` | Representative alpha ratio value for the file |
| `mean_alpha_ratio_db` | Mean of frame-wise alpha ratio values |
| `sd_alpha_ratio_db` | Standard deviation of frame-wise alpha ratio values |
| `min_alpha_ratio_db` | Minimum frame-wise alpha ratio value |
| `max_alpha_ratio_db` | Maximum frame-wise alpha ratio value |
| `framewise_csv` | Path to the frame-wise CSV file |

### 2. Frame-wise CSV

Example path:

```text
alpha_ratio_results/input_alpha_ratio_framewise.csv
```

Main columns:

| Column | Description |
|---|---|
| `frame_index` | Frame number |
| `start_sample` | Start sample of the frame |
| `end_sample` | End sample of the frame |
| `frame_length_samples` | Actual frame length |
| `start_time_s` | Frame start time in seconds |
| `end_time_s` | Frame end time in seconds |
| `rms` | RMS amplitude of the frame |
| `energy_50_1000` | Summed spectral power in the 50-1000 Hz band |
| `energy_1000_5000` | Summed spectral power in the 1000-5000 Hz band |
| `alpha_ratio_db_low_over_high` | Main alpha ratio value |
| `alpha_ratio_db_high_over_low_check` | Inverse ratio for sign checking |

---

## Recommended Method Description for a Manuscript

```text
Alpha ratio was calculated frame by frame using non-overlapping 2048-sample frames. The final frame was retained even when shorter than 2048 samples and was zero-padded for FFT computation. For each frame, a Hann window was applied, and the power spectrum was obtained by FFT. Band energies were computed by summing spectral power in the 50-1000 Hz and 1000-5000 Hz bands. Alpha ratio was defined as 10log10(E_50-1000/E_1000-5000). The median of the frame-wise alpha ratio values was used as the representative alpha ratio for each audio file.
```

---

## Notes on Interpretation

This implementation uses the low-over-high definition:

```text
10log10(E_50_1000 / E_1000_5000)
```

Therefore:

- A larger value indicates relatively stronger low-frequency energy.
- A smaller value indicates relatively stronger high-frequency energy.
- Negative values are possible when the 1000-5000 Hz band is stronger than the 50-1000 Hz band.

Because different software packages and publications may use different sign conventions, always report the exact formula used.

---

## Important Limitations

Alpha ratio is sensitive to:

- Recording level
- Microphone position
- Room acoustics
- Background noise
- Segment duration
- Presence of consonants or unvoiced frames
- Window size and overlap settings
- The selected frequency-band definition

For longitudinal comparison, recordings should be made under the same or comparable conditions whenever possible.

---

## Repository Structure Example

```text
alpha_ratio/
├── README.md
├── calculate_alpha_ratio_framewise_median.py
└── alpha_ratio_results/
    ├── alpha_ratio_summary.csv
    └── input_alpha_ratio_framewise.csv
```

---

## License

Add the license information for your repository here.
