# Y_Q_value: Frame-wise Yamashita-style LPC-Based Q Value

`Y_Q_value.py` calculates an LPC-based Q value for singing-voice WAV files using the procedure illustrated in Yamashita's Excel workbook `1_Cessa_mono_Q.xlsx`.

The original workbook demonstrates Q-value calculation for one 2048-sample frame. This script generalizes that workflow to the whole WAV file by calculating Q values for all complete 2048-sample frames and writing the frame-wise time series to CSV.

## Correction note: peak search range

The Yamashita Excel workbook uses a peak search range of **2400-4000 Hz** for this Q-value workflow. An earlier draft of this repository mistakenly inherited **2000-4000 Hz** from a different LPC-Q script. The present `Y_Q_value.py` uses **2400-4000 Hz** by default to match the workbook definition.

## Purpose

The program is intended for acoustic analysis of classical singing voice segments, such as `tanto` or `cessa`, where the researcher wants to quantify the sharpness of an LPC spectral-envelope peak in the singer's-formant-related region.

For each valid frame, the script estimates:

- `fa`: LPC-envelope peak frequency
- `fc`: lower -3 dB crossing frequency
- `fb`: upper -3 dB crossing frequency
- `Q = fa / (fb - fc)`

The representative Q value for each WAV file is the **median of valid frame-wise Q values**.

## Analysis settings

Default settings are:

| Parameter | Default |
|---|---:|
| Sampling rate | 44,100 Hz |
| Analysis region | Whole WAV file |
| Frame size | 2048 samples |
| Hop size | 2048 samples |
| Overlap | None |
| Window | Periodic Hann window |
| LPC order | 12 |
| LPC method | Autocorrelation / Yule-Walker |
| LPC-envelope frequency grid | 1024 one-sided bin-centered frequencies |
| Frequency grid formula | `(k + 0.5) * sample_rate / 2048` |
| Peak search range | 2400-4000 Hz |
| -3 dB crossing method | Linear interpolation |
| Representative file-level value | Median of valid frame-wise Q values |

At 44.1 kHz, one 2048-sample frame corresponds to approximately 46.44 ms.

## Definition

For each complete non-silent frame:

1. Multiply the 2048-sample frame by the periodic Hann window:

```text
w[n] = 0.5 * (1 - cos(2*pi*n/2048))
```

2. Estimate LPC coefficients by the autocorrelation / Yule-Walker method.

3. Evaluate the LPC spectral envelope on a 2048-point one-sided bin-centered frequency grid.

4. Search for candidate local peaks in the 2400-4000 Hz range.

5. For each candidate peak, find the lower and upper frequencies where the LPC envelope falls 3 dB below the peak.

6. Calculate:

```text
Q = fa / (fb - fc)
```

where:

- `fa` = selected LPC-envelope peak frequency
- `fc` = lower -3 dB crossing frequency
- `fb` = upper -3 dB crossing frequency
- `fb - fc` = -3 dB bandwidth

If more than one valid candidate peak is found, the default behavior is to select the candidate with the largest Q value. This can be changed to selecting the highest LPC-envelope peak by using `--peak-selection highest-peak`.

## Frame exclusion policy

The following frames are recorded in the frame-wise CSV but excluded from Q calculation and from the median:

1. Frames shorter than 2048 samples, usually the final incomplete frame.
2. Silent or near-silent frames.

By default, a frame is treated as silent or near-silent if:

```text
RMS <= -60 dBFS
```

Frames that are complete and non-silent are analyzed. However, some frames may still have no valid Q value if no valid 2400-4000 Hz LPC peak and two -3 dB crossings can be found. Such frames are also excluded from the median and marked as `no_valid_Q` in the CSV.

## Installation

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Required packages:

```text
numpy
scipy
soundfile
```

## Usage

### Single WAV file

```bash
python Y_Q_value.py input.wav
```

### Folder of WAV files

```bash
python Y_Q_value.py ./audio
```

### Include subfolders

```bash
python Y_Q_value.py ./audio --recursive
```

### Specify output folder

```bash
python Y_Q_value.py input.wav --output-dir Y_Q_value_results
```

### Select the highest LPC-envelope peak instead of the largest-Q candidate

```bash
python Y_Q_value.py input.wav --peak-selection highest-peak
```

### Allow sampling rates other than 44.1 kHz

```bash
python Y_Q_value.py input.wav --allow-non-44100
```

No resampling is performed. For reproducibility with the Yamashita Excel example, 44.1 kHz WAV files are recommended.

## Output files

The default output folder is:

```text
Y_Q_value_results/
```

For each WAV file, the script writes a frame-wise CSV:

```text
<wav_stem>_Y_Q_value_framewise.csv
```

It also writes a summary CSV:

```text
Y_Q_value_summary.csv
```

## Frame-wise CSV columns

Important columns include:

| Column | Meaning |
|---|---|
| `frame_index` | Frame number |
| `start_time_s` | Frame start time in seconds |
| `center_time_s` | Frame center time in seconds |
| `end_time_s` | Frame end time in seconds |
| `rms_dbfs` | Frame RMS level in dBFS |
| `is_short_frame` | Whether the frame has fewer than 2048 samples |
| `is_silent_frame` | Whether the frame is silent or near-silent |
| `include_in_median` | Whether the Q value is included in the file-level median |
| `exclusion_reason` | Reason for exclusion, if any |
| `fa_hz` | Selected LPC-envelope peak frequency |
| `fc_hz` | Lower -3 dB crossing frequency |
| `fb_hz` | Upper -3 dB crossing frequency |
| `bandwidth_hz` | `fb_hz - fc_hz` |
| `q_value` | Frame-wise Q value |
| `status` | `valid`, `no_valid_Q`, or exclusion/error status |

## Summary CSV columns

Important summary columns include:

| Column | Meaning |
|---|---|
| `n_valid_q_frames` | Number of frames included in the median |
| `n_short_frames_excluded` | Number of short frames excluded |
| `n_silent_frames_excluded` | Number of silent/near-silent frames excluded |
| `n_no_valid_q_frames` | Number of complete non-silent frames with no valid Q value |
| `median_Y_Q_value` | Main representative Q value for the WAV file |
| `mean_Y_Q_value` | Mean of valid frame-wise Q values |
| `sd_Y_Q_value` | Standard deviation of valid frame-wise Q values |
| `iqr_Y_Q_value` | Interquartile range of valid frame-wise Q values |

## Suggested Methods wording

```text
Y_Q_value was calculated frame by frame using 2048-sample non-overlapping frames. Each frame was multiplied by a periodic Hann window, and LPC coefficients were estimated by the autocorrelation/Yule-Walker method with LPC order 12. The LPC spectral envelope was evaluated on a 2048-point one-sided bin-centered frequency grid. Within the 2400-4000 Hz range, local LPC-envelope peaks were searched, and the lower and upper -3 dB crossing frequencies were estimated by linear interpolation. Q was defined as fa/(fb - fc), where fa is the selected peak frequency and fc and fb are the lower and upper -3 dB crossing frequencies. Short frames and silent or near-silent frames were excluded from Q calculation and from the file-level median. The median of valid frame-wise Q values was used as the representative Y_Q_value for each WAV file.
```

## Interpretation notes

`Y_Q_value` should be interpreted as a frame-wise estimate of the sharpness of a selected LPC spectral-envelope peak in the 2400-4000 Hz region. It is an operational acoustic index, not a direct physiological measurement.

A missing or blank median does not necessarily mean that the audio file is invalid. It means that, under the selected parameters, no complete non-silent frame yielded both a valid 2400-4000 Hz LPC peak and two -3 dB crossing frequencies.

## Reproducibility notes

For reproducibility, report at least:

- sampling rate
- frame size
- hop size
- window type
- LPC order
- LPC method
- LPC-envelope frequency grid
- peak search range
- peak selection policy
- silence threshold
- frame exclusion policy
- representative statistic, e.g., median of valid frame-wise values
