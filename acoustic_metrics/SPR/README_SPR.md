# Frame-wise Singing Power Ratio (SPR)

This repository contains a Python implementation for calculating **Singing Power Ratio (SPR)** from WAV files using a frame-wise approach.

The analysis settings are aligned with a 2048-sample frame design:

- **Frame size:** 2048 samples
- **FFT size:** 2048 points
- **Window:** Hann window
- **Hop size:** 2048 samples by default, i.e., non-overlapping frames
- **Final frame:** retained even when shorter than 2048 samples
- **Representative file-level value:** median of frame-wise SPR values

At a sampling rate of 44.1 kHz, 2048 samples correspond to approximately **46.4 ms**:

```text
2048 / 44100 = 0.04644 s = 46.4 ms
```

Therefore, this implementation should **not** be described as using a 92.9-ms segment unless the frame size is changed to 4096 samples. The value 92.9 ms comes from:

```text
4096 / 44100 = 0.09288 s = 92.9 ms
```

## Rationale

Earlier versions of some SPR scripts use a 4096-point FFT. If the actual frame length is also 4096 samples, the analysis window is approximately 92.9 ms at 44.1 kHz. However, if the study design defines each analysis frame as 2048 samples, then using a 2048-point FFT is clearer and more consistent with the frame duration.

For broad frequency-band analysis such as SPR, a 2048-point FFT at 44.1 kHz gives a frequency-bin spacing of approximately 21.5 Hz:

```text
44100 / 2048 = 21.53 Hz
```

This is sufficiently fine for comparing energy or peaks in the 0-2 kHz and 2-4 kHz frequency ranges.

## SPR definition

The default implementation follows the commonly used peak-ratio interpretation of SPR:

```text
SPR_dB = 10 * log10(P_peak_2_4kHz / P_peak_0_2kHz)
```

where:

- `P_peak_0_2kHz` is the greatest spectral peak power in the 0-2 kHz band.
- `P_peak_2_4kHz` is the greatest spectral peak power in the 2-4 kHz band.

Because power is calculated as squared FFT magnitude, this is equivalent to the dB difference between the corresponding peak amplitudes.

A higher SPR value indicates that the strongest component in the 2-4 kHz band is stronger relative to the strongest component in the 0-2 kHz band. In singing-voice research, this high-frequency band is often associated with resonant vocal quality and projection.

## Alternative band-power mode

The script also provides an optional band-power mode:

```bash
python calculate_spr_framewise_median.py input.wav --method bandpower
```

In this mode, SPR is calculated as:

```text
SPR_dB = 10 * log10(E_2_4kHz / E_0_2kHz)
```

where:

- `E_0_2kHz` is summed spectral power in the 0-2 kHz band.
- `E_2_4kHz` is summed spectral power in the 2-4 kHz band.

Use this mode only if the intended study definition is a broad-band power ratio rather than the original peak-ratio style.

## Requirements

Install the required Python packages:

```bash
pip install numpy scipy
```

The script was written for Python 3.10 or later, but it should also work with recent Python 3 versions that support type annotations used in the code.

## Basic usage

Analyze a single WAV file:

```bash
python calculate_spr_framewise_median.py input.wav
```

Analyze all WAV files in a folder:

```bash
python calculate_spr_framewise_median.py ./audio
```

Analyze all WAV files recursively:

```bash
python calculate_spr_framewise_median.py ./audio --recursive
```

Specify an output directory:

```bash
python calculate_spr_framewise_median.py input.wav --output-dir spr_results
```

## Default command used for this study

```bash
python calculate_spr_framewise_median.py input.wav \
    --frame-size 2048 \
    --hop-size 2048 \
    --n-fft 2048 \
    --method peak
```

These options are already the default settings, so the shorter command is equivalent:

```bash
python calculate_spr_framewise_median.py input.wav
```

## Output files

The script creates an output directory named `spr_results` by default.

For each WAV file, it writes a frame-wise CSV file:

```text
<filename>_spr_framewise.csv
```

It also writes a file-level summary CSV:

```text
spr_summary.csv
```

## Frame-wise CSV columns

| Column | Description |
|---|---|
| `file` | WAV filename |
| `frame_index` | Frame number |
| `start_sample` | Start sample of the frame |
| `end_sample` | End sample of the frame |
| `frame_length_samples` | Actual frame length; the final frame may be shorter than 2048 samples |
| `start_time_s` | Frame start time in seconds |
| `end_time_s` | Frame end time in seconds |
| `rms` | Frame RMS after frame-level DC removal and before windowing |
| `spr_db` | Frame-wise SPR value in dB |
| `low_band_value` | Low-band peak power or band power, depending on method |
| `high_band_value` | High-band peak power or band power, depending on method |
| `low_peak_freq_hz` | Frequency of the maximum spectral peak in the 0-2 kHz band |
| `high_peak_freq_hz` | Frequency of the maximum spectral peak in the 2-4 kHz band |
| `status` | `included` or `excluded_low_rms` |

## Summary CSV columns

| Column | Description |
|---|---|
| `file` | WAV filename |
| `path` | Path to the analyzed file |
| `sampling_rate_hz` | Sampling rate |
| `duration_s` | Duration of the file |
| `frame_size_samples` | Frame size in samples |
| `hop_size_samples` | Hop size in samples |
| `n_fft` | FFT size |
| `frame_duration_ms` | Frame duration in milliseconds |
| `method` | `peak` or `bandpower` |
| `low_band_hz` | Low-frequency band |
| `high_band_hz` | High-frequency band |
| `n_frames_total` | Total number of frames |
| `n_frames_included` | Number of frames included in median calculation |
| `median_spr_db` | Representative SPR value for the file |
| `mean_spr_db` | Mean of frame-wise SPR values |
| `sd_spr_db` | Standard deviation of frame-wise SPR values |
| `min_spr_db` | Minimum frame-wise SPR value |
| `max_spr_db` | Maximum frame-wise SPR value |
| `framewise_csv` | Path to the frame-wise CSV file |

## Treatment of the final frame

The final frame is included even if it contains fewer than 2048 samples. For FFT computation, the short frame is zero-padded to 2048 samples.

This prevents the final part of the audio token from being discarded.

## Optional low-RMS exclusion

By default, all frames are included:

```text
--min-rms 0.0
```

If very low-level frames should be excluded, set a threshold manually:

```bash
python calculate_spr_framewise_median.py input.wav --min-rms 0.001
```

For strict reproducibility, report the threshold if this option is used.

## Suggested Methods text

> SPR was calculated frame by frame using a 2048-sample Hann window and a 2048-point FFT. At a sampling rate of 44.1 kHz, each frame corresponded to approximately 46.4 ms. The final frame was retained even when shorter than 2048 samples and was zero-padded to 2048 samples for FFT computation. For each frame, SPR was defined as the dB ratio between the greatest spectral peak in the 2-4 kHz band and the greatest spectral peak in the 0-2 kHz band. The median of frame-wise SPR values was used as the representative SPR value for each audio file.

If using the band-power variant, replace the fourth sentence with:

> For each frame, SPR was defined as 10log10(E2-4kHz/E0-2kHz), where E2-4kHz and E0-2kHz denote summed spectral power in the 2-4 kHz and 0-2 kHz bands, respectively.

## Reproducibility checklist

When reporting results, document the following parameters:

- Sampling rate
- Mono conversion method
- Frame size
- Hop size
- FFT size
- Window function
- Frequency bands
- SPR definition: peak ratio or band-power ratio
- Treatment of the final frame
- Whether low-RMS frames were excluded
- Summary statistic used as the file-level representative value

## Notes

- This implementation uses **2048 samples** as both the frame size and FFT size by default.
- At 44.1 kHz, this corresponds to approximately **46.4 ms**, not 92.9 ms.
- A 4096-sample frame corresponds to approximately **92.9 ms** at 44.1 kHz.
- If a 4096-point FFT is used only for zero-padding a 2048-sample frame, the actual analysis window remains 46.4 ms; only the FFT interpolation grid changes.
