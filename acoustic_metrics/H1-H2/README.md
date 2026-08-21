# Frame-wise H1-H2 Calculation

This repository contains a Python script for calculating **H1-H2** from WAV files using a frame-wise workflow for singing-voice analysis.

The revised representative H1-H2 value for each WAV file is:

```text
median_h1_minus_h2_db
```

This value is the **median of valid frame-wise raw H1-H2 values** calculated over the **entire WAV file**.

---

## Purpose

Earlier versions of the analysis used a fixed analysis point or a short fixed segment, for example an analysis frame at the vowel midpoint or a segment such as 0.430–0.930 seconds. This revised script does **not** use a midpoint-only segment. Instead, it analyzes the entire WAV file using non-overlapping frames.

For each input file, the script:

1. reads the WAV file;
2. converts stereo audio to mono if necessary;
3. splits the entire signal into **2048-sample frames**;
4. retains the final frame even when it is shorter than 2048 samples;
5. zero-pads short final frames to the FFT size for spectral analysis;
6. estimates F0 in each frame using normalized autocorrelation;
7. estimates H1 and H2 amplitudes from the FFT spectrum near F0 and 2F0;
8. calculates raw H1-H2 for each valid frame;
9. writes frame-wise values to CSV;
10. writes a file-level summary CSV;
11. reports the **median of valid frame-wise raw H1-H2 values** as the representative H1-H2 value for each WAV file.

---

## Acoustic Definition

For each valid frame:

\[
\mathrm{H1-H2} = A(H1) - A(H2)
\]

where:

- **H1** is the amplitude of the first harmonic near F0;
- **H2** is the amplitude of the second harmonic near 2 × F0;
- harmonic amplitudes are expressed in dB;
- H1-H2 is expressed in dB.

A larger H1-H2 value means that H1 is stronger relative to H2. In voice-quality research, H1-H2 is commonly used as one acoustic correlate of phonation type or glottal-source characteristics. However, it is also affected by vowel, F0, formant frequencies, formant bandwidths, recording conditions, and F0-estimation accuracy.

---

## Default Analysis Parameters

| Parameter | Default value |
|---|---:|
| Frame size | 2048 samples |
| Hop size | 2048 samples |
| Overlap | none |
| FFT size | 2048 points |
| Window | Hann window |
| Final frame | retained even if shorter than 2048 samples |
| Final-frame FFT handling | zero-padded to 2048 samples |
| F0 estimation | normalized autocorrelation |
| Default F0 range | 70–1200 Hz |
| Minimum autocorrelation clarity | 0.30 |
| Minimum RMS level | -60 dBFS |
| Representative file-level value | median of valid frame-wise raw H1-H2 values |

At a sampling rate of 44.1 kHz, a 2048-sample frame corresponds to approximately:

\[
2048 / 44100 \approx 0.0464 \text{ s} = 46.4 \text{ ms}
\]

A 51.2 ms frame length corresponds to 2048 samples only when the sampling rate is 40 kHz:

\[
2048 / 40000 = 0.0512 \text{ s} = 51.2 \text{ ms}
\]

---

## Important Note on H1\*-H2\*

The primary output of this script is **raw, uncorrected H1-H2**:

```text
median_h1_minus_h2_db
```

The script can also report an approximate F1-corrected value, **H1\*-H2\***, when the option below is used:

```bash
--with-correction
```

The corrected value is written to:

```text
h1_star_minus_h2_star_db
median_h1_star_minus_h2_star_db
```

However, the corrected value should be treated as **diagnostic only** in this implementation. It should not be treated as a direct VoiceSauce-equivalent value unless the formant-correction procedure has been separately validated for the specific data.

### How to interpret H1\*-H2\* cautiously

The F1-corrected value attempts to compensate for the effect of the first formant on the measured amplitudes of H1 and H2. This is useful in principle because raw harmonic amplitudes are affected by vocal-tract filtering, not only by the glottal source.

However, the correction depends strongly on the estimated **F1 frequency** and **F1 bandwidth**. This is especially important when **H2 is above F1**.

When H2 is above F1, H2 lies on the descending side of the first-formant resonance. In that configuration, even a small error in estimated F1 or bandwidth can cause a large change in the corrected H2 amplitude. Consequently, H1\*-H2\* may change substantially, may exaggerate the apparent source difference, or may even suggest a misleading direction of change.

For this reason:

- use `median_h1_minus_h2_db` as the primary H1-H2 value;
- treat `median_h1_star_minus_h2_star_db` as a supplementary diagnostic value;
- inspect `formant_warning` in the frame-wise CSV;
- be especially cautious with frames marked:

```text
H2_above_F1_correction_sensitive_to_F1_estimate
```

- do not interpret corrected values from warning-marked frames as direct evidence of glottal closure or breathiness;
- do not compare corrected values across singers or vowels unless F1 estimation quality has been checked;
- for high-pitched singing, soprano passaggio/high notes, or short word-level tokens containing consonants, raw H1-H2 is usually more transparent than approximate frame-wise H1\*-H2\*.

In summary, H1\*-H2\* can be useful for checking how much the first formant may be influencing raw H1-H2, but the revised recommended output value for this workflow is still:

```text
median_h1_minus_h2_db
```

---

## Installation

The script requires Python 3 and the following packages:

```bash
pip install numpy scipy
```

---

## Basic Usage

### Single WAV file

```bash
python calculate_h1_h2_framewise_median.py input.wav
```

By default, output files are saved in:

```text
h1_h2_results/
```

### Folder of WAV files

```bash
python calculate_h1_h2_framewise_median.py ./audio
```

### Recursive folder processing

```bash
python calculate_h1_h2_framewise_median.py ./audio --recursive
```

### Optional F1-corrected diagnostic columns

```bash
python calculate_h1_h2_framewise_median.py input.wav --with-correction
```

This adds approximate H1\*, H2\*, and H1\*-H2\* values to the frame-wise and summary CSV files when F1 estimation succeeds. These values are diagnostic and should not replace raw H1-H2 unless justified in the method section.

---

## Optional Arguments

```bash
python calculate_h1_h2_framewise_median.py input.wav \
  --output-dir h1_h2_results \
  --frame-size 2048 \
  --hop-size 2048 \
  --fft-size 2048 \
  --min-f0 70 \
  --max-f0 1200 \
  --min-clarity 0.30
```

| Argument | Meaning | Default |
|---|---|---:|
| `--output-dir` | output folder | `h1_h2_results` |
| `--recursive` | process folders recursively | off |
| `--frame-size` | frame size in samples | 2048 |
| `--hop-size` | hop size in samples | 2048 |
| `--fft-size` | FFT size in samples | 2048 |
| `--min-f0` | minimum F0 for autocorrelation search | 70 Hz |
| `--max-f0` | maximum F0 for autocorrelation search | 1200 Hz |
| `--min-clarity` | minimum normalized autocorrelation peak for valid frames | 0.30 |
| `--min-rms-dbfs` | frames below this RMS level are marked invalid | -60 dBFS |
| `--harmonic-search-half-width-hz` | fixed half-width around expected H1/H2; if omitted, adaptive width is used | adaptive |
| `--with-correction` | also estimate approximate F1-corrected H1\*-H2\* | off |
| `--lpc-order` | LPC order for optional F1 estimation | 14 |
| `--preemphasis` | pre-emphasis coefficient for optional F1 estimation | 0.97 |
| `--min-f1` | minimum acceptable F1 candidate | 150 Hz |
| `--max-f1` | maximum acceptable F1 candidate | 1500 Hz |
| `--max-f1-bandwidth` | maximum acceptable F1 bandwidth | 1000 Hz |
| `--formant-proximity-warning-hz` | warning threshold when F1 is close to H1 or H2 | 100 Hz |

---

## Output Files

For each WAV file, the script creates a frame-wise CSV file:

```text
<filename>_h1_h2_framewise.csv
```

It also creates one summary file:

```text
h1_h2_summary.csv
```

---

## Frame-wise CSV Columns

| Column | Description |
|---|---|
| `frame_index` | frame number |
| `start_sample` | start sample of the frame |
| `end_sample` | end sample of the frame |
| `frame_length_samples` | actual frame length; the final frame may be shorter than 2048 samples |
| `start_time_s` | frame start time in seconds |
| `end_time_s` | frame end time in seconds |
| `rms_dbfs` | RMS amplitude in dBFS |
| `valid_h1h2` | whether H1-H2 could be calculated |
| `invalid_reason` | reason for invalid calculation, if any |
| `f0_hz` | estimated F0 in Hz |
| `autocorr_clarity` | normalized autocorrelation peak |
| `h1_frequency_hz` | spectral peak frequency used for H1 |
| `h2_frequency_hz` | spectral peak frequency used for H2 |
| `h1_db` | H1 amplitude in dB |
| `h2_db` | H2 amplitude in dB |
| `h1_minus_h2_db` | raw H1-H2 in dB |
| `f1_hz` | estimated F1, only when `--with-correction` is used |
| `b1_hz` | estimated F1 bandwidth, only when `--with-correction` is used |
| `h1_star_db` | approximate corrected H1\*, only when available |
| `h2_star_db` | approximate corrected H2\*, only when available |
| `h1_star_minus_h2_star_db` | approximate corrected H1\*-H2\*, only when available |
| `formant_warning` | warning or reason related to F1 correction |

Frames for which F0 or harmonic amplitudes cannot be estimated are retained in the CSV, but their H1-H2 value is left blank. The summary median is calculated from valid frame-wise H1-H2 values only.

---

## Summary CSV Columns

| Column | Description |
|---|---|
| `input_file` | input WAV file |
| `sampling_rate_hz` | sampling rate |
| `duration_s` | duration of the WAV file |
| `frame_size_samples` | frame size |
| `hop_size_samples` | hop size |
| `fft_size_samples` | FFT size |
| `n_frames_total` | total number of frames |
| `n_valid_h1h2_frames` | number of frames with valid raw H1-H2 values |
| `n_valid_corrected_frames` | number of frames with valid approximate H1\*-H2\* values |
| `n_corrected_warning_frames` | number of corrected frames with any formant warning |
| `n_h2_above_f1_warning_frames` | number of corrected frames where H2 was above F1 |
| `median_h1_minus_h2_db` | representative H1-H2 value for the file |
| `mean_h1_minus_h2_db` | mean of valid frame-wise H1-H2 values |
| `sd_h1_minus_h2_db` | standard deviation of valid frame-wise H1-H2 values |
| `median_h1_star_minus_h2_star_db` | diagnostic median H1\*-H2\*, only when correction is requested and valid |
| `mean_h1_star_minus_h2_star_db` | diagnostic mean H1\*-H2\*, only when correction is requested and valid |
| `sd_h1_star_minus_h2_star_db` | diagnostic SD H1\*-H2\*, only when correction is requested and valid |
| `framewise_csv` | path to frame-wise CSV |

---

## Recommended Output Value

For the revised analysis workflow, use this column as the H1-H2 value for each WAV file:

```text
median_h1_minus_h2_db
```

This corresponds to the median of all valid frame-wise raw H1-H2 values in the frame-wise CSV file.

---

## Suggested Methods Description

The following text may be adapted for a manuscript or supplementary material:

> H1-H2 was calculated frame by frame over the entire WAV file using non-overlapping 2048-sample frames. The final frame was retained even when shorter than 2048 samples and was zero-padded to 2048 samples for FFT computation. For each frame, F0 was estimated using normalized autocorrelation, and H1 and H2 amplitudes were estimated as the spectral peak amplitudes near F0 and 2F0, respectively. Raw H1-H2 was defined as the difference between the H1 and H2 amplitudes in dB. The median of valid frame-wise H1-H2 values was used as the representative H1-H2 value for each WAV file.

If the optional diagnostic correction is used, add:

> An approximate F1-corrected H1\*-H2\* value was also computed as a diagnostic measure when F1 estimation succeeded. Because the correction is sensitive to F1 and bandwidth estimation, especially when H2 is above F1, the primary reported value was the uncorrected median H1-H2.

---

## Notes for Singing-Voice Analysis

- H1-H2 requires a reliable F0 estimate.
- Frames containing silence, unvoiced consonants, or very low amplitude may not yield valid values.
- The script retains all frames in the CSV, including invalid frames, so that the analysis process remains transparent.
- The reported representative value is the median of valid frame-wise H1-H2 values, not a midpoint-only measurement and not an LTAS-based single value.
- For high-pitched singing, formant estimation and F1 correction can be unstable. Raw H1-H2 and formant-corrected H1\*-H2\* should not be conflated.
- Whole-word tokens such as `tanto` and `cessa` may include consonantal or transition frames. The median is used to reduce the influence of such local outliers, but invalid-frame counts should still be checked.

---

## Reproducibility Checklist

When reporting results, record the following:

- sampling rate;
- frame size;
- hop size;
- FFT size;
- window type;
- F0 estimation method;
- F0 search range;
- clarity threshold;
- treatment of final short frames;
- whether values are raw H1-H2 or formant-corrected H1\*-H2\*;
- aggregation statistic used for the file-level value.

For this script, the recommended file-level value is:

```text
median_h1_minus_h2_db
```

---

## References

- VoiceSauce documentation explains the convention that corrected harmonic measures are indicated with an asterisk, such as H1\*-H2\*.
- Iseli, Shue, and Alwan (2007) discuss acoustic measures related to the voice source and the need to compensate for the influence of formant frequencies and bandwidths.
- Recent work on H1-H2 emphasizes that both raw and corrected H1-H2 can involve practical issues, including error propagation from F0 and formant estimation.
