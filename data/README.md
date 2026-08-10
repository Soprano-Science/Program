# Data

This directory contains de-identified data files supporting the analyses reported in the revised manuscript:

**Toward a Cognitive Framework of Expert Singing-Voice Evaluation in Bel Canto Pedagogy: A Longitudinal Exploratory Study**

Participant IDs are anonymized and correspond to the IDs used in the manuscript and Supplementary Material.

## Files

### 01_rater_scores_and_terms.csv

Rater-level auditory-perceptual evaluation data for the 20 participants (E01–E07 and V01–V13).

Each participant was evaluated independently by four expert evaluators (T1–T4).

The file contains:

- `Subject`: anonymized participant ID
- `Rater`: evaluator code (T1–T4)
- `DevelopmentRating`: longitudinal development rating  
  - 0 = no development perceived  
  - 1 = some development perceived  
  - 2 = substantial development perceived
- `ParticipantMedian`: median of the four development ratings for that participant
- `SelectedTerms`: evaluative-term codes selected by the evaluator
- `TermCount`: number of evaluative terms selected
- `S01–S10`, `R01–R12`, `P01–P05`: binary indicators for the 27 evaluative terms  
  - 1 = selected  
  - 0 = not selected

The complete Italian evaluative vocabulary and English contextual glosses are reported in Supplementary Table S1 of the revised manuscript.

This file provides the source data for the rater-level rating and evaluative-term analyses, including Supplementary Table S2.

---

### 02_acoustic_features_Z1_Z2.csv

Participant-level values for the five supplementary acoustic measures at:

- `Z1`: university entry
- `Z2`: two years later

The file contains Z1, Z2, and change values for:

- CPPS (dB)
- Alpha ratio (dB)
- H1-H2 (dB)
- Spectral centroid (Hz)
- SPR (dB)

For all acoustic measures:

`Delta = Z2 - Z1`

The direction of a numerical change should not by itself be interpreted as improvement or deterioration. These measures are used descriptively to characterize participant-specific acoustic changes and are interpreted together with expert ratings, evaluative terms, SFR/Q trajectories, and pedagogical context.

---

### 03_M2025_acoustic_data.csv

Acoustic data for M2025, the panel-selected pedagogical reference recording.

M2025 is a recording of one professional female classical singer and is used as a pedagogical reference rather than as an independently validated acoustic gold standard.

The file contains:

- frame-by-frame SFR values
- frame-by-frame Q values
- frame index
- summary mean SFR
- summary mean Q
- CPPS
- Alpha ratio
- H1-H2
- Spectral centroid
- SPR

Rows labeled `timeseries` contain the SFR/Q time-series data.

The row labeled `summary` contains the aggregate acoustic values for M2025.

The M2025 SFR/Q series contains 138 frames.

---

### 04_normalized_DTW_results.csv

Participant-level Dynamic Time Warping (DTW) results comparing each participant's Z1 and Z2 SFR and Q time series with the M2025 reference.

For both SFR and Q, the file contains:

- Z1 cumulative DTW cost
- Z1 optimal path length
- Z1 path-normalized DTW distance
- Z2 cumulative DTW cost
- Z2 optimal path length
- Z2 path-normalized DTW distance
- Z1-to-Z2 change in normalized DTW distance
- number of frames in M2025, Z1, and Z2
- participant-level median expert rating
- descriptive convergence/divergence pattern

Path-normalized DTW distance was calculated as:

`normalized DTW distance = cumulative absolute local cost / optimal path length`

For the change scores:

`Delta_DTWnorm = Z2_DTWnorm - Z1_DTWnorm`

Therefore:

- negative Delta_DTWnorm = Z2 is acoustically closer to M2025 than Z1
- positive Delta_DTWnorm = Z2 is acoustically farther from M2025 than Z1

SFR and Q were analyzed separately and were not combined into a multivariate DTW score.

The original audio signals were not stretched, compressed, or resampled to equal duration. DTW was used only as an analytical alignment procedure between the original time series.

---

## Inter-rater reliability

The input matrix and reproducible Python analysis for the inter-rater reliability calculation are provided separately in:

`../ICC_Table4A/`

That directory contains the 20-participant × 4-rater rating matrix, Python script, package requirements, and calculation output for ICC(2,1) and ICC(2,k).

---

## Data interpretation

These files are provided to support transparency and reproducibility of the revised exploratory analysis.

The acoustic variables are treated as descriptive measures of limited aspects of the recorded signal. They are not interpreted as direct measurements of physiological states or as independent criteria of pedagogical quality.

Likewise, the evaluative terms are treated as pedagogical and perceptual annotations accompanying expert judgments rather than as independently validated physiological labels.
