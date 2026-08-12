# Singing Voice Acoustic Metrics in Python

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![Research code](https://img.shields.io/badge/status-research%20code-6f42c1.svg)](#research-use)
[![Metrics](https://img.shields.io/badge/acoustic%20metrics-7-2ea44f.svg)](#included-metrics)

Python implementations of seven acoustic measures used in singing-voice
research and classical singing pedagogy.

## Included metrics

| Metric | Main script | Method or reference basis |
|---|---|---|
| H1-H2 | [`acoustic_metrics/h1_h2/calculate_h1_h2_holmberg1995.py`](acoustic_metrics/h1_h2/calculate_h1_h2_holmberg1995.py) | Holmberg et al. (1995) |
| CPPS | [`acoustic_metrics/cpps/calculate_cpps_baker2024.py`](acoustic_metrics/cpps/calculate_cpps_baker2024.py) | Baker et al.; Praat PowerCepstrogram |
| SFR | [`acoustic_metrics/sfr/calculate_sfr_excel_definition.py`](acoustic_metrics/sfr/calculate_sfr_excel_definition.py) | SFR workbook definition |
| SPR | [`acoustic_metrics/spr/calculate_spr_omori1996.py`](acoustic_metrics/spr/calculate_spr_omori1996.py) | Omori et al. (1996) |
| Q value | [`acoustic_metrics/q_value/calculate_q_lpc_jstage2014.py`](acoustic_metrics/q_value/calculate_q_lpc_jstage2014.py) | LPC peak and -3 dB bandwidth |
| Alpha Ratio | [`acoustic_metrics/alpha_ratio/calculate_alpha_ratio_patel2010.py`](acoustic_metrics/alpha_ratio/calculate_alpha_ratio_patel2010.py) | Patel et al. (2010) |
| Spectral Centroid | [`acoustic_metrics/spectral_centroid/calculate_spectral_centroid_schubert_wolfe2006.py`](acoustic_metrics/spectral_centroid/calculate_spectral_centroid_schubert_wolfe2006.py) | Schubert and Wolfe (2006) |

Each metric folder contains its own `README.md` with the definition, analysis
settings, requirements, and usage instructions.

## Repository structure

```text
Program/
├── Audio Source_E01-E07 _tanto_ before and after/
│   └── Participant audio recordings for E01-E07
│
├── Audio Source_V01-V13_tanto_ before and after/
│   └── Participant audio recordings for V01-V13
│
├── ICC_Table4A/
│   ├── README_ICC.md
│   ├── table4a_ratings.csv
│   ├── icc_table4a.py
│   ├── icc_results.txt
│   └── requirements_icc.txt
│
├── acoustic_metrics/
│   ├── alpha_ratio/
│   │   ├── README.md
│   │   └── calculate_alpha_ratio_patel2010.py
│   │
│   ├── cpps/
│   │   ├── README.md
│   │   └── calculate_cpps_baker2024.py
│   │
│   ├── h1_h2/
│   │   ├── README.md
│   │   └── calculate_h1_h2_holmberg1995.py
│   │
│   ├── q_value/
│   │   ├── README.md
│   │   └── calculate_q_lpc_jstage2014.py
│   │
│   ├── sfr/
│   │   ├── README.md
│   │   └── calculate_sfr_excel_definition.py
│   │
│   ├── spectral_centroid/
│   │   ├── README.md
│   │   └── calculate_spectral_centroid_schubert_wolfe2006.py
│   │
│   └── spr/
│       ├── README.md
│       └── calculate_spr_omori1996.py
│
├── normalized_dtw/
│   ├── README.md
│   ├── calculate_path_normalized_dtw.py
│   └── example_time_series_long.csv
│
├── data/
│   ├── README.md
│   ├── 01_rater_scores_and_terms.csv
│   ├── 02_acoustic_features_Z1_Z2.csv
│   ├── 03_M2025_acoustic_data.csv
│   └── 04_normalized_DTW_results.csv
│
├── .gitignore
├── Audio Source_M2025_tanto_tanto_mono.wav
├── README.md
└── requirements.txt
```

The `acoustic_metrics/` directory contains a dedicated subdirectory for each of the seven acoustic measures, with the corresponding analysis script and measure-specific `README.md`. The `data/` directory contains the de-identified participant-level and rater-level data supporting the revised manuscript. The `ICC_Table4A/` directory contains the complete 20 × 4 rating matrix and a reproducible Python implementation of the inter-rater reliability analysis. Participant audio recordings and the M2025 pedagogical reference recording are also publicly available in this repository.

## Installation

Create a virtual environment and install the tested library versions:

```bash
python -m venv .venv
```

Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

macOS or Linux:

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

## Example commands

CPPS:

```bash
python acoustic_metrics/cpps/calculate_cpps_baker2024.py input.wav --start 0.430 --end 0.930
```

SPR:

```bash
python acoustic_metrics/spr/calculate_spr_omori1996.py input.wav --center-sec 0.680
```

Spectral Centroid:

```bash
python acoustic_metrics/spectral_centroid/calculate_spectral_centroid_schubert_wolfe2006.py input.wav
```

See the README in each metric folder for metric-specific instructions.

## Tested environment

The analysis scripts were tested with the following software and library versions:

- Python 3.13.5
- NumPy 2.3.5
- SciPy 1.17.0
- SoundFile 0.13.1
- Librosa 0.11.0
- Praat-Parselmouth 0.4.7 (embedded Praat 6.1.38)
- Praat 6.4.27 (used for manual audio inspection and mono conversion)

The package versions required for the acoustic-analysis scripts are also specified in `requirements.txt`. The ICC analysis has its own reproducibility environment documented in `ICC_Table4A/requirements_icc.txt`.

## Research use

These scripts provide transparent and reproducible research implementations of the acoustic measures used in the revised manuscript. The implemented definitions, analysis intervals, parameter settings, and relevant methodological choices are documented in the corresponding metric-specific `README.md` files and source-code comments.

The acoustic measures are intended as descriptive measures of specific properties of the recorded signal. They should not be interpreted as direct measurements of physiological states, clinical diagnostic criteria, or independent criteria of pedagogical quality or vocal development. Users should verify analysis intervals, recording conditions, sampling rates, and parameter settings before comparing results across recordings, participants, or studies.

## Data and privacy

Participant audio recordings used in this study are included in this repository. Public sharing of the audio recordings is permitted under the participants' consent and the approved ethical conditions.

The repository also provides de-identified rater-level ratings and evaluative-term annotations, participant-level acoustic-feature data, M2025 reference data, path-normalized DTW results, analysis scripts, package information, and reproducibility documentation.

The de-identified tabular data supporting the revised manuscript are available in the `data/` directory. Reproducible inter-rater reliability data and code are available in the `ICC_Table4A/` directory.

## Citation

When using the data or analysis scripts from this repository, please cite the associated manuscript and this repository:

**Soprano-Science. Program: Singing Voice Acoustic Metrics and Reproducibility Data. GitHub repository.**

https://github.com/Soprano-Science/Program

For individual acoustic measures, please also cite the original methodological publication identified in the corresponding metric-specific `README.md`.

After publication of the associated article, the final journal citation should be used together with the repository citation.
