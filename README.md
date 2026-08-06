# Singing Voice Acoustic Metrics in Python

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![Research code](https://img.shields.io/badge/status-research%20code-6f42c1.svg)](#research-use)
[![Metrics](https://img.shields.io/badge/acoustic%20metrics-7-2ea44f.svg)](#included-metrics)

Python implementations of seven acoustic measures used in singing-voice
research and classical singing pedagogy.

## Included metrics

| Metric | Main script | Method or reference basis |
|---|---|---|
| H1-H2 | [`h1_h2/calculate_h1_h2_holmberg1995.py`](h1_h2/calculate_h1_h2_holmberg1995.py) | Holmberg et al. (1995) |
| CPPS | [`cpps/calculate_cpps_baker2024.py`](cpps/calculate_cpps_baker2024.py) | Baker et al.; Praat PowerCepstrogram |
| SFR | [`sfr/calculate_sfr_excel_definition.py`](sfr/calculate_sfr_excel_definition.py) | SFR workbook definition |
| SPR | [`spr/calculate_spr_omori1996.py`](spr/calculate_spr_omori1996.py) | Omori et al. (1996) |
| Q value | [`q_value/calculate_q_lpc_jstage2014.py`](q_value/calculate_q_lpc_jstage2014.py) | LPC peak and -3 dB bandwidth |
| Alpha Ratio | [`alpha_ratio/calculate_alpha_ratio_patel2010.py`](alpha_ratio/calculate_alpha_ratio_patel2010.py) | Patel et al. (2010) |
| Spectral Centroid | [`spectral_centroid/calculate_spectral_centroid_schubert_wolfe2006.py`](spectral_centroid/calculate_spectral_centroid_schubert_wolfe2006.py) | Schubert and Wolfe (2006) |

Each metric folder contains its own `README.md` with the definition, analysis
settings, requirements, and usage instructions.

## Repository structure

```text
Program/
|-- README.md
|-- requirements.txt
|-- .gitignore
|-- h1_h2/
|   |-- README.md
|   `-- calculate_h1_h2_holmberg1995.py
|-- cpps/
|   |-- README.md
|   `-- calculate_cpps_baker2024.py
|-- sfr/
|   |-- README.md
|   `-- calculate_sfr_excel_definition.py
|-- spr/
|   |-- README.md
|   `-- calculate_spr_omori1996.py
|-- q_value/
|   |-- README.md
|   `-- calculate_q_lpc_jstage2014.py
|-- alpha_ratio/
|   |-- README.md
|   `-- calculate_alpha_ratio_patel2010.py
`-- spectral_centroid/
    |-- README.md
    `-- calculate_spectral_centroid_schubert_wolfe2006.py
```

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
python cpps/calculate_cpps_baker2024.py input.wav --start 0.430 --end 0.930
```

SPR:

```bash
python spr/calculate_spr_omori1996.py input.wav --center-sec 0.680
```

Spectral Centroid:

```bash
python spectral_centroid/calculate_spectral_centroid_schubert_wolfe2006.py input.wav
```

See the README in each metric folder for metric-specific instructions.

## Tested environment

- Python 3.13.5
- NumPy 2.3.5
- SciPy 1.17.0
- SoundFile 0.13.1
- Librosa 0.11.0
- Praat-Parselmouth 0.4.7
- Embedded Praat 6.1.38

## Research use

These scripts are transparent research implementations of methods described in
the cited literature or analysis definitions. Some source publications do not
report every internal software setting. Such implementation choices are
documented in the relevant metric README and source-code comments.

The scripts are not clinical diagnostic software. Users should verify analysis
intervals, recording conditions, sampling rates, and parameter settings before
comparing results across recordings or studies.

## Data and privacy

No participant audio recordings are included in this repository. Avoid
uploading identifiable or restricted research audio unless the relevant
consent, ethics approval, and data-sharing permissions explicitly allow it.

## Citation

When using a script, cite both this repository and the original methodological
publication identified in the corresponding metric README.
