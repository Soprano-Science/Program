# Path-Normalized Dynamic Time Warping (DTW)

This directory contains a reproducible Python implementation of the path-normalized Dynamic Time Warping (DTW) procedure used in the revised manuscript.

The analysis compares participant SFR and Q time series at Z1 and Z2 separately with the M2025 pedagogical reference.

## Files

* `calculate_path_normalized_dtw.py` — Python implementation of path-normalized one-dimensional DTW
* `example_time_series_long.csv` — example long-format input file

The participant-level DTW results used in the revised manuscript are provided separately in:

`../data/04_normalized_DTW_results.csv`

## DTW definition

For two one-dimensional time series, the local cost is defined as:

`local cost = abs(x_i - y_j)`

The permitted DTW steps are:

* diagonal
* vertical
* horizontal

No additional global warping window is imposed.

The cumulative DTW cost is the sum of the local absolute differences along the optimal alignment path.

The path-normalized DTW distance is defined as:

`DTWnorm = cumulative local cost / number of aligned point pairs in the optimal path`

A smaller normalized DTW distance therefore indicates greater acoustic proximity between the participant time series and the M2025 reference for the acoustic measure being analyzed.

## SFR and Q are analyzed separately

SFR and Q are treated as separate one-dimensional time series.

They are not combined into a multivariate DTW score because they describe different singer's-formant-related acoustic properties and have different numerical scales.

Accordingly, SFR-DTW and Q-DTW distances should be interpreted separately.

## No alteration of the original performance duration

This DTW procedure is used only as an analytical alignment between acoustic time series.

The script does **not**:

* stretch or compress the audio signal
* change the duration of the sung word
* resample recordings to an equal number of frames
* interpolate the time series
* perform temporal normalization
* scale or standardize the feature values
* apply smoothing inside the DTW calculation

Thus, the original temporal duration of each sung performance is preserved.

Any acoustic preprocessing used to generate the SFR and Q input series is performed before the DTW calculation and should be documented separately.

## Longitudinal sign convention

The revised analysis reports the longitudinal change as:

`ProximityChange = DTWnorm(Z1, M2025) - DTWnorm(Z2, M2025)`

Under this convention:

* **positive value** = Z2 is closer to M2025 than Z1 (**approach**)
* **negative value** = Z2 is farther from M2025 than Z1 (**divergence**)
* **zero** = no change in normalized distance

This sign convention was adopted to make the direction of acoustic proximity intuitive: positive values indicate movement toward the pedagogical reference.

Earlier working tables used the opposite distance-change convention:

`Z2 - Z1`

The two conventions contain the same underlying Z1 and Z2 DTW distances; only the sign used to summarize longitudinal change is reversed.

## Expected input format

The study-mode input is a UTF-8 CSV in long format with the following columns:

`Subject,Occasion,Frame,SFR,Q`

For the M2025 reference, use:

* `Subject = M2025`
* `Occasion = REF`

For each participant, provide both:

* `Occasion = Z1`
* `Occasion = Z2`

Example:

```text
Subject,Occasion,Frame,SFR,Q
M2025,REF,0,0,0
M2025,REF,1,1,1
M2025,REF,2,2,2
E01,Z1,0,0,0
E01,Z1,1,0,0
E01,Z1,2,0,0
E01,Z2,0,0,0
E01,Z2,1,1,1
E01,Z2,2,2,2
```

Frames are sorted numerically within each Subject/Occasion combination.

## Missing values

By default, missing or non-finite values produce an error.

This default behavior is intended to prevent silent changes to the analyzed time series.

The optional setting:

`--missing-policy drop`

is available only when removal of missing/non-finite rows is prespecified, justified, and reported.

## Reproducibility

The implementation uses deterministic tie handling.

If multiple predecessor costs are equal, the priority order is:

1. diagonal
2. vertical
3. horizontal

This makes the selected optimal path reproducible when exact ties occur.

## Usage

### Run the built-in self-test

```bash
python calculate_path_normalized_dtw.py self-test
```

A successful test prints:

```text
Self-test passed.
```

### Run the example study input

```bash
python calculate_path_normalized_dtw.py study example_time_series_long.csv example_dtw_results.csv
```

By default, the script analyzes the `SFR` and `Q` columns.

### Specify feature columns explicitly

```bash
python calculate_path_normalized_dtw.py study example_time_series_long.csv example_dtw_results.csv --features SFR,Q
```

## Output

For each participant and each acoustic feature, the output includes:

* Z1 cumulative DTW cost
* Z1 optimal path length
* Z1 path-normalized DTW distance
* Z2 cumulative DTW cost
* Z2 optimal path length
* Z2 path-normalized DTW distance
* longitudinal proximity change (`Z1 - Z2`)
* direction (`approach`, `divergence`, or `unchanged`)

## Interpretation

Path-normalized DTW quantifies acoustic proximity to the M2025 pedagogical reference while allowing analytical correspondence between time-series points of recordings with different durations.

It should not be interpreted as evidence that the original audio was physically time-warped or that longer or shorter performances were normalized to the same duration.

The measure reflects both local value differences and the correspondence selected along the optimal DTW path. It is therefore not a pure measure of contour shape alone.

SFR-DTW and Q-DTW are used descriptively as partial acoustic descriptors and are interpreted together with expert ratings, evaluative-term annotations, supplementary acoustic measures, and pedagogical context.
