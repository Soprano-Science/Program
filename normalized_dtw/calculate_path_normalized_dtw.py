#!/usr/bin/env python3
"""Compute path-normalized DTW distances for SFR and Q time series.

This script implements the method used in the revised manuscript:

* local cost: absolute difference, ``abs(x_i - y_j)``;
* permitted steps: diagonal, vertical, and horizontal;
* no additional global warping window;
* path-normalized distance: cumulative local cost divided by the number
  of aligned point pairs in the optimal path;
* SFR and Q are analysed separately;
* no scaling, standardization, smoothing, interpolation, resampling, or
  temporal normalization is performed inside this script.

The reported longitudinal quantity is an intuitively oriented proximity
change:

    proximity_change = DTWnorm(Z1, M2025) - DTWnorm(Z2, M2025)

Positive values therefore indicate approach toward M2025 at Z2, and
negative values indicate divergence.  This is the negative of the legacy
``Z2 - Z1`` distance-change convention; the underlying Z1 and Z2 distances
are unchanged.

Expected input format
---------------------
A UTF-8 CSV in long format with these columns:

    Subject,Occasion,Frame,SFR,Q

Use ``Subject=M2025`` and ``Occasion=REF`` for the reference series.  Each
participant must have ``Z1`` and ``Z2`` rows.  Frames are sorted numerically
within each Subject/Occasion combination.

Example:

    Subject,Occasion,Frame,SFR,Q
    M2025,REF,0,12.3,30.1
    M2025,REF,1,12.8,31.0
    E01,Z1,0,10.0,24.0
    E01,Z2,0,11.1,27.2

The script requires finite numeric series.  By default, missing or non-finite
values produce an error.  ``--missing-policy drop`` is available only for a
prespecified analysis in which removal of such rows is justified and reported.
"""

from __future__ import annotations

import argparse
import csv
import math
import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class DTWResult:
    """Result of a one-dimensional DTW comparison."""

    cumulative_cost: float
    path_length: int
    normalized_distance: float


@dataclass(frozen=True)
class SeriesKey:
    subject: str
    occasion: str


def _as_finite_1d(values: Sequence[float], *, label: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"{label} must be one-dimensional, got shape {array.shape}.")
    if array.size == 0:
        raise ValueError(f"{label} is empty.")
    if not np.all(np.isfinite(array)):
        bad = np.flatnonzero(~np.isfinite(array)).tolist()
        raise ValueError(f"{label} contains non-finite values at indices {bad[:10]}.")
    return array


def path_normalized_dtw(
    reference: Sequence[float],
    target: Sequence[float],
) -> DTWResult:
    """Return cumulative cost, optimal path length, and normalized distance.

    Dynamic programming uses absolute local cost and the standard three-step
    pattern.  Ties are resolved deterministically in the order diagonal,
    vertical, horizontal.  With continuous acoustic values exact ties are
    uncommon, but the explicit rule makes the implementation reproducible.
    """

    x = _as_finite_1d(reference, label="reference")
    y = _as_finite_1d(target, label="target")
    n, m = x.size, y.size

    costs = np.full((n + 1, m + 1), np.inf, dtype=np.float64)
    predecessor = np.full((n + 1, m + 1), -1, dtype=np.int8)
    costs[0, 0] = 0.0

    # Predecessor codes: 0=diagonal, 1=vertical, 2=horizontal.
    for i in range(1, n + 1):
        xi = x[i - 1]
        for j in range(1, m + 1):
            candidates = (
                costs[i - 1, j - 1],
                costs[i - 1, j],
                costs[i, j - 1],
            )
            step = int(np.argmin(candidates))
            costs[i, j] = abs(xi - y[j - 1]) + candidates[step]
            predecessor[i, j] = step

    i, j = n, m
    path_length = 0
    while i > 0 or j > 0:
        if i == 0:
            j -= 1
        elif j == 0:
            i -= 1
        else:
            step = int(predecessor[i, j])
            if step == 0:
                i -= 1
                j -= 1
            elif step == 1:
                i -= 1
            elif step == 2:
                j -= 1
            else:
                raise RuntimeError(f"Invalid predecessor at ({i}, {j}): {step}")
        path_length += 1

    cumulative_cost = float(costs[n, m])
    if path_length <= 0:
        raise RuntimeError("The optimal path has zero length.")
    normalized_distance = cumulative_cost / path_length
    return DTWResult(cumulative_cost, path_length, normalized_distance)


def proximity_change(z1_distance: float, z2_distance: float) -> float:
    """Return Z1 minus Z2 so positive values mean approach toward M2025."""

    return float(z1_distance - z2_distance)


def direction(value: float, *, tolerance: float = 0.0) -> str:
    """Classify a proximity change as approach, divergence, or unchanged."""

    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    if value > tolerance:
        return "approach"
    if value < -tolerance:
        return "divergence"
    return "unchanged"


def _parse_float(
    raw: str,
    *,
    label: str,
    missing_policy: str,
) -> float | None:
    text = raw.strip()
    if text == "":
        if missing_policy == "drop":
            return None
        raise ValueError(f"Missing value in {label}.")
    try:
        value = float(text)
    except ValueError as exc:
        raise ValueError(f"Non-numeric value {text!r} in {label}.") from exc
    if not math.isfinite(value):
        if missing_policy == "drop":
            return None
        raise ValueError(f"Non-finite value {text!r} in {label}.")
    return value


def read_long_csv(
    path: Path,
    *,
    subject_column: str,
    occasion_column: str,
    frame_column: str,
    feature_columns: Sequence[str],
    missing_policy: str,
) -> dict[SeriesKey, dict[str, list[float]]]:
    """Read and group a long-format time-series CSV."""

    if missing_policy not in {"error", "drop"}:
        raise ValueError("missing_policy must be 'error' or 'drop'")

    grouped_rows: dict[SeriesKey, list[tuple[float, dict[str, float | None]]]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header row.")
        required = {subject_column, occasion_column, frame_column, *feature_columns}
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")

        for line_number, row in enumerate(reader, start=2):
            subject = row[subject_column].strip()
            occasion = row[occasion_column].strip()
            if not subject or not occasion:
                raise ValueError(
                    f"Blank Subject/Occasion at line {line_number} in {path}."
                )
            frame = _parse_float(
                row[frame_column],
                label=f"{frame_column}, line {line_number}",
                missing_policy="error",
            )
            assert frame is not None
            features: dict[str, float | None] = {}
            for feature in feature_columns:
                features[feature] = _parse_float(
                    row[feature],
                    label=f"{feature}, line {line_number}",
                    missing_policy=missing_policy,
                )
            grouped_rows.setdefault(SeriesKey(subject, occasion), []).append(
                (frame, features)
            )

    output: dict[SeriesKey, dict[str, list[float]]] = {}
    for key, rows in grouped_rows.items():
        rows.sort(key=lambda item: item[0])
        output[key] = {feature: [] for feature in feature_columns}
        for _, features in rows:
            for feature, value in features.items():
                if value is not None:
                    output[key][feature].append(value)
        for feature in feature_columns:
            if not output[key][feature]:
                raise ValueError(f"No usable {feature} values for {key}.")
    return output


def run_study(args: argparse.Namespace) -> None:
    feature_columns = [item.strip() for item in args.features.split(",") if item.strip()]
    if not feature_columns:
        raise ValueError("At least one feature column is required.")

    series = read_long_csv(
        args.input_csv,
        subject_column=args.subject_column,
        occasion_column=args.occasion_column,
        frame_column=args.frame_column,
        feature_columns=feature_columns,
        missing_policy=args.missing_policy,
    )

    reference_key = SeriesKey(args.reference_subject, args.reference_occasion)
    if reference_key not in series:
        raise ValueError(f"Reference series {reference_key} was not found.")

    subjects = sorted(
        {
            key.subject
            for key in series
            if key.subject != args.reference_subject
            and key.occasion in {args.z1_label, args.z2_label}
        }
    )
    if not subjects:
        raise ValueError("No participant Z1/Z2 series were found.")

    fieldnames = ["Subject"]
    for feature in feature_columns:
        fieldnames.extend(
            [
                f"{feature}_Z1_cumulative_cost",
                f"{feature}_Z1_path_length",
                f"{feature}_Z1_DTWnorm",
                f"{feature}_Z2_cumulative_cost",
                f"{feature}_Z2_path_length",
                f"{feature}_Z2_DTWnorm",
                f"{feature}_ProximityChange_Z1_minus_Z2",
                f"{feature}_Direction",
            ]
        )

    rows: list[dict[str, str | int | float]] = []
    for subject in subjects:
        z1_key = SeriesKey(subject, args.z1_label)
        z2_key = SeriesKey(subject, args.z2_label)
        if z1_key not in series or z2_key not in series:
            raise ValueError(f"{subject} does not have both Z1 and Z2 series.")

        out: dict[str, str | int | float] = {"Subject": subject}
        for feature in feature_columns:
            reference = series[reference_key][feature]
            z1 = series[z1_key][feature]
            z2 = series[z2_key][feature]
            r1 = path_normalized_dtw(reference, z1)
            r2 = path_normalized_dtw(reference, z2)
            change = proximity_change(r1.normalized_distance, r2.normalized_distance)
            out.update(
                {
                    f"{feature}_Z1_cumulative_cost": r1.cumulative_cost,
                    f"{feature}_Z1_path_length": r1.path_length,
                    f"{feature}_Z1_DTWnorm": r1.normalized_distance,
                    f"{feature}_Z2_cumulative_cost": r2.cumulative_cost,
                    f"{feature}_Z2_path_length": r2.path_length,
                    f"{feature}_Z2_DTWnorm": r2.normalized_distance,
                    f"{feature}_ProximityChange_Z1_minus_Z2": change,
                    f"{feature}_Direction": direction(
                        change, tolerance=args.direction_tolerance
                    ),
                }
            )
        rows.append(out)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} participant rows to: {args.output_csv}")
    print("Sign convention: positive = approach; negative = divergence.")


def run_pair(args: argparse.Namespace) -> None:
    def read_column(path: Path, column: str) -> list[float]:
        values: list[float] = []
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(f"Column {column!r} was not found in {path}.")
            for line_number, row in enumerate(reader, start=2):
                value = _parse_float(
                    row[column],
                    label=f"{column}, line {line_number}",
                    missing_policy=args.missing_policy,
                )
                if value is not None:
                    values.append(value)
        if not values:
            raise ValueError(f"No usable values in {path}:{column}.")
        return values

    reference = read_column(args.reference_csv, args.reference_column)
    target = read_column(args.target_csv, args.target_column)
    result = path_normalized_dtw(reference, target)
    print(f"Cumulative cost:    {result.cumulative_cost:.12f}")
    print(f"Path length:        {result.path_length}")
    print(f"Normalized distance:{result.normalized_distance: .12f}")


def run_self_test() -> None:
    reference = [0.0, 1.0, 2.0]
    z1 = [0.0, 0.0, 0.0]
    z2 = [0.0, 1.0, 2.0]
    r1 = path_normalized_dtw(reference, z1)
    r2 = path_normalized_dtw(reference, z2)
    change = proximity_change(r1.normalized_distance, r2.normalized_distance)
    assert r2.normalized_distance == 0.0
    assert change > 0.0
    assert direction(change) == "approach"

    diverged = proximity_change(0.0, 1.0)
    assert diverged < 0.0
    assert direction(diverged) == "divergence"
    print("Self-test passed.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Path-normalized one-dimensional DTW with positive proximity "
            "change indicating approach toward M2025."
        )
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"Python {platform.python_version()}, NumPy {np.__version__}",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    study = subparsers.add_parser(
        "study", help="Process all participant Z1/Z2 series from a long-format CSV."
    )
    study.add_argument("input_csv", type=Path)
    study.add_argument("output_csv", type=Path)
    study.add_argument("--features", default="SFR,Q")
    study.add_argument("--subject-column", default="Subject")
    study.add_argument("--occasion-column", default="Occasion")
    study.add_argument("--frame-column", default="Frame")
    study.add_argument("--reference-subject", default="M2025")
    study.add_argument("--reference-occasion", default="REF")
    study.add_argument("--z1-label", default="Z1")
    study.add_argument("--z2-label", default="Z2")
    study.add_argument(
        "--missing-policy",
        choices=("error", "drop"),
        default="error",
        help="Default 'error' is recommended for reproducibility.",
    )
    study.add_argument(
        "--direction-tolerance",
        type=float,
        default=0.0,
        help="Optional non-negative equivalence band around zero.",
    )
    study.set_defaults(func=run_study)

    pair = subparsers.add_parser(
        "pair", help="Compare one reference column with one target column."
    )
    pair.add_argument("reference_csv", type=Path)
    pair.add_argument("target_csv", type=Path)
    pair.add_argument("--reference-column", required=True)
    pair.add_argument("--target-column", required=True)
    pair.add_argument(
        "--missing-policy", choices=("error", "drop"), default="error"
    )
    pair.set_defaults(func=run_pair)

    test = subparsers.add_parser("self-test", help="Run a built-in synthetic test.")
    test.set_defaults(func=lambda _args: run_self_test())
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if getattr(args, "direction_tolerance", 0.0) < 0:
            raise ValueError("--direction-tolerance must be non-negative.")
        args.func(args)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
