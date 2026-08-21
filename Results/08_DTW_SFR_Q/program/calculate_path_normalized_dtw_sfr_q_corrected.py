#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Path-normalized DTW analysis for frame-wise SFR and Q-value series.

This script is intended for the revised frame-wise acoustic-metric outputs in
which each WAV file has one CSV per metric. It reads the SFR and Q-value
frame-wise CSV files, keeps only valid frames that were included in the metric
median calculation, and computes path-normalized one-dimensional DTW distances
between the pedagogical reference (M2025) and each participant's Z1 and Z2
series.

Main outputs
------------
1. SFR_Q_DTW_distance_delta_summary.csv
   DTWnorm(M2025, Z1), DTWnorm(M2025, Z2), and
   Delta_Z2_minus_Z1 = DTWnorm(Z2, M2025) - DTWnorm(Z1, M2025).
   Negative delta means Z2 is closer to M2025 than Z1 is; positive delta means greater distance.

2. SFR_Q_valid_frame_series_normalized_position.csv
   Valid frame-wise SFR and Q-value series for M2025, Z1, and Z2.  The x-axis
   variable is normalized_frame_position (0=start, 1=end), not clock time.

3. SFR_Q_DTW_path_aligned_long.csv
   DTW optimal-path-aligned data for each Subject x Metric x Occasion comparison
   against M2025.  The x-axis variable is normalized_path_position.

4. PNG plots
   For each participant and each metric, a separate plot shows M2025, Z1, and
   Z2 against normalized frame position.  SFR and Q-value plots are separate.

DTW definition
--------------
- Local cost: absolute difference abs(x_i - y_j)
- Permitted steps: diagonal, vertical, horizontal
- No global warping window
- Path-normalized distance: cumulative local cost / optimal path length
- No scaling, standardization, smoothing, interpolation, or resampling is applied
  before DTW.

Expected input directory
------------------------
The input root should contain metric subdirectories named 'sfr' and 'q_value',
for example:

    acoustic_metric_results/
    ├── sfr/
    │   ├── E01before_..._sfr_framewise.csv
    │   ├── E01after_..._sfr_framewise.csv
    │   └── M2025_..._sfr_framewise.csv
    └── q_value/
        ├── E01before_..._Y_Q_value_framewise.csv
        ├── E01after_..._Y_Q_value_framewise.csv
        └── M2025_..._Y_Q_value_framewise.csv

The frame-wise CSVs must include:
    frame_index, start_time_s, end_time_s, include_in_median,
    and either sfr_percent or q_value.

"""

from __future__ import annotations

import argparse
import csv
import math
import platform
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:  # plotting can be disabled from command line
    plt = None


@dataclass(frozen=True)
class SeriesPoint:
    subject: str
    occasion: str
    metric: str
    source_file: str
    original_frame_index: int
    valid_frame_index: int
    n_valid_frames: int
    normalized_frame_position: float
    start_time_s: float | None
    end_time_s: float | None
    value: float


@dataclass(frozen=True)
class DTWResult:
    cumulative_cost: float
    path_length: int
    normalized_distance: float
    path: tuple[tuple[int, int], ...]


METRIC_CONFIG = {
    "SFR": {
        "directory": "sfr",
        "glob": "*_sfr_framewise.csv",
        "value_column": "sfr_percent",
        "y_label": "SFR (%)",
        "summary_prefix": "SFR",
        "plot_subdir": "sfr",
    },
    "Q_value": {
        "directory": "q_value",
        "glob": "*_Y_Q_value_framewise.csv",
        "value_column": "q_value",
        "y_label": "Q value",
        "summary_prefix": "Q_value",
        "plot_subdir": "q_value",
    },
}

SUBJECT_RE = re.compile(r"^(E\d{2}|V\d{2})", re.IGNORECASE)


def parse_subject_occasion(filename: str) -> tuple[str, str]:
    """Infer Subject and Occasion from a frame-wise CSV filename."""
    stem = Path(filename).stem
    lower = stem.lower()
    if lower.startswith("m2025"):
        return "M2025", "REF"

    match = SUBJECT_RE.match(stem)
    if not match:
        raise ValueError(f"Cannot infer subject from filename: {filename}")
    subject = match.group(1).upper()

    if "before" in lower:
        occasion = "Z1"
    elif "after" in lower:
        occasion = "Z2"
    else:
        raise ValueError(
            f"Cannot infer occasion from filename: {filename}. "
            "Expected 'before' or 'after', or M2025 reference."
        )
    return subject, occasion


def parse_bool(raw: str) -> bool:
    return raw.strip().lower() in {"true", "1", "yes", "y"}


def parse_optional_float(raw: str | None) -> float | None:
    if raw is None:
        return None
    text = raw.strip()
    if text == "":
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def parse_optional_int(raw: str | None) -> int | None:
    value = parse_optional_float(raw)
    if value is None:
        return None
    return int(round(value))


def discover_metric_files(input_root: Path, metric: str) -> list[Path]:
    cfg = METRIC_CONFIG[metric]
    metric_dir = input_root / cfg["directory"]
    if not metric_dir.exists():
        raise FileNotFoundError(f"Metric directory not found: {metric_dir}")
    files = sorted(metric_dir.glob(cfg["glob"]))
    if not files:
        raise FileNotFoundError(f"No {metric} frame-wise CSV files found in {metric_dir}")
    return files


def read_metric_series(input_root: Path, metric: str) -> dict[tuple[str, str], list[SeriesPoint]]:
    """Read valid frame-wise series for a metric.

    Only rows with include_in_median=True and finite value are retained.  Short
    frames, silent frames, and frames with invalid metric values have already
    been marked include_in_median=False by the metric-calculation scripts.
    """
    cfg = METRIC_CONFIG[metric]
    value_column = cfg["value_column"]
    output: dict[tuple[str, str], list[SeriesPoint]] = {}

    for path in discover_metric_files(input_root, metric):
        if path.name.endswith("summary.csv") or "summary" in path.name:
            continue
        subject, occasion = parse_subject_occasion(path.name)
        rows: list[tuple[int, float | None, float | None, float]] = []
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError(f"{path} has no header row.")
            required = {"frame_index", "include_in_median", value_column}
            missing = sorted(required.difference(reader.fieldnames))
            if missing:
                raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")
            for line_number, row in enumerate(reader, start=2):
                if not parse_bool(row.get("include_in_median", "")):
                    continue
                value = parse_optional_float(row.get(value_column))
                if value is None:
                    continue
                original_frame_index = parse_optional_int(row.get("frame_index"))
                if original_frame_index is None:
                    raise ValueError(f"Missing frame_index at {path}:{line_number}")
                start_time_s = parse_optional_float(row.get("start_time_s"))
                end_time_s = parse_optional_float(row.get("end_time_s"))
                rows.append((original_frame_index, start_time_s, end_time_s, value))

        rows.sort(key=lambda item: item[0])
        n_valid = len(rows)
        points: list[SeriesPoint] = []
        for valid_index, (original_frame_index, start_time_s, end_time_s, value) in enumerate(rows):
            if n_valid > 1:
                normalized = valid_index / (n_valid - 1)
            else:
                normalized = 0.0
            points.append(
                SeriesPoint(
                    subject=subject,
                    occasion=occasion,
                    metric=metric,
                    source_file=path.name,
                    original_frame_index=original_frame_index,
                    valid_frame_index=valid_index,
                    n_valid_frames=n_valid,
                    normalized_frame_position=float(normalized),
                    start_time_s=start_time_s,
                    end_time_s=end_time_s,
                    value=float(value),
                )
            )
        output[(subject, occasion)] = points

    return output


def finite_values(points: Sequence[SeriesPoint], *, label: str) -> np.ndarray:
    values = np.asarray([p.value for p in points], dtype=np.float64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError(f"{label} is empty.")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{label} contains non-finite values.")
    return values


def path_normalized_dtw(reference: Sequence[float], target: Sequence[float]) -> DTWResult:
    """Return cumulative cost, path length, normalized distance, and optimal path.

    Dynamic programming uses absolute local cost and the standard three-step
    pattern.  Ties are resolved deterministically in the order diagonal,
    vertical, horizontal for reproducibility.
    """
    x = np.asarray(reference, dtype=np.float64)
    y = np.asarray(target, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("DTW input series must be one-dimensional.")
    if x.size == 0 or y.size == 0:
        raise ValueError("DTW input series must not be empty.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("DTW input series must contain only finite values.")

    n, m = x.size, y.size
    costs = np.full((n + 1, m + 1), np.inf, dtype=np.float64)
    predecessor = np.full((n + 1, m + 1), -1, dtype=np.int8)
    costs[0, 0] = 0.0

    # Predecessor codes: 0=diagonal, 1=vertical, 2=horizontal.
    for i in range(1, n + 1):
        xi = x[i - 1]
        for j in range(1, m + 1):
            candidates = (costs[i - 1, j - 1], costs[i - 1, j], costs[i, j - 1])
            step = int(np.argmin(candidates))
            costs[i, j] = abs(xi - y[j - 1]) + candidates[step]
            predecessor[i, j] = step

    i, j = n, m
    reversed_path: list[tuple[int, int]] = []
    while i > 0 or j > 0:
        # i and j are one-based positions in the DP matrix.  Clamp index pairs
        # in the border cases; these cases should not normally occur with the
        # initialized DTW matrix, but this keeps the backtracking safe.
        if i == 0:
            j -= 1
        elif j == 0:
            i -= 1
        else:
            reversed_path.append((i - 1, j - 1))
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

    path = tuple(reversed(reversed_path))
    path_length = len(path)
    if path_length <= 0:
        raise RuntimeError("The optimal path has zero length.")
    cumulative_cost = float(costs[n, m])
    normalized_distance = cumulative_cost / path_length
    return DTWResult(cumulative_cost, path_length, normalized_distance, path)


def direction(delta: float | None, *, tolerance: float = 0.0) -> str:
    """Interpret Delta = DTWnorm(Z2) - DTWnorm(Z1).

    Negative delta means the distance to M2025 decreased (approach).
    Positive delta means the distance increased (divergence).
    """
    if delta is None or not math.isfinite(delta):
        return "not_computed"
    if delta < -tolerance:
        return "approach"
    if delta > tolerance:
        return "divergence"
    return "unchanged"


def write_csv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def format_float(value: float | None, digits: int = 12) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def build_valid_series_rows(all_series: dict[str, dict[tuple[str, str], list[SeriesPoint]]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for metric, series in all_series.items():
        for (subject, occasion), points in sorted(series.items()):
            for p in points:
                rows.append(
                    {
                        "Metric": metric,
                        "Subject": subject,
                        "Occasion": occasion,
                        "SourceFile": p.source_file,
                        "OriginalFrameIndex": str(p.original_frame_index),
                        "ValidFrameIndex": str(p.valid_frame_index),
                        "NValidFrames": str(p.n_valid_frames),
                        "NormalizedFramePosition": format_float(p.normalized_frame_position, 12),
                        "StartTimeSeconds_metadata_only": "" if p.start_time_s is None else format_float(p.start_time_s, 12),
                        "EndTimeSeconds_metadata_only": "" if p.end_time_s is None else format_float(p.end_time_s, 12),
                        "Value": format_float(p.value, 12),
                    }
                )
    return rows


def participant_subjects(series: dict[tuple[str, str], list[SeriesPoint]]) -> list[str]:
    subjects = sorted({key[0] for key in series if key[0] != "M2025"})
    return subjects


def run_dtw_analysis(
    all_series: dict[str, dict[tuple[str, str], list[SeriesPoint]]],
    *,
    tolerance: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    summary_rows: list[dict[str, str]] = []
    path_rows: list[dict[str, str]] = []

    # union of all participant subjects present in either metric
    all_subjects = sorted(
        {
            key[0]
            for series in all_series.values()
            for key in series
            if key[0] != "M2025"
        }
    )

    for subject in all_subjects:
        for metric, series in all_series.items():
            ref_points = series.get(("M2025", "REF"), [])
            z1_points = series.get((subject, "Z1"), [])
            z2_points = series.get((subject, "Z2"), [])
            row: dict[str, str] = {
                "Subject": subject,
                "Metric": metric,
                "ReferenceSubject": "M2025",
                "ReferenceOccasion": "REF",
                "N_reference_valid_frames": str(len(ref_points)),
                "N_Z1_valid_frames": str(len(z1_points)),
                "N_Z2_valid_frames": str(len(z2_points)),
                "Z1_cumulative_cost": "",
                "Z1_path_length": "",
                "Z1_DTWnorm": "",
                "Z2_cumulative_cost": "",
                "Z2_path_length": "",
                "Z2_DTWnorm": "",
                "Delta_Z2_minus_Z1": "",
                "Direction": "not_computed",
                "Status": "",
            }

            if not ref_points:
                row["Status"] = "missing_reference_series"
                summary_rows.append(row)
                continue
            if not z1_points or not z2_points:
                missing = []
                if not z1_points:
                    missing.append("Z1")
                if not z2_points:
                    missing.append("Z2")
                row["Status"] = "missing_or_no_valid_" + "_and_".join(missing) + "_series"
                # Still write path data for any available occasion.
                for occasion, target_points in (("Z1", z1_points), ("Z2", z2_points)):
                    if target_points:
                        add_path_rows(path_rows, subject, metric, occasion, ref_points, target_points)
                summary_rows.append(row)
                continue

            try:
                ref_values = finite_values(ref_points, label=f"{metric} M2025 REF")
                z1_values = finite_values(z1_points, label=f"{metric} {subject} Z1")
                z2_values = finite_values(z2_points, label=f"{metric} {subject} Z2")
                r1 = path_normalized_dtw(ref_values, z1_values)
                r2 = path_normalized_dtw(ref_values, z2_values)
                delta = r2.normalized_distance - r1.normalized_distance
                row.update(
                    {
                        "Z1_cumulative_cost": format_float(r1.cumulative_cost, 12),
                        "Z1_path_length": str(r1.path_length),
                        "Z1_DTWnorm": format_float(r1.normalized_distance, 12),
                        "Z2_cumulative_cost": format_float(r2.cumulative_cost, 12),
                        "Z2_path_length": str(r2.path_length),
                        "Z2_DTWnorm": format_float(r2.normalized_distance, 12),
                        "Delta_Z2_minus_Z1": format_float(delta, 12),
                        "Direction": direction(delta, tolerance=tolerance),
                        "Status": "computed",
                    }
                )
                add_path_rows(path_rows, subject, metric, "Z1", ref_points, z1_points, r1)
                add_path_rows(path_rows, subject, metric, "Z2", ref_points, z2_points, r2)
            except Exception as exc:
                row["Status"] = f"error: {exc}"
            summary_rows.append(row)

    return summary_rows, path_rows


def add_path_rows(
    path_rows: list[dict[str, str]],
    subject: str,
    metric: str,
    occasion: str,
    ref_points: Sequence[SeriesPoint],
    target_points: Sequence[SeriesPoint],
    dtw_result: DTWResult | None = None,
) -> None:
    if dtw_result is None:
        try:
            dtw_result = path_normalized_dtw(
                finite_values(ref_points, label=f"{metric} M2025 REF"),
                finite_values(target_points, label=f"{metric} {subject} {occasion}"),
            )
        except Exception:
            return
    path_length = dtw_result.path_length
    for path_index, (i_ref, i_target) in enumerate(dtw_result.path):
        ref = ref_points[i_ref]
        target = target_points[i_target]
        normalized_path_position = path_index / (path_length - 1) if path_length > 1 else 0.0
        local_cost = abs(ref.value - target.value)
        path_rows.append(
            {
                "Metric": metric,
                "Subject": subject,
                "Occasion": occasion,
                "ReferenceSubject": "M2025",
                "ReferenceOccasion": "REF",
                "PathIndex": str(path_index),
                "PathLength": str(path_length),
                "NormalizedPathPosition": format_float(normalized_path_position, 12),
                "ReferenceOriginalFrameIndex": str(ref.original_frame_index),
                "ReferenceValidFrameIndex": str(ref.valid_frame_index),
                "ReferenceNormalizedFramePosition": format_float(ref.normalized_frame_position, 12),
                "ReferenceValue": format_float(ref.value, 12),
                "TargetOriginalFrameIndex": str(target.original_frame_index),
                "TargetValidFrameIndex": str(target.valid_frame_index),
                "TargetNormalizedFramePosition": format_float(target.normalized_frame_position, 12),
                "TargetValue": format_float(target.value, 12),
                "LocalAbsoluteCost": format_float(local_cost, 12),
                "Comparison_DTWnorm": format_float(dtw_result.normalized_distance, 12),
            }
        )


def write_metric_specific_summaries(output_dir: Path, summary_rows: Sequence[dict[str, str]]) -> None:
    fieldnames = [
        "Subject",
        "Metric",
        "ReferenceSubject",
        "ReferenceOccasion",
        "N_reference_valid_frames",
        "N_Z1_valid_frames",
        "N_Z2_valid_frames",
        "Z1_cumulative_cost",
        "Z1_path_length",
        "Z1_DTWnorm",
        "Z2_cumulative_cost",
        "Z2_path_length",
        "Z2_DTWnorm",
        "Delta_Z2_minus_Z1",
        "Direction",
        "Status",
    ]
    write_csv(output_dir / "SFR_Q_DTW_distance_delta_summary.csv", summary_rows, fieldnames)

    for metric in METRIC_CONFIG:
        rows = [r for r in summary_rows if r["Metric"] == metric]
        filename = f"{metric}_DTW_distance_delta_summary.csv"
        write_csv(output_dir / filename, rows, fieldnames)

def write_delta_sign_table(output_dir: Path, summary_rows: Sequence[dict[str, str]]) -> None:
    """Write a compact participant-level delta/sign table for SFR and Q value."""
    subjects = sorted({row["Subject"] for row in summary_rows})
    by_key = {(row["Subject"], row["Metric"]): row for row in summary_rows}
    fieldnames = [
        "Subject",
        "SFR_Delta_Z2_minus_Z1",
        "SFR_Sign",
        "SFR_Direction",
        "SFR_Z1_DTWnorm",
        "SFR_Z2_DTWnorm",
        "Q_value_Delta_Z2_minus_Z1",
        "Q_value_Sign",
        "Q_value_Direction",
        "Q_value_Z1_DTWnorm",
        "Q_value_Z2_DTWnorm",
        "Q_value_Status",
    ]
    rows: list[dict[str, str]] = []
    for subject in subjects:
        sfr = by_key.get((subject, "SFR"), {})
        qv = by_key.get((subject, "Q_value"), {})
        def sign_for(row: dict[str, str]) -> str:
            direction_value = row.get("Direction", "")
            if direction_value == "approach":
                return "-"
            if direction_value == "divergence":
                return "+"
            if direction_value == "unchanged":
                return "0"
            return "NA"
        rows.append({
            "Subject": subject,
            "SFR_Delta_Z2_minus_Z1": sfr.get("Delta_Z2_minus_Z1", ""),
            "SFR_Sign": sign_for(sfr),
            "SFR_Direction": sfr.get("Direction", ""),
            "SFR_Z1_DTWnorm": sfr.get("Z1_DTWnorm", ""),
            "SFR_Z2_DTWnorm": sfr.get("Z2_DTWnorm", ""),
            "Q_value_Delta_Z2_minus_Z1": qv.get("Delta_Z2_minus_Z1", ""),
            "Q_value_Sign": sign_for(qv),
            "Q_value_Direction": qv.get("Direction", ""),
            "Q_value_Z1_DTWnorm": qv.get("Z1_DTWnorm", ""),
            "Q_value_Z2_DTWnorm": qv.get("Z2_DTWnorm", ""),
            "Q_value_Status": qv.get("Status", ""),
        })
    write_csv(output_dir / "SFR_Q_delta_sign_table.csv", rows, fieldnames)


def make_plots(
    all_series: dict[str, dict[tuple[str, str], list[SeriesPoint]]],
    output_dir: Path,
) -> list[dict[str, str]]:
    if plt is None:
        raise RuntimeError("matplotlib is not installed; cannot create plots.")
    index_rows: list[dict[str, str]] = []
    plot_root = output_dir / "plots"
    plot_root.mkdir(parents=True, exist_ok=True)

    all_subjects = sorted(
        {
            key[0]
            for series in all_series.values()
            for key in series
            if key[0] != "M2025"
        }
    )

    for metric, series in all_series.items():
        cfg = METRIC_CONFIG[metric]
        metric_plot_dir = plot_root / cfg["plot_subdir"]
        metric_plot_dir.mkdir(parents=True, exist_ok=True)
        ref_points = series.get(("M2025", "REF"), [])
        for subject in all_subjects:
            fig, ax = plt.subplots(figsize=(8.5, 5.0))
            plotted_any = False
            warnings: list[str] = []
            for key, label in [
                (("M2025", "REF"), "M2025 reference"),
                ((subject, "Z1"), "Z1 before"),
                ((subject, "Z2"), "Z2 after"),
            ]:
                points = series.get(key, [])
                if not points:
                    warnings.append(f"{label}: no valid frames")
                    continue
                x = [p.normalized_frame_position for p in points]
                y = [p.value for p in points]
                ax.plot(x, y, marker="o", markersize=2.5, linewidth=1.2, label=label)
                plotted_any = True
            ax.set_title(f"{subject}: {metric} frame-wise series")
            ax.set_xlabel("Normalized frame position (0=start, 1=end; not clock time)")
            ax.set_ylabel(cfg["y_label"])
            ax.grid(True, linewidth=0.5, alpha=0.4)
            ax.legend(loc="best")
            if warnings:
                ax.text(
                    0.02,
                    0.02,
                    "; ".join(warnings),
                    transform=ax.transAxes,
                    fontsize=8,
                    va="bottom",
                    ha="left",
                )
            fig.tight_layout()
            out_path = metric_plot_dir / f"{subject}_{metric}_M2025_Z1_Z2_normalized_frame_position.png"
            if plotted_any:
                fig.savefig(out_path, dpi=180)
            plt.close(fig)
            if plotted_any:
                index_rows.append(
                    {
                        "Subject": subject,
                        "Metric": metric,
                        "PlotFile": str(out_path.relative_to(output_dir)),
                        "X_axis": "normalized_frame_position",
                        "X_axis_note": "0=start, 1=end; not clock time; computed after excluding short/silent/invalid frames",
                    }
                )
    write_csv(
        output_dir / "plot_index.csv",
        index_rows,
        ["Subject", "Metric", "PlotFile", "X_axis", "X_axis_note"],
    )
    return index_rows


def copy_script_and_readme(output_dir: Path, script_path: Path, readme_path: Path | None = None) -> None:
    code_dir = output_dir / "program"
    code_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(script_path, code_dir / script_path.name)
    if readme_path and readme_path.exists():
        shutil.copy2(readme_path, code_dir / "README.md")


def run(args: argparse.Namespace) -> None:
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    for metric in metrics:
        if metric not in METRIC_CONFIG:
            raise ValueError(f"Unsupported metric: {metric}. Choose from {', '.join(METRIC_CONFIG)}")

    all_series = {metric: read_metric_series(input_root, metric) for metric in metrics}

    valid_rows = build_valid_series_rows(all_series)
    valid_fieldnames = [
        "Metric",
        "Subject",
        "Occasion",
        "SourceFile",
        "OriginalFrameIndex",
        "ValidFrameIndex",
        "NValidFrames",
        "NormalizedFramePosition",
        "StartTimeSeconds_metadata_only",
        "EndTimeSeconds_metadata_only",
        "Value",
    ]
    write_csv(output_dir / "SFR_Q_valid_frame_series_normalized_position.csv", valid_rows, valid_fieldnames)

    summary_rows, path_rows = run_dtw_analysis(all_series, tolerance=args.direction_tolerance)
    write_metric_specific_summaries(output_dir, summary_rows)
    write_delta_sign_table(output_dir, summary_rows)

    path_fieldnames = [
        "Metric",
        "Subject",
        "Occasion",
        "ReferenceSubject",
        "ReferenceOccasion",
        "PathIndex",
        "PathLength",
        "NormalizedPathPosition",
        "ReferenceOriginalFrameIndex",
        "ReferenceValidFrameIndex",
        "ReferenceNormalizedFramePosition",
        "ReferenceValue",
        "TargetOriginalFrameIndex",
        "TargetValidFrameIndex",
        "TargetNormalizedFramePosition",
        "TargetValue",
        "LocalAbsoluteCost",
        "Comparison_DTWnorm",
    ]
    write_csv(output_dir / "SFR_Q_DTW_path_aligned_long.csv", path_rows, path_fieldnames)

    if args.make_plots:
        make_plots(all_series, output_dir)

    metadata_rows = [
        {"Item": "Python", "Value": platform.python_version()},
        {"Item": "NumPy", "Value": np.__version__},
        {"Item": "Input root", "Value": str(input_root)},
        {"Item": "Output dir", "Value": str(output_dir)},
        {"Item": "Metrics", "Value": ",".join(metrics)},
        {"Item": "Local cost", "Value": "absolute difference abs(x_i - y_j)"},
        {"Item": "Permitted steps", "Value": "diagonal, vertical, horizontal"},
        {"Item": "Global warping window", "Value": "none"},
        {"Item": "Distance normalization", "Value": "cumulative cost divided by optimal path length"},
        {"Item": "Delta convention", "Value": "Delta_Z2_minus_Z1 = DTWnorm(M2025,Z2) - DTWnorm(M2025,Z1)"},
        {"Item": "Negative delta", "Value": "distance decreased: approach toward M2025 at Z2"},
        {"Item": "Positive delta", "Value": "distance increased: divergence from M2025 at Z2"},
        {"Item": "Frame inclusion", "Value": "only rows with include_in_median=True and finite metric values"},
        {"Item": "Plot x-axis", "Value": "normalized_frame_position, not clock time"},
    ]
    write_csv(output_dir / "analysis_metadata.csv", metadata_rows, ["Item", "Value"])

    print(f"Wrote DTW outputs to: {output_dir}")
    print(f"Summary rows: {len(summary_rows)}")
    print(f"Path-aligned rows: {len(path_rows)}")
    if args.make_plots:
        print(f"Plots written under: {output_dir / 'plots'}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute path-normalized DTW distances and plots for frame-wise SFR and Q-value series."
    )
    parser.add_argument(
        "input_root",
        type=Path,
        help="Root directory containing 'sfr' and 'q_value' frame-wise CSV subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("dtw_sfr_q_results"),
        help="Output directory. Default: dtw_sfr_q_results",
    )
    parser.add_argument(
        "--metrics",
        default="SFR,Q_value",
        help="Comma-separated metric list. Default: SFR,Q_value",
    )
    parser.add_argument(
        "--direction-tolerance",
        type=float,
        default=0.0,
        help="Optional non-negative equivalence band around zero for direction labels.",
    )
    parser.add_argument(
        "--no-plots",
        dest="make_plots",
        action="store_false",
        help="Do not create PNG plots.",
    )
    parser.set_defaults(make_plots=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.direction_tolerance < 0:
            raise ValueError("--direction-tolerance must be non-negative.")
        run(args)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
