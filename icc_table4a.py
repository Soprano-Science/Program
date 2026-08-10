#!/usr/bin/env python3
"""
Reproducible ICC analysis for Table 4A.

Design
------
20 targets (participants) x 4 raters, complete and balanced.

Primary coefficient
-------------------
ICC(2,k) in Shrout & Fleiss notation, equivalent to ICC(A,k) in
McGraw & Wong notation:
- two-way random-effects model
- absolute agreement
- average of k=4 raters

The script also reports ICC(2,1) / ICC(A,1) for transparency.

Point-estimate formulas are ANOVA mean-square formulas.
95% confidence intervals use the F-based McGraw-Wong procedure for
the two-way absolute-agreement ICC, with the average-measures CI
obtained by the standard k-rater transformation.

References
----------
Shrout PE, Fleiss JL. Intraclass correlations: uses in assessing
rater reliability. Psychological Bulletin. 1979;86(2):420-428.
doi:10.1037/0033-2909.86.2.420

McGraw KO, Wong SP. Forming inferences about some intraclass
correlation coefficients. Psychological Methods. 1996;1(1):30-46.
doi:10.1037/1082-989X.1.1.30
"""

from pathlib import Path
import csv
import numpy as np
from scipy.stats import f

DATA_FILE = Path(__file__).with_name("table4a_ratings.csv")
ALPHA = 0.05

def load_ratings(path):
    subjects = []
    rows = []
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        expected = ["Subject", "Rater1", "Rater2", "Rater3", "Rater4"]
        if reader.fieldnames != expected:
            raise ValueError(
                f"Expected columns {expected}, but found {reader.fieldnames}"
            )
        for rec in reader:
            subjects.append(rec["Subject"])
            rows.append([float(rec[c]) for c in expected[1:]])

    x = np.asarray(rows, dtype=float)

    if x.shape != (20, 4):
        raise ValueError(f"Expected a 20 x 4 matrix; found {x.shape}.")
    if not np.isfinite(x).all():
        raise ValueError("Missing or non-finite ratings detected.")
    if not np.isin(x, [0.0, 1.0, 2.0]).all():
        raise ValueError("Ratings must be coded 0, 1, or 2.")

    return subjects, x

def icc_two_way_absolute(x, alpha=0.05):
    """
    ICC(2,1) / ICC(A,1) and ICC(2,k) / ICC(A,k).

    Rows = targets; columns = raters.
    """
    n, k = x.shape
    grand = x.mean()
    row_mean = x.mean(axis=1)
    col_mean = x.mean(axis=0)

    ss_rows = k * np.sum((row_mean - grand) ** 2)
    ss_cols = n * np.sum((col_mean - grand) ** 2)
    ss_total = np.sum((x - grand) ** 2)
    ss_error = ss_total - ss_rows - ss_cols

    df_rows = n - 1
    df_cols = k - 1
    df_error = (n - 1) * (k - 1)

    ms_rows = ss_rows / df_rows
    ms_cols = ss_cols / df_cols
    ms_error = ss_error / df_error

    # Shrout-Fleiss ICC(2,1) = McGraw-Wong ICC(A,1)
    icc_2_1 = (
        (ms_rows - ms_error)
        / (
            ms_rows
            + (k - 1) * ms_error
            + k * (ms_cols - ms_error) / n
        )
    )

    # Shrout-Fleiss ICC(2,k) = McGraw-Wong ICC(A,k)
    icc_2_k = (
        (ms_rows - ms_error)
        / (ms_rows + (ms_cols - ms_error) / n)
    )

    # F test associated with the row/target effect
    F_value = ms_rows / ms_error
    p_value = f.sf(F_value, df_rows, df_error)

    # McGraw-Wong F-based CI for ICC(A,1)
    fj = ms_cols / ms_error
    v_num = df_error * (
        k * icc_2_1 * fj
        + n * (1 + (k - 1) * icc_2_1)
        - k * icc_2_1
    ) ** 2

    v_den = (
        df_rows * k**2 * icc_2_1**2 * fj**2
        + (
            n * (1 + (k - 1) * icc_2_1)
            - k * icc_2_1
        ) ** 2
    )
    v = v_num / v_den

    f_upper = f.ppf(1 - alpha / 2, n - 1, v)
    f_lower = f.ppf(1 - alpha / 2, v, n - 1)

    lower_2_1 = (
        n * (ms_rows - f_upper * ms_error)
        / (
            f_upper
            * (k * ms_cols + (k * n - k - n) * ms_error)
            + n * ms_rows
        )
    )

    upper_2_1 = (
        n * (f_lower * ms_rows - ms_error)
        / (
            k * ms_cols
            + (k * n - k - n) * ms_error
            + n * f_lower * ms_rows
        )
    )

    # Convert single-measure CI to average-measures CI.
    lower_2_k = k * lower_2_1 / (1 + (k - 1) * lower_2_1)
    upper_2_k = k * upper_2_1 / (1 + (k - 1) * upper_2_1)

    return {
        "n_targets": n,
        "k_raters": k,
        "MS_targets": ms_rows,
        "MS_raters": ms_cols,
        "MS_error": ms_error,
        "ICC_2_1": icc_2_1,
        "ICC_2_1_CI_low": lower_2_1,
        "ICC_2_1_CI_high": upper_2_1,
        "ICC_2_k": icc_2_k,
        "ICC_2_k_CI_low": lower_2_k,
        "ICC_2_k_CI_high": upper_2_k,
        "F": F_value,
        "df1": df_rows,
        "df2": df_error,
        "p": p_value,
    }

def main():
    subjects, x = load_ratings(DATA_FILE)
    r = icc_two_way_absolute(x, alpha=ALPHA)

    print("Table 4A inter-rater reliability")
    print(f"Targets (participants): {r['n_targets']}")
    print(f"Raters: {r['k_raters']}")
    print()
    print("Two-way random-effects, absolute agreement")
    print(
        "ICC(2,1) / ICC(A,1) = "
        f"{r['ICC_2_1']:.6f} "
        f"(95% CI {r['ICC_2_1_CI_low']:.6f} to {r['ICC_2_1_CI_high']:.6f})"
    )
    print(
        "ICC(2,k) / ICC(A,k) = "
        f"{r['ICC_2_k']:.6f} "
        f"(95% CI {r['ICC_2_k_CI_low']:.6f} to {r['ICC_2_k_CI_high']:.6f})"
    )
    print()
    print(
        f"F({r['df1']}, {r['df2']}) = {r['F']:.6f}, "
        f"p = {r['p']:.10g}"
    )
    print()
    print("ANOVA mean squares")
    print(f"MS_targets = {r['MS_targets']:.12f}")
    print(f"MS_raters  = {r['MS_raters']:.12f}")
    print(f"MS_error   = {r['MS_error']:.12f}")

if __name__ == "__main__":
    main()
