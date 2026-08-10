# ICC analysis for Table 4A

This folder contains the de-identified 20 × 4 development-rating matrix used in Table 4A and a reproducible Python script for inter-rater reliability.

## Files

- `table4a_ratings.csv` — 20 participants × 4 expert-rater development ratings (0, 1, 2).
- `icc_table4a.py` — calculates ICC(2,1) and ICC(2,k), equivalent to ICC(A,1) and ICC(A,k), using a two-way random-effects, absolute-agreement model.
- `icc_results.txt` — output obtained from the supplied Table 4A data.
- `requirements_icc.txt` — Python package versions used for the reproducibility check.

## Run

From this folder:

```bash
python -m pip install -r requirements_icc.txt
python icc_table4a.py
```

## Primary result

For the Table 4A matrix:

- ICC(2,k) / ICC(A,k) = **0.731**
- 95% CI = **0.436–0.884**
- k = 4 raters; n = 20 participants

For transparency, the corresponding single-measure ICC is:

- ICC(2,1) / ICC(A,1) = **0.404**
- 95% CI = **0.162–0.655**

The ICC is reported descriptively for this fixed four-evaluator panel and is not used as a basis for population-level generalization.

## Method

The point estimates use the balanced two-way ANOVA mean-square formulas for Shrout-Fleiss ICC(2,1) and ICC(2,k), corresponding to McGraw-Wong absolute-agreement ICC(A,1) and ICC(A,k). The 95% confidence interval uses the F-based McGraw-Wong procedure.

## References

Shrout PE, Fleiss JL. Intraclass correlations: uses in assessing rater reliability. *Psychological Bulletin*. 1979;86(2):420–428. doi:10.1037/0033-2909.86.2.420

McGraw KO, Wong SP. Forming inferences about some intraclass correlation coefficients. *Psychological Methods*. 1996;1(1):30–46. doi:10.1037/1082-989X.1.1.30
