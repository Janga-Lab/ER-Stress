# TRAIN Registry Age-at-Diagnosis Modeling: Comprehensive Interpretation Guide

## 1) Why this analysis exists

The central question is:

Can we estimate **age at official diagnosis** for phenotype outcomes from genotype information in the TRAIN registry?

This work is set up as a deliberately strict baseline so we can answer that question credibly before adding higher-complexity methods.

The current pipeline is designed to:

1. Use standardized and clinically interpretable labels.
2. Avoid obvious leakage from repeated genotype patterns.
3. Produce transparent models with understandable feature effects.
4. Create artifacts and visuals that support scientific discussion with collaborators.

## 2) Data and assumptions we enforced

Dataset used:

- `analysis/TRAIN_train_standardized.xlsx`
- 277 rows (patients)

Clinical interpretation assumptions (from collaborator guidance, enforced by preprocessing):

1. Numeric values are age at **official physician diagnosis** in years.
2. Decimal values are fractional years.
3. `Y` means present without diagnosis age.
4. `N` means absent at the time recorded (not lifetime guaranteed absence).
5. Blank/UNK/dont know are unknown, not negative.
6. Mixed free text was parsed, but baseline modeling only uses trustworthy exact numeric diagnosis ages.

Why this matters:

- If we mix "present without age", unknowns, and exact ages in the same target, the model would learn label noise rather than diagnosis timing.
- Keeping only exact numeric diagnosis ages gives cleaner supervision at the cost of sample size.

## 3) Why we restricted targets to `exact_numeric`

For each phenotype, we created the modeling set using:

- `P__age_type == "exact_numeric"`
- `P__age_dx_years` non-null

This choice was intentional:

1. It preserves label meaning and reduces ambiguity.
2. It gives a stable baseline benchmark.
3. It prevents implicit imputation assumptions from driving initial conclusions.

Tradeoff:

- We lose data volume and some potentially informative records (upper bounds, intervals, present-no-age).

Current target availability:

| Phenotype | Total rows | exact_numeric | present_no_age | absent | unknown-like | other parsed types |
|---|---:|---:|---:|---:|---:|---:|
| Diabetes Mellitus | 277 | 142 | 8 | 41 | 84 | 2 |
| Optic Atrophy | 277 | 136 | 17 | 22 | 90 | 12 |
| Hearing Loss | 277 | 63 | 18 | 84 | 96 | 16 |
| Diabetes Insipidus | 277 | 48 | 7 | 105 | 112 | 5 |

Interpretation:

- Diabetes Mellitus and Optic Atrophy are the strongest baseline candidates by label count.
- Hearing Loss is viable but lower-powered.
- Diabetes Insipidus is likely underpowered for this exact-only regression baseline.

## 4) Why these features were used

We focused on genotype-derived predictors that are biologically plausible and robustly available:

1. Mutation category abstraction:
   - `Mutation_1_primary`, `Mutation_2_primary`, pair combinations
   - truncating indicators and counts
   - Rationale: broad functional class can capture severity direction without requiring large n for exact variant strings.
2. Position features:
   - `AA_pos_1`, `AA_pos_2`, min/max/mean/span, missingness flags
   - Rationale: protein position effects are often informative for phenotype timing and severity.
3. Topology features:
   - `TM_1`, `TM_2`, `any_TM`, `both_TM`
   - Rationale: transmembrane region involvement may map to functional impact.
4. Allelic pattern features:
   - missing second allele, homozygous-known, compound-het indicator
   - Rationale: zygosity and completeness of allelic characterization can influence risk and uncertainty.

Optional covariates (`Sex_clean`, `Diagnosis_clean`) are available in the script but were not included in default runs to keep the baseline centered on genotype-only signal.

## 5) Why grouped cross-validation was used

Evaluation used `GroupKFold` with group key = sorted genotype signature (`Allele_1_clean|Allele_2_clean`).

Reason:

- Without grouping, the same or near-identical genotype pattern can appear in both train and test folds, inflating performance.
- Grouped splits force better estimation of generalization beyond previously seen genotype signatures.

Current grouped sample structure:

| Phenotype | Exact-n rows | Unique genotype signatures |
|---|---:|---:|
| Diabetes Mellitus | 142 | 112 |
| Optic Atrophy | 136 | 115 |
| Hearing Loss | 63 | 56 |

## 6) Why these two baseline models

We fit two complementary regressors:

1. Elastic Net:
   - linear, regularized, interpretable coefficients.
   - Good for checking whether simple additive signal exists.
2. Gradient Boosting Regressor:
   - non-linear, interaction-friendly.
   - Good for modest tabular datasets with mixed feature types.

Using both gives a practical diagnostic:

- If both models struggle, the limiting factor is likely data signal quality/quantity, not just model class.
- If tree model consistently wins, important non-linearities/interactions are likely present.

## 7) What each plot means and why it is in the workflow

Each phenotype now has these figures under:

- `analysis/model_outputs/<phenotype>/figures/`

### `target_distribution.png`

What it answers:

- Is the target heavily skewed?
- Are there tails that may dominate error metrics?

Why it matters:

- Age-at-diagnosis is not uniformly distributed. Strong right tails can make RMSE jump and can bias models toward middle ages.

Current outputs:

Diabetes Mellitus  
![Diabetes Mellitus target distribution](model_outputs/diabetes_mellitus/figures/target_distribution.png)

Optic Atrophy  
![Optic Atrophy target distribution](model_outputs/optic_atrophy/figures/target_distribution.png)

Hearing Loss  
![Hearing Loss target distribution](model_outputs/hearing_loss/figures/target_distribution.png)

### `cv_metrics.png`

What it answers:

- Which baseline model is better on average?
- How stable are fold-level errors (mean +/- SD)?

Why it matters:

- A lower mean with high SD indicates instability.
- Similar means but different SD can still inform model preference.

Current outputs:

Diabetes Mellitus  
![Diabetes Mellitus CV metrics](model_outputs/diabetes_mellitus/figures/cv_metrics.png)

Optic Atrophy  
![Optic Atrophy CV metrics](model_outputs/optic_atrophy/figures/cv_metrics.png)

Hearing Loss  
![Hearing Loss CV metrics](model_outputs/hearing_loss/figures/cv_metrics.png)

### `oof_observed_vs_predicted.png`

What it answers:

- Does the model track the identity line (perfect prediction)?
- Is there compression toward mean age (underpredict older cases, overpredict younger)?

Why it matters:

- This is the most intuitive check for calibration-like behavior in regression.

Current outputs:

Diabetes Mellitus  
![Diabetes Mellitus observed vs predicted](model_outputs/diabetes_mellitus/figures/oof_observed_vs_predicted.png)

Optic Atrophy  
![Optic Atrophy observed vs predicted](model_outputs/optic_atrophy/figures/oof_observed_vs_predicted.png)

Hearing Loss  
![Hearing Loss observed vs predicted](model_outputs/hearing_loss/figures/oof_observed_vs_predicted.png)

### `residual_distribution.png`

What it answers:

- Is there systematic directional bias?
- Are residuals wide or heavy-tailed?

Why it matters:

- Median/shape differences reveal whether errors are random noise or structurally biased.

Current outputs:

Diabetes Mellitus  
![Diabetes Mellitus residual distribution](model_outputs/diabetes_mellitus/figures/residual_distribution.png)

Optic Atrophy  
![Optic Atrophy residual distribution](model_outputs/optic_atrophy/figures/residual_distribution.png)

Hearing Loss  
![Hearing Loss residual distribution](model_outputs/hearing_loss/figures/residual_distribution.png)

### `feature_importance_top.png`

What it answers:

- Which encoded predictors dominate each model?
- Are results biologically plausible and consistent across model types?

Why it matters:

- It creates hypothesis seeds for mechanistic follow-up and feature redesign.

Current outputs:

Diabetes Mellitus  
![Diabetes Mellitus feature importance](model_outputs/diabetes_mellitus/figures/feature_importance_top.png)

Optic Atrophy  
![Optic Atrophy feature importance](model_outputs/optic_atrophy/figures/feature_importance_top.png)

Hearing Loss  
![Hearing Loss feature importance](model_outputs/hearing_loss/figures/feature_importance_top.png)

## 8) Current quantitative results and interpretation

Primary fold-averaged CV metrics (`cv_summary_metrics.csv`):

| Phenotype | Model | MAE mean | RMSE mean | Comment |
|---|---|---:|---:|---|
| Diabetes Mellitus | Gradient Boosting | 7.50 | 10.63 | Best MAE among tested models |
| Diabetes Mellitus | Elastic Net | 7.93 | 11.82 | Slightly worse overall |
| Optic Atrophy | Gradient Boosting | 9.92 | 13.65 | Best MAE |
| Optic Atrophy | Elastic Net | 10.54 | 13.57 | Similar RMSE, worse MAE |
| Hearing Loss | Gradient Boosting | 9.51 | 14.57 | Best MAE but high fold variability |
| Hearing Loss | Elastic Net | 9.79 | 14.09 | Slightly better fold-mean RMSE |

Global out-of-fold diagnostics (`model_insights_summary.csv`):

- R2 is negative for all phenotypes/models in current runs.
- This means predictions are still not explaining variance better than a naive mean predictor in OOF space.
- MAE remains useful as an absolute error scale, but current predictive strength is modest.

Interpretation:

1. There is some usable signal (MAE improvements and structured importances), but it is not yet strong enough for high-confidence individual prediction.
2. Model difficulty is highest in later diagnosis ages (older-tail cases).
3. Gradient boosting tends to capture the available non-linear structure better on MAE, but not uniformly on all summary metrics.

Age-quartile error pattern (best model per phenotype, from OOF):

| Phenotype | Lower-age bins MAE | Highest-age bin MAE | Pattern |
|---|---:|---:|---|
| Diabetes Mellitus | about 4.53 to 7.04 | about 14.01 | Large error inflation in oldest quartile |
| Optic Atrophy | about 4.89 to 9.25 | about 18.64 | Strong oldest-bin difficulty |
| Hearing Loss | about 5.35 to 7.20 | about 18.52 | Same tail-error inflation with smaller n |

Clinical story from this pattern:

- Current models are reasonably anchored for earlier-to-mid diagnosis ages.
- They are much less reliable for late diagnosis ages, likely due heterogeneity, sparser examples, and unmodeled covariates.

## 9) What the data is telling us right now

The most defensible current narrative:

1. **Genotype structure contains timing signal**, but it is incomplete.
2. **Amino acid position features repeatedly dominate**, suggesting positional biology is relevant to onset timing.
3. **Mutation class interactions matter**, but effects are less stable in smaller phenotypes.
4. **Late-diagnosis cases are the hardest**, indicating missing explanatory variables and/or heterogeneous pathways.
5. **Prediction quality is baseline-grade**, suitable for hypothesis generation, not clinical deployment.

## 10) What is still missing from the story

Key missing context that limits interpretation:

1. Follow-up time and censoring context for absent (`N`) cases.
2. Family structure and relatedness handling.
3. Richer variant-level biology (domain annotations, conservation, functional scores).
4. Treatment/environment covariates that influence diagnosis timing.
5. Temporal and site-related factors in clinical ascertainment.

Without these, the model can understate true signal or misattribute effects.

## 11) How to improve model quality (technical roadmap)

### Phase A: Data quality and target design (highest leverage)

1. Add follow-up age/date per patient.
2. Recode "present_no_age", upper-bound, and intervals into censor-aware targets.
3. Curate outlier text entries with expert review.
4. Add confidence flags per label and propagate into weighting.

Expected impact:

- Better target fidelity and larger usable sample.

### Phase B: Evaluation hardening

1. Keep grouped CV by genotype signature as default.
2. Add repeated grouped CV for confidence intervals.
3. Add permutation tests to quantify whether observed performance exceeds chance.
4. Reserve a final untouched holdout if sample allows.

Expected impact:

- More trustworthy claims and less optimistic bias.

### Phase C: Modeling upgrades

1. Tune hyperparameters via nested grouped CV.
2. Add robust losses (Huber, quantile regression) for tail behavior.
3. Model transformed target (`log1p(age)`) and compare tail calibration.
4. Add monotonic or GAM-style models for interpretable non-linearity.

Expected impact:

- Better robustness and clearer effect-shape understanding.

### Phase D: Censor-aware modeling

1. Build survival models once follow-up ages are available:
   - Cox PH / AFT baselines
   - Random survival forests
2. Treat `N` as right-censored at last follow-up instead of dropping.
3. Evaluate with concordance and integrated Brier score.

Expected impact:

- Uses much more data and aligns more closely with natural clinical time-to-event framing.

### Phase E: Biological feature enrichment

1. Protein domain membership and structural region features.
2. External variant scores (conservation, deleteriousness, splice impact).
3. Variant-level embeddings or grouped burden features.
4. Mutation interaction terms guided by prior biology.

Expected impact:

- Better biological specificity and potentially stronger predictive performance.

## 12) How to improve the visual and scientific story

A stronger story is not only lower error; it is clearer evidence logic.

Recommended narrative sequence for manuscripts/slides:

1. Cohort flow:
   - 277 rows -> parsed labels -> exact numeric modeling sets by phenotype.
2. Label semantics and uncertainty:
   - show why exact-only baseline is trustworthy.
3. Baseline performance:
   - MAE/RMSE with grouped CV and uncertainty.
4. Behavior diagnostics:
   - OOF scatter + residual distributions + age-bin error.
5. Mechanistic clues:
   - feature importance consistency across phenotypes/models.
6. Gap analysis and next study:
   - censored modeling plan with follow-up data.

This storyline is stronger because it separates:

- what is known now,
- what is uncertain,
- and what experiment closes each uncertainty.

## 13) Practical commands to reproduce

Train:

```bash
python3 analysis/train_age_dx_model.py \
  --data analysis/TRAIN_train_standardized.xlsx \
  --phenotype "Hearing Loss" \
  --output-dir analysis/model_outputs \
  --n-splits 5
```

Visualize:

```bash
MPLCONFIGDIR=/tmp/mpl XDG_CACHE_HOME=/tmp/.cache python3 analysis/visualize_age_dx_models.py \
  --data analysis/TRAIN_train_standardized.xlsx \
  --phenotype "Hearing Loss" \
  --output-dir analysis/model_outputs \
  --n-splits 5
```

Swap phenotype to `"Diabetes Mellitus"` or `"Optic Atrophy"` to regenerate corresponding outputs.

## 14) Bottom line

This baseline is doing what it should do:

1. It provides a transparent and defensible first estimate of genotype-to-diagnosis-age signal.
2. It identifies where the signal appears strongest (position and mutation pattern features).
3. It reveals key failure modes (late-age tails, limited variance capture).
4. It defines a concrete path to improve both predictive performance and scientific confidence.

Use current outputs for hypothesis generation and collaborator discussion, not for clinical decision support.
