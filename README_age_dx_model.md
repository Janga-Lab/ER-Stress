# Age-at-Diagnosis Baseline Modeling (TRAIN standardized data)

This folder now includes a reproducible training script:

- `train_age_dx_model.py`
- `visualize_age_dx_models.py`

Comprehensive interpretation and roadmap:

- `README_modeling_story.md`

## What it does

1. Loads standardized TRAIN data (`.xlsx` or `.csv`).
2. Filters to rows with `target__age_type == "exact_numeric"` and non-null `target__age_dx_years`.
3. Builds genotype-derived features from mutation classes, truncation flags, AA positions, TM flags, and allele patterns.
4. Uses grouped cross-validation (`GroupKFold`) by sorted genotype signature (`Allele_1_clean|Allele_2_clean`).
5. Trains and compares:
   - Elastic Net regression
   - Gradient Boosting regressor
6. Saves metrics, feature importances, model artifacts, and metadata.

## Run

From repo root (`/Users/kvand/Documents/Diabeties`):

```bash
python3 train_age_dx_model.py \
  --data TRAIN_train_standardized.xlsx \
  --phenotype "Diabetes Mellitus" \
  --output-dir model_outputs \
  --n-splits 5
```

Optional covariates:

```bash
python3 train_age_dx_model.py \
  --data TRAIN_train_standardized.xlsx \
  --phenotype "Optic Atrophy" \
  --include-sex \
  --include-diagnosis
```

## Outputs

For each phenotype, outputs are under:

- `model_outputs/<phenotype_slug>/`

Key files:

- `cv_fold_metrics.csv`
- `cv_summary_metrics.csv`
- `analysis_dataset_snapshot.csv`
- `elastic_net_pipeline.joblib`
- `gradient_boosting_pipeline.joblib`
- `elastic_net_feature_importance.csv`
- `gradient_boosting_feature_importance.csv`
- `run_metadata.json`

## Visualizations and model insights

Generate insight plots and out-of-fold diagnostics:

```bash
MPLCONFIGDIR=/tmp/mpl XDG_CACHE_HOME=/tmp/.cache python3 visualize_age_dx_models.py \
  --data TRAIN_train_standardized.xlsx \
  --phenotype "Hearing Loss" \
  --output-dir model_outputs \
  --n-splits 5
```

This writes:

- `cv_oof_predictions.csv`
- `model_insights_summary.csv`
- `MODEL_INSIGHTS.md`
- `figures/target_distribution.png`
- `figures/cv_metrics.png`
- `figures/oof_observed_vs_predicted.png`
- `figures/residual_distribution.png`
- `figures/feature_importance_top.png`

## Prediction function (new standardized rows)

The script exposes:

- `predict_age_at_diagnosis(standardized_rows, model_path, include_sex=False, include_diagnosis=False)`

Example:

```python
from pathlib import Path
import pandas as pd
from train_age_dx_model import predict_age_at_diagnosis

rows = pd.read_excel("TRAIN_train_standardized.xlsx").head(3)
pred = predict_age_at_diagnosis(
    standardized_rows=rows,
    model_path=Path("model_outputs/diabetes_mellitus/gradient_boosting_pipeline.joblib"),
    include_sex=False,
    include_diagnosis=False,
)
print(pred)
```
