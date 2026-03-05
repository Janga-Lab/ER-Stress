#!/usr/bin/env python3
"""Train baseline age-at-diagnosis models from standardized TRAIN registry data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import joblib
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


RANDOM_STATE = 42


def make_one_hot_encoder() -> OneHotEncoder:
    """Build a OneHotEncoder compatible across sklearn versions."""
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def load_standardized_data(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported input format: {path.suffix}")


def normalize_text(series: pd.Series) -> pd.Series:
    values = series.astype("string").str.strip()
    values = values.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "NA": pd.NA})
    return values


def combine_pair_sorted(left: pd.Series, right: pd.Series, missing_label: str = "missing") -> pd.Series:
    left_values = normalize_text(left).fillna(missing_label)
    right_values = normalize_text(right).fillna(missing_label)
    combined = []
    for lv, rv in zip(left_values.tolist(), right_values.tolist()):
        combined.append("|".join(sorted([str(lv), str(rv)])))
    return pd.Series(combined, index=left.index, dtype="string")


def build_genotype_signature(df: pd.DataFrame) -> pd.Series:
    allele_1 = normalize_text(df["Allele_1_clean"])
    allele_2 = normalize_text(df["Allele_2_clean"])
    signatures: List[str] = []
    for a1, a2 in zip(allele_1.tolist(), allele_2.tolist()):
        present = [a for a in [a1, a2] if isinstance(a, str) and a]
        if not present:
            signatures.append("MISSING_ALLELES")
            continue
        signatures.append("|".join(sorted(present)))
    return pd.Series(signatures, index=df.index, dtype="string")


def build_feature_frame(
    df: pd.DataFrame,
    include_sex: bool = False,
    include_diagnosis: bool = False,
) -> pd.DataFrame:
    features = pd.DataFrame(index=df.index)

    features["Mutation_1_primary"] = normalize_text(df["Mutation_1_primary"]).fillna("unknown")
    features["Mutation_2_primary"] = normalize_text(df["Mutation_2_primary"]).fillna("unknown")
    features["Mutation_pair_primary"] = combine_pair_sorted(
        df["Mutation_1_primary"], df["Mutation_2_primary"], missing_label="unknown"
    )

    m1_trunc = pd.to_numeric(df["Mutation_1_is_truncating"], errors="coerce")
    m2_trunc = pd.to_numeric(df["Mutation_2_is_truncating"], errors="coerce")
    features["Mutation_1_is_truncating"] = m1_trunc
    features["Mutation_2_is_truncating"] = m2_trunc
    features["truncating_count"] = m1_trunc.fillna(0) + m2_trunc.fillna(0)
    features["any_truncating"] = ((m1_trunc.fillna(0) > 0) | (m2_trunc.fillna(0) > 0)).astype(int)

    aa1 = pd.to_numeric(df["AA_pos_1"], errors="coerce")
    aa2 = pd.to_numeric(df["AA_pos_2"], errors="coerce")
    aa_stack = pd.concat([aa1, aa2], axis=1)
    features["AA_pos_1"] = aa1
    features["AA_pos_2"] = aa2
    features["AA_pos_min"] = aa_stack.min(axis=1, skipna=True)
    features["AA_pos_max"] = aa_stack.max(axis=1, skipna=True)
    features["AA_pos_mean"] = aa_stack.mean(axis=1, skipna=True)
    features["AA_pos_span"] = features["AA_pos_max"] - features["AA_pos_min"]
    features["AA_pos_missing_count"] = aa_stack.isna().sum(axis=1)
    features["AA_pos_any_missing"] = (features["AA_pos_missing_count"] > 0).astype(int)

    tm1 = pd.to_numeric(df["TM_1"], errors="coerce")
    tm2 = pd.to_numeric(df["TM_2"], errors="coerce")
    features["TM_1"] = tm1
    features["TM_2"] = tm2
    features["any_TM"] = ((tm1.fillna(0) > 0) | (tm2.fillna(0) > 0)).astype(int)
    features["both_TM"] = ((tm1.fillna(0) > 0) & (tm2.fillna(0) > 0)).astype(int)

    allele_1 = normalize_text(df["Allele_1_clean"])
    allele_2 = normalize_text(df["Allele_2_clean"])
    features["allele1_missing"] = allele_1.isna().astype(int)
    features["allele2_missing"] = allele_2.isna().astype(int)
    features["allele_count_known"] = (~allele_1.isna()).astype(int) + (~allele_2.isna()).astype(int)
    features["compound_het"] = ((~allele_1.isna()) & (~allele_2.isna()) & (allele_1 != allele_2)).astype(int)
    features["homozygous_known"] = ((~allele_1.isna()) & (~allele_2.isna()) & (allele_1 == allele_2)).astype(int)

    if include_sex:
        features["Sex_clean"] = normalize_text(df["Sex_clean"]).fillna("unknown")
    if include_diagnosis:
        features["Diagnosis_clean"] = normalize_text(df["Diagnosis_clean"]).fillna("unknown")

    return features


def split_columns(features: pd.DataFrame) -> Tuple[List[str], List[str]]:
    categorical_cols = [col for col in features.columns if features[col].dtype.name in {"object", "string"}]
    numeric_cols = [col for col in features.columns if col not in categorical_cols]
    return categorical_cols, numeric_cols


def make_preprocessor(categorical_cols: Iterable[str], numeric_cols: Iterable[str]) -> ColumnTransformer:
    categorical_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", make_one_hot_encoder()),
        ]
    )
    numeric_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )
    return ColumnTransformer(
        transformers=[
            ("categorical", categorical_pipeline, list(categorical_cols)),
            ("numeric", numeric_pipeline, list(numeric_cols)),
        ]
    )


def make_cv_splits(
    X: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    n_splits: int,
) -> Tuple[List[Tuple[np.ndarray, np.ndarray]], str]:
    n_samples = len(X)
    if n_samples < 2:
        raise ValueError("Need at least 2 samples for cross-validation.")

    n_splits_effective = min(n_splits, n_samples)
    unique_groups = int(pd.Series(groups).nunique(dropna=False))

    if unique_groups >= n_splits_effective and n_splits_effective >= 2:
        splitter = GroupKFold(n_splits=n_splits_effective)
        splits = list(splitter.split(X, y, groups=groups))
        return splits, f"GroupKFold(n_splits={n_splits_effective})"

    if n_splits_effective >= 2:
        splitter = KFold(n_splits=n_splits_effective, shuffle=True, random_state=RANDOM_STATE)
        splits = list(splitter.split(X, y))
        return splits, f"KFold(n_splits={n_splits_effective}, shuffle=True)"

    raise ValueError("Unable to create cross-validation splits.")


def evaluate_model_cv(
    model_name: str,
    estimator,
    preprocessor: ColumnTransformer,
    X: pd.DataFrame,
    y: pd.Series,
    cv_splits: List[Tuple[np.ndarray, np.ndarray]],
) -> pd.DataFrame:
    records = []
    for fold_idx, (train_idx, test_idx) in enumerate(cv_splits, start=1):
        pipeline = Pipeline(
            steps=[
                ("preprocessor", clone(preprocessor)),
                ("model", clone(estimator)),
            ]
        )
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        pipeline.fit(X_train, y_train)
        preds = pipeline.predict(X_test)

        mae = mean_absolute_error(y_test, preds)
        rmse = float(np.sqrt(mean_squared_error(y_test, preds)))
        records.append(
            {
                "model": model_name,
                "fold": fold_idx,
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "mae": float(mae),
                "rmse": rmse,
            }
        )

    return pd.DataFrame.from_records(records)


def feature_importance_frame(pipeline: Pipeline) -> pd.DataFrame:
    preprocessor: ColumnTransformer = pipeline.named_steps["preprocessor"]
    model = pipeline.named_steps["model"]
    feature_names = preprocessor.get_feature_names_out()

    if hasattr(model, "coef_"):
        importance = np.ravel(model.coef_)
    elif hasattr(model, "feature_importances_"):
        importance = np.ravel(model.feature_importances_)
    else:
        raise ValueError("Model does not expose coef_ or feature_importances_.")

    frame = pd.DataFrame({"feature": feature_names, "importance": importance})
    frame["abs_importance"] = frame["importance"].abs()
    frame = frame.sort_values("abs_importance", ascending=False).reset_index(drop=True)
    return frame


def safe_slug(text: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in text).strip("_")


def predict_age_at_diagnosis(
    standardized_rows: pd.DataFrame,
    model_path: Path,
    include_sex: bool = False,
    include_diagnosis: bool = False,
) -> np.ndarray:
    """Predict age-at-diagnosis from standardized patient rows."""
    pipeline: Pipeline = joblib.load(model_path)
    feature_frame = build_feature_frame(
        standardized_rows, include_sex=include_sex, include_diagnosis=include_diagnosis
    )
    return pipeline.predict(feature_frame)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("analysis/TRAIN_train_standardized.xlsx"),
        help="Path to standardized TRAIN file (xlsx/csv).",
    )
    parser.add_argument(
        "--phenotype",
        type=str,
        default="Diabetes Mellitus",
        help="Target phenotype name prefix (e.g., 'Diabetes Mellitus').",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/model_outputs"),
        help="Directory to write model artifacts and reports.",
    )
    parser.add_argument(
        "--n-splits",
        type=int,
        default=5,
        help="Cross-validation folds.",
    )
    parser.add_argument(
        "--include-sex",
        action="store_true",
        help="Include Sex_clean as an additional categorical covariate.",
    )
    parser.add_argument(
        "--include-diagnosis",
        action="store_true",
        help="Include Diagnosis_clean as an additional categorical covariate.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = load_standardized_data(args.data)

    age_type_col = f"{args.phenotype}__age_type"
    age_col = f"{args.phenotype}__age_dx_years"
    for col in [age_type_col, age_col]:
        if col not in df.columns:
            raise KeyError(f"Missing required phenotype column: {col}")

    mask = (df[age_type_col] == "exact_numeric") & df[age_col].notna()
    analysis_df = df.loc[mask].copy()
    analysis_df[age_col] = pd.to_numeric(analysis_df[age_col], errors="coerce")
    analysis_df = analysis_df.loc[analysis_df[age_col].notna()].copy()

    if analysis_df.empty:
        raise ValueError("No exact_numeric diagnosis ages found for selected phenotype.")

    X = build_feature_frame(
        analysis_df, include_sex=args.include_sex, include_diagnosis=args.include_diagnosis
    )
    y = analysis_df[age_col].astype(float)
    groups = build_genotype_signature(analysis_df)
    categorical_cols, numeric_cols = split_columns(X)
    preprocessor = make_preprocessor(categorical_cols=categorical_cols, numeric_cols=numeric_cols)

    model_specs: Dict[str, object] = {
        "elastic_net": ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=20000),
        "gradient_boosting": GradientBoostingRegressor(
            random_state=RANDOM_STATE,
            n_estimators=350,
            learning_rate=0.03,
            max_depth=2,
            min_samples_leaf=3,
            subsample=0.9,
        ),
    }

    cv_splits, cv_strategy = make_cv_splits(X, y, groups=groups, n_splits=args.n_splits)

    fold_metrics_frames = []
    for model_name, estimator in model_specs.items():
        fold_metrics = evaluate_model_cv(
            model_name=model_name,
            estimator=estimator,
            preprocessor=preprocessor,
            X=X,
            y=y,
            cv_splits=cv_splits,
        )
        fold_metrics_frames.append(fold_metrics)

    cv_fold_metrics = pd.concat(fold_metrics_frames, ignore_index=True)
    cv_summary = (
        cv_fold_metrics.groupby("model", as_index=False)
        .agg(
            mae_mean=("mae", "mean"),
            mae_std=("mae", "std"),
            rmse_mean=("rmse", "mean"),
            rmse_std=("rmse", "std"),
        )
        .sort_values(["mae_mean", "rmse_mean"], ascending=[True, True])
        .reset_index(drop=True)
    )

    phenotype_slug = safe_slug(args.phenotype)
    run_output_dir = args.output_dir / phenotype_slug
    run_output_dir.mkdir(parents=True, exist_ok=True)

    cv_fold_metrics.to_csv(run_output_dir / "cv_fold_metrics.csv", index=False)
    cv_summary.to_csv(run_output_dir / "cv_summary_metrics.csv", index=False)

    analysis_snapshot = pd.DataFrame(index=analysis_df.index)
    analysis_snapshot["Patient ID"] = analysis_df.get("Patient ID")
    analysis_snapshot["target_age_dx_years"] = y
    analysis_snapshot["genotype_signature"] = groups
    analysis_snapshot.to_csv(run_output_dir / "analysis_dataset_snapshot.csv", index=False)

    trained_model_paths = {}
    feature_importance_paths = {}
    for model_name, estimator in model_specs.items():
        pipeline = Pipeline(
            steps=[
                ("preprocessor", clone(preprocessor)),
                ("model", clone(estimator)),
            ]
        )
        pipeline.fit(X, y)

        model_path = run_output_dir / f"{model_name}_pipeline.joblib"
        joblib.dump(pipeline, model_path)
        trained_model_paths[model_name] = str(model_path)

        importance = feature_importance_frame(pipeline)
        importance_path = run_output_dir / f"{model_name}_feature_importance.csv"
        importance.to_csv(importance_path, index=False)
        feature_importance_paths[model_name] = str(importance_path)

    metadata = {
        "data_path": str(args.data),
        "phenotype": args.phenotype,
        "target_columns": {"age_type_col": age_type_col, "age_col": age_col},
        "n_samples": int(len(analysis_df)),
        "n_features": int(X.shape[1]),
        "cv_strategy": cv_strategy,
        "n_unique_genotype_signatures": int(groups.nunique(dropna=False)),
        "include_sex": bool(args.include_sex),
        "include_diagnosis": bool(args.include_diagnosis),
        "model_paths": trained_model_paths,
        "feature_importance_paths": feature_importance_paths,
        "assumptions": {
            "ages_are_fractional_years": True,
            "used_only_exact_numeric_target_rows": True,
            "excluded_present_no_age_absent_unknown": True,
        },
    }
    with open(run_output_dir / "run_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    best_row = cv_summary.iloc[0]
    print(f"Phenotype: {args.phenotype}")
    print(f"Training rows (exact numeric only): {len(analysis_df)}")
    print(f"CV strategy: {cv_strategy}")
    print("\nCV summary metrics (lower is better):")
    print(cv_summary.to_string(index=False))
    print(
        f"\nBest baseline by mean MAE: {best_row['model']} "
        f"(MAE={best_row['mae_mean']:.3f}, RMSE={best_row['rmse_mean']:.3f})"
    )
    print(f"Artifacts saved under: {run_output_dir}")


if __name__ == "__main__":
    main()
