#!/usr/bin/env python3
"""Generate model-insight visualizations for age-at-diagnosis baselines."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.base import clone
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import Pipeline

from train_age_dx_model import (
    RANDOM_STATE,
    build_feature_frame,
    build_genotype_signature,
    feature_importance_frame,
    load_standardized_data,
    make_cv_splits,
    make_preprocessor,
    safe_slug,
    split_columns,
)


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
        default="Hearing Loss",
        help="Target phenotype name prefix.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/model_outputs"),
        help="Base output directory used by training artifacts.",
    )
    parser.add_argument(
        "--n-splits",
        type=int,
        default=5,
        help="Cross-validation folds.",
    )
    parser.add_argument(
        "--top-n-features",
        type=int,
        default=15,
        help="Number of top features to show per model in importance plots.",
    )
    parser.add_argument(
        "--include-sex",
        action="store_true",
        help="Include Sex_clean as covariate in feature frame.",
    )
    parser.add_argument(
        "--include-diagnosis",
        action="store_true",
        help="Include Diagnosis_clean as covariate in feature frame.",
    )
    return parser.parse_args()


def ensure_theme() -> None:
    sns.set_theme(context="talk", style="whitegrid")
    plt.rcParams.update(
        {
            "axes.titlesize": 16,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 11,
            "figure.titlesize": 18,
            "axes.titleweight": "bold",
            "axes.edgecolor": "#1f2937",
            "axes.linewidth": 1.0,
            "savefig.dpi": 180,
            "savefig.bbox": "tight",
        }
    )


def model_specs() -> Dict[str, object]:
    return {
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


def compute_oof_predictions(
    X: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    preprocessor,
    models: Dict[str, object],
    cv_splits,
) -> pd.DataFrame:
    rows: List[dict] = []
    for model_name, estimator in models.items():
        oof_pred = pd.Series(index=X.index, dtype=float)
        fold_ids = pd.Series(index=X.index, dtype="Int64")

        for fold_idx, (train_idx, test_idx) in enumerate(cv_splits, start=1):
            pipeline = Pipeline(
                steps=[
                    ("preprocessor", clone(preprocessor)),
                    ("model", clone(estimator)),
                ]
            )
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train = y.iloc[train_idx]
            pipeline.fit(X_train, y_train)

            fold_pred = pipeline.predict(X_test)
            oof_pred.iloc[test_idx] = fold_pred
            fold_ids.iloc[test_idx] = fold_idx

        for idx in X.index:
            rows.append(
                {
                    "row_index": int(idx),
                    "model": model_name,
                    "fold": int(fold_ids.loc[idx]),
                    "y_true": float(y.loc[idx]),
                    "y_pred": float(oof_pred.loc[idx]),
                    "residual": float(y.loc[idx] - oof_pred.loc[idx]),
                    "abs_error": float(abs(y.loc[idx] - oof_pred.loc[idx])),
                    "genotype_signature": str(groups.loc[idx]),
                }
            )

    return pd.DataFrame(rows)


def fold_metrics_from_oof(oof_df: pd.DataFrame) -> pd.DataFrame:
    metrics = []
    for (model_name, fold_id), chunk in oof_df.groupby(["model", "fold"]):
        mae = mean_absolute_error(chunk["y_true"], chunk["y_pred"])
        rmse = float(np.sqrt(mean_squared_error(chunk["y_true"], chunk["y_pred"])))
        metrics.append(
            {
                "model": model_name,
                "fold": int(fold_id),
                "n_test": int(len(chunk)),
                "mae": float(mae),
                "rmse": rmse,
            }
        )
    return pd.DataFrame(metrics).sort_values(["model", "fold"]).reset_index(drop=True)


def summary_metrics(fold_df: pd.DataFrame) -> pd.DataFrame:
    return (
        fold_df.groupby("model", as_index=False)
        .agg(
            mae_mean=("mae", "mean"),
            mae_std=("mae", "std"),
            rmse_mean=("rmse", "mean"),
            rmse_std=("rmse", "std"),
        )
        .sort_values(["mae_mean", "rmse_mean"], ascending=[True, True])
        .reset_index(drop=True)
    )


def plot_target_distribution(y: pd.Series, phenotype: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(y, kde=True, bins=14, color="#0f766e", edgecolor="white", alpha=0.9, ax=ax)
    ax.set_title(f"{phenotype}: Diagnosis Age Distribution (exact numeric)")
    ax.set_xlabel("Age at diagnosis (years)")
    ax.set_ylabel("Count")
    fig.savefig(out_path)
    plt.close(fig)


def plot_cv_metrics(summary_df: pd.DataFrame, phenotype: str, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = {"gradient_boosting": "#0f766e", "elastic_net": "#d97706"}

    axes[0].bar(
        summary_df["model"],
        summary_df["mae_mean"],
        yerr=summary_df["mae_std"].fillna(0),
        color=[colors.get(m, "#374151") for m in summary_df["model"]],
        capsize=8,
        alpha=0.92,
    )
    axes[0].set_title("MAE across CV folds")
    axes[0].set_ylabel("MAE (years)")
    axes[0].set_xlabel("Model")
    axes[0].tick_params(axis="x", rotation=12)

    axes[1].bar(
        summary_df["model"],
        summary_df["rmse_mean"],
        yerr=summary_df["rmse_std"].fillna(0),
        color=[colors.get(m, "#374151") for m in summary_df["model"]],
        capsize=8,
        alpha=0.92,
    )
    axes[1].set_title("RMSE across CV folds")
    axes[1].set_ylabel("RMSE (years)")
    axes[1].set_xlabel("Model")
    axes[1].tick_params(axis="x", rotation=12)

    fig.suptitle(f"{phenotype}: Baseline CV Performance (mean ± SD)")
    fig.savefig(out_path)
    plt.close(fig)


def plot_oof_scatter(oof_df: pd.DataFrame, phenotype: str, out_path: Path) -> None:
    models = sorted(oof_df["model"].unique().tolist())
    fig, axes = plt.subplots(1, len(models), figsize=(7 * len(models), 6), sharex=True, sharey=True)
    if len(models) == 1:
        axes = [axes]

    x_min = float(np.nanmin(oof_df["y_true"]))
    x_max = float(np.nanmax(oof_df["y_true"]))
    padding = max(1.0, 0.06 * (x_max - x_min))
    lo, hi = x_min - padding, x_max + padding
    palette = {"gradient_boosting": "#0f766e", "elastic_net": "#d97706"}

    for ax, model_name in zip(axes, models):
        chunk = oof_df.loc[oof_df["model"] == model_name].copy()
        mae = mean_absolute_error(chunk["y_true"], chunk["y_pred"])
        rmse = float(np.sqrt(mean_squared_error(chunk["y_true"], chunk["y_pred"])))
        r2 = float(r2_score(chunk["y_true"], chunk["y_pred"]))

        sns.scatterplot(
            data=chunk,
            x="y_true",
            y="y_pred",
            alpha=0.85,
            s=70,
            color=palette.get(model_name, "#374151"),
            edgecolor="white",
            linewidth=0.5,
            ax=ax,
        )
        ax.plot([lo, hi], [lo, hi], linestyle="--", color="#111827", linewidth=1.3)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_title(model_name.replace("_", " ").title())
        ax.set_xlabel("Observed age at diagnosis")
        ax.set_ylabel("Out-of-fold predicted age")
        ax.text(
            0.03,
            0.97,
            f"MAE: {mae:.2f}\nRMSE: {rmse:.2f}\nR²: {r2:.2f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=11,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "#cbd5e1", "boxstyle": "round,pad=0.3"},
        )

    fig.suptitle(f"{phenotype}: Observed vs Predicted (Out-of-Fold)")
    fig.savefig(out_path)
    plt.close(fig)


def plot_residuals(oof_df: pd.DataFrame, phenotype: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    palette = {"gradient_boosting": "#0f766e", "elastic_net": "#d97706"}
    sns.violinplot(
        data=oof_df,
        x="model",
        y="residual",
        hue="model",
        palette=palette,
        legend=False,
        inner="quartile",
        cut=0,
        linewidth=1.1,
        ax=ax,
    )
    ax.axhline(0, linestyle="--", color="#111827", linewidth=1.2)
    ax.set_title(f"{phenotype}: Residual Distribution by Model")
    ax.set_xlabel("Model")
    ax.set_ylabel("Residual (observed - predicted years)")
    ax.tick_params(axis="x", rotation=12)
    fig.savefig(out_path)
    plt.close(fig)


def plot_top_features(
    importances: Dict[str, pd.DataFrame], phenotype: str, out_path: Path, top_n: int
) -> None:
    model_names = sorted(importances.keys())
    fig, axes = plt.subplots(1, len(model_names), figsize=(8 * len(model_names), 7), sharex=False)
    if len(model_names) == 1:
        axes = [axes]
    palette = {"gradient_boosting": "#0f766e", "elastic_net": "#d97706"}

    for ax, model_name in zip(axes, model_names):
        top = importances[model_name].head(top_n).copy()
        top = top.iloc[::-1]
        ax.barh(
            top["feature"],
            top["abs_importance"],
            color=palette.get(model_name, "#374151"),
            alpha=0.92,
        )
        ax.set_title(f"{model_name.replace('_', ' ').title()}: Top {top_n} Features")
        ax.set_xlabel("Absolute importance")
        ax.set_ylabel("Encoded feature")
        ax.grid(axis="x", linestyle=":", alpha=0.5)

    fig.suptitle(f"{phenotype}: Feature Importance")
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    ensure_theme()

    phenotype_slug = safe_slug(args.phenotype)
    run_output_dir = args.output_dir / phenotype_slug
    figures_dir = run_output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

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
        raise ValueError("No exact_numeric rows available for this phenotype.")

    X = build_feature_frame(
        analysis_df, include_sex=args.include_sex, include_diagnosis=args.include_diagnosis
    )
    y = analysis_df[age_col].astype(float)
    groups = build_genotype_signature(analysis_df)
    categorical_cols, numeric_cols = split_columns(X)
    preprocessor = make_preprocessor(categorical_cols=categorical_cols, numeric_cols=numeric_cols)
    cv_splits, cv_strategy = make_cv_splits(X, y, groups=groups, n_splits=args.n_splits)

    specs = model_specs()
    oof_df = compute_oof_predictions(
        X=X,
        y=y,
        groups=groups,
        preprocessor=preprocessor,
        models=specs,
        cv_splits=cv_splits,
    )
    fold_df = fold_metrics_from_oof(oof_df)
    summary_df = summary_metrics(fold_df)

    oof_df.to_csv(run_output_dir / "cv_oof_predictions.csv", index=False)
    fold_df.to_csv(run_output_dir / "cv_fold_metrics_recomputed.csv", index=False)
    summary_df.to_csv(run_output_dir / "cv_summary_metrics_recomputed.csv", index=False)

    importances: Dict[str, pd.DataFrame] = {}
    for model_name, estimator in specs.items():
        pipeline = Pipeline(
            steps=[
                ("preprocessor", clone(preprocessor)),
                ("model", clone(estimator)),
            ]
        )
        pipeline.fit(X, y)
        importance = feature_importance_frame(pipeline)
        importance.to_csv(run_output_dir / f"{model_name}_feature_importance_recomputed.csv", index=False)
        importances[model_name] = importance

    plot_target_distribution(y, args.phenotype, figures_dir / "target_distribution.png")
    plot_cv_metrics(summary_df, args.phenotype, figures_dir / "cv_metrics.png")
    plot_oof_scatter(oof_df, args.phenotype, figures_dir / "oof_observed_vs_predicted.png")
    plot_residuals(oof_df, args.phenotype, figures_dir / "residual_distribution.png")
    plot_top_features(
        importances=importances,
        phenotype=args.phenotype,
        out_path=figures_dir / "feature_importance_top.png",
        top_n=args.top_n_features,
    )

    insights = []
    for model_name in summary_df["model"].tolist():
        chunk = oof_df.loc[oof_df["model"] == model_name]
        signed_bias = float(chunk["residual"].mean())
        insights.append(
            {
                "model": model_name,
                "mae": float(mean_absolute_error(chunk["y_true"], chunk["y_pred"])),
                "rmse": float(np.sqrt(mean_squared_error(chunk["y_true"], chunk["y_pred"]))),
                "r2": float(r2_score(chunk["y_true"], chunk["y_pred"])),
                "mean_signed_residual": signed_bias,
            }
        )
    insights_df = pd.DataFrame(insights).sort_values("mae").reset_index(drop=True)
    insights_df.to_csv(run_output_dir / "model_insights_summary.csv", index=False)
    summary_table = insights_df.round(3).to_string(index=False)

    summary_lines = [
        f"# {args.phenotype} model insights",
        "",
        f"- Samples used: {len(analysis_df)} (exact_numeric only)",
        f"- CV strategy: {cv_strategy}",
        "",
        "## Model summary",
        "```text",
        summary_table,
        "```",
        "",
        "## Generated figures",
        "- figures/target_distribution.png",
        "- figures/cv_metrics.png",
        "- figures/oof_observed_vs_predicted.png",
        "- figures/residual_distribution.png",
        "- figures/feature_importance_top.png",
    ]
    (run_output_dir / "MODEL_INSIGHTS.md").write_text("\n".join(summary_lines), encoding="utf-8")

    print(f"Phenotype: {args.phenotype}")
    print(f"Rows used: {len(analysis_df)}")
    print(f"CV strategy: {cv_strategy}")
    print("Saved:")
    print(f"- {run_output_dir / 'cv_oof_predictions.csv'}")
    print(f"- {run_output_dir / 'model_insights_summary.csv'}")
    print(f"- {run_output_dir / 'MODEL_INSIGHTS.md'}")
    print(f"- {figures_dir}")


if __name__ == "__main__":
    main()
