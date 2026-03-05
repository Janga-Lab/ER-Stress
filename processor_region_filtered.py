#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from typing import Any, Dict, List, Optional, Tuple

from config import (
    CURATED_DATA_DIR,
    DEFAULT_CLASS_WEIGHTS,
    DEFAULT_CONFIDENCE_HIGH_THRESHOLD,
    DEFAULT_CONFIDENCE_LOW_THRESHOLD,
    DEFAULT_INHERITANCE_MAP,
    DEFAULT_SCORE_BIN_THRESHOLDS,
    DEFAULT_UNMAPPED_WARNING_THRESHOLD,
    DEFAULT_VARIANT_PREDICTION_THRESHOLD,
    PROCESSED_DATA_DIR,
    PipelineConfig,
    ensure_data_directories,
)
from processor import (
    build_decision_tables,
    configure_runtime_env,
    evaluate_variant_prediction,
    process_variants,
    read_csv_rows,
    score_bin,
    write_csv,
)

DEFAULT_INPUT_CURATED_VARIANTS = CURATED_DATA_DIR / "curated_variants_clinical_association.csv"
DEFAULT_OUT_ROOT = CURATED_DATA_DIR / "new"
DEFAULT_SPLICEAI_WIDE_CSV = PROCESSED_DATA_DIR / "spliceai_other_scores_wide.csv"

DEFAULT_KEEP_REGIONS = [
    "tm",
    "luminal",
    "kinase_domain",
    "cytosolic",
    "disulfide_cysteine",
    "signal_peptide",
    "b_chain",
]

SCORE_BIN_ORDER = ["low", "mild", "moderate", "high", "very_high"]

PREDICTOR_SCORE_COLUMNS = [
    "alphamissense_score",
    "primateai3d_score",
    "revel_max_score",
    "sift_max_score",
    "promoterai_score",
]

PREDICTOR_POINTS_COLUMNS = [
    "alphamissense_points",
    "cadd_points",
    "phylop_points",
    "polyphen_max_points",
    "revel_max_points",
    "sift_max_points",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter curated clinical-association variants by useful region buckets, split by gene, "
            "and run per-gene scoring/plots into data/curated_data/new."
        )
    )
    parser.add_argument(
        "--curated-variants",
        default=str(DEFAULT_INPUT_CURATED_VARIANTS),
        help="Input curated variants CSV.",
    )
    parser.add_argument(
        "--out-root",
        default=str(DEFAULT_OUT_ROOT),
        help="Output root directory (per-gene folders will be created here).",
    )
    parser.add_argument(
        "--keep-regions",
        default=",".join(DEFAULT_KEEP_REGIONS),
        help="Comma-separated region_bucket values to keep.",
    )
    parser.add_argument(
        "--spliceai-wide",
        default=str(DEFAULT_SPLICEAI_WIDE_CSV),
        help="SpliceAI other-scores wide CSV to merge as secondary predictor evidence.",
    )
    parser.add_argument(
        "--score-mode",
        choices=["base", "augmented"],
        default="augmented",
        help="Use base score only or score augmented with predictor bonus.",
    )
    parser.add_argument("--variant-pred-threshold", type=float, default=DEFAULT_VARIANT_PREDICTION_THRESHOLD)
    return parser.parse_args()


def normalize_keep_regions(text: str) -> List[str]:
    return [x.strip().lower() for x in text.split(",") if x.strip()]


def grouped_rows_by_gene(rows: List[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
    by_gene: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        gene = (row.get("gene_symbol") or "").strip().upper()
        if not gene:
            continue
        row["gene_symbol"] = gene
        by_gene[gene].append(row)
    return dict(sorted(by_gene.items()))


def parse_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        if isinstance(value, float) and math.isnan(value):
            return None
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None
    return None


def parse_int(value: Any) -> Optional[int]:
    f = parse_float(value)
    return int(f) if f is not None else None


def _predictor_vote_summary(row: Dict[str, Any]) -> Dict[str, Any]:
    alpha = parse_float(row.get("alphamissense_score"))
    primate = parse_float(row.get("primateai3d_score"))
    revel = parse_float(row.get("revel_max_score"))
    sift = parse_float(row.get("sift_max_score"))

    damaging_votes = 0
    benign_votes = 0

    if alpha is not None:
        if alpha >= 0.8:
            damaging_votes += 1
        elif alpha <= 0.2:
            benign_votes += 1

    if primate is not None:
        if primate >= 0.8:
            damaging_votes += 1
        elif primate <= 0.2:
            benign_votes += 1

    if revel is not None:
        if revel >= 0.7:
            damaging_votes += 1
        elif revel <= 0.3:
            benign_votes += 1

    if sift is not None:
        if sift <= 0.05:
            damaging_votes += 1
        elif sift >= 0.2:
            benign_votes += 1

    total_votes = damaging_votes + benign_votes
    if total_votes == 0:
        signal = "missing"
        bonus = 0.0
    elif damaging_votes >= 2 and benign_votes == 0:
        signal = "strong_damaging"
        bonus = 1.0
    elif damaging_votes >= 1 and benign_votes == 0:
        signal = "moderate_damaging"
        bonus = 0.5
    elif benign_votes >= 2 and damaging_votes == 0:
        signal = "supportive_benign"
        bonus = -0.5
    elif benign_votes >= 1 and damaging_votes == 0:
        signal = "weak_benign"
        bonus = -0.25
    else:
        signal = "mixed"
        bonus = 0.0

    return {
        "predictor_votes_damaging": damaging_votes,
        "predictor_votes_benign": benign_votes,
        "predictor_votes_total": total_votes,
        "predictor_signal": signal,
        "predictor_bonus": round(bonus, 4),
    }


def load_predictor_map(spliceai_wide_csv: Path) -> Dict[Tuple[str, str], Dict[str, Any]]:
    if not spliceai_wide_csv.exists():
        return {}

    rows = read_csv_rows(spliceai_wide_csv)
    best: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for row in rows:
        gene = (row.get("gene_symbol") or "").strip().upper()
        var_id = str(row.get("clinvar_variation_id") or "").strip()
        if not gene or not var_id:
            continue

        status = (row.get("status") or "").strip().lower()
        n_scores = sum(1 for c in PREDICTOR_SCORE_COLUMNS if parse_float(row.get(c)) is not None)
        rank = (1 if status == "ok" else 0, n_scores)

        key = (gene, var_id)
        prev = best.get(key)
        if prev is not None:
            prev_rank = prev.get("_rank") or (0, 0)
            if rank <= prev_rank:
                continue

        rec: Dict[str, Any] = {
            "_rank": rank,
            "predictor_status": status or "",
            "predictor_missing_predictors": (row.get("missing_predictors") or "").strip(),
            "predictor_error": (row.get("error") or "").strip(),
            "predictor_table_rows": parse_int(row.get("table_rows")) or 0,
            "predictor_lookup_url": (row.get("lookup_url") or "").strip(),
            "predictor_n_scores": n_scores,
        }
        for col in PREDICTOR_SCORE_COLUMNS + PREDICTOR_POINTS_COLUMNS:
            rec[col] = parse_float(row.get(col))
        rec.update(_predictor_vote_summary(rec))
        best[key] = rec

    for rec in best.values():
        rec.pop("_rank", None)
    return best


def augment_variant_rows_with_predictors(
    variant_rows: List[Dict[str, Any]],
    predictor_map: Dict[Tuple[str, str], Dict[str, Any]],
    cfg: PipelineConfig,
    score_mode: str,
) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    has_map = bool(predictor_map)
    augmented_rows = 0
    rows_with_any_predictor = 0

    for row in variant_rows:
        gene = (row.get("gene_symbol") or "").strip().upper()
        var_id = str(row.get("clinvar_variation_id") or "").strip()
        key = (gene, var_id)
        pred = predictor_map.get(key)

        base_score = parse_float(row.get("allele_severity_points")) or 0.0
        base_bin = (row.get("score_bin") or "").strip().lower() or score_bin(base_score, cfg.score_bin_thresholds)
        row["allele_severity_points_base"] = round(base_score, 4)
        row["score_bin_base"] = base_bin

        if pred is None:
            row["predictor_status"] = "missing_key" if has_map else "not_loaded"
            row["predictor_missing_predictors"] = ""
            row["predictor_error"] = ""
            row["predictor_table_rows"] = ""
            row["predictor_lookup_url"] = ""
            row["predictor_n_scores"] = 0
            row["predictor_votes_damaging"] = 0
            row["predictor_votes_benign"] = 0
            row["predictor_votes_total"] = 0
            row["predictor_signal"] = "missing"
            row["predictor_bonus"] = 0.0
            for col in PREDICTOR_SCORE_COLUMNS + PREDICTOR_POINTS_COLUMNS:
                row[col] = ""
            bonus = 0.0
        else:
            for col in (
                [
                    "predictor_status",
                    "predictor_missing_predictors",
                    "predictor_error",
                    "predictor_table_rows",
                    "predictor_lookup_url",
                    "predictor_n_scores",
                    "predictor_votes_damaging",
                    "predictor_votes_benign",
                    "predictor_votes_total",
                    "predictor_signal",
                    "predictor_bonus",
                ]
                + PREDICTOR_SCORE_COLUMNS
                + PREDICTOR_POINTS_COLUMNS
            ):
                v = pred.get(col)
                row[col] = "" if v is None else v
            bonus = parse_float(pred.get("predictor_bonus")) or 0.0
            if (parse_int(pred.get("predictor_n_scores")) or 0) > 0:
                rows_with_any_predictor += 1

        augmented_score = max(0.0, base_score + bonus)
        row["predictor_bonus_applied"] = round(bonus, 4)
        row["allele_severity_points_augmented"] = round(augmented_score, 4)
        row["score_bin_augmented"] = score_bin(augmented_score, cfg.score_bin_thresholds)

        if score_mode == "augmented":
            row["allele_severity_points"] = row["allele_severity_points_augmented"]
            row["score_bin"] = row["score_bin_augmented"]
            if abs(augmented_score - base_score) > 1e-9:
                augmented_rows += 1

    return variant_rows, {
        "rows_total": len(variant_rows),
        "rows_with_any_predictor": rows_with_any_predictor,
        "rows_with_adjusted_score": augmented_rows,
    }


def rank_values(values: List[float]) -> List[float]:
    indexed = sorted(enumerate(values), key=lambda x: x[1])
    ranks = [0.0] * len(values)
    i = 0
    rank = 1.0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = (rank + rank + (j - i)) / 2.0
        for k in range(i, j + 1):
            orig_idx = indexed[k][0]
            ranks[orig_idx] = avg_rank
        rank += (j - i + 1)
        i = j + 1
    return ranks


def pearson_corr(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    mx = mean(xs)
    my = mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = sum((x - mx) ** 2 for x in xs)
    den_y = sum((y - my) ** 2 for y in ys)
    den = math.sqrt(den_x * den_y)
    if den <= 0:
        return None
    return num / den


def spearman_corr(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    rx = rank_values(xs)
    ry = rank_values(ys)
    return pearson_corr(rx, ry)


def linear_slope(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    mx = mean(xs)
    my = mean(ys)
    den = sum((x - mx) ** 2 for x in xs)
    if den <= 0:
        return None
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return num / den


def _plot_info_panel(title: str, message: str, out_path: Path) -> Optional[Path]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.2, 4.5))
    ax.axis("off")
    ax.set_title(title, fontsize=12)
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        fontsize=11,
        color="#37474f",
        transform=ax.transAxes,
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close()
    return out_path


def summarize_variant_enrichment(
    variant_rows: List[Dict[str, Any]],
    score_col: str = "allele_severity_points",
    score_bin_col: str = "score_bin",
) -> Dict[str, Any]:
    labeled = [r for r in variant_rows if parse_int(r.get("pathogenic_binary")) in {0, 1}]
    by_bin_counts: Dict[str, int] = {b: 0 for b in SCORE_BIN_ORDER}
    by_bin_path: Dict[str, int] = {b: 0 for b in SCORE_BIN_ORDER}

    xs: List[float] = []
    ys: List[float] = []
    for row in labeled:
        b = (row.get(score_bin_col) or "low").lower()
        if b in by_bin_counts:
            by_bin_counts[b] += 1
            by_bin_path[b] += int(row.get("pathogenic_binary") or 0)
        score = parse_float(row.get(score_col))
        y = parse_int(row.get("pathogenic_binary"))
        if score is not None and y is not None:
            xs.append(score)
            ys.append(float(y))

    by_bin_rate: Dict[str, Optional[float]] = {}
    for b in SCORE_BIN_ORDER:
        n = by_bin_counts[b]
        by_bin_rate[b] = (by_bin_path[b] / n) if n else None

    split_method = "bin_split"
    low_n = by_bin_counts["low"] + by_bin_counts["mild"]
    high_n = by_bin_counts["high"] + by_bin_counts["very_high"]
    low_path = by_bin_path["low"] + by_bin_path["mild"]
    high_path = by_bin_path["high"] + by_bin_path["very_high"]

    low_rate = (low_path / low_n) if low_n else None
    high_rate = (high_path / high_n) if high_n else None

    # If strict low/high bins are sparse after filtering, use a score-median split.
    if (low_rate is None or high_rate is None) and len(xs) >= 6:
        split_method = "score_median_split"
        pairs = sorted(zip(xs, ys), key=lambda t: t[0])
        mid = len(pairs) // 2
        lower = pairs[:mid]
        upper = pairs[mid:]
        if lower and upper:
            low_n = len(lower)
            high_n = len(upper)
            low_path = int(sum(y for _, y in lower))
            high_path = int(sum(y for _, y in upper))
            low_rate = low_path / low_n
            high_rate = high_path / high_n

    rho = spearman_corr(xs, ys) if len(xs) >= 3 else None

    if len(labeled) < 10 or low_rate is None or high_rate is None:
        status = "insufficient_data"
        note = "Not enough labeled variants or missing low/high bin representation."
    elif high_rate > low_rate and (rho is None or rho > 0):
        status = "supports"
        note = "Higher-score bins show stronger pathogenic enrichment."
    else:
        status = "does_not_support"
        note = "Enrichment pattern is weak or non-monotonic."

    return {
        "status": status,
        "note": note,
        "n_total": len(variant_rows),
        "n_labeled": len(labeled),
        "low_bin_n": low_n,
        "high_bin_n": high_n,
        "low_bin_pathogenic_n": low_path,
        "high_bin_pathogenic_n": high_path,
        "low_bin_pathogenic_rate": low_rate,
        "high_bin_pathogenic_rate": high_rate,
        "spearman_score_vs_pathogenic": rho,
        "split_method": split_method,
        "by_bin_counts": by_bin_counts,
        "by_bin_pathogenic_counts": by_bin_path,
        "by_bin_rates": by_bin_rate,
    }


def plot_variant_enrichment(
    variant_rows: List[Dict[str, Any]],
    out_path: Path,
    score_col: str = "allele_severity_points",
    score_bin_col: str = "score_bin",
    title_suffix: str = "",
) -> Optional[Path]:
    summary = summarize_variant_enrichment(variant_rows, score_col=score_col, score_bin_col=score_bin_col)
    if int(summary["n_labeled"]) == 0:
        return _plot_info_panel(
            "Pathogenic Signal vs Score",
            "No benign/pathogenic labels are available for this gene, so enrichment cannot be evaluated.",
            out_path,
        )

    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.2, 5.6))

    bin_counts: Dict[str, int] = summary["by_bin_counts"]  # type: ignore[assignment]
    bin_paths: Dict[str, int] = summary["by_bin_pathogenic_counts"]  # type: ignore[assignment]
    bin_rates: Dict[str, Optional[float]] = summary["by_bin_rates"]  # type: ignore[assignment]

    bins = list(SCORE_BIN_ORDER)
    x = list(range(len(bins)))
    rates: List[float] = []
    colors: List[str] = []
    palette = {
        "low": "#90a4ae",
        "mild": "#4fc3f7",
        "moderate": "#29b6f6",
        "high": "#ef6c00",
        "very_high": "#d84315",
    }

    for b in bins:
        r = bin_rates.get(b)
        rates.append(float(r) if r is not None else 0.0)
        colors.append(palette.get(b, "#90a4ae"))

    bars = ax.bar(x, rates, color=colors, edgecolor="white", linewidth=1.0, width=0.66, alpha=0.94)
    ax.set_xticks(x, bins)
    ax.set_ylabel("Pathogenic-like fraction (0 to 1)")
    ax.set_ylim(0, 1.08)

    title = "Pathogenic Signal by Severity Bin"
    if title_suffix:
        title = f"{title} ({title_suffix})"
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.2)

    for i, (bar, b) in enumerate(zip(bars, bins)):
        n = int(bin_counts.get(b, 0))
        p = int(bin_paths.get(b, 0))
        r = bin_rates.get(b)
        if n == 0 or r is None:
            y = 0.03
            label = "n=0"
            color = "#455a64"
            va = "bottom"
            bbox = None
        else:
            r_float = float(r)
            pct = int(round(r_float * 100))
            label = f"{p}/{n}\n{pct}%"
            if r_float >= 0.18:
                # Put the label inside bars when possible for consistent readability.
                y = r_float / 2.0
                va = "center"
            else:
                y = min(1.02, r_float + 0.04)
                va = "bottom"
            color = "#102027"
            bbox = {"boxstyle": "round,pad=0.22", "facecolor": "white", "alpha": 0.78, "edgecolor": "none"}
        ax.text(i, y, label, ha="center", va=va, fontsize=8.5, color=color, fontweight=700, bbox=bbox)

    # Add a simple trend line across non-empty bins for visual directionality.
    trend_x = []
    trend_y = []
    for i, b in enumerate(bins):
        r = bin_rates.get(b)
        n = int(bin_counts.get(b, 0))
        if n > 0 and r is not None:
            trend_x.append(i)
            trend_y.append(float(r))
    if len(trend_x) >= 2:
        ax.plot(trend_x, trend_y, color="#37474f", linewidth=1.5, linestyle="--", marker="o", markersize=3.5, alpha=0.8)

    rho = summary["spearman_score_vs_pathogenic"]
    split_method = summary["split_method"]
    non_empty_bins = [b for b in bins if int(bin_counts.get(b, 0)) > 0]
    if non_empty_bins:
        bin_profile = ", ".join([f"{b}:{int(bin_counts.get(b, 0))}" for b in non_empty_bins])
        interpretation = (
            "Five-bin profile shown directly on bars (pathogenic-like count/total + percent). "
            f"Non-empty bins in this gene: {bin_profile}."
        )
    else:
        interpretation = "No score bins contain labeled variants in this subset."
    meta_line = (
        f"Status={summary['status']} | labeled={summary['n_labeled']} | split={split_method} | "
        f"spearman={'' if rho is None else f'{rho:.3f}'}"
    )
    fig.text(0.5, 0.048, interpretation, ha="center", va="center", fontsize=8.6, color="#455a64")
    fig.text(0.5, 0.023, meta_line, ha="center", va="center", fontsize=8.1, color="#607d8b")

    plt.tight_layout(rect=[0, 0.11, 1, 1])
    plt.savefig(out_path, dpi=180)
    plt.close()
    return out_path


def summarize_decision_clarity(decision_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    total = len(decision_rows)
    wet = sum(1 for r in decision_rows if parse_int(r.get("needs_wetlab_validation")) == 1)
    depri = total - wet
    tier_counts = Counter((r.get("decision_tier") or "").strip() for r in decision_rows)

    if total == 0:
        status = "insufficient_data"
        note = "No decision rows."
    elif wet > 0 and depri > 0:
        status = "supports"
        note = "Decision output distinguishes wet-lab priority vs deprioritized cases."
    else:
        status = "supports_with_limitations"
        note = "Decision output is explicit but one-sided for this filtered subset."

    return {
        "status": status,
        "note": note,
        "n_decision_rows": total,
        "n_wetlab_priority": wet,
        "n_deprioritized": depri,
        "tier_counts": dict(tier_counts),
    }


def plot_decision_clarity(decision_rows: List[Dict[str, Any]], out_path: Path) -> Optional[Path]:
    summary = summarize_decision_clarity(decision_rows)
    if int(summary["n_decision_rows"]) == 0:
        return _plot_info_panel("Decision Clarity", "No decision rows available.", out_path)

    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.8, 4.9))

    tiers = ["Tier_A", "Tier_B", "Tier_C", "Tier_D"]
    tier_counts = [int(summary["tier_counts"].get(t, 0)) for t in tiers]  # type: ignore[index]
    wet = int(summary["n_wetlab_priority"])
    dep = int(summary["n_deprioritized"])
    total = int(summary["n_decision_rows"]) or 1

    categories = ["Wet-lab priority", "Deprioritized"]
    values = [wet, dep]
    bar_colors = ["#ef6c00", "#1e88e5"]
    bars = ax.barh(categories, values, color=bar_colors, alpha=0.92, edgecolor="white", linewidth=1.0)
    ax.set_xlim(0, max(total * 1.12, 1))
    ax.set_xlabel("Number of cases")
    ax.set_title("What should we do next with these cases?")
    ax.grid(axis="x", alpha=0.2)

    for bar, count in zip(bars, values):
        pct = (100.0 * count / total) if total else 0.0
        label = f"{count} ({pct:.1f}%)"
        width = bar.get_width()
        y = bar.get_y() + bar.get_height() / 2.0
        if width > total * 0.15:
            x = width * 0.5
            color = "white"
            ha = "center"
        else:
            x = width + total * 0.015
            color = "#263238"
            ha = "left"
        ax.text(x, y, label, ha=ha, va="center", fontsize=9, color=color, fontweight="bold")

    status_line = f"Status={summary['status']} | {summary['note']}"
    tier_line = f"Tier counts: A={tier_counts[0]}, B={tier_counts[1]}, C={tier_counts[2]}, D={tier_counts[3]}"
    fig.text(0.5, 0.045, status_line, ha="center", va="center", fontsize=8.5, color="#455a64")
    fig.text(0.5, 0.02, tier_line, ha="center", va="center", fontsize=8, color="#607d8b")

    plt.tight_layout(rect=[0, 0.10, 1, 1])
    plt.savefig(out_path, dpi=180)
    plt.close()
    return out_path


def plot_success_overview(success_rows: List[Dict[str, Any]], out_path: Path) -> Optional[Path]:
    if not success_rows:
        return _plot_info_panel("Success Overview", "No genes processed.", out_path)
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    status_order = ["supports", "supports_with_limitations", "insufficient_data", "does_not_support"]
    status_to_color = {
        "supports": "#2e7d32",
        "supports_with_limitations": "#f9a825",
        "insufficient_data": "#90a4ae",
        "does_not_support": "#c62828",
    }

    genes = [str(r["gene_symbol"]) for r in success_rows]
    criteria_cols = [
        "criterion_1_variant_enrichment_status",
        "criterion_3_decision_clarity_status",
    ]
    criteria_labels = [
        "Variant Enrichment",
        "Decision Clarity",
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.4, max(3.6, 0.7 * len(genes) + 1.8)))
    ax.set_xlim(0, len(criteria_cols))
    ax.set_ylim(0, len(genes))
    ax.invert_yaxis()
    ax.set_xticks([i + 0.5 for i in range(len(criteria_cols))], criteria_labels)
    ax.set_yticks([i + 0.5 for i in range(len(genes))], genes)
    ax.set_title("Per-Gene Success Criteria Overview")

    for i, row in enumerate(success_rows):
        for j, col in enumerate(criteria_cols):
            status = str(row.get(col, "insufficient_data"))
            color = status_to_color.get(status, "#90a4ae")
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=color, edgecolor="white", linewidth=1.5))
            short = {
                "supports": "S",
                "supports_with_limitations": "SL",
                "insufficient_data": "I",
                "does_not_support": "N",
            }.get(status, "?")
            ax.text(j + 0.5, i + 0.5, short, ha="center", va="center", color="white", fontsize=9, fontweight="bold")

    from matplotlib.lines import Line2D

    legend_items = [
        Line2D([0], [0], marker="s", color="w", label=s.replace("_", " "), markerfacecolor=status_to_color[s], markersize=9)
        for s in status_order
    ]
    ax.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=False, fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close()
    return out_path


def main() -> int:
    args = parse_args()
    ensure_data_directories()

    curated_path = Path(args.curated_variants)
    if not curated_path.exists():
        raise SystemExit(f"Curated variants file not found: {curated_path}")

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    configure_runtime_env(out_root)
    keep_regions = normalize_keep_regions(args.keep_regions)
    keep_region_set = set(keep_regions)
    if not keep_region_set:
        raise SystemExit("No keep regions provided.")

    cfg = PipelineConfig(
        class_weights=dict(DEFAULT_CLASS_WEIGHTS),
        inheritance_map=dict(DEFAULT_INHERITANCE_MAP),
        score_bin_thresholds=dict(DEFAULT_SCORE_BIN_THRESHOLDS),
        variant_prediction_threshold=float(args.variant_pred_threshold),
        confidence_low_threshold=DEFAULT_CONFIDENCE_LOW_THRESHOLD,
        confidence_high_threshold=DEFAULT_CONFIDENCE_HIGH_THRESHOLD,
        unmapped_warning_threshold=DEFAULT_UNMAPPED_WARNING_THRESHOLD,
    )
    spliceai_wide_path = Path(args.spliceai_wide)
    predictor_map = load_predictor_map(spliceai_wide_path)
    use_augmented = args.score_mode == "augmented" and bool(predictor_map)
    if args.score_mode == "augmented" and not predictor_map:
        print(
            f"[WARN] score-mode=augmented requested but predictor file not available/usable: {spliceai_wide_path}. "
            "Falling back to base score."
        )

    curated_rows_all = read_csv_rows(curated_path)
    filtered_rows = [
        row
        for row in curated_rows_all
        if (row.get("region_bucket") or "").strip().lower() in keep_region_set
    ]
    by_gene = grouped_rows_by_gene(filtered_rows)

    # Write top-level filtered snapshot and summary.
    write_csv(out_root / "curated_variants_region_filtered.csv", filtered_rows)

    original_by_gene = Counter((row.get("gene_symbol") or "").strip().upper() for row in curated_rows_all)
    filtered_by_gene = Counter((row.get("gene_symbol") or "").strip().upper() for row in filtered_rows)

    summary_rows: List[Dict[str, object]] = []
    for gene in sorted(set(original_by_gene) | set(filtered_by_gene)):
        total = int(original_by_gene.get(gene, 0))
        kept = int(filtered_by_gene.get(gene, 0))
        dropped = total - kept
        keep_frac = (kept / total) if total else 0.0
        summary_rows.append(
            {
                "gene_symbol": gene,
                "input_rows": total,
                "kept_rows": kept,
                "dropped_rows": dropped,
                "keep_fraction": round(keep_frac, 4),
            }
        )
    write_csv(out_root / "region_filter_summary.csv", summary_rows)

    run_summaries: List[Dict[str, object]] = []
    success_rows: List[Dict[str, object]] = []

    for gene, gene_curated_rows in by_gene.items():
        gene_dir = out_root / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        for path in gene_dir.iterdir():
            if path.is_file():
                path.unlink()

        gene_curated_csv = gene_dir / f"curated_variants_{gene}.csv"
        write_csv(gene_curated_csv, gene_curated_rows)

        variant_rows = process_variants(gene_curated_rows, cfg)
        variant_rows, predictor_stats = augment_variant_rows_with_predictors(
            variant_rows=variant_rows,
            predictor_map=predictor_map,
            cfg=cfg,
            score_mode="augmented" if use_augmented else "base",
        )
        prediction_summary = evaluate_variant_prediction(variant_rows, cfg.variant_prediction_threshold)
        decision_rows, priority_ranked = build_decision_tables(variant_rows, prediction_summary, cfg)

        variant_processed_csv = gene_dir / "variant_processed.csv"
        decision_csv = gene_dir / "decision_tiers.csv"
        priority_csv = gene_dir / "priority_ranked.csv"

        write_csv(variant_processed_csv, variant_rows)
        write_csv(decision_csv, decision_rows)
        write_csv(priority_csv, priority_ranked)

        # Success-criteria evaluators and dedicated figures
        base_variant_signal = summarize_variant_enrichment(
            variant_rows,
            score_col="allele_severity_points_base",
            score_bin_col="score_bin_base",
        )
        if predictor_map:
            augmented_variant_signal = summarize_variant_enrichment(
                variant_rows,
                score_col="allele_severity_points_augmented",
                score_bin_col="score_bin_augmented",
            )
        else:
            augmented_variant_signal = base_variant_signal

        variant_signal = augmented_variant_signal if use_augmented else base_variant_signal
        decision_signal = summarize_decision_clarity(decision_rows)

        plot_variant_enrichment(
            variant_rows,
            gene_dir / "criterion1_variant_enrichment_base.png",
            score_col="allele_severity_points_base",
            score_bin_col="score_bin_base",
            title_suffix="base score",
        )
        plot_variant_enrichment(
            variant_rows,
            gene_dir / "criterion1_variant_enrichment_augmented.png",
            score_col="allele_severity_points_augmented",
            score_bin_col="score_bin_augmented",
            title_suffix=("augmented score" if predictor_map else "augmented score (same as base; predictors unavailable)"),
        )
        plot_variant_enrichment(
            variant_rows,
            gene_dir / "criterion1_variant_enrichment.png",
            score_col=("allele_severity_points_augmented" if use_augmented else "allele_severity_points_base"),
            score_bin_col=("score_bin_augmented" if use_augmented else "score_bin_base"),
            title_suffix=("augmented score" if use_augmented else "base score"),
        )
        plot_decision_clarity(decision_rows, gene_dir / "criterion3_decision_clarity.png")

        write_csv(
            gene_dir / "criterion1_enrichment_compare.csv",
            [
                {
                    "mode": "base",
                    "status": base_variant_signal["status"],
                    "n_labeled": base_variant_signal["n_labeled"],
                    "split_method": base_variant_signal["split_method"],
                    "low_bin_n": base_variant_signal["low_bin_n"],
                    "high_bin_n": base_variant_signal["high_bin_n"],
                    "low_bin_pathogenic_rate": (
                        "" if base_variant_signal["low_bin_pathogenic_rate"] is None else round(float(base_variant_signal["low_bin_pathogenic_rate"]), 4)
                    ),
                    "high_bin_pathogenic_rate": (
                        "" if base_variant_signal["high_bin_pathogenic_rate"] is None else round(float(base_variant_signal["high_bin_pathogenic_rate"]), 4)
                    ),
                    "spearman_score_vs_pathogenic": (
                        "" if base_variant_signal["spearman_score_vs_pathogenic"] is None else round(float(base_variant_signal["spearman_score_vs_pathogenic"]), 4)
                    ),
                },
                {
                    "mode": "augmented",
                    "status": augmented_variant_signal["status"],
                    "n_labeled": augmented_variant_signal["n_labeled"],
                    "split_method": augmented_variant_signal["split_method"],
                    "low_bin_n": augmented_variant_signal["low_bin_n"],
                    "high_bin_n": augmented_variant_signal["high_bin_n"],
                    "low_bin_pathogenic_rate": (
                        "" if augmented_variant_signal["low_bin_pathogenic_rate"] is None else round(float(augmented_variant_signal["low_bin_pathogenic_rate"]), 4)
                    ),
                    "high_bin_pathogenic_rate": (
                        "" if augmented_variant_signal["high_bin_pathogenic_rate"] is None else round(float(augmented_variant_signal["high_bin_pathogenic_rate"]), 4)
                    ),
                    "spearman_score_vs_pathogenic": (
                        "" if augmented_variant_signal["spearman_score_vs_pathogenic"] is None else round(float(augmented_variant_signal["spearman_score_vs_pathogenic"]), 4)
                    ),
                },
            ],
        )

        criterion1_rows: List[Dict[str, object]] = []
        for bin_name in SCORE_BIN_ORDER:
            n_bin = int(variant_signal["by_bin_counts"].get(bin_name, 0))  # type: ignore[index]
            rate = variant_signal["by_bin_rates"].get(bin_name)  # type: ignore[index]
            criterion1_rows.append(
                {
                    "score_bin": bin_name,
                    "n_labeled": n_bin,
                    "pathogenic_rate": "" if rate is None else round(float(rate), 4),
                }
            )
        write_csv(gene_dir / "criterion1_score_bin_enrichment.csv", criterion1_rows)

        success_metric_row = {
            "gene_symbol": gene,
            "criterion_1_variant_enrichment_status": variant_signal["status"],
            "criterion_1_variant_enrichment_note": variant_signal["note"],
            "criterion_1_n_labeled": variant_signal["n_labeled"],
            "criterion_1_split_method": variant_signal["split_method"],
            "criterion_1_low_bin_n": variant_signal["low_bin_n"],
            "criterion_1_high_bin_n": variant_signal["high_bin_n"],
            "criterion_1_low_bin_pathogenic_rate": (
                "" if variant_signal["low_bin_pathogenic_rate"] is None else round(float(variant_signal["low_bin_pathogenic_rate"]), 4)
            ),
            "criterion_1_high_bin_pathogenic_rate": (
                "" if variant_signal["high_bin_pathogenic_rate"] is None else round(float(variant_signal["high_bin_pathogenic_rate"]), 4)
            ),
            "criterion_1_spearman_score_vs_pathogenic": (
                "" if variant_signal["spearman_score_vs_pathogenic"] is None else round(float(variant_signal["spearman_score_vs_pathogenic"]), 4)
            ),
            "criterion_3_decision_clarity_status": decision_signal["status"],
            "criterion_3_decision_clarity_note": decision_signal["note"],
            "criterion_3_n_decision_rows": decision_signal["n_decision_rows"],
            "criterion_3_n_wetlab_priority": decision_signal["n_wetlab_priority"],
            "criterion_3_n_deprioritized": decision_signal["n_deprioritized"],
            "scoring_mode": "augmented" if use_augmented else "base",
            "predictor_rows_with_any_score": predictor_stats["rows_with_any_predictor"],
            "predictor_rows_with_adjusted_score": predictor_stats["rows_with_adjusted_score"],
        }
        write_csv(gene_dir / "success_metrics.csv", [success_metric_row])

        n_labeled = 0
        if prediction_summary:
            n_labeled = int((prediction_summary[0].get("n_labeled") or 0))
        run_summaries.append(
            {
                "gene_symbol": gene,
                "filtered_curated_rows": len(gene_curated_rows),
                "processed_rows": len(variant_rows),
                "prediction_labeled_rows": n_labeled,
                "decision_rows": len(decision_rows),
                "scoring_mode": "augmented" if use_augmented else "base",
                "predictor_rows_with_any_score": predictor_stats["rows_with_any_predictor"],
                "predictor_rows_with_adjusted_score": predictor_stats["rows_with_adjusted_score"],
                "criterion_1_variant_enrichment_status": variant_signal["status"],
                "criterion_3_decision_clarity_status": decision_signal["status"],
            }
        )

        success_rows.append(
            {
                "gene_symbol": gene,
                "criterion_1_variant_enrichment_status": variant_signal["status"],
                "criterion_1_variant_enrichment_note": variant_signal["note"],
                "criterion_1_n_labeled": variant_signal["n_labeled"],
                "criterion_1_split_method": variant_signal["split_method"],
                "criterion_1_low_bin_n": variant_signal["low_bin_n"],
                "criterion_1_high_bin_n": variant_signal["high_bin_n"],
                "criterion_1_low_bin_pathogenic_rate": (
                    "" if variant_signal["low_bin_pathogenic_rate"] is None else round(float(variant_signal["low_bin_pathogenic_rate"]), 4)
                ),
                "criterion_1_high_bin_pathogenic_rate": (
                    "" if variant_signal["high_bin_pathogenic_rate"] is None else round(float(variant_signal["high_bin_pathogenic_rate"]), 4)
                ),
                "criterion_1_spearman_score_vs_pathogenic": (
                    "" if variant_signal["spearman_score_vs_pathogenic"] is None else round(float(variant_signal["spearman_score_vs_pathogenic"]), 4)
                ),
                "criterion_3_decision_clarity_status": decision_signal["status"],
                "criterion_3_decision_clarity_note": decision_signal["note"],
                "criterion_3_n_decision_rows": decision_signal["n_decision_rows"],
                "criterion_3_n_wetlab_priority": decision_signal["n_wetlab_priority"],
                "criterion_3_n_deprioritized": decision_signal["n_deprioritized"],
                "scoring_mode": "augmented" if use_augmented else "base",
                "predictor_rows_with_any_score": predictor_stats["rows_with_any_predictor"],
                "predictor_rows_with_adjusted_score": predictor_stats["rows_with_adjusted_score"],
            }
        )

        print(f"[OK] {gene}: wrote per-gene outputs to {gene_dir}")

    write_csv(out_root / "per_gene_run_summary.csv", run_summaries)
    write_csv(out_root / "success_criteria_summary.csv", success_rows)
    (out_root / "keep_regions.json").write_text(json.dumps(keep_regions, indent=2), encoding="utf-8")
    predictor_manifest = {
        "spliceai_wide_csv": str(spliceai_wide_path),
        "predictor_map_variants": len(predictor_map),
        "score_mode_requested": args.score_mode,
        "score_mode_effective": "augmented" if use_augmented else "base",
    }
    (out_root / "predictor_integration.json").write_text(json.dumps(predictor_manifest, indent=2), encoding="utf-8")
    stale_overview = out_root / "success_criteria_overview.png"
    if stale_overview.exists():
        stale_overview.unlink()

    print(f"[OK] Wrote {out_root / 'curated_variants_region_filtered.csv'} ({len(filtered_rows)} rows)")
    print(f"[OK] Wrote {out_root / 'region_filter_summary.csv'} ({len(summary_rows)} rows)")
    print(f"[OK] Wrote {out_root / 'per_gene_run_summary.csv'} ({len(run_summaries)} rows)")
    print(f"[OK] Wrote {out_root / 'success_criteria_summary.csv'} ({len(success_rows)} rows)")
    print(f"[OK] Wrote {out_root / 'predictor_integration.json'}")
    print(f"[OK] Kept regions: {keep_regions}")
    print(f"[OK] Score mode: {'augmented' if use_augmented else 'base'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
