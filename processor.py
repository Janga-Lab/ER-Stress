#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional, Tuple

from config import (
    CURATED_VARIANTS_CSV_PATH,
    DEFAULT_CLASS_WEIGHTS,
    DEFAULT_CONFIDENCE_HIGH_THRESHOLD,
    DEFAULT_CONFIDENCE_LOW_THRESHOLD,
    DEFAULT_INHERITANCE_MAP,
    DEFAULT_SCORE_BIN_THRESHOLDS,
    DEFAULT_UNMAPPED_WARNING_THRESHOLD,
    DEFAULT_VARIANT_PREDICTION_THRESHOLD,
    PROCESSED_DATA_DIR,
    PROCESSED_DECISION_TIERS_CSV_PATH,
    PROCESSED_PRIORITY_RANKED_CSV_PATH,
    PROCESSED_RUN_MANIFEST_PATH,
    PROCESSED_VARIANT_CSV_PATH,
    PROCESSED_VARIANT_PREDICTION_CSV_PATH,
    PROCESSED_VARIANT_SUMMARY_CSV_PATH,
    PipelineConfig,
    ensure_data_directories,
)

PATHOGENIC_SIG_VALUES = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "likely pathogenic/pathogenic",
}

BENIGN_SIG_VALUES = {
    "benign",
    "likely benign",
    "benign/likely benign",
    "likely benign/benign",
}

HIGH_IMPACT_CLASS_SET = {"lof_truncating", "splice", "cnv_sv"}
SCORE_BIN_ORDER = ["low", "mild", "moderate", "high", "very_high"]


def clean_optional_text(value: Any) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip()
        return text or None
    return str(value)


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
    if f is None:
        return None
    return int(f)


def canonical_variation_id(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def read_csv_rows(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    columns = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(65536)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def score_bin(score: float, thresholds: Dict[str, float]) -> str:
    if score >= thresholds["very_high"]:
        return "very_high"
    if score >= thresholds["high"]:
        return "high"
    if score >= thresholds["moderate"]:
        return "moderate"
    if score >= thresholds["mild"]:
        return "mild"
    return "low"


def clinvar_binary_label(sig: Any) -> Optional[int]:
    text = clean_optional_text(sig)
    if not text:
        return None
    s = text.lower()
    if s in PATHOGENIC_SIG_VALUES:
        return 1
    if s in BENIGN_SIG_VALUES:
        return 0
    return None


def auc_rank(y_true: List[int], y_score: List[float]) -> Optional[float]:
    if len(y_true) != len(y_score) or not y_true:
        return None

    n_pos = sum(1 for x in y_true if x == 1)
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None

    pairs = sorted(zip(y_score, y_true), key=lambda x: x[0])
    ranks = [0.0] * len(pairs)
    i = 0
    rank = 1
    while i < len(pairs):
        j = i
        while j + 1 < len(pairs) and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        avg_rank = (rank + rank + (j - i)) / 2.0
        for k in range(i, j + 1):
            ranks[k] = avg_rank
        rank += (j - i + 1)
        i = j + 1

    sum_ranks_pos = sum(r for r, (_, y) in zip(ranks, pairs) if y == 1)
    auc = (sum_ranks_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def threshold_metrics(y_true: List[int], y_score: List[float], threshold: float) -> Dict[str, float]:
    if not y_true:
        return {
            "accuracy": 0.0,
            "precision": 0.0,
            "recall": 0.0,
            "f1": 0.0,
        }

    y_pred = [1 if score >= threshold else 0 for score in y_score]
    tp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 1)
    tn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 0)
    fp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 1)
    fn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 0)

    n = len(y_true)
    accuracy = (tp + tn) / n if n else 0.0
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def process_variants(curated_rows: List[Dict[str, Any]], cfg: PipelineConfig) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in curated_rows:
        variant_bucket = clean_optional_text(row.get("variant_bucket")) or "other"
        critical = int(parse_int(row.get("is_critical_region")) or 0) == 1

        base = cfg.class_weights.get(variant_bucket, 0.0)
        if base <= 0:
            points = 0.0
        else:
            points = base + (cfg.critical_bonus if critical else 0.0)

        sig = clean_optional_text(row.get("clinvar_clinical_significance")) or ""
        path_binary = clinvar_binary_label(sig)
        if path_binary is None:
            path_binary = parse_int(row.get("pathogenic_binary"))

        rec = dict(row)
        rec["gene_symbol"] = (clean_optional_text(rec.get("gene_symbol")) or "").upper()
        rec["clinvar_variation_id"] = canonical_variation_id(rec.get("clinvar_variation_id"))
        rec["allele_severity_points"] = round(points, 4)
        rec["score_bin"] = score_bin(points, cfg.score_bin_thresholds)
        rec["pathogenic_binary"] = "" if path_binary is None else int(path_binary)
        rec["is_critical_region"] = int(critical)
        out.append(rec)

    return out


def summarize_variants(variant_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    counts: Counter[Tuple[str, str, str]] = Counter()
    points_store: Dict[Tuple[str, str, str], List[float]] = defaultdict(list)

    for row in variant_rows:
        gene = row.get("gene_symbol") or ""
        kind = row.get("variant_bucket") or ""
        region = row.get("region_bucket") or ""
        key = (gene, kind, region)
        counts[key] += 1
        points_store[key].append(float(row.get("allele_severity_points") or 0.0))

    summary: List[Dict[str, Any]] = []
    for key in sorted(counts):
        gene, kind, region = key
        pts = points_store[key]
        summary.append(
            {
                "gene_symbol": gene,
                "variant_bucket": kind,
                "region_bucket": region,
                "n": counts[key],
                "mean_points": round(sum(pts) / len(pts), 4) if pts else 0.0,
            }
        )
    return summary


def evaluate_variant_prediction(
    variant_rows: List[Dict[str, Any]],
    threshold: float,
) -> List[Dict[str, Any]]:
    by_gene: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in variant_rows:
        label = parse_int(row.get("pathogenic_binary"))
        if label is None:
            continue
        by_gene[row.get("gene_symbol") or "UNKNOWN"].append(row)

    summary_rows: List[Dict[str, Any]] = []

    all_rows = [r for rows in by_gene.values() for r in rows]
    if all_rows:
        by_gene = {"ALL": all_rows, **dict(sorted(by_gene.items()))}
    else:
        by_gene = dict(sorted(by_gene.items()))

    for gene, rows in by_gene.items():
        y = [int(r["pathogenic_binary"]) for r in rows]
        score = [float(r.get("allele_severity_points") or 0.0) for r in rows]
        lof_binary = [1.0 if (r.get("variant_bucket") == "lof_truncating") else 0.0 for r in rows]
        critical_binary = [float(int(r.get("is_critical_region") or 0)) for r in rows]

        score_auc = auc_rank(y, score)
        lof_auc = auc_rank(y, lof_binary)
        critical_auc = auc_rank(y, critical_binary)
        metrics = threshold_metrics(y, score, threshold)

        summary_rows.append(
            {
                "gene_symbol": gene,
                "n_labeled": len(rows),
                "n_pathogenic_like": sum(y),
                "n_benign_like": len(rows) - sum(y),
                "auc_score": round(score_auc, 4) if score_auc is not None else "",
                "auc_baseline_lof_binary": round(lof_auc, 4) if lof_auc is not None else "",
                "auc_baseline_critical_binary": round(critical_auc, 4) if critical_auc is not None else "",
                "threshold": threshold,
                "accuracy_at_threshold": round(metrics["accuracy"], 4),
                "precision_at_threshold": round(metrics["precision"], 4),
                "recall_at_threshold": round(metrics["recall"], 4),
                "f1_at_threshold": round(metrics["f1"], 4),
            }
        )

    return summary_rows


def plot_auc_comparison_dumbbell(summary_rows: List[Dict[str, Any]], out_path: Path) -> Optional[Path]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    filtered = []
    for row in summary_rows:
        gene = clean_optional_text(row.get("gene_symbol")) or ""
        if not gene:
            continue
        auc_score = parse_float(row.get("auc_score"))
        auc_lof = parse_float(row.get("auc_baseline_lof_binary"))
        auc_crit = parse_float(row.get("auc_baseline_critical_binary"))
        if auc_score is None and auc_lof is None and auc_crit is None:
            continue
        filtered.append((gene, auc_score, auc_lof, auc_crit))

    if not filtered:
        return None

    filtered.sort(key=lambda x: (0, x[0]) if x[0] == "ALL" else (1, x[0]))
    y = list(range(len(filtered)))
    score_vals = [x[1] for x in filtered]
    lof_vals = [x[2] for x in filtered]
    crit_vals = [x[3] for x in filtered]
    labels = [x[0] for x in filtered]

    plt.figure(figsize=(8.6, max(4.2, 0.44 * len(filtered) + 1.6)))
    for yi, vals in enumerate(zip(score_vals, lof_vals, crit_vals)):
        present = [v for v in vals if v is not None]
        if present:
            plt.hlines(yi, min(present), max(present), color="#cfd8dc", linewidth=2.0, zorder=1)

    def scatter(values: List[Optional[float]], color: str, label: str, marker: str) -> None:
        xs, ys = [], []
        for yi, val in enumerate(values):
            if val is None:
                continue
            xs.append(val)
            ys.append(yi)
        if xs:
            plt.scatter(xs, ys, s=66, color=color, label=label, marker=marker, edgecolor="white", linewidth=0.6, zorder=3)

    scatter(score_vals, "#1e88e5", "Score model", "o")
    scatter(lof_vals, "#ef6c00", "LOF baseline", "s")
    scatter(crit_vals, "#6a1b9a", "Critical baseline", "D")

    plt.yticks(y, labels)
    plt.xlim(0.45, 1.02)
    plt.xlabel("AUC")
    plt.title("Prediction Comparison Across Genes")
    plt.grid(axis="x", alpha=0.2, linewidth=0.8)
    plt.legend(loc="lower right", fontsize=8)
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=180)
    plt.close()
    return out_path


def plot_gene_scorebin_heatmap(variant_rows: List[Dict[str, Any]], out_path: Path) -> Optional[Path]:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        return None

    labeled_rows = [r for r in variant_rows if parse_int(r.get("pathogenic_binary")) in {0, 1}]
    if not labeled_rows:
        return None

    genes = sorted({clean_optional_text(r.get("gene_symbol")) or "UNKNOWN" for r in labeled_rows})
    genes = ["ALL"] + [g for g in genes if g != "ALL"]

    counts: Dict[Tuple[str, str], int] = defaultdict(int)
    pathogenic: Dict[Tuple[str, str], int] = defaultdict(int)

    for row in labeled_rows:
        g = clean_optional_text(row.get("gene_symbol")) or "UNKNOWN"
        b = (clean_optional_text(row.get("score_bin")) or "low").lower()
        if b not in SCORE_BIN_ORDER:
            continue
        y = int(row.get("pathogenic_binary") or 0)
        counts[(g, b)] += 1
        pathogenic[(g, b)] += y
        counts[("ALL", b)] += 1
        pathogenic[("ALL", b)] += y

    matrix = np.full((len(genes), len(SCORE_BIN_ORDER)), np.nan)
    labels = [["" for _ in SCORE_BIN_ORDER] for _ in genes]

    for i, gene in enumerate(genes):
        for j, bin_name in enumerate(SCORE_BIN_ORDER):
            n = counts.get((gene, bin_name), 0)
            if n == 0:
                labels[i][j] = "NA"
                continue
            rate = pathogenic.get((gene, bin_name), 0) / n
            matrix[i, j] = rate
            labels[i][j] = f"{int(round(rate * 100))}%\n(n={n})"

    plt.figure(figsize=(1.75 * len(SCORE_BIN_ORDER) + 3.1, 0.55 * len(genes) + 2.0))
    im = plt.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0.0, vmax=1.0)
    plt.xticks(range(len(SCORE_BIN_ORDER)), SCORE_BIN_ORDER)
    plt.yticks(range(len(genes)), genes)
    plt.colorbar(im, fraction=0.032, pad=0.02, label="Pathogenic-like fraction")

    for i in range(len(genes)):
        for j in range(len(SCORE_BIN_ORDER)):
            text = labels[i][j]
            val = matrix[i, j]
            color = "white" if (not np.isnan(val) and val > 0.6) else "#2f3e46"
            plt.text(j, i, text, ha="center", va="center", fontsize=8, color=color)

    plt.title("Pathogenic Enrichment by Score Bin and Gene")
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=180)
    plt.close()
    return out_path


def plot_impact_region_profile(variant_rows: List[Dict[str, Any]], out_path: Path) -> Optional[Path]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    labeled_rows = [r for r in variant_rows if parse_int(r.get("pathogenic_binary")) in {0, 1}]
    if not labeled_rows:
        return None

    categories = ["high+critical", "high_only", "critical_only", "other"]
    counts = Counter()
    pathogenic = Counter()

    for row in labeled_rows:
        high_impact = (row.get("variant_bucket") or "") in HIGH_IMPACT_CLASS_SET
        critical = int(row.get("is_critical_region") or 0) == 1
        if high_impact and critical:
            cat = "high+critical"
        elif high_impact and not critical:
            cat = "high_only"
        elif critical and not high_impact:
            cat = "critical_only"
        else:
            cat = "other"
        counts[cat] += 1
        pathogenic[cat] += int(row.get("pathogenic_binary") or 0)

    rates = [(pathogenic[c] / counts[c]) if counts[c] else 0.0 for c in categories]
    ns = [counts[c] for c in categories]

    plt.figure(figsize=(8.2, 4.5))
    bars = plt.bar(categories, rates, color=["#ef6c00", "#8d6e63", "#6a1b9a", "#546e7a"], alpha=0.85)
    plt.ylim(0, 1.0)
    plt.ylabel("Pathogenic-like fraction")
    plt.title("Pathogenic Signal by Impact/Region Profile")
    plt.grid(axis="y", alpha=0.2)

    for bar, n, rate in zip(bars, ns, rates):
        x = bar.get_x() + bar.get_width() / 2
        plt.text(x, rate + 0.03, f"n={n}\n{int(round(rate * 100))}%", ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=180)
    plt.close()
    return out_path


def build_decision_tables(
    variant_rows: List[Dict[str, Any]],
    prediction_summary: List[Dict[str, Any]],
    cfg: PipelineConfig,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    gene_stats = {
        row["gene_symbol"]: {
            "n_labeled": parse_int(row.get("n_labeled")) or 0,
            "auc_score": parse_float(row.get("auc_score")),
        }
        for row in prediction_summary
        if clean_optional_text(row.get("gene_symbol"))
    }

    decision_rows: List[Dict[str, Any]] = []
    for row in variant_rows:
        gene = row.get("gene_symbol") or ""
        stats = gene_stats.get(gene, {"n_labeled": 0, "auc_score": None})

        n_labeled = stats["n_labeled"]
        auc_score = stats["auc_score"]

        evidence_component = min(1.0, n_labeled / 80.0)
        model_component = 0.0
        if auc_score is not None:
            model_component = max(0.0, min(1.0, (auc_score - 0.5) / 0.4))
        mapping_component = 1.0 if canonical_variation_id(row.get("clinvar_variation_id")) else 0.0

        confidence_score = round(0.45 * evidence_component + 0.45 * model_component + 0.10 * mapping_component, 4)
        if confidence_score >= cfg.confidence_high_threshold:
            confidence_level = "high"
        elif confidence_score >= cfg.confidence_low_threshold:
            confidence_level = "medium"
        else:
            confidence_level = "low"

        bin_name = row.get("score_bin") or "low"
        sig_text = (clean_optional_text(row.get("clinvar_clinical_significance")) or "").lower()
        benign_like = sig_text in BENIGN_SIG_VALUES
        pathogenic_like = sig_text in PATHOGENIC_SIG_VALUES

        predicted_severe = bin_name in {"high", "very_high"}
        predicted_low = bin_name in {"low", "mild"}

        major_discordance = bool(
            (predicted_severe and benign_like)
            or (predicted_low and pathogenic_like)
        )

        needs_wetlab = bool(confidence_level == "low" or major_discordance)

        if predicted_severe and confidence_level == "high" and not major_discordance:
            tier = "Tier_A"
        elif predicted_severe and confidence_level in {"high", "medium"} and not major_discordance:
            tier = "Tier_B"
        elif needs_wetlab:
            tier = "Tier_C"
        else:
            tier = "Tier_D"

        reason_codes = []
        if predicted_severe:
            reason_codes.append("predicted_severe")
        if confidence_level == "low":
            reason_codes.append("low_confidence")
        if major_discordance:
            reason_codes.append("major_discordance")
        if int(row.get("is_critical_region") or 0) == 1:
            reason_codes.append("critical_region")

        decision_rows.append(
            {
                "gene": gene,
                "variant_or_genotype_id": canonical_variation_id(row.get("clinvar_variation_id")),
                "predicted_severity": float(row.get("allele_severity_points") or 0.0),
                "score_bin": bin_name,
                "confidence_score": confidence_score,
                "confidence_level": confidence_level,
                "major_discordance": int(major_discordance),
                "needs_wetlab_validation": int(needs_wetlab),
                "decision_tier": tier,
                "reason_codes": ";".join(reason_codes),
            }
        )

    ranked = sorted(
        decision_rows,
        key=lambda r: (
            -int(r["needs_wetlab_validation"]),
            float(r["confidence_score"]),
            -float(r["predicted_severity"]),
            -int(r["major_discordance"]),
            r["gene"],
            r["variant_or_genotype_id"],
        ),
    )

    return decision_rows, ranked


def plot_decision_landscape(decision_rows: List[Dict[str, Any]], out_path: Path) -> Optional[Path]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    if not decision_rows:
        return None

    x = [float(r.get("predicted_severity") or 0.0) for r in decision_rows]
    y = [float(r.get("confidence_score") or 0.0) for r in decision_rows]
    c = ["#ef6c00" if int(r.get("needs_wetlab_validation") or 0) == 1 else "#1e88e5" for r in decision_rows]

    plt.figure(figsize=(7.5, 5.2))
    plt.scatter(x, y, c=c, alpha=0.75, edgecolors="white", linewidths=0.5)
    plt.xlabel("Predicted Severity Points")
    plt.ylabel("Confidence Score")
    plt.title("Decision Landscape: Severity vs Confidence")
    plt.grid(alpha=0.2)

    from matplotlib.lines import Line2D

    legend_items = [
        Line2D([0], [0], marker="o", color="w", label="Needs wet-lab", markerfacecolor="#ef6c00", markersize=8),
        Line2D([0], [0], marker="o", color="w", label="No immediate wet-lab", markerfacecolor="#1e88e5", markersize=8),
    ]
    plt.legend(handles=legend_items, loc="lower right")
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=180)
    plt.close()
    return out_path


def build_run_manifest(
    input_files: Iterable[Path],
    output_files: Iterable[Path],
    cfg: PipelineConfig,
    runtime_env_overrides: Dict[str, str],
) -> Dict[str, Any]:
    return {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_env_overrides": runtime_env_overrides,
        "config": {
            "variant_prediction_threshold": cfg.variant_prediction_threshold,
            "critical_bonus": cfg.critical_bonus,
            "class_weights": cfg.class_weights,
            "inheritance_map": cfg.inheritance_map,
            "confidence_low_threshold": cfg.confidence_low_threshold,
            "confidence_high_threshold": cfg.confidence_high_threshold,
            "unmapped_warning_threshold": cfg.unmapped_warning_threshold,
        },
        "input_files": [
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "exists": path.exists(),
            }
            for path in input_files
        ],
        "output_files": [
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "exists": path.exists(),
            }
            for path in output_files
        ],
    }


def configure_runtime_env(out_dir: Path) -> Dict[str, str]:
    """
    Ensure plotting libraries use writable cache/config directories.
    """
    applied: Dict[str, str] = {}

    if not os.getenv("MPLBACKEND"):
        os.environ["MPLBACKEND"] = "Agg"
        applied["MPLBACKEND"] = "Agg"

    if not os.getenv("MPLCONFIGDIR"):
        mpl_dir = Path(tempfile.gettempdir()) / "erscore_mplconfig"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)
        applied["MPLCONFIGDIR"] = str(mpl_dir)

    if not os.getenv("XDG_CACHE_HOME"):
        cache_dir = Path(tempfile.gettempdir()) / "erscore_cache"
        (cache_dir / "fontconfig").mkdir(parents=True, exist_ok=True)
        os.environ["XDG_CACHE_HOME"] = str(cache_dir)
        applied["XDG_CACHE_HOME"] = str(cache_dir)

    return applied


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Process curated data into scores, predictions, triage, and plots")
    parser.add_argument("--curated-variants", default=str(CURATED_VARIANTS_CSV_PATH))
    parser.add_argument("--out-dir", default=str(PROCESSED_DATA_DIR))
    parser.add_argument("--variant-pred-threshold", type=float, default=DEFAULT_VARIANT_PREDICTION_THRESHOLD)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_data_directories()

    curated_variants_path = Path(args.curated_variants)
    if not curated_variants_path.exists():
        raise SystemExit(f"Curated variants file not found: {curated_variants_path}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    runtime_env_overrides = configure_runtime_env(out_dir)

    cfg = PipelineConfig(
        class_weights=dict(DEFAULT_CLASS_WEIGHTS),
        inheritance_map=dict(DEFAULT_INHERITANCE_MAP),
        score_bin_thresholds=dict(DEFAULT_SCORE_BIN_THRESHOLDS),
        variant_prediction_threshold=float(args.variant_pred_threshold),
        confidence_low_threshold=DEFAULT_CONFIDENCE_LOW_THRESHOLD,
        confidence_high_threshold=DEFAULT_CONFIDENCE_HIGH_THRESHOLD,
        unmapped_warning_threshold=DEFAULT_UNMAPPED_WARNING_THRESHOLD,
    )

    curated_rows = read_csv_rows(curated_variants_path)
    variant_rows = process_variants(curated_rows, cfg)
    variant_summary = summarize_variants(variant_rows)
    prediction_summary = evaluate_variant_prediction(variant_rows, cfg.variant_prediction_threshold)
    decision_rows, priority_ranked = build_decision_tables(variant_rows, prediction_summary, cfg)

    variant_processed_csv = out_dir / PROCESSED_VARIANT_CSV_PATH.name
    variant_summary_csv = out_dir / PROCESSED_VARIANT_SUMMARY_CSV_PATH.name
    prediction_csv = out_dir / PROCESSED_VARIANT_PREDICTION_CSV_PATH.name
    decision_csv = out_dir / PROCESSED_DECISION_TIERS_CSV_PATH.name
    priority_csv = out_dir / PROCESSED_PRIORITY_RANKED_CSV_PATH.name

    write_csv(variant_processed_csv, variant_rows)
    write_csv(variant_summary_csv, variant_summary)
    write_csv(prediction_csv, prediction_summary)
    write_csv(decision_csv, decision_rows)
    write_csv(priority_csv, priority_ranked)

    # Variant-level plots
    plot_auc_comparison_dumbbell(prediction_summary, out_dir / "prediction_auc_dumbbell.png")
    plot_gene_scorebin_heatmap(variant_rows, out_dir / "gene_scorebin_pathogenic_heatmap.png")
    plot_impact_region_profile(variant_rows, out_dir / "impact_region_pathogenic_profile.png")
    plot_decision_landscape(decision_rows, out_dir / "decision_landscape.png")

    run_manifest = build_run_manifest(
        input_files=[curated_variants_path],
        output_files=[
            variant_processed_csv,
            variant_summary_csv,
            prediction_csv,
            decision_csv,
            priority_csv,
        ],
        cfg=cfg,
        runtime_env_overrides=runtime_env_overrides,
    )

    manifest_path = out_dir / PROCESSED_RUN_MANIFEST_PATH.name
    manifest_path.write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")

    print(f"[OK] Wrote {variant_processed_csv} ({len(variant_rows)} rows)")
    print(f"[OK] Wrote {variant_summary_csv} ({len(variant_summary)} rows)")
    print(f"[OK] Wrote {prediction_csv} ({len(prediction_summary)} rows)")
    print(f"[OK] Wrote {decision_csv} ({len(decision_rows)} rows)")
    print(f"[OK] Wrote {priority_csv} ({len(priority_ranked)} rows)")
    print(f"[OK] Wrote {manifest_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
