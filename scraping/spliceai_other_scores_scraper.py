#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import time
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
from urllib.parse import quote

try:
    from config import PROCESSED_DATA_DIR, RAW_CLINICAL_ASSOC_JSONL_PATH, SPLICEAI_LOOKUP_BASE_URL
except Exception:
    PROCESSED_DATA_DIR = Path("data/processed_data")
    RAW_CLINICAL_ASSOC_JSONL_PATH = Path("data/raw_data/harvest_raw_clinical_association.jsonl")
    SPLICEAI_LOOKUP_BASE_URL = "https://spliceailookup.broadinstitute.org"

PREDICTOR_LABEL_TO_KEY = {
    "AlphaMissense": "alphamissense",
    "CADD": "cadd",
    "PhyloP": "phylop",
    "PolyPhen (max)": "polyphen_max",
    "PrimateAI-3D": "primateai3d",
    "PromoterAI": "promoterai",
    "REVEL": "revel_max",
    "SIFT (max)": "sift_max",
}

PREDICTOR_ORDER = [
    "alphamissense",
    "cadd",
    "phylop",
    "polyphen_max",
    "primateai3d",
    "promoterai",
    "revel_max",
    "sift_max",
]

POINT_BASED_PREDICTORS = {
    "alphamissense",
    "cadd",
    "phylop",
    "polyphen_max",
    "revel_max",
    "sift_max",
}

RETRYABLE_STATUSES = {
    "lookup_error",
}


def log_line(message: str, log_handle=None) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {message}"
    print(line, flush=True)
    if log_handle:
        log_handle.write(line + "\n")
        log_handle.flush()


def normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def normalize_label(text: str) -> str:
    return normalize_space(text).lower()


LABEL_NORMALIZED_TO_KEY = {
    normalize_label(label): key for label, key in PREDICTOR_LABEL_TO_KEY.items()
}


def resolve_predictor_key(label: str) -> Optional[str]:
    norm = normalize_label(label)
    direct = LABEL_NORMALIZED_TO_KEY.get(norm)
    if direct:
        return direct

    if "alphamissense" in norm:
        return "alphamissense"
    if "polyphen" in norm:
        return "polyphen_max"
    if norm.startswith("sift"):
        return "sift_max"
    if norm.startswith("revel"):
        return "revel_max"
    if "primateai" in norm:
        return "primateai3d"
    if "promoterai" in norm:
        return "promoterai"
    if "phylop" in norm:
        return "phylop"
    if norm.startswith("cadd"):
        return "cadd"
    return None


def parse_float_from_text(text: str) -> Optional[float]:
    if not text:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?", text)
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def parse_point_from_text(text: str) -> str:
    if not text:
        return ""
    match = re.search(r"[+-]\d+", text)
    return match.group(0) if match else ""


def is_rate_limit_error_text(text: str) -> bool:
    norm = normalize_space(text).lower()
    if not norm:
        return False
    return ("rate limit exceeded" in norm) or ("only supports interactive use" in norm)


def build_lookup_url(
    variant: str,
    assembly: str = "37",
    bc: str = "basic",
    distance: int = 500,
    mask: int = 0,
    ra: int = 0,
) -> str:
    variant_q = quote(variant, safe="")
    return (
        f"{SPLICEAI_LOOKUP_BASE_URL}/"
        f"#variant={variant_q}&hg={assembly}&bc={bc}&distance={distance}&mask={mask}&ra={ra}"
    )


def load_variants(input_csv: Path, limit: Optional[int] = None) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    seen: set[Tuple[str, str, str]] = set()

    def append_row(gene: str, variation_id: str, hgvs_c: str, hgvs_p: str) -> None:
        query_variant = hgvs_c or hgvs_p
        key = (gene, variation_id, query_variant)
        if key in seen:
            return
        seen.add(key)
        rows.append(
            {
                "gene_symbol": gene,
                "clinvar_variation_id": variation_id,
                "query_variant": query_variant,
            }
        )

    if input_csv.suffix.lower() == ".jsonl":
        with input_csv.open("r", encoding="utf-8") as handle:
            for line in handle:
                text = line.strip()
                if not text:
                    continue
                try:
                    row = json.loads(text)
                except json.JSONDecodeError:
                    continue
                if row.get("record_type") != "variant":
                    continue
                gene = (row.get("gene_symbol") or "").strip()
                variation_id = (row.get("clinvar_variation_id") or "").strip()
                hgvs_c = (row.get("hgvs_c") or row.get("HGVS_c") or "").strip()
                hgvs_p = (row.get("hgvs_p") or row.get("HGVS_p") or "").strip()
                append_row(gene, variation_id, hgvs_c, hgvs_p)
                if limit is not None and len(rows) >= limit:
                    break
        return rows

    with input_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            gene = (row.get("gene_symbol") or "").strip()
            variation_id = (row.get("clinvar_variation_id") or "").strip()
            hgvs_c = (row.get("hgvs_c") or "").strip()
            hgvs_p = (row.get("hgvs_p") or "").strip()
            append_row(gene, variation_id, hgvs_c, hgvs_p)
            if limit is not None and len(rows) >= limit:
                break
    return rows


def build_variant_key(row: Dict[str, Any]) -> Tuple[str, str, str]:
    return (
        (row.get("gene_symbol") or "").strip(),
        (row.get("clinvar_variation_id") or "").strip(),
        (row.get("query_variant") or "").strip(),
    )


def read_existing_status_map(path: Path) -> Dict[Tuple[str, str, str], str]:
    if not path.exists():
        return {}
    status_map: Dict[Tuple[str, str, str], str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            status_map[build_variant_key(row)] = (row.get("status") or "").strip().lower()
    return status_map


def read_keys_with_status(path: Path, target_status: str) -> set[Tuple[str, str, str]]:
    target = (target_status or "").strip().lower()
    if not path.exists() or not target:
        return set()
    keys: set[Tuple[str, str, str]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            status = (row.get("status") or "").strip().lower()
            if status == target:
                keys.add(build_variant_key(row))
    return keys


def should_retry_status(status: str) -> bool:
    s = (status or "").strip().lower()
    if not s:
        return True
    return s in RETRYABLE_STATUSES


def load_latest_wide_rows(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    latest_by_key: Dict[Tuple[str, str, str], Tuple[int, Dict[str, str]]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for idx, row in enumerate(reader):
            latest_by_key[build_variant_key(row)] = (idx, row)
    return [row for _, row in sorted(latest_by_key.values(), key=lambda x: x[0])]


def create_driver(browser: str, headless: bool):
    try:
        from selenium import webdriver
    except Exception as exc:
        raise SystemExit(
            "Selenium is required. Install with: pip install selenium"
        ) from exc

    if browser == "chrome":
        options = webdriver.ChromeOptions()
        if headless:
            options.add_argument("--headless=new")
        options.add_argument("--disable-gpu")
        options.add_argument("--window-size=1920,1080")
        options.add_argument("--disable-dev-shm-usage")
        options.add_argument("--no-sandbox")
        return webdriver.Chrome(options=options)

    if browser == "firefox":
        options = webdriver.FirefoxOptions()
        if headless:
            options.add_argument("-headless")
        return webdriver.Firefox(options=options)

    raise SystemExit(f"Unsupported browser: {browser}")


def load_lookup_page(driver, timeout_seconds: int) -> None:
    from selenium.webdriver.support.ui import WebDriverWait

    driver.get(f"{SPLICEAI_LOOKUP_BASE_URL}/")
    js_ready = """
        return !!document.getElementById("search-box")
            && !!document.getElementById("submit-button")
            && !!document.getElementById("other-predictors-header");
    """
    WebDriverWait(driver, timeout_seconds).until(lambda d: bool(d.execute_script(js_ready)))


def submit_lookup_request(
    driver,
    *,
    query_variant: str,
    assembly: str,
    bc: str,
    distance: int,
    mask: str,
) -> Dict[str, Any]:
    js_submit = """
        const variant = arguments[0];
        const hg = String(arguments[1]);
        const bc = String(arguments[2]);
        const distance = String(arguments[3]);
        const mask = String(arguments[4]);

        const out = {ok: false, error: "", hashBefore: window.location.hash || ""};

        const searchBox = document.getElementById("search-box");
        const submitButton = document.getElementById("submit-button");
        if (!searchBox) {
            out.error = "missing_search_box";
            return out;
        }
        if (!submitButton) {
            out.error = "missing_submit_button";
            return out;
        }

        // Clear previous table rows to avoid stale-read bugs.
        document.querySelectorAll("#other-predictors-table tr").forEach((tr) => {
            if (tr.querySelectorAll("td").length > 0) tr.remove();
        });

        const errorBox = document.getElementById("error-box");
        if (errorBox) {
            errorBox.innerHTML = "";
            errorBox.style.display = "none";
        }

        searchBox.value = variant;
        searchBox.dispatchEvent(new Event("input", { bubbles: true }));
        searchBox.dispatchEvent(new Event("change", { bubbles: true }));

        const hgInput = document.querySelector(`input[name='hg'][value='${hg}']`);
        if (hgInput) {
            hgInput.checked = true;
            hgInput.dispatchEvent(new Event("change", { bubbles: true }));
        }

        const bcInput = document.querySelector(`input[name='gencode-gene-set'][value='${bc}']`);
        if (bcInput) {
            bcInput.checked = true;
            bcInput.dispatchEvent(new Event("change", { bubbles: true }));
        }

        const maxDistanceInput = document.getElementById("max-distance-input");
        if (maxDistanceInput) {
            maxDistanceInput.value = distance;
            maxDistanceInput.dispatchEvent(new Event("input", { bubbles: true }));
            maxDistanceInput.dispatchEvent(new Event("change", { bubbles: true }));
        }

        const maskInput = document.querySelector("input[name='mask']");
        if (maskInput) {
            maskInput.checked = mask === "1";
            maskInput.dispatchEvent(new Event("change", { bubbles: true }));
        }

        submitButton.click();
        out.ok = true;
        return out;
    """
    result = driver.execute_script(
        js_submit,
        query_variant,
        assembly,
        bc,
        int(distance),
        str(mask),
    )
    return result or {"ok": False, "error": "submit_js_failed"}


def wait_for_submit_cycle(driver, timeout_seconds: int) -> None:
    from selenium.common.exceptions import TimeoutException
    from selenium.webdriver.support.ui import WebDriverWait

    js_loading = """
        const b = document.getElementById("submit-button");
        return !!(b && b.classList.contains("loading"));
    """
    js_done = """
        const submitButton = document.getElementById("submit-button");
        const isLoading = submitButton ? submitButton.classList.contains("loading") : false;
        const errorText = ((document.getElementById("error-box") || {}).innerText || "").trim();
        const rows = Array.from(document.querySelectorAll("#other-predictors-table tr"))
            .filter((tr) => tr.querySelectorAll("td").length > 0);
        const rowCount = rows.length;
        return (!isLoading) && (rowCount > 0 || errorText.length > 0);
    """

    try:
        WebDriverWait(driver, min(15, timeout_seconds)).until(lambda d: bool(d.execute_script(js_loading)))
    except TimeoutException:
        # Some runs can be very fast or fail before entering loading state.
        pass

    WebDriverWait(driver, timeout_seconds).until(lambda d: bool(d.execute_script(js_done)))


def get_other_scores_rows_and_error(driver) -> Tuple[List[List[str]], str, str, str]:
    js_extract = """
        const out = {rows: [], error: "", firstRowText: "", currentHash: window.location.hash || ""};
        out.error = ((document.getElementById("error-box") || {}).innerText || "").replace(/\\s+/g, " ").trim();
        const header = document.getElementById("other-predictors-header");
        if (!header) return out;
        const table = header.closest("table");
        if (!table) return out;
        const dataRows = Array.from(table.querySelectorAll("tr")).filter(
            (tr) => tr.querySelectorAll("td").length > 0
        );
        for (const row of dataRows) {
            const cells = Array.from(row.querySelectorAll("td")).map(td => td.innerText.replace(/\\s+/g, " ").trim());
            out.rows.push(cells);
        }
        if (out.rows.length > 0 && out.rows[0].length > 0) {
            out.firstRowText = out.rows[0][0];
        }
        return out;
    """
    payload = driver.execute_script(js_extract) or {}
    rows = payload.get("rows") or []
    error = (payload.get("error") or "").strip()
    first_row_text = (payload.get("firstRowText") or "").strip()
    current_hash = (payload.get("currentHash") or "").strip()
    return rows, error, first_row_text, current_hash


@dataclass
class ParsedTable:
    score_values: Dict[str, Optional[float]]
    point_values: Dict[str, str]
    primateai3d_gene_threshold: Optional[float]
    missing_predictors: List[str]
    raw_rows_json: str
    n_rows: int
    first_row_text: str
    page_hash: str


def parse_other_scores_rows(rows: List[List[str]], first_row_text: str = "", page_hash: str = "") -> ParsedTable:
    score_values: Dict[str, Optional[float]] = {k: None for k in PREDICTOR_ORDER}
    point_values: Dict[str, str] = {k: "" for k in POINT_BASED_PREDICTORS}
    primateai3d_gene_threshold: Optional[float] = None

    for cells in rows:
        key: Optional[str] = None
        label_idx = -1
        for idx, cell in enumerate(cells):
            mapped = resolve_predictor_key(cell)
            if mapped:
                key = mapped
                label_idx = idx
                break
        if not key:
            continue

        score_text = cells[label_idx + 1] if label_idx + 1 < len(cells) else ""
        points_text = cells[label_idx + 2] if label_idx + 2 < len(cells) else ""

        if key == "primateai3d":
            match = re.search(
                r"([-+]?\d+(?:\.\d+)?)\s*\(gene-specific threshold:\s*([-+]?\d+(?:\.\d+)?)\)",
                score_text,
                flags=re.IGNORECASE,
            )
            if match:
                score_values[key] = float(match.group(1))
                primateai3d_gene_threshold = float(match.group(2))
            else:
                score_values[key] = parse_float_from_text(score_text)
            continue

        score_values[key] = parse_float_from_text(score_text)
        if key in point_values:
            point_values[key] = parse_point_from_text(points_text)

    missing_predictors = [k for k in PREDICTOR_ORDER if score_values.get(k) is None]
    return ParsedTable(
        score_values=score_values,
        point_values=point_values,
        primateai3d_gene_threshold=primateai3d_gene_threshold,
        missing_predictors=missing_predictors,
        raw_rows_json=json.dumps(rows, ensure_ascii=False),
        n_rows=len(rows),
        first_row_text=first_row_text,
        page_hash=page_hash,
    )


def build_wide_row(
    base: Dict[str, str],
    lookup_url: str,
    assembly: str,
    status: str,
    error: str,
    elapsed_seconds: float,
    parsed: Optional[ParsedTable],
) -> Dict[str, Any]:
    row: Dict[str, Any] = {
        "gene_symbol": base.get("gene_symbol", ""),
        "clinvar_variation_id": base.get("clinvar_variation_id", ""),
        "query_variant": base.get("query_variant", ""),
        "assembly": assembly,
        "lookup_url": lookup_url,
        "status": status,
        "error": error,
        "elapsed_seconds": round(elapsed_seconds, 3),
        "table_rows": parsed.n_rows if parsed else 0,
        "table_first_row": parsed.first_row_text if parsed else "",
        "page_hash": parsed.page_hash if parsed else "",
        "missing_predictors": ";".join(parsed.missing_predictors) if parsed else ";".join(PREDICTOR_ORDER),
        "raw_rows_json": parsed.raw_rows_json if parsed else "[]",
        "alphamissense_score": "",
        "alphamissense_points": "",
        "cadd_score": "",
        "cadd_points": "",
        "phylop_score": "",
        "phylop_points": "",
        "polyphen_max_score": "",
        "polyphen_max_points": "",
        "primateai3d_score": "",
        "primateai3d_gene_threshold": "",
        "promoterai_score": "",
        "revel_max_score": "",
        "revel_max_points": "",
        "sift_max_score": "",
        "sift_max_points": "",
    }
    if not parsed:
        return row

    for key in PREDICTOR_ORDER:
        val = parsed.score_values.get(key)
        if val is not None:
            row[f"{key}_score"] = val

    row["primateai3d_gene_threshold"] = (
        parsed.primateai3d_gene_threshold
        if parsed.primateai3d_gene_threshold is not None
        else ""
    )

    for key in POINT_BASED_PREDICTORS:
        points = parsed.point_values.get(key, "")
        row[f"{key}_points"] = points
    return row


def wide_row_to_long_rows(wide_row: Dict[str, Any]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for key in PREDICTOR_ORDER:
        score_val = wide_row.get(f"{key}_score", "")
        points_val = wide_row.get(f"{key}_points", "") if key in POINT_BASED_PREDICTORS else ""
        if score_val in ("", None):
            continue
        extra: Dict[str, Any] = {}
        if key == "primateai3d":
            extra["gene_threshold"] = wide_row.get("primateai3d_gene_threshold", "")

        out.append(
            {
                "gene_symbol": wide_row.get("gene_symbol", ""),
                "clinvar_variation_id": wide_row.get("clinvar_variation_id", ""),
                "query_variant": wide_row.get("query_variant", ""),
                "assembly": wide_row.get("assembly", ""),
                "predictor": key,
                "score": score_val,
                "points": points_val,
                "status": wide_row.get("status", ""),
                "error": wide_row.get("error", ""),
                "extra_json": json.dumps(extra, ensure_ascii=False) if extra else "",
            }
        )
    return out


def write_header_if_needed(path: Path, fieldnames: List[str]) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()


def append_csv_rows(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, Any]]) -> int:
    count = 0
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def build_predictor_coverage_summary(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    total = len(rows)
    if total == 0:
        return []
    out: List[Dict[str, Any]] = []
    for key in PREDICTOR_ORDER:
        col = f"{key}_score"
        n_scored = sum(1 for row in rows if str(row.get(col, "")).strip() != "")
        out.append(
            {
                "predictor": key,
                "n_variants": total,
                "n_scored": n_scored,
                "percent_scored": round(100.0 * n_scored / total, 4),
            }
        )
    return out


def build_gene_predictor_coverage_summary(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_gene: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        gene = (row.get("gene_symbol") or "").strip() or "UNKNOWN"
        by_gene.setdefault(gene, []).append(row)

    out: List[Dict[str, Any]] = []
    for gene in sorted(by_gene):
        subset = by_gene[gene]
        n = len(subset)
        n_ok = sum(1 for r in subset if (r.get("status") or "") == "ok")
        n_lookup_error = sum(1 for r in subset if (r.get("status") or "") == "lookup_error")
        n_exception = sum(1 for r in subset if (r.get("status") or "") == "exception")
        n_no_rows = sum(1 for r in subset if (r.get("status") or "") == "no_rows")
        n_variant_mismatch = sum(1 for r in subset if (r.get("status") or "") == "variant_mismatch")
        for key in PREDICTOR_ORDER:
            col = f"{key}_score"
            n_scored = sum(1 for row in subset if str(row.get(col, "")).strip() != "")
            out.append(
                {
                    "gene_symbol": gene,
                    "predictor": key,
                    "n_variants": n,
                    "n_ok": n_ok,
                    "n_lookup_error": n_lookup_error,
                    "n_exception": n_exception,
                    "n_no_rows": n_no_rows,
                    "n_variant_mismatch": n_variant_mismatch,
                    "n_scored": n_scored,
                    "percent_scored": round(100.0 * n_scored / n, 4) if n else 0.0,
                }
            )
    return out


def write_correlation_outputs(rows: List[Dict[str, Any]], out_dir: Path) -> None:
    try:
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        return

    score_cols = [f"{k}_score" for k in PREDICTOR_ORDER]
    df = pd.DataFrame(rows)
    if df.empty:
        return

    numeric = df[score_cols].apply(pd.to_numeric, errors="coerce")
    numeric = numeric.dropna(how="all")
    if numeric.empty:
        return

    pearson = numeric.corr(method="pearson", min_periods=15)
    spearman = numeric.corr(method="spearman", min_periods=15)

    out_dir.mkdir(parents=True, exist_ok=True)
    pearson.to_csv(out_dir / "spliceai_other_scores_correlation_pearson.csv")
    spearman.to_csv(out_dir / "spliceai_other_scores_correlation_spearman.csv")

    for matrix_name, matrix in [("pearson", pearson), ("spearman", spearman)]:
        arr = matrix.to_numpy(dtype=float)
        plt.figure(figsize=(8.2, 6.8))
        im = plt.imshow(arr, cmap="coolwarm", vmin=-1.0, vmax=1.0)
        plt.colorbar(im, fraction=0.046, pad=0.04, label=f"{matrix_name.title()} correlation")
        plt.xticks(range(len(score_cols)), score_cols, rotation=45, ha="right")
        plt.yticks(range(len(score_cols)), score_cols)

        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                val = arr[i, j]
                if np.isnan(val):
                    text = "NA"
                    color = "#333333"
                else:
                    text = f"{val:.2f}"
                    color = "white" if abs(val) >= 0.5 else "#333333"
                plt.text(j, i, text, ha="center", va="center", fontsize=8, color=color)

        plt.title(f"SpliceAI Other Scores Correlation ({matrix_name.title()})")
        plt.tight_layout()
        plt.savefig(out_dir / f"spliceai_other_scores_correlation_{matrix_name}.png", dpi=180)
        plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-scrape SpliceAI Lookup 'Other scores' for project variants."
    )
    parser.add_argument(
        "--input-csv",
        default=str(RAW_CLINICAL_ASSOC_JSONL_PATH),
        help="Input variants file (.csv or .jsonl) with gene_symbol, clinvar_variation_id, and hgvs_c/hgvs_p.",
    )
    parser.add_argument(
        "--out-csv",
        default=str(PROCESSED_DATA_DIR / "spliceai_other_scores_wide.csv"),
        help="Output wide CSV path.",
    )
    parser.add_argument(
        "--out-long-csv",
        default=str(PROCESSED_DATA_DIR / "spliceai_other_scores_long.csv"),
        help="Output long CSV path.",
    )
    parser.add_argument(
        "--summary-csv",
        default=str(PROCESSED_DATA_DIR / "spliceai_other_scores_summary.csv"),
        help="Coverage summary CSV path.",
    )
    parser.add_argument(
        "--gene-summary-csv",
        default=str(PROCESSED_DATA_DIR / "spliceai_other_scores_gene_summary.csv"),
        help="Per-gene predictor coverage CSV path.",
    )
    parser.add_argument(
        "--log-file",
        default=str(PROCESSED_DATA_DIR / "spliceai_other_scores_run.log"),
        help="Log file path. Logs are always printed to stdout and also appended here.",
    )
    parser.add_argument(
        "--log-every",
        type=int,
        default=1,
        help="Emit a progress row log every N processed variants (default: 1).",
    )
    parser.add_argument("--assembly", choices=["37", "38"], default="37")
    parser.add_argument("--distance", type=int, default=500)
    parser.add_argument("--mask", choices=["0", "1"], default="0")
    parser.add_argument("--bc", choices=["basic", "comprehensive"], default="basic")
    parser.add_argument("--ra", choices=["0", "1"], default="0")
    parser.add_argument("--browser", choices=["chrome", "firefox"], default="chrome")
    parser.add_argument("--headless", action="store_true", default=False)
    parser.add_argument("--no-headless", dest="headless", action="store_false")
    parser.add_argument("--timeout-seconds", type=int, default=120)
    parser.add_argument("--sleep-seconds", type=float, default=3.0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--overwrite", action="store_true", default=False)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_every = max(1, int(args.log_every))

    input_path = Path(args.input_csv)
    out_csv = Path(args.out_csv)
    out_long_csv = Path(args.out_long_csv)
    summary_csv = Path(args.summary_csv)
    gene_summary_csv = Path(args.gene_summary_csv)
    log_file = Path(args.log_file)
    out_dir = out_csv.parent

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_handle = log_file.open("a", encoding="utf-8")
    interrupted = False
    started = time.time()

    try:
        log_line("-" * 72, log_handle)
        log_line("Starting SpliceAI Other Scores scrape run.", log_handle)
        log_line(f"Input file: {input_path}", log_handle)
        log_line(f"Output (wide): {out_csv}", log_handle)
        log_line(f"Output (long): {out_long_csv}", log_handle)
        log_line(f"Summary CSVs: {summary_csv}, {gene_summary_csv}", log_handle)
        log_line(
            f"Settings: assembly={args.assembly}, browser={args.browser}, headless={args.headless}, "
            f"timeout={args.timeout_seconds}s, sleep={args.sleep_seconds}s, log_every={log_every}",
            log_handle,
        )

        if not input_path.exists():
            log_line(f"[ERROR] Input file not found: {input_path}", log_handle)
            return 2

        if args.overwrite:
            for p in (out_csv, out_long_csv, summary_csv, gene_summary_csv):
                if p.exists():
                    p.unlink()
            log_line("[INFO] Overwrite mode: existing output files removed.", log_handle)

        variants = load_variants(input_path, limit=args.limit)
        if not variants:
            log_line(f"[ERROR] No rows found in input file: {input_path}", log_handle)
            return 2

        existing_status_map = read_existing_status_map(out_csv) if args.resume and out_csv.exists() else {}
        historical_lookup_error_keys = (
            read_keys_with_status(out_csv, "lookup_error") if args.resume and out_csv.exists() else set()
        )
        resume_enabled = bool(existing_status_map)
        if resume_enabled:
            completed_count = sum(1 for s in existing_status_map.values() if not should_retry_status(s))
            retryable_count = len(existing_status_map) - completed_count
            log_line(
                f"[INFO] Resume mode: found {len(existing_status_map)} prior rows "
                f"(completed={completed_count}, retryable={retryable_count}) in {out_csv}",
                log_handle,
            )
            if historical_lookup_error_keys:
                log_line(
                    f"[INFO] Resume mode: found {len(historical_lookup_error_keys)} unique keys with "
                    f"historical status=lookup_error. These keys will be retried each run.",
                    log_handle,
                )

        def should_process_key(key: Tuple[str, str, str]) -> bool:
            if not args.resume:
                return True
            if historical_lookup_error_keys:
                # When lookup_error history exists, always retry those keys on every run.
                # Also process keys not seen before so newly added variants are included.
                return key in historical_lookup_error_keys or key not in existing_status_map
            prev_status = existing_status_map.get(key, "")
            return should_retry_status(prev_status)

        wide_fields = list(
            build_wide_row(
                base={"gene_symbol": "", "clinvar_variation_id": "", "query_variant": ""},
                lookup_url="",
                assembly=args.assembly,
                status="",
                error="",
                elapsed_seconds=0.0,
                parsed=None,
            ).keys()
        )
        long_fields = [
            "gene_symbol",
            "clinvar_variation_id",
            "query_variant",
            "assembly",
            "predictor",
            "score",
            "points",
            "status",
            "error",
            "extra_json",
        ]

        write_header_if_needed(out_csv, wide_fields)
        write_header_if_needed(out_long_csv, long_fields)

        remaining = 0
        for base_row in variants:
            key = build_variant_key(base_row)
            if should_process_key(key):
                remaining += 1

        written_wide = 0
        written_long = 0
        processed = 0
        status_counts: Dict[str, int] = {}

        if remaining == 0:
            log_line("[INFO] No remaining variants to scrape. Rebuilding summaries from existing output.", log_handle)
        else:
            log_line(f"[INFO] Variants loaded: total={len(variants)}, to_process={remaining}", log_handle)
            driver = create_driver(args.browser, args.headless)
            driver.set_page_load_timeout(args.timeout_seconds)
            try:
                load_lookup_page(driver, timeout_seconds=args.timeout_seconds)
                for base_row in variants:
                    key = build_variant_key(base_row)
                    if not should_process_key(key):
                        continue

                    processed += 1
                    query_variant = (base_row.get("query_variant") or "").strip()
                    if not query_variant or query_variant == "-":
                        status = "skipped_no_variant"
                        wide_row = build_wide_row(
                            base=base_row,
                            lookup_url="",
                            assembly=args.assembly,
                            status=status,
                            error="No HGVS query string available in input row.",
                            elapsed_seconds=0.0,
                            parsed=None,
                        )
                        append_csv_rows(out_csv, wide_fields, [wide_row])
                        written_wide += 1
                        status_counts[status] = status_counts.get(status, 0) + 1
                        if processed == 1 or processed % log_every == 0:
                            log_line(
                                f"[{processed}/{remaining}] {base_row.get('gene_symbol','')}:{base_row.get('clinvar_variation_id','')} "
                                f"status={status}",
                                log_handle,
                            )
                        continue

                    lookup_url = build_lookup_url(
                        variant=query_variant,
                        assembly=args.assembly,
                        bc=args.bc,
                        distance=args.distance,
                        mask=int(args.mask),
                        ra=int(args.ra),
                    )

                    t0 = time.time()
                    status = "ok"
                    error = ""
                    parsed: Optional[ParsedTable] = None
                    try:
                        submit_result = submit_lookup_request(
                            driver,
                            query_variant=query_variant,
                            assembly=args.assembly,
                            bc=args.bc,
                            distance=args.distance,
                            mask=args.mask,
                        )
                        if not submit_result.get("ok"):
                            raise RuntimeError(f"submit_failed:{submit_result.get('error')}")

                        wait_for_submit_cycle(driver, timeout_seconds=args.timeout_seconds)
                        rows, error_text, first_row_text, current_hash = get_other_scores_rows_and_error(driver)

                        # Some page errors can appear before the "Other scores" table finishes rendering.
                        if error_text and not rows:
                            settle_seconds = 3.0 if is_rate_limit_error_text(error_text) else 1.0
                            settle_deadline = time.time() + settle_seconds
                            while time.time() < settle_deadline:
                                time.sleep(0.25)
                                retry_rows, retry_error_text, retry_first_row_text, retry_hash = (
                                    get_other_scores_rows_and_error(driver)
                                )
                                if retry_rows:
                                    rows = retry_rows
                                    error_text = retry_error_text
                                    first_row_text = retry_first_row_text
                                    current_hash = retry_hash
                                    break

                        parsed = parse_other_scores_rows(
                            rows,
                            first_row_text=first_row_text,
                            page_hash=current_hash,
                        )

                        if parsed.n_rows == 0:
                            if error_text:
                                status = "lookup_error"
                                error = error_text
                            else:
                                status = "no_rows"
                                error = "Other scores table had no rows."
                        else:
                            if error_text:
                                normalized_error_text = normalize_space(error_text)
                                err_preview = normalized_error_text[:220]
                                suffix = "..." if len(normalized_error_text) > 220 else ""
                                if is_rate_limit_error_text(error_text):
                                    log_line(
                                        f"[WARN] Non-fatal rate-limit/error text with parsed Other scores rows: "
                                        f"{err_preview}{suffix}",
                                        log_handle,
                                    )
                                else:
                                    log_line(
                                        f"[WARN] Non-fatal error text with parsed Other scores rows: "
                                        f"{err_preview}{suffix}",
                                        log_handle,
                                    )

                        if parsed.n_rows > 0 and parsed.first_row_text:
                            # Guard against stale table parsing: first row should reference the current query variant.
                            if query_variant.lower() not in parsed.first_row_text.lower():
                                status = "variant_mismatch"
                                error = (
                                    f"Parsed table appears to reference a different variant. "
                                    f"query={query_variant}; first_row={parsed.first_row_text}"
                                )
                    except Exception as exc:
                        status = "exception"
                        error = str(exc)

                    elapsed = time.time() - t0
                    wide_row = build_wide_row(
                        base=base_row,
                        lookup_url=lookup_url,
                        assembly=args.assembly,
                        status=status,
                        error=error,
                        elapsed_seconds=elapsed,
                        parsed=parsed,
                    )
                    append_csv_rows(out_csv, wide_fields, [wide_row])
                    written_wide += 1

                    long_rows = wide_row_to_long_rows(wide_row)
                    if long_rows:
                        written_long += append_csv_rows(out_long_csv, long_fields, long_rows)

                    status_counts[status] = status_counts.get(status, 0) + 1
                    if processed == 1 or processed % log_every == 0 or status != "ok":
                        log_line(
                            f"[{processed}/{remaining}] {base_row.get('gene_symbol','')}:{base_row.get('clinvar_variation_id','')} "
                            f"status={status} rows={wide_row.get('table_rows', 0)} elapsed={elapsed:.1f}s",
                            log_handle,
                        )
                    if args.sleep_seconds > 0:
                        time.sleep(args.sleep_seconds)
            except KeyboardInterrupt:
                interrupted = True
                log_line("[WARN] Interrupted by user. Writing partial summaries from collected rows.", log_handle)
            finally:
                driver.quit()

        all_rows = load_latest_wide_rows(out_csv)

        summary_rows = build_predictor_coverage_summary(all_rows)
        gene_summary_rows = build_gene_predictor_coverage_summary(all_rows)
        if summary_rows:
            summary_fields = list(summary_rows[0].keys())
            with summary_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=summary_fields)
                writer.writeheader()
                writer.writerows(summary_rows)
        if gene_summary_rows:
            gene_summary_fields = list(gene_summary_rows[0].keys())
            with gene_summary_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=gene_summary_fields)
                writer.writeheader()
                writer.writerows(gene_summary_rows)

        write_correlation_outputs(all_rows, out_dir)

        total_elapsed = time.time() - started
        if status_counts:
            counts_text = ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items()))
            log_line(f"[OK] Status counts this run: {counts_text}", log_handle)
        log_line(f"[OK] Wrote {out_csv} (+{written_wide} rows this run)", log_handle)
        log_line(f"[OK] Wrote {out_long_csv} (+{written_long} rows this run)", log_handle)
        if summary_rows:
            log_line(f"[OK] Wrote {summary_csv}", log_handle)
        if gene_summary_rows:
            log_line(f"[OK] Wrote {gene_summary_csv}", log_handle)
        if summary_rows or gene_summary_rows:
            log_line(f"[OK] Wrote correlation tables/plots under {out_dir}", log_handle)
        log_line(f"[OK] Finished in {total_elapsed:.1f}s", log_handle)
        return 130 if interrupted else 0
    finally:
        log_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
