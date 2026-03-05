#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from config import (
    CURATED_GENE_FEATURES_JSON_PATH,
    CURATED_VARIANTS_CSV_PATH,
    CURATION_VALIDATION_REPORT_PATH,
    PREPARED_FOR_PARSER_JSON_PATH,
    RAW_JSONL_PATH,
    ensure_data_directories,
)
from prepare_for_parser import load_or_prepare_payload

PLACEHOLDER_TEXT_VALUES = {
    "",
    "not provided",
    "none provided",
    "not specified",
    "unknown",
    "n/a",
    "na",
}

LOF_TERMS = {
    "frameshift_variant",
    "stop_gained",
    "start_lost",
    "transcript_ablation",
    "feature_truncation",
    "exon_loss_variant",
}

SPLICE_TERMS = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}

INFRAME_TERMS = {"inframe_insertion", "inframe_deletion", "protein_altering_variant"}
MISSENSE_TERMS = {"missense_variant"}
SYNONYMOUS_TERMS = {"synonymous_variant", "stop_retained_variant"}

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


@dataclass
class GeneFeatureContext:
    transmembrane: List[Tuple[int, int, str]]
    topological: List[Tuple[int, int, str]]
    domains_regions: List[Tuple[int, int, str]]
    signal_peptides: List[Tuple[int, int, str]]
    chains: List[Tuple[int, int, str]]
    propeptides: List[Tuple[int, int, str]]
    disulfide_bonds: List[Tuple[int, int, str]]


def clean_optional_text(value: Any) -> Optional[str]:
    if not isinstance(value, str):
        return None
    text = value.strip()
    if not text:
        return None
    if text.lower() in PLACEHOLDER_TEXT_VALUES:
        return None
    return text


def safe_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if math.isnan(value):
            return None
        return int(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            return int(float(text))
        except ValueError:
            return None
    return None


def canonical_variation_id(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def parse_consequence_terms(rec: Dict[str, Any]) -> List[str]:
    terms = ((rec.get("vep") or {}).get("consequence_terms") or [])
    if isinstance(terms, str):
        terms = [terms]
    if not isinstance(terms, list):
        return []
    return [t for t in terms if isinstance(t, str) and t.strip()]


def feature_ranges(uniprot_feature_list: Any) -> List[Tuple[int, int, str]]:
    out: List[Tuple[int, int, str]] = []
    if not isinstance(uniprot_feature_list, list):
        return out

    for item in uniprot_feature_list:
        if not isinstance(item, dict):
            continue
        start = safe_int(item.get("start"))
        end = safe_int(item.get("end"))
        desc = clean_optional_text(item.get("description")) or ""
        if start is None and end is None:
            continue
        if start is None:
            start = end
        if end is None:
            end = start
        if start is None or end is None:
            continue
        if end < start:
            start, end = end, start
        out.append((start, end, desc))

    return out


def build_gene_feature_context(gene_record: Dict[str, Any]) -> GeneFeatureContext:
    up = gene_record.get("uniprot") or {}
    return GeneFeatureContext(
        transmembrane=feature_ranges(up.get("transmembrane")),
        topological=feature_ranges(up.get("topological_domains")),
        domains_regions=feature_ranges(up.get("domains_regions")),
        signal_peptides=feature_ranges(up.get("signal_peptides")),
        chains=feature_ranges(up.get("chains")),
        propeptides=feature_ranges(up.get("propeptides")),
        disulfide_bonds=feature_ranges(up.get("disulfide_bonds")),
    )


def find_overlapping_feature(position: int, ranges: List[Tuple[int, int, str]]) -> Optional[Tuple[int, int, str]]:
    for start, end, desc in ranges:
        if start <= position <= end:
            return (start, end, desc)
    return None


def extract_protein_position(rec: Dict[str, Any]) -> Optional[int]:
    vep = rec.get("vep") or {}
    start = safe_int(vep.get("protein_start"))
    if start is not None:
        return start
    return safe_int(vep.get("protein_end"))


def classify_variant_bucket(rec: Dict[str, Any]) -> str:
    terms = set(parse_consequence_terms(rec))
    variant_type = (clean_optional_text(rec.get("type")) or "").lower()

    if terms & LOF_TERMS:
        return "lof_truncating"
    if terms & SPLICE_TERMS:
        return "splice"
    if terms & INFRAME_TERMS:
        return "inframe_indel"
    if terms & MISSENSE_TERMS:
        return "missense"
    if terms & SYNONYMOUS_TERMS:
        return "synonymous"

    structural_tokens = (
        "copy number",
        "deletion",
        "duplication",
        "inversion",
        "translocation",
        "structural",
        "cnv",
    )
    if any(tok in variant_type for tok in structural_tokens):
        return "cnv_sv"

    return "other"


def assign_ins_region_bucket(position: Optional[int], ctx: GeneFeatureContext) -> str:
    if position is None:
        return "non_protein_or_unknown"

    if find_overlapping_feature(position, ctx.signal_peptides):
        return "signal_peptide"

    if find_overlapping_feature(position, ctx.disulfide_bonds):
        return "disulfide_cysteine"

    for start, end, desc in ctx.chains:
        if start <= position <= end:
            d = desc.lower()
            if "b chain" in d:
                return "b_chain"
            if "a chain" in d:
                return "a_chain"

    for start, end, desc in ctx.propeptides:
        if start <= position <= end:
            d = desc.lower()
            if "c-peptide" in d or "connecting peptide" in d:
                return "c_peptide"

    cleavage_positions: set[int] = set()
    for start, end, _ in (ctx.chains + ctx.propeptides):
        for pos in (start - 1, start, end, end + 1):
            if pos >= 1:
                cleavage_positions.add(pos)
    if position in cleavage_positions:
        return "cleavage_site"

    return "other_protein_region"


def assign_region_bucket(gene: str, rec: Dict[str, Any], ctx: GeneFeatureContext) -> str:
    position = extract_protein_position(rec)
    if gene == "INS":
        return assign_ins_region_bucket(position, ctx)

    if position is None:
        return "non_protein_or_unknown"

    tm_hit = find_overlapping_feature(position, ctx.transmembrane)
    topo_hit = find_overlapping_feature(position, ctx.topological)
    dom_hit = find_overlapping_feature(position, ctx.domains_regions)

    if gene == "WFS1":
        if tm_hit:
            return "tm"
        if topo_hit:
            desc = topo_hit[2].lower()
            if "luminal" in desc or "lumen" in desc or "extracellular" in desc:
                return "luminal"
            if "cyto" in desc:
                return "cytosolic"
            return "topological_non_tm"
        return "non_tm"

    if gene == "EIF2AK3":
        if dom_hit and any(x in dom_hit[2].lower() for x in ("kinase", "protein kinase")):
            return "kinase_domain"
        if tm_hit:
            return "tm"
        if topo_hit:
            desc = topo_hit[2].lower()
            if "luminal" in desc or "lumen" in desc:
                return "luminal"
            if "cyto" in desc:
                return "cytosolic"
            return "topological_other"
        if dom_hit:
            return "domain_or_region"
        return "other_protein_region"

    if gene == "CISD2":
        if dom_hit and any(x in dom_hit[2].lower() for x in ("cdgsh", "fe-s", "2fe-2s", "iron-sulfur")):
            return "fes_binding_domain"
        if tm_hit:
            return "tm"
        if dom_hit:
            return "domain_or_region"
        return "other_protein_region"

    if gene == "IER3IP1":
        if tm_hit:
            return "tm"
        if dom_hit:
            return "domain_or_region"
        return "non_tm"

    if dom_hit:
        return "domain_or_region"
    if tm_hit:
        return "tm"
    if topo_hit:
        return "topological_other"
    return "other_protein_region"


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


def load_raw_jsonl(path: Path) -> Tuple[Dict[str, Dict[str, Any]], List[Dict[str, Any]]]:
    gene_records: Dict[str, Dict[str, Any]] = {}
    variants: List[Dict[str, Any]] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, 1):
            text = line.strip()
            if not text:
                continue
            try:
                obj = json.loads(text)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON at line {line_num}: {exc}") from exc

            record_type = obj.get("record_type")
            if record_type == "gene":
                gene_symbol = (clean_optional_text(obj.get("gene_symbol")) or "").upper()
                if gene_symbol:
                    gene_records[gene_symbol] = obj
            elif record_type == "variant":
                variants.append(obj)

    return gene_records, variants


def build_curated_rows(
    gene_records: Dict[str, Dict[str, Any]],
    variants: List[Dict[str, Any]],
    critical_region_map: Dict[str, set[str]],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    feature_context = {
        gene: build_gene_feature_context(gene_rec)
        for gene, gene_rec in gene_records.items()
    }

    seen_keys = set()
    curated_rows: List[Dict[str, Any]] = []

    validation = {
        "input_variant_rows": len(variants),
        "output_curated_rows": 0,
        "duplicates_skipped": 0,
        "missing_gene_symbol": 0,
        "missing_variation_id": 0,
        "rows_with_missing_hgvs": 0,
        "rows_without_vep_terms": 0,
    }

    for rec in variants:
        gene = (clean_optional_text(rec.get("gene_symbol")) or "").upper()
        if not gene:
            validation["missing_gene_symbol"] += 1
            continue

        variation_id = canonical_variation_id(rec.get("clinvar_variation_id"))
        if not variation_id:
            validation["missing_variation_id"] += 1
            continue

        key = (gene, variation_id)
        if key in seen_keys:
            validation["duplicates_skipped"] += 1
            continue
        seen_keys.add(key)

        hgvs_c = clean_optional_text(rec.get("hgvs_c")) or ""
        hgvs_p = clean_optional_text(rec.get("hgvs_p")) or ""
        if not hgvs_c and not hgvs_p:
            validation["rows_with_missing_hgvs"] += 1

        terms = parse_consequence_terms(rec)
        if not terms:
            validation["rows_without_vep_terms"] += 1

        variant_bucket = classify_variant_bucket(rec)
        ctx = feature_context.get(gene, GeneFeatureContext([], [], [], [], [], [], []))
        region_bucket = assign_region_bucket(gene, rec, ctx)
        is_critical = int(region_bucket in critical_region_map.get(gene, set()))

        vep = rec.get("vep") or {}
        clinvar = rec.get("clinvar") or {}

        curated_rows.append(
            {
                "gene_symbol": gene,
                "clinvar_variation_id": variation_id,
                "hgvs_c": hgvs_c,
                "hgvs_p": hgvs_p,
                "dbsnp_rsid": clean_optional_text(rec.get("dbsnp_rsid")) or "",
                "variant_type": clean_optional_text(rec.get("type")) or "",
                "consequence_terms": ";".join(terms),
                "impact": clean_optional_text(vep.get("impact")) or "",
                "transcript_id": clean_optional_text(vep.get("transcript_id")) or "",
                "protein_start": safe_int(vep.get("protein_start")) or "",
                "protein_end": safe_int(vep.get("protein_end")) or "",
                "variant_bucket": variant_bucket,
                "region_bucket": region_bucket,
                "is_critical_region": is_critical,
                "clinvar_clinical_significance": clean_optional_text(rec.get("clinvar_clinical_significance")) or "",
                "pathogenic_binary": (
                    clinvar_binary_label(rec.get("clinvar_clinical_significance"))
                    if clinvar_binary_label(rec.get("clinvar_clinical_significance")) is not None
                    else ""
                ),
                "clinvar_review_status": clean_optional_text(clinvar.get("germline_review_status")) or "",
                "clinvar_last_evaluated": clean_optional_text(clinvar.get("germline_last_evaluated")) or "",
                "clinvar_conditions": ";".join([x for x in (rec.get("clinvar_conditions") or []) if clean_optional_text(x)]),
                "clinvar_pubmed_count": len(rec.get("clinvar_pubmed_ids") or []),
                "pubmed_count": len(rec.get("pubmed") or []),
                "uniprot_accession": clean_optional_text(rec.get("uniprot_accession")) or "",
            }
        )

    validation["output_curated_rows"] = len(curated_rows)
    return curated_rows, validation


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    columns = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_and_curate(
    raw_jsonl: Path,
    prepared_for_parser_json: Path,
    curated_variants_csv: Path,
    curated_gene_features_json: Path,
    validation_report_json: Path,
) -> None:
    gene_records, variants = load_raw_jsonl(raw_jsonl)
    prepared_existed = prepared_for_parser_json.exists()
    prepared_payload = load_or_prepare_payload(
        raw_jsonl=raw_jsonl,
        prepared_json=prepared_for_parser_json,
        gene_records=gene_records,
        force_refresh=False,
    )
    critical_map_raw = prepared_payload.get("critical_region_map") or {}
    critical_map: Dict[str, set[str]] = {}
    if isinstance(critical_map_raw, dict):
        for gene, labels in critical_map_raw.items():
            if not isinstance(gene, str):
                continue
            if not isinstance(labels, list):
                labels = []
            critical_map[gene.upper()] = {
                label.strip().lower()
                for label in labels
                if isinstance(label, str) and label.strip()
            }
    for gene in gene_records:
        critical_map.setdefault(gene, set())

    curated_rows, validation = build_curated_rows(gene_records, variants, critical_map)
    write_csv(curated_variants_csv, curated_rows)
    curated_gene_features_json.write_text(json.dumps(gene_records, indent=2, ensure_ascii=False), encoding="utf-8")
    validation_report_json.write_text(json.dumps(validation, indent=2), encoding="utf-8")

    print(f"[OK] Wrote {curated_variants_csv} ({len(curated_rows)} rows)")
    print(f"[OK] Wrote {curated_gene_features_json} ({len(gene_records)} genes)")
    action = "Loaded existing" if prepared_existed else "Created"
    print(f"[OK] {action} {prepared_for_parser_json}")
    print(f"[OK] UniProt dynamic critical map: {json.dumps({k: sorted(v) for k, v in critical_map.items()}, ensure_ascii=False)}")
    print(f"[OK] Wrote {validation_report_json}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parse, clean, validate, and curate raw imported data")
    parser.add_argument("--raw-jsonl", default=str(RAW_JSONL_PATH))
    parser.add_argument("--prepared-for-parser-json", default=str(PREPARED_FOR_PARSER_JSON_PATH))
    parser.add_argument("--out-curated-variants", default=str(CURATED_VARIANTS_CSV_PATH))
    parser.add_argument("--out-gene-features", default=str(CURATED_GENE_FEATURES_JSON_PATH))
    parser.add_argument("--out-validation-report", default=str(CURATION_VALIDATION_REPORT_PATH))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_data_directories()

    raw_jsonl = Path(args.raw_jsonl)
    if not raw_jsonl.exists():
        raise SystemExit(f"Raw JSONL not found: {raw_jsonl}")

    parse_and_curate(
        raw_jsonl=raw_jsonl,
        prepared_for_parser_json=Path(args.prepared_for_parser_json),
        curated_variants_csv=Path(args.out_curated_variants),
        curated_gene_features_json=Path(args.out_gene_features),
        validation_report_json=Path(args.out_validation_report),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
