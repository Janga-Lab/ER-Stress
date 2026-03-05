#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from config import (
    PREPARED_FOR_PARSER_JSON_PATH,
    RAW_JSONL_PATH,
    ensure_data_directories,
)


def clean_optional_text(value: Any) -> Optional[str]:
    if not isinstance(value, str):
        return None
    text = value.strip()
    return text or None


def safe_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
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


def derive_labels_from_uniprot_features(gene: str, uniprot: Dict[str, Any]) -> List[str]:
    labels: set[str] = set()

    transmembrane = uniprot.get("transmembrane") or []
    topological = uniprot.get("topological_domains") or []
    domains = uniprot.get("domains_regions") or []
    signal_peptides = uniprot.get("signal_peptides") or []
    chains = uniprot.get("chains") or []
    propeptides = uniprot.get("propeptides") or []
    disulfide_bonds = uniprot.get("disulfide_bonds") or []

    if isinstance(transmembrane, list) and transmembrane:
        labels.add("tm")

    for topo in topological if isinstance(topological, list) else []:
        if not isinstance(topo, dict):
            continue
        desc = (clean_optional_text(topo.get("description")) or "").lower()
        if any(x in desc for x in ("luminal", "lumen", "extracellular")):
            labels.add("luminal")
        if "cyto" in desc:
            labels.add("cytosolic")

    for dom in domains if isinstance(domains, list) else []:
        if not isinstance(dom, dict):
            continue
        desc = (clean_optional_text(dom.get("description")) or "").lower()
        if "kinase" in desc:
            labels.add("kinase_domain")
        if any(x in desc for x in ("cdgsh", "fe-s", "2fe-2s", "iron-sulfur")):
            labels.add("fes_binding_domain")

    if gene == "INS":
        if isinstance(signal_peptides, list) and signal_peptides:
            labels.add("signal_peptide")
        if isinstance(disulfide_bonds, list) and disulfide_bonds:
            labels.add("disulfide_cysteine")
        if (isinstance(propeptides, list) and propeptides) or (isinstance(chains, list) and chains):
            labels.add("cleavage_site")
        for chain in chains if isinstance(chains, list) else []:
            if not isinstance(chain, dict):
                continue
            desc = (clean_optional_text(chain.get("description")) or "").lower()
            if "a chain" in desc:
                labels.add("a_chain")
            if "b chain" in desc:
                labels.add("b_chain")

    return sorted(labels)


def load_gene_records_from_raw_jsonl(raw_jsonl: Path) -> Dict[str, Dict[str, Any]]:
    gene_records: Dict[str, Dict[str, Any]] = {}
    with raw_jsonl.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            obj = json.loads(text)
            if obj.get("record_type") != "gene":
                continue
            gene = (clean_optional_text(obj.get("gene_symbol")) or "").upper()
            if gene:
                gene_records[gene] = obj
    return gene_records


def build_prepared_payload(
    gene_records: Dict[str, Dict[str, Any]],
    raw_jsonl: Path,
) -> Dict[str, Any]:
    critical_map: Dict[str, List[str]] = {}
    accession_map: Dict[str, str] = {}
    feature_counts: Dict[str, Dict[str, int]] = {}

    for gene in sorted(gene_records):
        record = gene_records[gene]
        uniprot = record.get("uniprot") or {}
        critical_map[gene] = derive_labels_from_uniprot_features(gene, uniprot)
        accession_map[gene] = clean_optional_text(record.get("uniprot_accession")) or ""
        feature_counts[gene] = {
            "transmembrane": len(feature_ranges(uniprot.get("transmembrane"))),
            "topological_domains": len(feature_ranges(uniprot.get("topological_domains"))),
            "domains_regions": len(feature_ranges(uniprot.get("domains_regions"))),
            "signal_peptides": len(feature_ranges(uniprot.get("signal_peptides"))),
            "chains": len(feature_ranges(uniprot.get("chains"))),
            "propeptides": len(feature_ranges(uniprot.get("propeptides"))),
            "disulfide_bonds": len(feature_ranges(uniprot.get("disulfide_bonds"))),
        }

    return {
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_raw_jsonl": str(raw_jsonl),
        "n_genes": len(gene_records),
        "genes": sorted(gene_records.keys()),
        "critical_region_map": critical_map,
        "uniprot_accessions": accession_map,
        "uniprot_feature_counts": feature_counts,
    }


def write_prepared_payload(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def load_prepared_payload(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(obj, dict):
        return None
    critical = obj.get("critical_region_map")
    if not isinstance(critical, dict):
        return None
    return obj


def load_or_prepare_payload(
    raw_jsonl: Path,
    prepared_json: Path,
    *,
    gene_records: Optional[Dict[str, Dict[str, Any]]] = None,
    force_refresh: bool = False,
) -> Dict[str, Any]:
    if not force_refresh:
        existing = load_prepared_payload(prepared_json)
        if existing is not None:
            return existing

    if gene_records is None:
        gene_records = load_gene_records_from_raw_jsonl(raw_jsonl)
    payload = build_prepared_payload(gene_records, raw_jsonl)
    write_prepared_payload(prepared_json, payload)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare parser inputs (dynamic UniProt critical-region map)")
    parser.add_argument("--raw-jsonl", default=str(RAW_JSONL_PATH))
    parser.add_argument("--out-json", default=str(PREPARED_FOR_PARSER_JSON_PATH))
    parser.add_argument("--force-refresh", action="store_true", default=False)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_data_directories()

    raw_jsonl = Path(args.raw_jsonl)
    if not raw_jsonl.exists():
        raise SystemExit(f"Raw JSONL not found: {raw_jsonl}")

    out_json = Path(args.out_json)
    payload = load_or_prepare_payload(
        raw_jsonl=raw_jsonl,
        prepared_json=out_json,
        gene_records=None,
        force_refresh=args.force_refresh,
    )
    print(f"[OK] Wrote {out_json}")
    print(f"[OK] critical_region_map: {json.dumps(payload.get('critical_region_map') or {}, ensure_ascii=False)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

