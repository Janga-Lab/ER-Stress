#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from config import (
    CLINICALTABLES_VARIANTS_SEARCH_URL,
    ENSEMBL_VEP_GRCH37_BASE_URL,
    ENSEMBL_VEP_GRCH38_BASE_URL,
    NCBI_EUTILS_BASE_URL,
    NCBI_VARIATION_BASE_URL,
    UNIPROT_ENTRY_BASE_URL,
    UNIPROT_SEARCH_URL,
    DEFAULT_GENES,
    PipelineConfig,
    RAW_CLINICAL_ASSOC_JSONL_PATH,
    RAW_JSONL_PATH,
    ensure_data_directories,
)

PLACEHOLDER_TEXT_VALUES = {
    "",
    "not provided",
    "none provided",
    "not specified",
    "unknown",
    "n/a",
    "na",
}

RELEVANT_CLINVAR_TERMS = {
    "missense_variant",
    "frameshift_variant",
    "stop_gained",
    "start_lost",
    "stop_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
}

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


def clean_optional_text(value: Any) -> Optional[str]:
    if not isinstance(value, str):
        return None
    text = value.strip()
    if not text:
        return None
    if text.lower() in PLACEHOLDER_TEXT_VALUES:
        return None
    return text


def value_has_nonempty_text(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return clean_optional_text(value) is not None
    if isinstance(value, list):
        return any(value_has_nonempty_text(item) for item in value)
    if isinstance(value, dict):
        return any(value_has_nonempty_text(item) for item in value.values())
    return bool(str(value).strip())


def pick_best_hgvs(record: Dict[str, Any]) -> Optional[str]:
    hgvs_c = record.get("HGVS_c")
    if isinstance(hgvs_c, str) and hgvs_c.strip():
        return hgvs_c.strip()
    hgvs_p = record.get("HGVS_p")
    if isinstance(hgvs_p, str) and hgvs_p.strip():
        return hgvs_p.strip()
    return None


def pick_transcript_consequence_for_input(vep_record: Dict[str, Any], input_hgvs: str) -> Optional[Dict[str, Any]]:
    for consequence in vep_record.get("transcript_consequences", []) or []:
        if consequence.get("hgvsc") == input_hgvs or consequence.get("hgvsp") == input_hgvs:
            return consequence
    consequences = vep_record.get("transcript_consequences") or []
    return consequences[0] if consequences else None


def summarize_vep_for_input(vep_record: Dict[str, Any], input_hgvs: str) -> Dict[str, Any]:
    consequence = pick_transcript_consequence_for_input(vep_record, input_hgvs)
    if not consequence:
        return {"input": vep_record.get("input")}

    return {
        "input": vep_record.get("input"),
        "consequence_terms": consequence.get("consequence_terms"),
        "impact": consequence.get("impact"),
        "gene_symbol": consequence.get("gene_symbol"),
        "transcript_id": consequence.get("transcript_id"),
        "hgvsc": consequence.get("hgvsc"),
        "hgvsp": consequence.get("hgvsp"),
        "protein_start": consequence.get("protein_start"),
        "protein_end": consequence.get("protein_end"),
        "amino_acids": consequence.get("amino_acids"),
        "codons": consequence.get("codons"),
    }


def should_fetch_clinvar(vep_summary: Dict[str, Any]) -> bool:
    terms = vep_summary.get("consequence_terms") or []
    return any(term in RELEVANT_CLINVAR_TERMS for term in terms)


class DatabaseImporter:
    def __init__(self, cfg: PipelineConfig):
        self.cfg = cfg
        self.session = self._make_retry_session()
        self._clinvar_cache: Dict[str, Dict[str, Any]] = {}

    def _make_retry_session(self) -> requests.Session:
        session = requests.Session()
        retry = Retry(
            total=5,
            backoff_factor=0.6,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(["GET", "POST"]),
            raise_on_status=False,
        )
        session.mount("https://", HTTPAdapter(max_retries=retry))
        session.headers.update({"User-Agent": "er-stress-variant-harvester/2.0"})
        return session

    def _http_get_json(
        self,
        url: str,
        *,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: int = 30,
    ) -> Any:
        response = self.session.get(url, params=params, headers=headers, timeout=timeout)
        response.raise_for_status()
        return response.json()

    def _ncbi_common_params(self) -> Dict[str, str]:
        params: Dict[str, str] = {"tool": self.cfg.ncbi_tool}
        if self.cfg.ncbi_email:
            params["email"] = self.cfg.ncbi_email
        if self.cfg.ncbi_api_key:
            params["api_key"] = self.cfg.ncbi_api_key
        return params

    def fetch_clinicaltables_variants_for_gene(self, gene_symbol: str) -> List[Dict[str, Any]]:
        ef_fields = [
            "VariationID",
            "GeneSymbol",
            "HGVS_c",
            "HGVS_p",
            "AminoAcidChange",
            "dbSNP",
            "phenotypes",
            "phenotype",
            "Chromosome",
            "Start",
            "Stop",
            "Type",
            "ReferenceAllele",
            "AlternateAllele",
            "ChromosomeAccession",
            "GenomicLocation",
        ]
        results: List[Dict[str, Any]] = []
        count = min(int(self.cfg.max_variants_per_gene), 500)
        offset = 0

        while True:
            params = {
                "terms": gene_symbol,
                "q": f"GeneSymbol:{gene_symbol}",
                "count": count,
                "offset": offset,
                "ef": ",".join(ef_fields),
                "df": "VariationID,Name",
            }
            payload = self._http_get_json(CLINICALTABLES_VARIANTS_SEARCH_URL, params=params, timeout=30)
            if not isinstance(payload, list) or len(payload) < 3:
                break

            total = int(payload[0] or 0)
            codes = payload[1] or []
            extra = payload[2] or {}
            if not codes:
                break

            for idx, variation_id in enumerate(codes):
                rec = {"gene_symbol": gene_symbol, "clinvar_variation_id": str(variation_id)}
                for field in ef_fields:
                    values = extra.get(field)
                    if isinstance(values, list) and idx < len(values):
                        rec[field] = values[idx]
                results.append(rec)

            offset += len(codes)
            if offset >= total or offset >= self.cfg.max_variants_per_gene:
                break
            time.sleep(self.cfg.sleep_seconds)

        return results

    def fetch_vep_for_hgvs(self, hgvs: str) -> Optional[Dict[str, Any]]:
        base_url = ENSEMBL_VEP_GRCH37_BASE_URL if self.cfg.assembly.lower() == "grch37" else ENSEMBL_VEP_GRCH38_BASE_URL
        url = f"{base_url}/{quote(hgvs, safe='')}"
        try:
            data = self._http_get_json(url, headers={"Accept": "application/json"}, timeout=45)
        except requests.HTTPError:
            return None

        if isinstance(data, list) and data:
            return data[0]
        if isinstance(data, dict):
            return data
        return None

    def fetch_uniprot_primary_accession(self, gene_symbol: str, taxon_id: int = 9606) -> Optional[str]:
        for reviewed in ("true", "false"):
            params = {
                "query": f"(gene_exact:{gene_symbol}) AND (organism_id:{taxon_id}) AND (reviewed:{reviewed})",
                "format": "json",
                "size": 1,
                "fields": "accession",
            }
            data = self._http_get_json(UNIPROT_SEARCH_URL, params=params, timeout=30)
            results = data.get("results") or []
            if results:
                return results[0].get("primaryAccession")
        return None

    def fetch_uniprot_features(self, accession: str) -> Dict[str, Any]:
        url = f"{UNIPROT_ENTRY_BASE_URL}/{accession}.json"
        try:
            data = self._http_get_json(url, timeout=30)
        except requests.HTTPError:
            return {}

        tm: List[Dict[str, Any]] = []
        topo: List[Dict[str, Any]] = []
        domains: List[Dict[str, Any]] = []
        signal_peptides: List[Dict[str, Any]] = []
        chains: List[Dict[str, Any]] = []
        propeptides: List[Dict[str, Any]] = []
        disulfide_bonds: List[Dict[str, Any]] = []
        for feature in (data.get("features") or []):
            ftype = feature.get("type")
            loc = feature.get("location") or {}
            start = (loc.get("start") or {}).get("value")
            end = (loc.get("end") or {}).get("value")
            desc = feature.get("description")
            if ftype == "Transmembrane":
                tm.append({"start": start, "end": end, "description": desc})
            elif ftype == "Topological domain":
                topo.append({"start": start, "end": end, "description": desc})
            elif ftype in ("Domain", "Region"):
                domains.append({"start": start, "end": end, "description": desc, "type": ftype})
            elif ftype in ("Signal peptide", "Signal"):
                signal_peptides.append({"start": start, "end": end, "description": desc})
            elif ftype in ("Chain", "Peptide"):
                chains.append({"start": start, "end": end, "description": desc})
            elif ftype == "Propeptide":
                propeptides.append({"start": start, "end": end, "description": desc})
            elif ftype == "Disulfide bond":
                disulfide_bonds.append({"start": start, "end": end, "description": desc})

        protein_name = (
            (((data.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName") or {}).get("value")
        )
        return {
            "accession": accession,
            "protein_name": protein_name,
            "transmembrane": tm,
            "topological_domains": topo,
            "domains_regions": domains,
            "signal_peptides": signal_peptides,
            "chains": chains,
            "propeptides": propeptides,
            "disulfide_bonds": disulfide_bonds,
        }

    def fetch_spdis_for_hgvs(self, hgvs: str) -> List[str]:
        url = f"{NCBI_VARIATION_BASE_URL}/hgvs/{quote(hgvs, safe='')}/contextuals"
        try:
            data = self._http_get_json(url, timeout=30)
        except requests.HTTPError:
            return []

        spdis: List[str] = []
        if isinstance(data, dict):
            for key in ("spdis", "data"):
                if key not in data:
                    continue
                val = data[key]
                if isinstance(val, list):
                    for item in val:
                        if isinstance(item, str):
                            spdis.append(item)
                        elif isinstance(item, dict) and "spdi" in item:
                            spdis.append(item["spdi"])
                elif isinstance(val, dict) and isinstance(val.get("spdis"), list):
                    for item in val["spdis"]:
                        if isinstance(item, str):
                            spdis.append(item)
        elif isinstance(data, list):
            for item in data:
                if isinstance(item, str):
                    spdis.append(item)

        seen = set()
        out: List[str] = []
        for spdi in spdis:
            if spdi not in seen:
                seen.add(spdi)
                out.append(spdi)
        return out

    def fetch_rsids_for_spdi(self, spdi: str) -> List[str]:
        url = f"{NCBI_VARIATION_BASE_URL}/spdi/{quote(spdi, safe='')}/rsids"
        try:
            data = self._http_get_json(url, timeout=30)
        except requests.HTTPError:
            return []

        rsids: List[Any] = []
        if isinstance(data, dict):
            for key in ("rsids", "data"):
                if key not in data:
                    continue
                val = data[key]
                if isinstance(val, list):
                    rsids.extend(val)
                elif isinstance(val, dict) and isinstance(val.get("rsids"), list):
                    rsids.extend(val["rsids"])

        out: List[str] = []
        for item in rsids:
            if isinstance(item, int):
                out.append(f"rs{item}")
            elif isinstance(item, str) and item.startswith("rs"):
                out.append(item)
            elif isinstance(item, str) and item.isdigit():
                out.append(f"rs{item}")
        return sorted(set(out))

    def ensure_rsid(self, record: Dict[str, Any]) -> Optional[str]:
        rsid = record.get("dbSNP")
        if isinstance(rsid, str) and rsid.startswith("rs"):
            return rsid

        hgvs = pick_best_hgvs(record)
        if not hgvs:
            return None

        spdis = self.fetch_spdis_for_hgvs(hgvs)
        if not spdis:
            return None

        rsids = self.fetch_rsids_for_spdi(spdis[0])
        return rsids[0] if rsids else None

    def fetch_clinvar_esummary(self, variation_id: str) -> Dict[str, Any]:
        params: Dict[str, Any] = {
            "db": "clinvar",
            "id": str(variation_id),
            "retmode": "json",
            **self._ncbi_common_params(),
        }
        data = self._http_get_json(f"{NCBI_EUTILS_BASE_URL}/esummary.fcgi", params=params, timeout=30)
        return ((data.get("result") or {}).get(str(variation_id)) or {})

    def fetch_clinvar_elink_pubmed(self, variation_id: str) -> List[str]:
        params: Dict[str, Any] = {
            "dbfrom": "clinvar",
            "db": "pubmed",
            "id": str(variation_id),
            "retmode": "json",
            **self._ncbi_common_params(),
        }
        data = self._http_get_json(f"{NCBI_EUTILS_BASE_URL}/elink.fcgi", params=params, timeout=30)

        pmids: List[str] = []
        for linkset in (data.get("linksets") or []):
            if not isinstance(linkset, dict):
                continue
            for db in (linkset.get("linksetdbs") or []):
                if not isinstance(db, dict):
                    continue
                if db.get("dbto") == "pubmed":
                    pmids.extend([str(x) for x in (db.get("links") or [])])

        seen = set()
        out: List[str] = []
        for pmid in pmids:
            if pmid not in seen:
                seen.add(pmid)
                out.append(pmid)
        return out

    def extract_clinvar_germline(self, record: Dict[str, Any]) -> Dict[str, Any]:
        out: Dict[str, Any] = {
            "germline_description": None,
            "germline_review_status": None,
            "germline_last_evaluated": None,
            "traits": [],
        }

        gc = record.get("germline_classification") or {}
        if isinstance(gc, list):
            gc = next((x for x in gc if isinstance(x, dict)), {})
        if not isinstance(gc, dict):
            gc = {}

        legacy_cs = record.get("clinical_significance") or {}
        if isinstance(legacy_cs, list):
            legacy_cs = next((x for x in legacy_cs if isinstance(x, dict)), {})
        if not isinstance(legacy_cs, dict):
            legacy_cs = {}

        out["germline_description"] = clean_optional_text(gc.get("description") or legacy_cs.get("description"))
        out["germline_review_status"] = clean_optional_text(gc.get("review_status") or legacy_cs.get("review_status"))
        out["germline_last_evaluated"] = clean_optional_text(gc.get("last_evaluated"))

        raw_traits = gc.get("trait_set") or record.get("trait_set") or []
        if isinstance(raw_traits, dict):
            raw_traits = [raw_traits]
        if not isinstance(raw_traits, list):
            raw_traits = []

        traits: List[Dict[str, Any]] = []
        for trait in raw_traits:
            if not isinstance(trait, dict):
                continue
            name = clean_optional_text(trait.get("trait_name") or trait.get("name"))
            if not name:
                continue
            xrefs = trait.get("trait_xrefs") or []
            if isinstance(xrefs, dict):
                xrefs = [xrefs]
            if not isinstance(xrefs, list):
                xrefs = []
            traits.append({"name": name, "xrefs": xrefs})

        out["traits"] = traits
        return out

    def fetch_pubmed_esearch(self, query: str) -> List[str]:
        params: Dict[str, Any] = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": self.cfg.pubmed_max_pmids,
            **self._ncbi_common_params(),
        }
        try:
            data = self._http_get_json(f"{NCBI_EUTILS_BASE_URL}/esearch.fcgi", params=params, timeout=30)
        except requests.HTTPError:
            return []
        ids = (((data or {}).get("esearchresult") or {}).get("idlist") or [])
        return [str(x) for x in ids]

    def fetch_pubmed_esummary(self, pmids: List[str]) -> List[Dict[str, Any]]:
        if not pmids:
            return []

        params: Dict[str, Any] = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "json",
            **self._ncbi_common_params(),
        }
        try:
            data = self._http_get_json(f"{NCBI_EUTILS_BASE_URL}/esummary.fcgi", params=params, timeout=30)
        except requests.HTTPError:
            return []

        result = data.get("result") or {}
        out: List[Dict[str, Any]] = []
        for pmid in pmids:
            item = result.get(pmid)
            if not item:
                continue
            out.append(
                {
                    "pmid": pmid,
                    "title": item.get("title"),
                    "pubdate": item.get("pubdate"),
                    "source": item.get("source"),
                    "authors": [a.get("name") for a in (item.get("authors") or []) if isinstance(a, dict)],
                }
            )
        return out

    def harvest_gene(self, gene_symbol: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
        accession = self.fetch_uniprot_primary_accession(gene_symbol)
        gene_info: Dict[str, Any] = {
            "gene_symbol": gene_symbol,
            "uniprot_accession": accession,
        }
        if accession:
            uniprot = self.fetch_uniprot_features(accession)
            gene_info["uniprot"] = uniprot

        variants = self.fetch_clinicaltables_variants_for_gene(gene_symbol)

        enriched: List[Dict[str, Any]] = []
        clinvar_calls = 0
        clinvar_call_limit = self.cfg.max_variants_per_gene

        for variant in variants:
            hgvs = pick_best_hgvs(variant)
            rsid = self.ensure_rsid(variant)

            vep_raw = self.fetch_vep_for_hgvs(hgvs) if hgvs else None
            vep = summarize_vep_for_input(vep_raw, hgvs) if (vep_raw and hgvs) else {}

            variation_id = variant.get("clinvar_variation_id")
            clinvar_extra: Dict[str, Any] = {
                "germline_description": None,
                "germline_review_status": None,
                "germline_last_evaluated": None,
                "traits": [],
                "pubmed_ids": [],
            }

            fetch_ok = bool(variation_id) and should_fetch_clinvar(vep)
            if fetch_ok:
                cache_key = str(variation_id)
                if cache_key in self._clinvar_cache:
                    clinvar_extra = self._clinvar_cache[cache_key]
                elif clinvar_calls < clinvar_call_limit:
                    try:
                        esummary_record = self.fetch_clinvar_esummary(cache_key)
                        germline = self.extract_clinvar_germline(esummary_record)
                        pmids = self.fetch_clinvar_elink_pubmed(cache_key)
                        clinvar_extra = {
                            "germline_description": germline.get("germline_description"),
                            "germline_review_status": germline.get("germline_review_status"),
                            "germline_last_evaluated": germline.get("germline_last_evaluated"),
                            "traits": germline.get("traits") or [],
                            "pubmed_ids": pmids[: self.cfg.pubmed_max_pmids],
                        }
                        self._clinvar_cache[cache_key] = clinvar_extra
                        clinvar_calls += 1
                        time.sleep(self.cfg.sleep_seconds)
                    except Exception:
                        clinvar_extra = {
                            "germline_description": None,
                            "germline_review_status": None,
                            "germline_last_evaluated": None,
                            "traits": [],
                            "pubmed_ids": [],
                        }

            pmids = list(clinvar_extra.get("pubmed_ids") or [])
            if not pmids:
                if rsid:
                    pmids = self.fetch_pubmed_esearch(f"{rsid}[All Fields] AND {gene_symbol}[All Fields]")
                elif hgvs:
                    pmids = self.fetch_pubmed_esearch(f"{gene_symbol}[All Fields] AND \"{hgvs}\"[All Fields]")

            pubs = self.fetch_pubmed_esummary(pmids[: self.cfg.pubmed_max_pmids]) if pmids else []

            out = {
                "gene_symbol": gene_symbol,
                "clinvar_variation_id": str(variation_id) if variation_id else None,
                "hgvs_c": variant.get("HGVS_c"),
                "hgvs_p": variant.get("HGVS_p"),
                "amino_acid_change": variant.get("AminoAcidChange"),
                "dbsnp_rsid": rsid,
                "phenotype": variant.get("phenotype"),
                "phenotypes": variant.get("phenotypes"),
                "chromosome": variant.get("Chromosome"),
                "start": variant.get("Start"),
                "stop": variant.get("Stop"),
                "type": variant.get("Type"),
                "ref": variant.get("ReferenceAllele"),
                "alt": variant.get("AlternateAllele"),
                "vep": vep,
                "clinvar": clinvar_extra,
                "clinvar_clinical_significance": clinvar_extra.get("germline_description"),
                "clinvar_conditions": list(
                    dict.fromkeys(
                        [
                            clean_optional_text(t.get("name"))
                            for t in (clinvar_extra.get("traits") or [])
                            if isinstance(t, dict) and clean_optional_text(t.get("name"))
                        ]
                    )
                ),
                "clinvar_pubmed_ids": clinvar_extra.get("pubmed_ids"),
                "pubmed": pubs,
                "uniprot_accession": accession,
            }
            enriched.append(out)
            time.sleep(self.cfg.sleep_seconds)

        return gene_info, enriched

    def run(self, genes: List[str], out_jsonl: Path) -> None:
        with out_jsonl.open("w", encoding="utf-8") as out_f:
            for gene_symbol in genes:
                gene_info, variants = self.harvest_gene(gene_symbol)
                out_f.write(json.dumps({"record_type": "gene", **gene_info}, ensure_ascii=False) + "\n")
                for variant in variants:
                    out_f.write(json.dumps({"record_type": "variant", **variant}, ensure_ascii=False) + "\n")
                print(f"[OK] {gene_symbol}: {len(variants)} variants written")


def variant_has_clinical_association(record: Dict[str, Any], mode: str = "clinvar") -> bool:
    has_clinvar_sig = value_has_nonempty_text(record.get("clinvar_clinical_significance"))
    has_clinvar_conditions = value_has_nonempty_text(record.get("clinvar_conditions"))
    if mode == "clinvar":
        return has_clinvar_sig or has_clinvar_conditions

    has_clinicaltables_phenotype = value_has_nonempty_text(record.get("phenotype")) or value_has_nonempty_text(
        record.get("phenotypes")
    )
    if mode == "broad":
        return has_clinvar_sig or has_clinvar_conditions or has_clinicaltables_phenotype

    raise ValueError(f"Unsupported clinical association mode: {mode}")


def variant_has_clinical_label(record: Dict[str, Any], mode: str = "pathogenic_only") -> bool:
    sig = (clean_optional_text(record.get("clinvar_clinical_significance")) or "").lower()
    if not sig:
        return False
    if mode == "pathogenic_only":
        return sig in PATHOGENIC_SIG_VALUES
    if mode == "labeled":
        return sig in PATHOGENIC_SIG_VALUES or sig in BENIGN_SIG_VALUES
    raise ValueError(f"Unsupported clinical label mode: {mode}")


def write_clinical_association_subset(
    raw_jsonl: Path,
    out_jsonl: Path,
    mode: str = "clinvar",
) -> Dict[str, Any]:
    if raw_jsonl.resolve() == out_jsonl.resolve():
        raise ValueError("Clinical-association output path must be different from full raw JSONL path.")

    gene_count = 0
    variant_total = 0
    variant_kept = 0
    total_by_gene: Counter[str] = Counter()
    kept_by_gene: Counter[str] = Counter()

    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with raw_jsonl.open("r", encoding="utf-8") as in_f, out_jsonl.open("w", encoding="utf-8") as out_f:
        for line in in_f:
            text = line.strip()
            if not text:
                continue
            try:
                obj = json.loads(text)
            except json.JSONDecodeError:
                continue

            record_type = obj.get("record_type")
            if record_type == "gene":
                gene_count += 1
                out_f.write(json.dumps(obj, ensure_ascii=False) + "\n")
                continue

            if record_type != "variant":
                continue

            variant_total += 1
            gene = (clean_optional_text(obj.get("gene_symbol")) or "UNKNOWN").upper()
            total_by_gene[gene] += 1

            if variant_has_clinical_association(obj, mode=mode):
                variant_kept += 1
                kept_by_gene[gene] += 1
                out_f.write(json.dumps(obj, ensure_ascii=False) + "\n")

    return {
        "mode": mode,
        "gene_records_written": gene_count,
        "variant_total": variant_total,
        "variant_kept": variant_kept,
        "variant_keep_fraction": round((variant_kept / variant_total), 4) if variant_total else 0.0,
        "total_by_gene": dict(sorted(total_by_gene.items())),
        "kept_by_gene": dict(sorted(kept_by_gene.items())),
    }


def write_clinical_label_subset(
    raw_jsonl: Path,
    out_jsonl: Path,
    mode: str = "pathogenic_only",
) -> Dict[str, Any]:
    if raw_jsonl.resolve() == out_jsonl.resolve():
        raise ValueError("Clinical-label subset output path must be different from full raw JSONL path.")

    gene_count = 0
    variant_total = 0
    variant_kept = 0
    total_by_gene: Counter[str] = Counter()
    kept_by_gene: Counter[str] = Counter()

    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with raw_jsonl.open("r", encoding="utf-8") as in_f, out_jsonl.open("w", encoding="utf-8") as out_f:
        for line in in_f:
            text = line.strip()
            if not text:
                continue
            try:
                obj = json.loads(text)
            except json.JSONDecodeError:
                continue

            record_type = obj.get("record_type")
            if record_type == "gene":
                gene_count += 1
                out_f.write(json.dumps(obj, ensure_ascii=False) + "\n")
                continue

            if record_type != "variant":
                continue

            variant_total += 1
            gene = (clean_optional_text(obj.get("gene_symbol")) or "UNKNOWN").upper()
            total_by_gene[gene] += 1

            if variant_has_clinical_label(obj, mode=mode):
                variant_kept += 1
                kept_by_gene[gene] += 1
                out_f.write(json.dumps(obj, ensure_ascii=False) + "\n")

    return {
        "mode": mode,
        "gene_records_written": gene_count,
        "variant_total": variant_total,
        "variant_kept": variant_kept,
        "variant_keep_fraction": round((variant_kept / variant_total), 4) if variant_total else 0.0,
        "total_by_gene": dict(sorted(total_by_gene.items())),
        "kept_by_gene": dict(sorted(kept_by_gene.items())),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Import raw variant data from external databases")
    parser.add_argument("--genes", nargs="+", default=list(DEFAULT_GENES), help="Gene symbols to harvest")
    parser.add_argument("--out-jsonl", default=str(RAW_JSONL_PATH), help="Raw JSONL output path")
    parser.add_argument(
        "--out-clinical-jsonl",
        default=str(RAW_CLINICAL_ASSOC_JSONL_PATH),
        help="Clinical-association-only raw JSONL output path.",
    )
    parser.add_argument(
        "--clinical-association-mode",
        choices=["clinvar", "broad"],
        default="clinvar",
        help="Subset rule for clinical-association raw output.",
    )
    parser.add_argument(
        "--skip-clinical-subset",
        action="store_true",
        default=False,
        help="Skip writing the clinical-association-only raw subset file.",
    )
    parser.add_argument(
        "--out-pathogenic-jsonl",
        default="",
        help="Optional output path for clinically labeled/pathogenic subset raw JSONL. Disabled when empty.",
    )
    parser.add_argument(
        "--pathogenic-subset-mode",
        choices=["pathogenic_only", "labeled"],
        default="pathogenic_only",
        help="Subset rule for --out-pathogenic-jsonl.",
    )
    parser.add_argument("--assembly", choices=["grch37", "grch38"], default="grch37")
    parser.add_argument("--max-variants-per-gene", type=int, default=500)
    parser.add_argument("--pubmed-max", type=int, default=10)
    parser.add_argument("--sleep-seconds", type=float, default=0.12)
    parser.add_argument("--ncbi-email", default=os.getenv("NCBI_EMAIL"))
    parser.add_argument("--ncbi-api-key", default=os.getenv("NCBI_API_KEY"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_data_directories()

    genes = sorted({g.strip().upper() for g in args.genes if g and g.strip()})
    if not genes:
        raise SystemExit("No valid genes provided")

    cfg = PipelineConfig(
        genes=genes,
        assembly=args.assembly,
        max_variants_per_gene=args.max_variants_per_gene,
        pubmed_max_pmids=args.pubmed_max,
        sleep_seconds=args.sleep_seconds,
        ncbi_email=args.ncbi_email,
        ncbi_api_key=args.ncbi_api_key,
    )

    importer = DatabaseImporter(cfg)
    out_jsonl = Path(args.out_jsonl)
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    importer.run(cfg.genes, out_jsonl)
    print(f"Done. Wrote {out_jsonl}")

    if not args.skip_clinical_subset:
        out_clinical_jsonl = Path(args.out_clinical_jsonl)
        summary = write_clinical_association_subset(
            raw_jsonl=out_jsonl,
            out_jsonl=out_clinical_jsonl,
            mode=args.clinical_association_mode,
        )
        print(f"Done. Wrote {out_clinical_jsonl}")
        print(
            "[OK] Clinical-association subset "
            f"(mode={summary['mode']}): kept {summary['variant_kept']}/{summary['variant_total']} "
            f"({summary['variant_keep_fraction']:.1%}) variants"
        )
        print(f"[OK] Kept-by-gene: {summary['kept_by_gene']}")

    if args.out_pathogenic_jsonl:
        out_pathogenic_jsonl = Path(args.out_pathogenic_jsonl)
        path_summary = write_clinical_label_subset(
            raw_jsonl=out_jsonl,
            out_jsonl=out_pathogenic_jsonl,
            mode=args.pathogenic_subset_mode,
        )
        print(f"Done. Wrote {out_pathogenic_jsonl}")
        print(
            "[OK] Clinical-label subset "
            f"(mode={path_summary['mode']}): kept {path_summary['variant_kept']}/{path_summary['variant_total']} "
            f"({path_summary['variant_keep_fraction']:.1%}) variants"
        )
        print(f"[OK] Kept-by-gene: {path_summary['kept_by_gene']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
