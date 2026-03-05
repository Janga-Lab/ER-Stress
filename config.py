from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

# ----------------------------
# Project paths
# ----------------------------

ROOT_DIR = Path(__file__).resolve().parent
DATA_DIR = ROOT_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw_data"
CURATED_DATA_DIR = DATA_DIR / "curated_data"
PROCESSED_DATA_DIR = DATA_DIR / "processed_data"

RAW_JSONL_PATH = RAW_DATA_DIR / "harvest_raw.jsonl"
RAW_CLINICAL_ASSOC_JSONL_PATH = RAW_DATA_DIR / "harvest_raw_clinical_association.jsonl"
PREPARED_FOR_PARSER_JSON_PATH = CURATED_DATA_DIR / "prepared_for_parser.json"
CURATED_VARIANTS_CSV_PATH = CURATED_DATA_DIR / "curated_variants.csv"
CURATED_GENE_FEATURES_JSON_PATH = CURATED_DATA_DIR / "curated_gene_features.json"
CURATION_VALIDATION_REPORT_PATH = CURATED_DATA_DIR / "curation_validation_report.json"

PROCESSED_VARIANT_CSV_PATH = PROCESSED_DATA_DIR / "variant_processed.csv"
PROCESSED_VARIANT_SUMMARY_CSV_PATH = PROCESSED_DATA_DIR / "variant_summary.csv"
PROCESSED_VARIANT_PREDICTION_CSV_PATH = PROCESSED_DATA_DIR / "variant_prediction_summary.csv"
PROCESSED_DECISION_TIERS_CSV_PATH = PROCESSED_DATA_DIR / "decision_tiers.csv"
PROCESSED_PRIORITY_RANKED_CSV_PATH = PROCESSED_DATA_DIR / "priority_ranked.csv"
PROCESSED_RUN_MANIFEST_PATH = PROCESSED_DATA_DIR / "run_manifest.json"

# ----------------------------
# External APIs (single source of truth)
# ----------------------------

CLINICALTABLES_VARIANTS_SEARCH_URL = "https://clinicaltables.nlm.nih.gov/api/variants/v4/search"
ENSEMBL_VEP_GRCH37_BASE_URL = "https://grch37.rest.ensembl.org/vep/human/hgvs"
ENSEMBL_VEP_GRCH38_BASE_URL = "https://rest.ensembl.org/vep/human/hgvs"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_BASE_URL = "https://rest.uniprot.org/uniprotkb"
NCBI_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_VARIATION_BASE_URL = "https://api.ncbi.nlm.nih.gov/variation/v0"
PUBMED_WEB_BASE_URL = "https://pubmed.ncbi.nlm.nih.gov"
SPLICEAI_LOOKUP_BASE_URL = "https://spliceailookup.broadinstitute.org"

# ----------------------------
# Default scientific settings
# ----------------------------

DEFAULT_GENES: List[str] = ["WFS1", "CISD2", "EIF2AK3", "IER3IP1", "INS"]
DEFAULT_ASSEMBLY = "grch37"

DEFAULT_MAX_VARIANTS_PER_GENE = 500
DEFAULT_PUBMED_MAX_PMIDS = 10
DEFAULT_SLEEP_SECONDS = 0.12
DEFAULT_NCBI_TOOL = "er_stress_variant_harvester"

DEFAULT_CLASS_WEIGHTS: Dict[str, float] = {
    "lof_truncating": 4.0,
    "splice": 3.0, # Splicing defects 
    "cnv_sv": 3.0, 
    "inframe_indel": 2.5,
    "missense": 2.0,
    "synonymous": 0.0,
    "other": 0.0,
}

DEFAULT_CRITICAL_BONUS = 1.0

DEFAULT_INHERITANCE_MAP: Dict[str, str] = {
    "WFS1": "recessive",
    "CISD2": "recessive",
    "EIF2AK3": "recessive",
    "IER3IP1": "recessive",
    "INS": "dominant",
}

DEFAULT_SCORE_BIN_THRESHOLDS: Dict[str, float] = {
    "very_high": 8.0,
    "high": 6.0,
    "moderate": 4.0,
    "mild": 2.0,
}

DEFAULT_VARIANT_PREDICTION_THRESHOLD = 4.0
DEFAULT_CONFIDENCE_LOW_THRESHOLD = 0.50
DEFAULT_CONFIDENCE_HIGH_THRESHOLD = 0.75
DEFAULT_UNMAPPED_WARNING_THRESHOLD = 0.25


@dataclass
class PipelineConfig:
    genes: List[str] = field(default_factory=lambda: list(DEFAULT_GENES))
    assembly: str = DEFAULT_ASSEMBLY
    max_variants_per_gene: int = DEFAULT_MAX_VARIANTS_PER_GENE
    pubmed_max_pmids: int = DEFAULT_PUBMED_MAX_PMIDS
    sleep_seconds: float = DEFAULT_SLEEP_SECONDS
    ncbi_tool: str = DEFAULT_NCBI_TOOL
    ncbi_email: str | None = None
    ncbi_api_key: str | None = None

    class_weights: Dict[str, float] = field(default_factory=lambda: dict(DEFAULT_CLASS_WEIGHTS))
    critical_bonus: float = DEFAULT_CRITICAL_BONUS
    inheritance_map: Dict[str, str] = field(default_factory=lambda: dict(DEFAULT_INHERITANCE_MAP))
    score_bin_thresholds: Dict[str, float] = field(default_factory=lambda: dict(DEFAULT_SCORE_BIN_THRESHOLDS))

    variant_prediction_threshold: float = DEFAULT_VARIANT_PREDICTION_THRESHOLD
    confidence_low_threshold: float = DEFAULT_CONFIDENCE_LOW_THRESHOLD
    confidence_high_threshold: float = DEFAULT_CONFIDENCE_HIGH_THRESHOLD
    unmapped_warning_threshold: float = DEFAULT_UNMAPPED_WARNING_THRESHOLD


def ensure_data_directories() -> None:
    RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)
    CURATED_DATA_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)
