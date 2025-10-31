#!/usr/bin/env python3
"""
Cohort Builder: discovery ? QC ? metadata retrieval ? merge

Core functionality for assembling a scientifically defensible cohort with
explicit inclusion/exclusion reasoning and robust logging.

Exit codes:
 0: success (n_after_metadata_qc > 0)
 2: biologically empty (ran fine but n_after_metadata_qc == 0)
 1: internal failure (system/HTTP failures across sources)
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import requests

# Optional pandas import for convenience (not required for core path)
try:
    import pandas as pd  # noqa: F401
except Exception:  # pragma: no cover
    pd = None


# ------------------------- Constants and Config -------------------------
ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FILEREPORT_SAMPLE_RESULT = "sample"
EUTILS_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EUTILS_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

USER_AGENT = "Cohort-Builder/1.0"

# ENA filereport field whitelist (explicit, valid fields for read_run)
ENA_FIELDS = [
    # Accessions
    "run_accession", "study_accession", "secondary_study_accession",
    "sample_accession", "secondary_sample_accession", "experiment_accession",
    "submission_accession",
    # Taxon / organism
    "tax_id", "scientific_name",
    # Platform / instrument
    "instrument_platform", "instrument_model",
    # Library
    "library_name", "library_strategy", "library_source", "library_selection", "library_layout",
    # Counts / sizes
    "read_count", "base_count",
    # Dates / center
    "first_public", "last_updated", "center_name",
    # Titles / descriptions
    "study_title", "experiment_title", "sample_title", "sample_description",
    # Sample-level (present on read_run where available)
    "collection_date",
    # FASTQ file pointers
    "fastq_ftp", "fastq_md5", "fastq_bytes",
    # Submitted file pointers
    "submitted_ftp", "submitted_md5", "submitted_bytes"
]

# Discovery exclude keywords (case-insensitive)
DISCOVERY_EXCLUDE_KEYWORDS = [
    "microbiome", "metagenome", "stool", "fecal", "oral", "16s", "amplicon", "shotgun metagenomics"
]

# Allowed strategies for RNA-Seq-like data (case-insensitive containment)
LIBRARY_STRATEGY_ALLOWLIST = [
    # common bulk RNA terms
    "rna-seq", "rna seq", "mrna-seq", "mrna seq", "mrna sequencing",
    "total rna", "totalrna", "total rna sequencing", "cdna", "whole transcriptome",
    "transcriptome", "transcriptomic", "stranded rna", "poly a", "polya", "ribo-zero", "ribodeplete", "ribodepleted",
    # single-cell / nucleus terms
    "scrna", "sc rna", "single cell", "single-cell", "snrna", "sn rna",
    # common protocols
    "smart-seq", "smartseq", "10x", "3' rna", "5' rna", "3prime", "5prime"
]

HUMAN_TAX_IDS = {"9606"}
HUMAN_NAMES = {"homo sapiens", "human"}

# ---------------------------- Utilities ----------------------------

def ensure_dirs(out_dir: Path) -> Dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = out_dir / "logs"
    logs_dir.mkdir(exist_ok=True)
    return {
        "out": out_dir,
        "logs": logs_dir,
        "metadata_tsv": out_dir / "comprehensive_metadata.tsv",
        "discard_tsv": out_dir / "discarded_accessions.tsv",
        "manifest_json": out_dir / "cohort_manifest.json",
        "ena_log": logs_dir / "ena_requests.log",
        "efetch_log": logs_dir / "efetch_runinfo.log",
        "builder_log": logs_dir / "cohort_builder.log",
    }


def log_append(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(line.rstrip("\n") + "\n")


def retry_request(method: str, url: str, session: requests.Session, *,
                  params: Optional[dict] = None,
                  headers: Optional[dict] = None,
                  data: Optional[dict] = None,
                  retries: int = 3,
                  backoff_base: float = 1.0) -> requests.Response:
    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            resp = session.request(method, url, params=params, headers=headers, data=data, timeout=30)
            # Basic validation
            if resp.status_code == 200:
                return resp
            last_exc = RuntimeError(f"HTTP {resp.status_code}: {resp.text[:200]}")
        except Exception as e:  # pragma: no cover
            last_exc = e
        time.sleep(backoff_base * (2 ** (attempt - 1)))
    raise last_exc if last_exc else RuntimeError("Unknown HTTP error")


# ----------------------- YAML Term Resolution -----------------------

def load_yaml_minimal(yaml_path: Path) -> Dict[str, dict]:
    """
    Minimal YAML loader for a very simple schema (string keys + lists/strings).
    If PyYAML is available, use it; otherwise parse a constrained subset.
    """
    try:
        import yaml  # type: ignore
        with open(yaml_path, "r", encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}
    except Exception:
        # Fallback: constrained parser for:
        # key:\n  field: "..."\n  list_field:\n    - "..."
        text = yaml_path.read_text(encoding="utf-8")
        data: Dict[str, dict] = {}
        current_key: Optional[str] = None
        current_list_field: Optional[str] = None
        for line in text.splitlines():
            if re.match(r"^\s*#", line) or not line.strip():
                continue
            # Top-level key
            m_key = re.match(r"^(\w[\w_]*):\s*$", line)
            if m_key:
                current_key = m_key.group(1)
                data[current_key] = {}
                current_list_field = None
                continue
            if current_key is None:
                continue
            # List field header (must be matched BEFORE scalar field)
            m_list = re.match(r"^\s{2}(\w[\w_]*):\s*$", line)
            if m_list:
                field = m_list.group(1)
                data[current_key][field] = []
                current_list_field = field
                continue
            # Scalar field "  name: value"
            m_scalar = re.match(r"^\s{2}(\w[\w_]*):\s*\"?(.*?)\"?$", line)
            if m_scalar:
                field, val = m_scalar.group(1), m_scalar.group(2)
                # Only treat as scalar if a value is present on the same line
                if val != "":
                    data[current_key][field] = val
                    current_list_field = None
                continue
            # List item
            m_item = re.match(r"^\s{4}-\s*\"?(.*?)\"?$", line)
            if m_item and current_list_field:
                data[current_key][current_list_field].append(m_item.group(1))
        return data


def normalize_cancer_key(label: str) -> str:
    k = label.strip().lower()
    k = re.sub(r"[^a-z0-9]+", "_", k)
    k = re.sub(r"_+", "_", k).strip("_")
    # Specific mapping for common variants
    mapping = {
        "anus": "anus_anal_canal_anorectum",
        "anus_anal_canal_anorectum": "anus_anal_canal_anorectum",
        "anus_anal_canal_&_anorectum": "anus_anal_canal_anorectum",
        "anus_anal_canal_and_anorectum": "anus_anal_canal_anorectum",
    }
    return mapping.get(k, k)


def resolve_cancer_terms(cancer_label: str, yaml_path: Path) -> dict:
    db = load_yaml_minimal(yaml_path) if yaml_path.exists() else {}
    key = normalize_cancer_key(cancer_label)
    entry = db.get(key)
    if not entry:
        # Fallback minimal definition to allow immediate use
        entry = {
            "display_name": cancer_label,
            "anatomy_terms": ["anal", "anal canal", "anorectal", "rectal"],
            "malignancy_terms": [
                "cancer", "carcinoma", "adenocarcinoma", "squamous cell carcinoma", "malignant"
            ],
            "search_terms": [
                "anal squamous cell carcinoma",
                "anorectal carcinoma",
                "rectal adenocarcinoma",
                "anal canal carcinoma",
                "rectal cancer",
            ],
            "notes": "exclude stool/oral microbiome; focus on human epithelial malignancy",
        }
    # Coerce types and ensure lists are lists; if missing or empty, backfill sensible defaults
    def as_list(v, default: List[str]) -> List[str]:
        if isinstance(v, list):
            return v
        return default

    entry_anatomy = as_list(entry.get("anatomy_terms"), ["anal", "anal canal", "anorectal", "rectal"])
    # Ensure common synonyms are present
    for syn in ["anus", "rectum", "anorectum"]:
        if syn not in entry_anatomy:
            entry_anatomy.append(syn)
    entry_malign = as_list(entry.get("malignancy_terms"), [
        "cancer", "carcinoma", "adenocarcinoma", "squamous cell carcinoma", "malignant"
    ])
    entry_search = as_list(entry.get("search_terms"), [
        "anal squamous cell carcinoma", "anorectal carcinoma", "rectal adenocarcinoma", "anal canal carcinoma", "rectal cancer"
    ])

    return {
        "normalized_key": key,
        "display_name": entry.get("display_name", cancer_label),
        "anatomy_terms": entry_anatomy,
        "malignancy_terms": entry_malign,
        "search_terms": entry_search,
        "notes": entry.get("notes", ""),
    }


# ------------------------ Discovery and Filters ------------------------

def parse_sample_attribute(attr_str: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not attr_str:
        return out
    # ENA style: key1=value1; key2=value2; ...
    parts = [p.strip() for p in attr_str.split(";") if p.strip()]
    for part in parts:
        if "=" in part:
            k, v = part.split("=", 1)
        elif ":" in part:
            k, v = part.split(":", 1)
        else:
            # Unkeyed attribute, skip
            continue
        key = re.sub(r"[^a-z0-9]+", "_", k.strip().lower()).strip("_")
        out[key] = v.strip()
    return out


def looks_human(scientific_name: Optional[str], tax_id: Optional[str]) -> bool:
    sname = (scientific_name or "").strip().lower()
    tid = (str(tax_id or "").strip())
    if tid in HUMAN_TAX_IDS:
        return True
    if any(h in sname for h in HUMAN_NAMES):
        return True
    return False


def contains_any(text: str, terms: List[str]) -> bool:
    t = text.lower()
    return any(term.lower() in t for term in terms)


def candidate_text_blurb(rec: dict) -> str:
    fields = [
        rec.get("study_title", ""),
        rec.get("experiment_title", ""),
        rec.get("sample_title", ""),
        rec.get("sample_description", ""),
    ]
    # include parsed attributes textual values
    attrs = parse_sample_attribute(rec.get("sample_attribute", ""))
    fields.extend(attrs.values())
    # include biosample-derived attributes if present
    for k, v in rec.items():
        if k.startswith("biosample_") and isinstance(v, str):
            fields.append(v)
    # include bioproject-derived attributes if present
    for k, v in rec.items():
        if k.startswith("bioproject_") and isinstance(v, str):
            fields.append(v)
    return " \n ".join([f for f in fields if f])


def apply_discovery_filters(candidates: List[dict], anatomy_terms: List[str], malignancy_terms: List[str]) -> Tuple[List[dict], List[dict]]:
    kept: List[dict] = []
    dropped: List[dict] = []
    for rec in candidates:
        run = rec.get("run_accession") or rec.get("Run")
        sname = rec.get("scientific_name") or rec.get("ScientificName")
        tid = str(rec.get("tax_id") or rec.get("TaxID") or "").strip()
        text = candidate_text_blurb(rec)
        if not looks_human(sname, tid):
            dropped.append({"accession": run, "stage_dropped": "discovery", "reason": "non_human"})
            continue
        if contains_any(text, DISCOVERY_EXCLUDE_KEYWORDS):
            dropped.append({"accession": run, "stage_dropped": "discovery", "reason": "exclude_keyword"})
            continue
        if not (contains_any(text, anatomy_terms) and contains_any(text, malignancy_terms)):
            # Allow pass-through if the matched query term itself carries both signals
            qterm = (rec.get("query_term") or "").lower()
            if not (contains_any(qterm, anatomy_terms) and contains_any(qterm, malignancy_terms)):
                dropped.append({"accession": run, "stage_dropped": "discovery", "reason": "missing_anatomy_or_malignancy_terms"})
                continue
        kept.append(rec)
    return kept, dropped


def search_ena(term: str, limit: int, session: requests.Session, logs: Dict[str, Path]) -> List[dict]:
    params = {
        "dataPortal": "ena",
        "result": "read_run",
        "format": "json",
        "fields": ",".join([
            "run_accession","study_accession","sample_accession","experiment_accession",
            "scientific_name","tax_id","library_strategy","library_source","library_selection","library_layout",
            "study_title","experiment_title","sample_title","sample_description"
        ]),
        # Use ENA query syntax; exact phrase search for the term
        "query": f'tax_eq(9606) AND (study_title="{term}" OR sample_title="{term}" OR experiment_title="{term}")',
        "limit": str(limit),
    }
    try:
        resp = retry_request("GET", ENA_SEARCH_URL, session, params=params, headers={"User-Agent": USER_AGENT})
        data = resp.json() if resp.headers.get("Content-Type", "").startswith("application/json") else []
        if isinstance(data, dict):
            entries = data.get("entries") or data.get("data") or []
        elif isinstance(data, list):
            entries = data
        else:
            entries = []
        for e in entries:
            e.setdefault("source", "ENA")
            e.setdefault("query_term", term)
        return entries
    except Exception as e:  # pragma: no cover
        log_append(logs["ena_log"], f"search_ena error for term '{term}': {e}")
        return []


def search_sra(term: str, limit: int, session: requests.Session, logs: Dict[str, Path]) -> List[dict]:
    # E-utilities esearch to get SRA UIDs
    try:
        es_params = {
            "db": "sra",
            "term": f"{term} AND Homo sapiens[Organism]",
            "retmax": str(limit),
            "retmode": "json",
        }
        es_resp = retry_request("GET", EUTILS_ESEARCH_URL, session, params=es_params, headers={"User-Agent": USER_AGENT})
        es_json = es_resp.json()
        idlist = es_json.get("esearchresult", {}).get("idlist", [])
        if not idlist:
            return []
        # efetch runinfo for the list of UIDs
        ef_params = {
            "db": "sra",
            "id": ",".join(idlist[:200]),  # cap batch size to avoid URL limit
            "rettype": "runinfo",
            "retmode": "text",
        }
        ef_resp = retry_request("GET", EUTILS_EFETCH_URL, session, params=ef_params, headers={"User-Agent": USER_AGENT})
        text = ef_resp.text
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if len(lines) <= 1:
            return []
        reader = csv.DictReader(lines)
        out: List[dict] = []
        for row in reader:
            if row.get("Run"):
                out.append({
                    "run_accession": row.get("Run"),
                    "scientific_name": row.get("ScientificName"),
                    "tax_id": row.get("TaxID"),
                    "study_title": "",  # not present here
                    "sample_title": row.get("SampleName", ""),
                    "sample_description": "",
                    "library_strategy": row.get("LibraryStrategy", ""),
                    "library_source": row.get("LibrarySource", ""),
                    "library_selection": row.get("LibrarySelection", ""),
                    "library_layout": row.get("LibraryLayout", ""),
                    "source": "SRA",
                    "query_term": term
                })
        return out
    except Exception as e:  # pragma: no cover
        log_append(logs["ena_log"], f"search_sra error for term '{term}': {e}")
        return []


# ----------------------- Metadata Retrieval/QC -----------------------

def fetch_ena_metadata(accession: str, session: requests.Session, logs: Dict[str, Path], retries: int = 3) -> Tuple[Optional[dict], Optional[str]]:
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "download": "true",
    }
    try:
        resp = retry_request("GET", ENA_FILEREPORT_URL, session, params=params, headers={"User-Agent": USER_AGENT}, retries=retries)
        text = resp.text
        if not text or text.lstrip().lower().startswith("please add column config") or "<html" in text.lower():
            reason = "ena_filereport_invalid_response"
            log_append(logs["ena_log"], f"{accession}\t{reason}")
            return None, reason
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if len(lines) < 2:
            reason = "ena_filereport_empty"
            log_append(logs["ena_log"], f"{accession}\t{reason}")
            return None, reason
        reader = csv.DictReader(lines, delimiter='\t')
        rows = list(reader)
        if not rows:
            reason = "ena_filereport_no_rows"
            log_append(logs["ena_log"], f"{accession}\t{reason}")
            return None, reason
        row = rows[0]
        row.setdefault("run_accession", accession)
        return row, None
    except Exception as e:  # pragma: no cover
        reason = f"ena_filereport_exception:{e}"[:200]
        log_append(logs["ena_log"], f"{accession}\t{reason}")
        return None, reason


def fetch_ncbi_runinfo(accession: str, session: requests.Session, logs: Dict[str, Path], retries: int = 3) -> Tuple[Optional[dict], Optional[str]]:
    params = {
        "db": "sra",
        "id": accession,
        "rettype": "runinfo",
        "retmode": "text",
    }
    try:
        resp = retry_request("GET", EUTILS_EFETCH_URL, session, params=params, headers={"User-Agent": USER_AGENT}, retries=retries)
        text = resp.text
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if len(lines) <= 1:
            reason = "runinfo_empty"
            log_append(logs["efetch_log"], f"{accession}\t{reason}")
            return None, reason
        reader = csv.DictReader(lines)
        for row in reader:
            if row.get("Run") == accession:
                out = {
                    "run_accession": accession,
                    "ScientificName": row.get("ScientificName"),
                    "TaxID": row.get("TaxID"),
                    "LibraryStrategy": row.get("LibraryStrategy"),
                    "LibrarySource": row.get("LibrarySource"),
                    "LibrarySelection": row.get("LibrarySelection"),
                    "LibraryLayout": row.get("LibraryLayout"),
                    "BioSample": row.get("BioSample"),
                    "BioProject": row.get("BioProject"),
                }
                return out, None
        reason = "runinfo_no_matching_row"
        log_append(logs["efetch_log"], f"{accession}\t{reason}")
        return None, reason
    except Exception as e:  # pragma: no cover
        reason = f"runinfo_exception:{e}"[:200]
        log_append(logs["efetch_log"], f"{accession}\t{reason}")
        return None, reason


def merge_metadata_row(ena_row: Optional[dict], runinfo_row: Optional[dict]) -> dict:
    merged: dict = {}
    if ena_row:
        merged.update(ena_row)
    if runinfo_row:
        for k, v in runinfo_row.items():
            if k not in merged or not merged[k]:
                merged[k] = v
    return merged


def fetch_biosample(biosample_id: str, session: requests.Session, logs: Dict[str, Path], retries: int = 3) -> Tuple[Optional[dict], Optional[str]]:
    if not biosample_id:
        return None, "no_biosample_id"
    params = {
        "db": "biosample",
        "id": biosample_id,
        "retmode": "json",
    }
    try:
        resp = retry_request("GET", EUTILS_EFETCH_URL, session, params=params, headers={"User-Agent": USER_AGENT}, retries=retries)
        data = resp.json()
        # Parse common structures
        bios = []
        if isinstance(data, dict):
            if "BioSampleSet" in data:
                bs = data["BioSampleSet"].get("BioSample")
                if isinstance(bs, list):
                    bios = bs
                elif isinstance(bs, dict):
                    bios = [bs]
        if not bios:
            return None, "biosample_empty"
        bs = bios[0]
        out: dict = {"biosample_accession": bs.get("accession", biosample_id)}
        # Title/Description/Organism
        for fld in ("Title", "Description"):
            if fld in bs:
                out[f"biosample_{fld.lower()}"] = bs[fld]
        org = bs.get("Organism") or {}
        if isinstance(org, dict):
            if "#text" in org:
                out["biosample_organism"] = org["#text"]
            if "@taxonomy_id" in org:
                out["biosample_taxonomy_id"] = org["@taxonomy_id"]
        # Attributes
        attrs = bs.get("Attributes", {}).get("Attribute")
        if isinstance(attrs, dict):
            attrs = [attrs]
        if isinstance(attrs, list):
            for a in attrs:
                name = a.get("@attribute_name") or a.get("harmonized_name") or a.get("attribute_name")
                val = a.get("#text") or a.get("value") or a.get("text")
                if name and val:
                    key = re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")
                    out[f"biosample_{key}"] = str(val)
        return out, None
    except Exception as e:
        reason = f"biosample_exception:{e}"[:200]
        log_append(logs["efetch_log"], f"{biosample_id}\t{reason}")
        return None, reason


def fetch_bioproject(bioproject_id: str, session: requests.Session, logs: Dict[str, Path], retries: int = 3) -> Tuple[Optional[dict], Optional[str]]:
    if not bioproject_id:
        return None, "no_bioproject_id"
    params = {
        "db": "bioproject",
        "id": bioproject_id,
        "retmode": "json",
    }
    try:
        resp = retry_request("GET", EUTILS_EFETCH_URL, session, params=params, headers={"User-Agent": USER_AGENT}, retries=retries)
        data = resp.json()
        project = None
        if isinstance(data, dict) and "Project" in data:
            project = data["Project"]
        if not project:
            return None, "bioproject_empty"
        out: dict = {
            "bioproject_accession": bioproject_id,
        }
        pid = project.get("ProjectID", {})
        arch = pid.get("ArchiveID", {})
        if isinstance(arch, dict):
            out["bioproject_archive_id"] = arch.get("#text", "")
        descr = project.get("ProjectDescr", {})
        if isinstance(descr, dict):
            out["bioproject_title"] = descr.get("Title", "")
            out["bioproject_description"] = descr.get("Description", "")
        # Study descriptor
        study = project.get("Study", {})
        if isinstance(study, dict):
            sdesc = study.get("Descriptor", {})
            if isinstance(sdesc, dict):
                out["bioproject_study_title"] = sdesc.get("StudyTitle", "")
                out["bioproject_study_description"] = sdesc.get("StudyAbstract", "")
                out["bioproject_study_type"] = sdesc.get("StudyType", "")
        # Target organism
        ptype = project.get("ProjectType", {})
        target = ptype.get("Target", {}) if isinstance(ptype, dict) else {}
        org = target.get("Organism", {}) if isinstance(target, dict) else {}
        if isinstance(org, dict):
            out["bioproject_target_organism"] = org.get("#text", "")
            out["bioproject_target_taxonomy_id"] = org.get("@taxonomy_id", "")
        return out, None
    except Exception as e:
        reason = f"bioproject_exception:{e}"[:200]
        log_append(logs["efetch_log"], f"{bioproject_id}\t{reason}")
        return None, reason


def fetch_ena_sample(sample_accession: str, session: requests.Session, logs: Dict[str, Path], retries: int = 3) -> Tuple[Optional[dict], Optional[str]]:
    """Fetch ENA 'sample' result for a sample_accession and parse into biosample_* fields.
    Tries to request sample_attributes; on 400 fallback to minimal fields.
    """
    if not sample_accession:
        return None, "no_sample_accession"
    def do_request(fields: List[str]) -> Tuple[Optional[dict], Optional[str]]:
        params = {
            "accession": sample_accession,
            "result": ENA_FILEREPORT_SAMPLE_RESULT,
            "format": "tsv",
            "download": "true",
            "fields": ",".join(fields),
        }
        try:
            resp = retry_request("GET", ENA_FILEREPORT_URL, session, params=params, headers={"User-Agent": USER_AGENT}, retries=retries)
            text = resp.text
            if not text or "<html" in text.lower() or text.lower().startswith("please add column config"):
                return None, "ena_sample_invalid_response"
            lines = [ln for ln in text.splitlines() if ln.strip()]
            if len(lines) < 2:
                return None, "ena_sample_empty"
            reader = csv.DictReader(lines, delimiter='\t')
            rows = list(reader)
            if not rows:
                return None, "ena_sample_no_rows"
            row = rows[0]
            out: dict = {
                "biosample_accession": row.get("sample_accession", sample_accession)
            }
            # map core sample fields
            mapping = {
                "scientific_name": "biosample_scientific_name",
                "tax_id": "biosample_tax_id",
                "sample_title": "biosample_title",
                "sample_description": "biosample_description",
                "first_public": "biosample_first_public",
                "last_updated": "biosample_last_updated",
                "center_name": "biosample_center_name",
                "collection_date": "biosample_collection_date",
            }
            for src, dst in mapping.items():
                if row.get(src):
                    out[dst] = row.get(src)
            # parse sample_attributes if present
            attr_val = row.get("sample_attributes") or row.get("sample_attribute")
            if isinstance(attr_val, str) and attr_val.strip():
                attrs = parse_sample_attribute(attr_val)
                for k, v in attrs.items():
                    out[f"biosample_{k}"] = v
            return out, None
        except Exception as e:
            return None, f"ena_sample_exception:{e}"[:200]
    # Try with attributes first
    fields_with_attr = [
        "sample_accession","scientific_name","tax_id","sample_title","sample_description",
        "first_public","last_updated","center_name","collection_date","sample_attributes"
    ]
    row, err = do_request(fields_with_attr)
    if row is None and err and "invalid" in err:
        # Fallback minimal
        fields_min = ["sample_accession","scientific_name","tax_id","sample_title","first_public","last_updated","center_name"]
        row, err = do_request(fields_min)
    if row is None and err:
        log_append(logs["ena_log"], f"{sample_accession}\t{err}")
    return row, err


def strategy_allowed(strategy: str, library_source: str) -> bool:
    s = (strategy or "").lower()
    src = (library_source or "").lower()
    # hard excludes
    if any(kw in s for kw in ["amplicon", "16s", "metagenomic", "metatranscriptomic"]):
        return False
    # accept if explicitly transcriptomic in source
    if "transcriptomic" in src:
        return True
    # accept if strategy contains any known rna terms
    if any(allow in s for allow in LIBRARY_STRATEGY_ALLOWLIST):
        return True
    # conservative default
    return False


def apply_metadata_qc(row: dict, anatomy_terms: List[str], malignancy_terms: List[str]) -> Tuple[bool, Optional[str]]:
    # Species/metagenome checks
    if not looks_human(row.get("scientific_name") or row.get("ScientificName"), str(row.get("tax_id") or row.get("TaxID") or "")):
        return False, "non_human"
    text = candidate_text_blurb(row)
    if contains_any(text, DISCOVERY_EXCLUDE_KEYWORDS):
        return False, "exclude_keyword"
    # Anatomy + malignancy requirement
    if not (contains_any(text, anatomy_terms) and contains_any(text, malignancy_terms)):
        return False, "missing_anatomy_or_malignancy_terms"
    # Library strategy
    if not strategy_allowed(row.get("library_strategy", ""), row.get("library_source", "")):
        return False, "library_strategy_not_allowed"
    # Parsed sample attributes can indicate metagenome
    attrs = parse_sample_attribute(row.get("sample_attribute", ""))
    if any(contains_any(str(v), ["metagenome", "microbiome"]) for v in attrs.values()):
        return False, "attributes_indicate_microbiome"
    return True, None


# -------------------------- Writing outputs --------------------------

def write_tsv(rows: List[dict], path: Path) -> None:
    if not rows:
        # Write empty file with no header
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("")
        return
    # Determine columns by union of keys, stable order with run_accession first
    cols = ["run_accession"]
    seen = set(cols)
    for row in rows:
        for k in row.keys():
            if k not in seen:
                cols.append(k)
                seen.add(k)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def write_manifest(manifest: dict, path: Path) -> None:
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


# ------------------------------ Orchestrator ------------------------------

def build_cohort(cancer_label: str, out_dir: Path, yaml_path: Path, limit: int = 2000) -> int:
    paths = ensure_dirs(out_dir)
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    # Resolve terms
    terms = resolve_cancer_terms(cancer_label, yaml_path)
    anatomy_terms = terms["anatomy_terms"]
    malignancy_terms = terms["malignancy_terms"]
    search_terms = terms["search_terms"]

    # Discovery: ENA + SRA
    raw_candidates: List[dict] = []
    for term in search_terms:
        ena_res = search_ena(term, limit, session, paths)
        sra_res = search_sra(term, limit, session, paths)
        raw_candidates.extend(ena_res)
        raw_candidates.extend(sra_res)

    # Deduplicate by run_accession (canonical key)
    dedup: Dict[str, dict] = {}
    for rec in raw_candidates:
        acc = rec.get("run_accession") or rec.get("Run")
        if not acc:
            continue
        if acc not in dedup:
            dedup[acc] = rec
    candidates = list(dedup.values())

    n_candidates_raw = len(candidates)

    # Discovery filters
    kept_candidates, dropped_discovery = apply_discovery_filters(candidates, anatomy_terms, malignancy_terms)

    # Metadata retrieval + merge + QC
    merged_rows: List[dict] = []
    dropped_qc: List[dict] = []
    internal_failures = 0

    biosample_cache: Dict[str, dict] = {}
    for rec in kept_candidates:
        acc = rec.get("run_accession") or rec.get("Run")
        ena_row, ena_err = fetch_ena_metadata(acc, session, paths)
        runinfo_row, runinfo_err = fetch_ncbi_runinfo(acc, session, paths)
        if ena_err and runinfo_err:
            internal_failures += 1
            dropped_qc.append({"accession": acc, "stage_dropped": "metadata_qc", "reason": "no_metadata_sources"})
            continue
        merged = merge_metadata_row(ena_row, runinfo_row)
        # Fetch and merge BioSample attributes if available (prefer ENA sample endpoint; fallback NCBI)
        biosample_id = merged.get("BioSample") or merged.get("biosample_accession") or merged.get("sample_accession")
        if isinstance(biosample_id, str) and biosample_id.strip():
            bsid = biosample_id.strip()
            bs_row = biosample_cache.get(bsid)
            if bs_row is None:
                # Try ENA sample first for attributes
                bs_row, bs_err = fetch_ena_sample(bsid, session, paths)
                if bs_row is None:
                    # Fallback to NCBI BioSample
                    bs_row, bs_err = fetch_biosample(bsid, session, paths)
                if bs_row is not None:
                    biosample_cache[bsid] = bs_row
            if bs_row:
                for k, v in bs_row.items():
                    if k not in merged or not merged[k]:
                        merged[k] = v
        # Fetch and merge BioProject metadata if available
        bioproject_id = merged.get("BioProject") or merged.get("bioproject_accession")
        if isinstance(bioproject_id, str) and bioproject_id.strip():
            bpid = bioproject_id.strip()
            # simple in-memory cache via biosample_cache as well
            cache_key = f"__BP__{bpid}"
            bp_row = biosample_cache.get(cache_key)
            if bp_row is None:
                bp_row, bp_err = fetch_bioproject(bpid, session, paths)
                if bp_row is not None:
                    biosample_cache[cache_key] = bp_row
            if bp_row:
                for k, v in bp_row.items():
                    if k not in merged or not merged[k]:
                        merged[k] = v
        keep, reason = apply_metadata_qc(merged, anatomy_terms, malignancy_terms)
        if keep:
            merged_rows.append(merged)
        else:
            dropped_qc.append({"accession": acc, "stage_dropped": "metadata_qc", "reason": reason or "failed_qc"})

    # Outputs
    write_tsv(merged_rows, paths["metadata_tsv"])
    # Normalize discarded entries to always include run_accession
    normalized_discards: List[dict] = []
    for d in (dropped_discovery + dropped_qc):
        nd = dict(d)
        if nd.get("run_accession") in (None, "") and nd.get("accession"):
            nd["run_accession"] = nd.get("accession")
        normalized_discards.append(nd)
    write_tsv(normalized_discards, paths["discard_tsv"])  # simple 2+ column TSV

    # Aggregate reasons
    reason_counts = collections.Counter([d.get("stage_dropped", "") + ": " + d.get("reason", "") for d in (dropped_discovery + dropped_qc)])

    # Diagnostics: counts by source and query term (on deduped candidates)
    src_counts = collections.Counter([str(rec.get("source") or "").upper() or "UNKNOWN" for rec in candidates])
    term_counts = collections.Counter([rec.get("query_term") or "" for rec in candidates])

    # Field coverage across merged rows for key fields
    key_fields = [
        "run_accession", "study_accession", "sample_accession", "experiment_accession",
        "library_strategy", "library_source", "instrument_model", "read_count", "base_count",
        "scientific_name", "tax_id", "biosample_accession", "biosample_title", "biosample_description",
        "bioproject_accession", "center_name", "fastq_ftp"
    ]
    coverage: Dict[str, dict] = {}
    total_kept = max(1, len(merged_rows))
    for f in key_fields:
        present = sum(1 for r in merged_rows if str(r.get(f, "")).strip() != "")
        coverage[f] = {"count": present, "fraction": round(present / total_kept, 3)}

    manifest = {
        "normalized_cancer_key": terms["normalized_key"],
        "display_name": terms["display_name"],
        "search_terms": search_terms,
        "anatomy_terms": anatomy_terms,
        "malignancy_terms": malignancy_terms,
        "n_candidates_raw": n_candidates_raw,
        "n_after_discovery_filter": len(kept_candidates),
        "n_after_metadata_qc": len(merged_rows),
        "top_discard_reasons": dict(reason_counts.most_common()),
        "source_counts": dict(src_counts),
        "term_counts": dict(term_counts),
        "field_coverage": coverage,
        "qc_criteria": {
            "species_filter": "Homo sapiens / TaxID 9606",
            "anatomy_terms_required": True,
            "malignancy_terms_required": True,
            "discovery_exclude_keywords": DISCOVERY_EXCLUDE_KEYWORDS,
            "library_strategy_allowlist": LIBRARY_STRATEGY_ALLOWLIST,
        },
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    write_manifest(manifest, paths["manifest_json"])

    # Exit codes
    if len(merged_rows) > 0:
        return 0
    # If all metadata sources failed heavily, treat as internal failure
    if internal_failures >= max(1, len(kept_candidates)):
        return 1
    return 2


def main():
    ap = argparse.ArgumentParser(description="Cohort builder: discovery?QC?metadata?merge")
    ap.add_argument("-c", "--cancer-label", required=True, help="Cancer label (e.g., 'Anus')")
    ap.add_argument("-o", "--out-dir", required=True, help="Output metadata directory (cancer-specific)")
    ap.add_argument("-y", "--yaml", default="data/cancer_terms.yml", help="YAML mapping for cancer terms")
    ap.add_argument("-m", "--max-results", type=int, default=1000, help="Max results per source per term")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    yaml_path = Path(args.yaml)

    try:
        code = build_cohort(args.cancer_label, out_dir, yaml_path, limit=args.max_results)
        sys.exit(code)
    except Exception as e:  # pragma: no cover
        # Best-effort fallback logging
        try:
            paths = ensure_dirs(out_dir)
            log_append(paths["builder_log"], f"fatal_exception\t{e}")
        except Exception:
            pass
        print(f"ERROR: Cohort builder failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
