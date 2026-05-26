#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Script: PYSCRIPT_Standardize_annotations_of_bacterial_genome.py
#
# Description:
# This script standardizes gene and CDS annotations in a GenBank file.
# It synchronizes /gene and /product qualifiers between paired gene and
# CDS features, giving priority to CDS annotations. It uses a local mapping
# table first, and when needed queries UniProtKB restricted to cyanobacteria
# to infer missing gene abbreviations from product descriptions or missing
# product descriptions from gene abbreviations. Annotations whose standard_name
# indicates hypothetical proteins are excluded from logging and external lookup,
# but their product qualifier is standardized to "hypothetical protein".
# Screen logging is written both to the terminal and to a plain-text log file.
# -----------------------------------------------------------------------------

from __future__ import annotations

import logging
import os
import re
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

import requests
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


# -----------------------------------------------------------------------------
# User mapping table
# -----------------------------------------------------------------------------
GENE_PRODUCT_MAP: Dict[str, Dict[str, str]] = {
    # Example entries:
    # "psbA": {"gene": "psbA", "product": "photosystem II D1 protein"},
    # "photosystem II D1 protein": {"gene": "psbA", "product": "photosystem II D1 protein"},
}


# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_TIMEOUT = int(os.environ.get("UNIPROT_REQUEST_TIMEOUT", "20"))
UNIPROT_MAX_RESULTS = int(os.environ.get("UNIPROT_MAX_RESULTS", "200"))

# UniProt taxonomy: Cyanobacteriota
CYANOBACTERIA_TAXON_ID = 1117

# Common bacterial/cyanobacterial-style gene symbol pattern
GENE_SYMBOL_RE = re.compile(r"^[A-Za-z]{3,4}[A-Za-z0-9]?[A-Za-z]?$")

# ANSI terminal colors
COLOR_WHITE = "\033[97m"
COLOR_YELLOW = "\033[93m"
COLOR_RED = "\033[91m"
COLOR_RESET = "\033[0m"


# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
class ColorConsoleFormatter(logging.Formatter):
    """
    Colorize only annotation status lines in the console.
    """

    def format(self, record: logging.LogRecord) -> str:
        message = super().format(record)
        color = getattr(record, "color", None)
        if color:
            return f"\n{color}{message}{COLOR_RESET}"
        return message


def setup_logging(log_file_path: str) -> logging.Logger:
    """
    Configure console logging with colors and file logging without colors.
    """
    logger = logging.getLogger("annotation_standardizer")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(ColorConsoleFormatter("%(message)s"))
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file_path, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(file_handler)

    return logger


# -----------------------------------------------------------------------------
# Text helpers
# -----------------------------------------------------------------------------
def clean_text(x: Optional[str]) -> Optional[str]:
    if x is None:
        return None
    return re.sub(r"\s+", " ", x.strip().strip('"'))


def normalize_text(x: Optional[str]) -> Optional[str]:
    if x is None:
        return None
    x = clean_text(x)
    if not x:
        return None
    x = x.lower().replace("-", " ")
    x = re.sub(r"[^\w\s/]", " ", x)
    x = re.sub(r"\s+", " ", x).strip()
    return x


def is_hypothetical(product: Optional[str]) -> bool:
    p = normalize_text(product)
    return bool(p and "hypothetical protein" in p)


def is_hypothetical_standard_name(standard_name: Optional[str]) -> bool:
    s = normalize_text(standard_name)
    if not s:
        return False
    return s in {"hypothetical protein cds", "hypothetical protein gene"}


def looks_like_gene(x: Optional[str]) -> bool:
    x = clean_text(x)
    return bool(x and " " not in x and GENE_SYMBOL_RE.match(x))


def first_q(f: Optional[SeqFeature], key: str) -> Optional[str]:
    if f is None:
        return None
    vals = f.qualifiers.get(key)
    return clean_text(vals[0]) if vals else None


def set_q(f: SeqFeature, key: str, val: str) -> None:
    f.qualifiers[key] = [val]


def choose(cds_val: Optional[str], gene_val: Optional[str]) -> Optional[str]:
    return clean_text(cds_val) or clean_text(gene_val)


def location_key(f: SeqFeature) -> Tuple[int, int, int]:
    return (
        int(f.location.start),
        int(f.location.end),
        int(f.location.strand or 0),
    )


def location_string(f: Optional[SeqFeature]) -> str:
    if f is None:
        return "NA"
    start = int(f.location.start)
    end = int(f.location.end)
    strand = int(f.location.strand or 0)
    strand_symbol = "+" if strand == 1 else "-" if strand == -1 else "?"
    return f"{start}-{end}({strand_symbol})"


def display_name(gene_f: Optional[SeqFeature], cds_f: Optional[SeqFeature]) -> str:
    """
    Display standard_name rather than locus_tag.
    Prefer CDS standard_name, then gene standard_name, then gene, then NA.
    """
    std_name = choose(first_q(cds_f, "standard_name"), first_q(gene_f, "standard_name"))
    if std_name:
        return std_name
    gene_name = choose(first_q(cds_f, "gene"), first_q(gene_f, "gene"))
    if gene_name:
        return gene_name
    return "NA"


def values_changed(
    gene_f: Optional[SeqFeature],
    cds_f: Optional[SeqFeature],
    final_gene: str,
    final_product: str,
) -> bool:
    """
    Determine whether writing final values would change any existing feature qualifier.
    """
    original_pairs = []

    if gene_f is not None:
        original_pairs.append((first_q(gene_f, "gene"), first_q(gene_f, "product")))

    if cds_f is not None:
        original_pairs.append((first_q(cds_f, "gene"), first_q(cds_f, "product")))

    for orig_gene, orig_prod in original_pairs:
        if clean_text(orig_gene) != clean_text(final_gene) or clean_text(orig_prod) != clean_text(final_product):
            return True

    return False


def is_successful_resolution(final_gene: str, final_product: str) -> bool:
    """
    A successful resolution must avoid unresolved fallback values.
    """
    if final_gene == "unknown_gene":
        return False
    if final_product == "hypothetical protein":
        return False
    return True


# -----------------------------------------------------------------------------
# Screen and file logging for annotation events
# -----------------------------------------------------------------------------
def log_processing_result(
    logger: logging.Logger,
    gene_f: Optional[SeqFeature],
    cds_f: Optional[SeqFeature],
    final_gene: str,
    final_product: str,
    resolution_steps: List[str],
    changed: bool,
    success: bool,
) -> None:
    """
    White  = no change
    Yellow = successful change
    Red    = unsuccessful resolution
    """
    disp_name = display_name(gene_f, cds_f)
    loc = location_string(cds_f or gene_f)
    resolution = ";".join(resolution_steps) if resolution_steps else "no_change"

    if not changed:
        color = COLOR_WHITE
        status = "UNCHANGED"
    elif success:
        color = COLOR_YELLOW
        status = "CHANGED"
    else:
        color = COLOR_RED
        status = "UNRESOLVED"

    message = (
        f"[{status}] standard_name={disp_name} | location={loc} | "
        f"final_gene={final_gene} | final_product={final_product} | "
        f"resolution={resolution}"
    )
    logger.info(message, extra={"color": color})


# -----------------------------------------------------------------------------
# UniProt helpers
# -----------------------------------------------------------------------------
def fetch_uniprot(query: str, fields: str) -> List[dict]:
    params = {
        "query": query,
        "format": "json",
        "size": UNIPROT_MAX_RESULTS,
        "fields": fields,
    }
    try:
        response = requests.get(UNIPROT_URL, params=params, timeout=UNIPROT_TIMEOUT)
        response.raise_for_status()
        return response.json().get("results", [])
    except Exception:
        return []


def extract_gene(entry: dict) -> Optional[str]:
    genes = entry.get("genes") or []
    for g in genes:
        gene_name = g.get("geneName")
        if isinstance(gene_name, dict):
            value = gene_name.get("value")
            if value:
                return clean_text(value)
    return None


def extract_product(entry: dict) -> Optional[str]:
    desc = entry.get("proteinDescription") or {}

    rec = desc.get("recommendedName") or {}
    full = rec.get("fullName")
    if isinstance(full, dict) and full.get("value"):
        return clean_text(full["value"])

    for section in ("submissionNames", "alternativeNames"):
        vals = desc.get(section) or []
        for v in vals:
            full = v.get("fullName")
            if isinstance(full, dict) and full.get("value"):
                return clean_text(full["value"])

    return None


def product_match_score(candidate: Optional[str], target: Optional[str]) -> Tuple[int, int]:
    c = normalize_text(candidate) or ""
    t = normalize_text(target) or ""

    if not c or not t:
        return (0, 0)
    if c == t:
        return (3, len(c))
    if c in t or t in c:
        return (2, min(len(c), len(t)))

    c_words = set(c.split())
    t_words = set(t.split())
    overlap = len(c_words & t_words)
    return (1 if overlap else 0, overlap)


def infer_gene_from_product_uniprot(product: str) -> Tuple[Optional[str], List[Tuple[str, int]]]:
    """
    Given a product description, query UniProt restricted to cyanobacteria
    and return the most common plausible gene abbreviation.
    """
    if not product or is_hypothetical(product):
        return None, []

    fields = "gene_primary,protein_name"
    query = f'(protein_name:"{product}") AND (taxonomy_id:{CYANOBACTERIA_TAXON_ID})'
    results = fetch_uniprot(query, fields)

    if not results:
        broad_query = f'("{product}") AND (taxonomy_id:{CYANOBACTERIA_TAXON_ID})'
        results = fetch_uniprot(broad_query, fields)

    weighted = Counter()
    exact = Counter()
    seen = Counter()

    for entry in results:
        g = extract_gene(entry)
        p = extract_product(entry)
        if not g or not looks_like_gene(g):
            continue

        score_class, _ = product_match_score(p, product)
        if score_class == 0:
            continue

        seen[g] += 1
        if score_class == 3:
            exact[g] += 1
            weighted[g] += 5
        elif score_class == 2:
            weighted[g] += 3
        else:
            weighted[g] += 1

    if not weighted:
        return None, []

    ranked = sorted(
        weighted.items(),
        key=lambda x: (-x[1], -exact.get(x[0], 0), -seen.get(x[0], 0), x[0]),
    )
    return ranked[0][0], ranked


def infer_product_from_gene_uniprot(gene: str) -> Tuple[Optional[str], List[Tuple[str, int]]]:
    """
    Given a gene symbol, query UniProt restricted to cyanobacteria
    and return the most common product description.
    """
    if not gene or not looks_like_gene(gene):
        return None, []

    fields = "gene_primary,protein_name"
    query = f'(gene:{gene}) AND (taxonomy_id:{CYANOBACTERIA_TAXON_ID})'
    results = fetch_uniprot(query, fields)

    counts = Counter()
    for entry in results:
        g = extract_gene(entry)
        p = extract_product(entry)
        if not g or not p:
            continue
        if clean_text(g) != clean_text(gene):
            continue
        if is_hypothetical(p):
            continue
        counts[p] += 1

    if not counts:
        return None, []

    ranked = counts.most_common()
    return ranked[0][0], ranked


# -----------------------------------------------------------------------------
# Mapping helpers
# -----------------------------------------------------------------------------
def apply_local_mapping(
    gene: Optional[str],
    product: Optional[str],
    standard_name: Optional[str],
) -> Tuple[Optional[str], Optional[str]]:
    gene = clean_text(gene)
    product = clean_text(product)
    standard_name = clean_text(standard_name)

    final_gene = gene
    final_product = product

    if final_gene and final_gene in GENE_PRODUCT_MAP:
        mapped = GENE_PRODUCT_MAP[final_gene]
        final_gene = mapped["gene"]
        final_product = final_product or mapped["product"]

    if final_product and final_product in GENE_PRODUCT_MAP:
        mapped = GENE_PRODUCT_MAP[final_product]
        final_gene = final_gene or mapped["gene"]
        final_product = mapped["product"]

    if standard_name:
        if standard_name in GENE_PRODUCT_MAP:
            mapped = GENE_PRODUCT_MAP[standard_name]
            final_gene = final_gene or mapped["gene"]
            final_product = final_product or mapped["product"]
        elif looks_like_gene(standard_name):
            final_gene = final_gene or standard_name
        elif not final_product:
            final_product = standard_name

    return final_gene, final_product


# -----------------------------------------------------------------------------
# Conflict reporting
# -----------------------------------------------------------------------------
def add_conflict(
    report: List[str],
    standard_name_value: str,
    issue: str,
    gene_gene: Optional[str],
    gene_prod: Optional[str],
    cds_gene: Optional[str],
    cds_prod: Optional[str],
    final_gene: Optional[str],
    final_prod: Optional[str],
    resolution: str,
) -> None:
    report.append(
        "\t".join(
            [
                standard_name_value or "NA",
                issue,
                gene_gene or "",
                gene_prod or "",
                cds_gene or "",
                cds_prod or "",
                final_gene or "",
                final_prod or "",
                resolution,
            ]
        )
    )


# -----------------------------------------------------------------------------
# Core logic
# -----------------------------------------------------------------------------
def process_pair(
    gene_f: Optional[SeqFeature],
    cds_f: Optional[SeqFeature],
    report: List[str],
    logger: logging.Logger,
) -> None:
    gene_gene = first_q(gene_f, "gene")
    gene_prod = first_q(gene_f, "product")
    gene_std = first_q(gene_f, "standard_name")

    cds_gene = first_q(cds_f, "gene")
    cds_prod = first_q(cds_f, "product")
    cds_std = first_q(cds_f, "standard_name")

    standard_name_value = display_name(gene_f, cds_f)

    merged_gene = choose(cds_gene, gene_gene)
    merged_prod = choose(cds_prod, gene_prod)
    merged_std = choose(cds_std, gene_std)

    # Special rule for annotations explicitly labeled as hypothetical protein
    # by standard_name: do not log, do not query UniProt, but standardize
    # product to "hypothetical protein".
    if is_hypothetical_standard_name(merged_std):
        final_gene = merged_gene or "unknown_gene"
        final_prod = "hypothetical protein"

        if gene_f is not None:
            set_q(gene_f, "gene", final_gene)
            set_q(gene_f, "product", final_prod)

        if cds_f is not None:
            set_q(cds_f, "gene", final_gene)
            set_q(cds_f, "product", final_prod)

        return

    resolution_steps: List[str] = []

    if gene_gene and cds_gene and gene_gene != cds_gene:
        add_conflict(
            report,
            standard_name_value,
            "gene_conflict",
            gene_gene,
            gene_prod,
            cds_gene,
            cds_prod,
            merged_gene,
            merged_prod,
            "CDS_priority_for_gene",
        )
        resolution_steps.append("gene_conflict_CDS_priority")

    if gene_prod and cds_prod and gene_prod != cds_prod:
        add_conflict(
            report,
            standard_name_value,
            "product_conflict",
            gene_gene,
            gene_prod,
            cds_gene,
            cds_prod,
            merged_gene,
            merged_prod,
            "CDS_priority_for_product",
        )
        resolution_steps.append("product_conflict_CDS_priority")

    final_gene, final_prod = apply_local_mapping(merged_gene, merged_prod, merged_std)

    if final_gene and not final_prod and final_gene in GENE_PRODUCT_MAP:
        final_prod = GENE_PRODUCT_MAP[final_gene]["product"]
        resolution_steps.append("local_map_gene_to_product")

    if final_prod and not final_gene and final_prod in GENE_PRODUCT_MAP:
        final_gene = GENE_PRODUCT_MAP[final_prod]["gene"]
        resolution_steps.append("local_map_product_to_gene")

    if not is_hypothetical(final_prod):
        if (not final_gene or not looks_like_gene(final_gene)) and final_prod:
            inferred_gene, ranked = infer_gene_from_product_uniprot(final_prod)
            if inferred_gene:
                final_gene = inferred_gene
                resolution_steps.append(
                    "uniprot_product_to_gene:" + ",".join(f"{k}:{v}" for k, v in ranked[:5])
                )

        if final_gene and looks_like_gene(final_gene) and not final_prod:
            inferred_prod, ranked = infer_product_from_gene_uniprot(final_gene)
            if inferred_prod:
                final_prod = inferred_prod
                resolution_steps.append(
                    "uniprot_gene_to_product:" + ",".join(f"{k}:{v}" for k, v in ranked[:5])
                )

    if not final_gene:
        final_gene = "unknown_gene"
        resolution_steps.append("fallback_unknown_gene")

    if not final_prod:
        final_prod = "hypothetical protein"
        resolution_steps.append("fallback_hypothetical_product")

    changed = values_changed(gene_f, cds_f, final_gene, final_prod)
    success = is_successful_resolution(final_gene, final_prod)

    if gene_f is not None:
        set_q(gene_f, "gene", final_gene)
        set_q(gene_f, "product", final_prod)

    if cds_f is not None:
        set_q(cds_f, "gene", final_gene)
        set_q(cds_f, "product", final_prod)

    if resolution_steps:
        add_conflict(
            report,
            standard_name_value,
            "resolution",
            gene_gene,
            gene_prod,
            cds_gene,
            cds_prod,
            final_gene,
            final_prod,
            ";".join(resolution_steps),
        )

    log_processing_result(
        logger=logger,
        gene_f=gene_f,
        cds_f=cds_f,
        final_gene=final_gene,
        final_product=final_prod,
        resolution_steps=resolution_steps,
        changed=changed,
        success=success,
    )


def standardize_record(record: SeqRecord, report: List[str], logger: logging.Logger) -> None:
    by_locus = defaultdict(list)
    by_coord = defaultdict(list)

    for f in record.features:
        if f.type not in ("gene", "CDS"):
            continue
        locus = first_q(f, "locus_tag")
        if locus:
            by_locus[locus].append(f)
        by_coord[location_key(f)].append(f)

    processed = set()

    for feats in by_locus.values():
        gene_f = next((f for f in feats if f.type == "gene"), None)
        cds_f = next((f for f in feats if f.type == "CDS"), None)
        process_pair(gene_f, cds_f, report, logger)
        if gene_f is not None:
            processed.add(id(gene_f))
        if cds_f is not None:
            processed.add(id(cds_f))

    for feats in by_coord.values():
        remaining = [f for f in feats if id(f) not in processed]
        if not remaining:
            continue
        gene_f = next((f for f in remaining if f.type == "gene"), None)
        cds_f = next((f for f in remaining if f.type == "CDS"), None)
        process_pair(gene_f, cds_f, report, logger)
        if gene_f is not None:
            processed.add(id(gene_f))
        if cds_f is not None:
            processed.add(id(cds_f))

    for f in record.features:
        if f.type not in ("gene", "CDS"):
            continue
        if id(f) in processed:
            continue
        if f.type == "gene":
            process_pair(f, None, report, logger)
        else:
            process_pair(None, f, report, logger)
        processed.add(id(f))


def main() -> None:
    if len(sys.argv) != 3:
        print(
            "Usage: python PYSCRIPT_Standardize_annotations_of_bacterial_genome.py input.gb output.gb",
            file=sys.stderr,
        )
        sys.exit(1)

    input_gb = sys.argv[1]
    output_gb = sys.argv[2]
    report_file = output_gb + ".conflicts.tsv"
    screen_log_file = output_gb + ".screen.log"

    logger = setup_logging(screen_log_file)

    logger.info("[START] Reading GenBank input...")
    records = list(SeqIO.parse(input_gb, "genbank"))
    if not records:
        print("No GenBank records found in input file.", file=sys.stderr)
        sys.exit(1)

    logger.info(f"[START] Loaded {len(records)} record(s).")
    report: List[str] = [
        "\t".join(
            [
                "standard_name",
                "issue",
                "gene_feature_gene",
                "gene_feature_product",
                "cds_feature_gene",
                "cds_feature_product",
                "final_gene",
                "final_product",
                "resolution",
            ]
        )
    ]

    for i, rec in enumerate(records, start=1):
        logger.info(f"[RECORD] Processing record {i}/{len(records)}: {rec.id}")
        standardize_record(rec, report, logger)

    logger.info("[WRITE] Writing standardized GenBank output...")
    SeqIO.write(records, output_gb, "genbank")

    logger.info("[WRITE] Writing conflict report...")
    with open(report_file, "w", encoding="utf-8") as out:
        for line in report:
            out.write(line + "\n")

    logger.info(f"[DONE] Written standardized GenBank file to: {output_gb}")
    logger.info(f"[DONE] Written conflict report to: {report_file}")
    logger.info(f"[DONE] Written screen log to: {screen_log_file}")


if __name__ == "__main__":
    main()
