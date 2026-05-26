#!/usr/bin/env python3
"""
Evaluates if annotated CDS reading frames are intact in a bacterial genome
(available in GenBank flatfile format).

Typical usage:
    python PYSCRIPT_Evaluate_reading_frames_of_genes.py -i genome.gb -o genome_summary

Outputs:
    <prefix>.replicons.tsv           annotation-summary metrics per record/replicon
    <prefix>.genome.tsv              whole-genome annotation-summary metrics
    <prefix>.cds_frame_integrity.tsv  CDS-level reading-frame validation results

Notes
-----
- Works with one or many GenBank records in the same file (chromosome + plasmids).
- Uses annotated features already present in the GenBank file.
- The CDS table validates annotated CDS reading-frame integrity; it is not
    de novo gene prediction.
- "Complete vs partial rRNA" is inferred heuristically from GenBank feature
  locations/qualifiers:
    * partial if location is fuzzy (< or >) or qualifier "partial" is present
    * otherwise counted as complete
- "Intergenic spacers" are computed from merged annotated intervals formed from
  gene/CDS/RNA/pseudogene features, so the result is annotation-centric.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import statistics
import warnings
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


# ----------------------------- helpers --------------------------------- #

def safe_mean(values: Sequence[float]) -> float:
    return float(statistics.mean(values)) if values else float("nan")


def safe_median(values: Sequence[float]) -> float:
    return float(statistics.median(values)) if values else float("nan")


def gc_percent(seq: str) -> float:
    seq = seq.upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    atgc = a + c + g + t
    if atgc == 0:
        return float("nan")
    return 100.0 * (g + c) / atgc


def format_value(v):
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return "nan"
        return f"{v:.1f}"
    return v


def get_feature_span(feature: SeqFeature) -> Tuple[int, int]:
    """
    Return zero-based half-open [start, end) span across the full feature.
    For joined/compound features, this is the min start to max end span.
    """
    try:
        loc = feature.location
        if isinstance(loc, CompoundLocation):
            starts = [int(part.start) for part in loc.parts]
            ends = [int(part.end) for part in loc.parts]
            return min(starts), max(ends)
        return int(loc.start), int(loc.end)
    except Exception:
        return 0, 0


def get_feature_length(feature: SeqFeature) -> int:
    try:
        return int(len(feature.location))
    except Exception:
        start, end = get_feature_span(feature)
        return max(0, end - start)


def get_feature_boundaries(feature: SeqFeature) -> Tuple[int, int]:
    """
    Return one-based inclusive start/end coordinates for reporting.
    """
    start0, end0 = get_feature_span(feature)
    return start0 + 1, end0


def has_partial_location(feature: SeqFeature) -> bool:
    """
    Heuristic detection of partial/fuzzy locations.
    In Biopython, '<' and '>' usually become BeforePosition / AfterPosition.
    """
    try:
        loc = feature.location
    except Exception:
        return False

    def _pos_is_partial(pos) -> bool:
        cls = pos.__class__.__name__
        return cls in {"BeforePosition", "AfterPosition", "WithinPosition", "OneOfPosition"}

    if isinstance(loc, CompoundLocation):
        for part in loc.parts:
            if _pos_is_partial(part.start) or _pos_is_partial(part.end):
                return True
        return False

    return _pos_is_partial(loc.start) or _pos_is_partial(loc.end)


def qualifier_contains(feature: SeqFeature, key: str, pattern: str) -> bool:
    vals = feature.qualifiers.get(key, [])
    regex = re.compile(pattern, flags=re.IGNORECASE)
    return any(regex.search(v) for v in vals)


def get_first_qualifier(feature: SeqFeature, key: str) -> str:
    vals = feature.qualifiers.get(key, [])
    return vals[0] if vals else ""


def get_joined_qualifier(feature: SeqFeature, key: str) -> str:
    vals = feature.qualifiers.get(key, [])
    return " | ".join(vals) if vals else ""


def is_pseudogene_feature(feature: SeqFeature) -> bool:
    if feature.type == "pseudogene":
        return True
    if "pseudo" in feature.qualifiers:
        return True
    if "pseudogene" in feature.qualifiers:
        return True
    return False


def infer_cds_partiality(feature: SeqFeature) -> bool:
    """
    Detect expected partial CDSs from fuzzy locations and explicit qualifiers.
    """
    if has_partial_location(feature):
        return True
    if "partial" in feature.qualifiers:
        return True
    if qualifier_contains(feature, "note", r"\bpartial\b|\bfragment\b|\bincomplete\b|\btruncated\b"):
        return True
    if qualifier_contains(feature, "product", r"\bpartial\b|\bfragment\b|\bincomplete\b|\btruncated\b"):
        return True
    if qualifier_contains(feature, "exception", r"\bpartial\b|\bfragment\b|\bincomplete\b|\btruncated\b"):
        return True
    return False


def parse_translation_table(feature: SeqFeature) -> Tuple[int, bool]:
    """
    Return the selected translation table and whether the qualifier had to be defaulted.
    """
    vals = feature.qualifiers.get("transl_table", [])
    if not vals:
        return 11, False
    try:
        table_id = int(str(vals[0]).strip())
    except Exception:
        return 11, True
    if table_id not in CodonTable.unambiguous_dna_by_id:
        return 11, True
    return table_id, False


def normalize_translation_text(text: str) -> str:
    return re.sub(r"\s+", "", text).upper()


def manual_translate_cds(cds_seq: str, table_id: int) -> Tuple[str, int, int, str, str, int]:
    """
    Translate in-frame codons directly so reporting can continue even if
    Seq.translate raises on malformed features.
    """
    table = CodonTable.unambiguous_dna_by_id[table_id]
    full_length = len(cds_seq) - (len(cds_seq) % 3)
    codons = [cds_seq[i : i + 3] for i in range(0, full_length, 3)]
    amino_acids: List[str] = []
    ambiguous_codon_count = 0
    internal_stop_count = 0

    for index, codon in enumerate(codons):
        if re.fullmatch(r"[ACGT]{3}", codon) is None:
            ambiguous_codon_count += 1
            amino_acids.append("X")
            continue
        if codon in table.stop_codons:
            amino_acids.append("*")
            if index < len(codons) - 1:
                internal_stop_count += 1
            continue
        amino_acids.append(table.forward_table.get(codon, "X"))

    start_codon = codons[0] if codons else ""
    stop_codon = codons[-1] if codons else ""
    return (
        "".join(amino_acids),
        ambiguous_codon_count,
        internal_stop_count,
        start_codon,
        stop_codon,
        len(cds_seq) % 3,
    )


def infer_rrna_subtype(feature: SeqFeature) -> str:
    """
    Return one of: 5S, 16S, 23S, other, unknown
    """
    texts = []
    for key in ("product", "gene", "note", "standard_name"):
        texts.extend(feature.qualifiers.get(key, []))
    blob = " ".join(texts).lower()

    if "16s" in blob or "small subunit ribosomal rna" in blob or "ssu ribosomal rna" in blob:
        return "16S"
    if "23s" in blob or "large subunit ribosomal rna" in blob or "lsu ribosomal rna" in blob:
        return "23S"
    if "5s" in blob:
        return "5S"
    if "rrna" in blob or "ribosomal rna" in blob:
        return "other"
    return "unknown"


def is_nc_rna(feature: SeqFeature) -> bool:
    if feature.type in {"ncRNA", "misc_RNA"}:
        return True
    if feature.type == "RNA":
        return True
    return False


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    ints = sorted((s, e) for s, e in intervals if e > s)
    if not ints:
        return []
    merged = [ints[0]]
    for s, e in ints[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e:
            merged[-1] = (last_s, max(last_e, e))
        else:
            merged.append((s, e))
    return merged


def interval_gaps_and_overlaps(intervals: Iterable[Tuple[int, int]]) -> Tuple[List[int], List[int]]:
    """
    For sorted merged or unmerged intervals, compute positive gaps and positive overlaps
    between adjacent intervals after sorting by start,end.
    """
    ints = sorted((s, e) for s, e in intervals if e > s)
    gaps: List[int] = []
    overlaps: List[int] = []
    if len(ints) < 2:
        return gaps, overlaps

    prev_s, prev_e = ints[0]
    for s, e in ints[1:]:
        diff = s - prev_e
        if diff > 0:
            gaps.append(diff)
        elif diff < 0:
            overlaps.append(-diff)
        if e > prev_e:
            prev_s, prev_e = s, e
        else:
            prev_s, prev_e = prev_s, prev_e
    return gaps, overlaps


@dataclass
class RepliconSummary:
    replicon_id: str
    replicon_name: str
    record_type: str
    topology: str
    length_bp: int
    gc_percent: float

    genes_total: int
    cds_total: int
    cds_with_translation: int
    pseudogenes_total: int

    trna_total: int
    rrna_total: int
    rrna_5s: int
    rrna_16s: int
    rrna_23s: int
    rrna_other: int
    rrna_complete: int
    rrna_partial: int

    ncrna_total: int
    tmrna_total: int
    other_rna_total: int

    avg_gene_length_bp: float
    median_gene_length_bp: float
    avg_cds_length_bp: float
    median_cds_length_bp: float

    annotated_bases_bp: int
    coding_bases_bp: int
    coding_density_percent: float

    intergenic_spacers_n: int
    intergenic_spacers_total_bp: int
    intergenic_spacers_mean_bp: float
    intergenic_spacers_median_bp: float
    intergenic_spacers_min_bp: float
    intergenic_spacers_max_bp: float

    annotated_overlaps_n: int
    annotated_overlaps_total_bp: int
    annotated_overlaps_mean_bp: float

    cds_frame_total: int
    cds_frame_intact_total: int
    cds_frame_partial_total: int
    cds_frame_pseudogene_total: int
    cds_frame_problematic_total: int
    cds_frame_error_total: int
    cds_frame_internal_stop_total: int
    cds_frame_length_error_total: int
    cds_frame_missing_start_total: int
    cds_frame_missing_stop_total: int
    cds_frame_ambiguous_codon_total: int
    cds_frame_translation_mismatch_total: int
    cds_frame_non_partial_non_pseudogene_intact_percent: float


@dataclass
class CDSFrameIntegrity:
    replicon_id: str
    locus_tag: str
    gene: str
    product: str
    feature_start: int
    feature_end: int
    strand: int
    cds_length_bp: int
    transl_table: int
    is_partial: bool
    is_pseudogene: bool
    length_mod_3: int
    start_codon: str
    stop_codon: str
    has_valid_start: bool
    has_valid_stop: bool
    internal_stop_count: int
    ambiguous_codon_count: int
    translation_qualifier_present: bool
    translation_matches_qualifier: bool | None
    integrity_status: str
    issues: str


def evaluate_cds_feature(record: SeqRecord, feature: SeqFeature) -> CDSFrameIntegrity:
    issues: List[str] = []
    is_partial = infer_cds_partiality(feature)
    is_pseudogene = is_pseudogene_feature(feature)
    transl_table, table_defaulted = parse_translation_table(feature)
    if table_defaulted and "transl_table" in feature.qualifiers:
        issues.append("translation_table_defaulted_to_11")

    feature_start, feature_end = get_feature_boundaries(feature)
    strand = feature.location.strand if feature.location is not None else None
    locus_tag = get_first_qualifier(feature, "locus_tag")
    gene = get_first_qualifier(feature, "gene")
    product = get_joined_qualifier(feature, "product")
    translation_qualifier_present = "translation" in feature.qualifiers and bool(feature.qualifiers.get("translation"))

    try:
        cds_seq = str(feature.extract(record.seq)).upper().replace("U", "T")
    except Exception as exc:
        issues.append(f"extraction_error:{exc}")
        return CDSFrameIntegrity(
            replicon_id=record.id,
            locus_tag=locus_tag,
            gene=gene,
            product=product,
            feature_start=feature_start,
            feature_end=feature_end,
            strand=strand,
            cds_length_bp=get_feature_length(feature),
            transl_table=transl_table,
            is_partial=is_partial,
            is_pseudogene=is_pseudogene,
            length_mod_3=-1,
            start_codon="",
            stop_codon="",
            has_valid_start=False,
            has_valid_stop=False,
            internal_stop_count=-1,
            ambiguous_codon_count=-1,
            translation_qualifier_present=translation_qualifier_present,
            translation_matches_qualifier=None,
            integrity_status="error",
            issues=";".join(issues),
        )

    cds_length_bp = len(cds_seq)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            script_translation = str(Seq(cds_seq).translate(table=transl_table, to_stop=False))
        translation_failed = False
    except Exception as exc:
        issues.append(f"translation_error:{exc}")
        script_translation = ""
        translation_failed = True

    manual_translation, ambiguous_codon_count, internal_stop_count, start_codon, stop_codon, length_mod_3 = manual_translate_cds(cds_seq, transl_table)
    table = CodonTable.unambiguous_dna_by_id[transl_table]
    has_valid_start = bool(start_codon) and start_codon in table.start_codons
    has_valid_stop = bool(stop_codon) and stop_codon in table.stop_codons

    translation_matches_qualifier: bool | None = None
    if translation_qualifier_present and not translation_failed:
        qualifier_translation = normalize_translation_text(get_joined_qualifier(feature, "translation"))
        candidate_translation = normalize_translation_text(script_translation or manual_translation)
        translation_matches_qualifier = candidate_translation.rstrip("*") == qualifier_translation.rstrip("*")
        if not translation_matches_qualifier:
            issues.append("translation_mismatch")

    if not is_partial and not is_pseudogene:
        if length_mod_3 != 0:
            issues.append("length_not_divisible_by_three")
        if not has_valid_start:
            issues.append("missing_start_codon")
        if not has_valid_stop:
            issues.append("missing_stop_codon")
        if internal_stop_count > 0:
            issues.append("internal_stop_codon")
        if ambiguous_codon_count > 0:
            issues.append("ambiguous_codon")

    if translation_failed:
        integrity_status = "error"
    elif is_pseudogene:
        integrity_status = "pseudogene"
        issues.insert(0, "pseudogene_feature")
    elif is_partial:
        integrity_status = "partial"
        issues.insert(0, "partial_feature")
    elif any(tag in issues for tag in {
        "length_not_divisible_by_three",
        "missing_start_codon",
        "missing_stop_codon",
        "internal_stop_codon",
        "ambiguous_codon",
        "translation_mismatch",
    }):
        integrity_status = "problematic"
    else:
        integrity_status = "intact"

    if integrity_status == "error" and not issues:
        issues.append("translation_or_extraction_error")

    return CDSFrameIntegrity(
        replicon_id=record.id,
        locus_tag=locus_tag,
        gene=gene,
        product=product,
        feature_start=feature_start,
        feature_end=feature_end,
        strand=strand,
        cds_length_bp=cds_length_bp,
        transl_table=transl_table,
        is_partial=is_partial,
        is_pseudogene=is_pseudogene,
        length_mod_3=length_mod_3,
        start_codon=start_codon,
        stop_codon=stop_codon,
        has_valid_start=has_valid_start,
        has_valid_stop=has_valid_stop,
        internal_stop_count=internal_stop_count,
        ambiguous_codon_count=ambiguous_codon_count,
        translation_qualifier_present=translation_qualifier_present,
        translation_matches_qualifier=translation_matches_qualifier,
        integrity_status=integrity_status,
        issues=";".join(dict.fromkeys(issues)),
    )


# -------------------------- core summarization -------------------------- #

def summarize_record(record: SeqRecord) -> Tuple[RepliconSummary, List[CDSFrameIntegrity]]:
    seq = str(record.seq)
    length_bp = len(seq)
    gc = gc_percent(seq)

    genes_total = 0
    cds_total = 0
    cds_with_translation = 0
    pseudogenes_total = 0

    trna_total = 0
    rrna_total = 0
    rrna_5s = 0
    rrna_16s = 0
    rrna_23s = 0
    rrna_other = 0
    rrna_complete = 0
    rrna_partial = 0

    ncrna_total = 0
    tmrna_total = 0
    other_rna_total = 0

    gene_lengths: List[int] = []
    cds_lengths: List[int] = []

    annotation_intervals: List[Tuple[int, int]] = []
    coding_intervals: List[Tuple[int, int]] = []
    cds_rows: List[CDSFrameIntegrity] = []

    cds_frame_total = 0
    cds_frame_intact_total = 0
    cds_frame_partial_total = 0
    cds_frame_pseudogene_total = 0
    cds_frame_problematic_total = 0
    cds_frame_error_total = 0
    cds_frame_internal_stop_total = 0
    cds_frame_length_error_total = 0
    cds_frame_missing_start_total = 0
    cds_frame_missing_stop_total = 0
    cds_frame_ambiguous_codon_total = 0
    cds_frame_translation_mismatch_total = 0

    for feat in record.features:
        ftype = feat.type

        if is_pseudogene_feature(feat):
            pseudogenes_total += 1

        if ftype == "gene":
            genes_total += 1
            gene_lengths.append(get_feature_length(feat))
            annotation_intervals.append(get_feature_span(feat))

        elif ftype == "CDS":
            cds_total += 1
            cds_len = get_feature_length(feat)
            cds_lengths.append(cds_len)
            annotation_intervals.append(get_feature_span(feat))
            coding_intervals.append(get_feature_span(feat))
            if "translation" in feat.qualifiers:
                cds_with_translation += 1

            cds_frame_total += 1
            cds_row = evaluate_cds_feature(record, feat)
            cds_rows.append(cds_row)

            if cds_row.integrity_status == "intact":
                cds_frame_intact_total += 1
            elif cds_row.integrity_status == "partial":
                cds_frame_partial_total += 1
            elif cds_row.integrity_status == "pseudogene":
                cds_frame_pseudogene_total += 1
            elif cds_row.integrity_status == "problematic":
                cds_frame_problematic_total += 1
            elif cds_row.integrity_status == "error":
                cds_frame_error_total += 1

            if cds_row.internal_stop_count > 0:
                cds_frame_internal_stop_total += 1
            if cds_row.length_mod_3 != -1 and cds_row.length_mod_3 != 0:
                cds_frame_length_error_total += 1
            if cds_row.integrity_status == "problematic" and not cds_row.has_valid_start:
                cds_frame_missing_start_total += 1
            if cds_row.integrity_status == "problematic" and not cds_row.has_valid_stop:
                cds_frame_missing_stop_total += 1
            if cds_row.ambiguous_codon_count > 0:
                cds_frame_ambiguous_codon_total += 1
            if cds_row.translation_qualifier_present and cds_row.translation_matches_qualifier is False:
                cds_frame_translation_mismatch_total += 1

        elif ftype == "tRNA":
            trna_total += 1
            annotation_intervals.append(get_feature_span(feat))

        elif ftype == "rRNA":
            rrna_total += 1
            annotation_intervals.append(get_feature_span(feat))
            subtype = infer_rrna_subtype(feat)
            if subtype == "5S":
                rrna_5s += 1
            elif subtype == "16S":
                rrna_16s += 1
            elif subtype == "23S":
                rrna_23s += 1
            else:
                rrna_other += 1

            partial = (
                has_partial_location(feat)
                or qualifier_contains(feat, "note", r"\bpartial\b")
                or qualifier_contains(feat, "product", r"\bpartial\b")
            )
            if partial:
                rrna_partial += 1
            else:
                rrna_complete += 1

        elif ftype == "tmRNA":
            tmrna_total += 1
            annotation_intervals.append(get_feature_span(feat))

        elif is_nc_rna(feat):
            ncrna_total += 1
            annotation_intervals.append(get_feature_span(feat))

        elif ftype.endswith("RNA"):
            other_rna_total += 1
            annotation_intervals.append(get_feature_span(feat))

    merged_annotation = merge_intervals(annotation_intervals)
    merged_coding = merge_intervals(coding_intervals)

    annotated_bases_bp = sum(e - s for s, e in merged_annotation)
    coding_bases_bp = sum(e - s for s, e in merged_coding)
    coding_density_percent = (100.0 * coding_bases_bp / length_bp) if length_bp else float("nan")

    gaps, overlaps = interval_gaps_and_overlaps(merged_annotation)

    _, raw_overlaps = interval_gaps_and_overlaps(annotation_intervals)

    source = None
    for feat in record.features:
        if feat.type == "source":
            source = feat
            break

    topology = "unknown"
    if source is not None:
        topology_vals = source.qualifiers.get("topology", [])
        if topology_vals:
            topology = topology_vals[0]

    record_type = "replicon"
    descr = (record.description or "").lower()
    if "plasmid" in descr:
        record_type = "plasmid"
    else:
        source_plasmid = source.qualifiers.get("plasmid", []) if source else []
        if source_plasmid:
            record_type = "plasmid"
        else:
            record_type = "chromosome_or_other"

    return RepliconSummary(
        replicon_id=record.id,
        replicon_name=record.name,
        record_type=record_type,
        topology=topology,
        length_bp=length_bp,
        gc_percent=gc,

        genes_total=genes_total,
        cds_total=cds_total,
        cds_with_translation=cds_with_translation,
        pseudogenes_total=pseudogenes_total,

        trna_total=trna_total,
        rrna_total=rrna_total,
        rrna_5s=rrna_5s,
        rrna_16s=rrna_16s,
        rrna_23s=rrna_23s,
        rrna_other=rrna_other,
        rrna_complete=rrna_complete,
        rrna_partial=rrna_partial,

        ncrna_total=ncrna_total,
        tmrna_total=tmrna_total,
        other_rna_total=other_rna_total,

        avg_gene_length_bp=safe_mean(gene_lengths),
        median_gene_length_bp=safe_median(gene_lengths),
        avg_cds_length_bp=safe_mean(cds_lengths),
        median_cds_length_bp=safe_median(cds_lengths),

        annotated_bases_bp=annotated_bases_bp,
        coding_bases_bp=coding_bases_bp,
        coding_density_percent=coding_density_percent,

        intergenic_spacers_n=len(gaps),
        intergenic_spacers_total_bp=sum(gaps),
        intergenic_spacers_mean_bp=safe_mean(gaps),
        intergenic_spacers_median_bp=safe_median(gaps),
        intergenic_spacers_min_bp=min(gaps) if gaps else float("nan"),
        intergenic_spacers_max_bp=max(gaps) if gaps else float("nan"),

        annotated_overlaps_n=len(raw_overlaps),
        annotated_overlaps_total_bp=sum(raw_overlaps),
        annotated_overlaps_mean_bp=safe_mean(raw_overlaps),

        cds_frame_total=cds_frame_total,
        cds_frame_intact_total=cds_frame_intact_total,
        cds_frame_partial_total=cds_frame_partial_total,
        cds_frame_pseudogene_total=cds_frame_pseudogene_total,
        cds_frame_problematic_total=cds_frame_problematic_total,
        cds_frame_error_total=cds_frame_error_total,
        cds_frame_internal_stop_total=cds_frame_internal_stop_total,
        cds_frame_length_error_total=cds_frame_length_error_total,
        cds_frame_missing_start_total=cds_frame_missing_start_total,
        cds_frame_missing_stop_total=cds_frame_missing_stop_total,
        cds_frame_ambiguous_codon_total=cds_frame_ambiguous_codon_total,
        cds_frame_translation_mismatch_total=cds_frame_translation_mismatch_total,
        cds_frame_non_partial_non_pseudogene_intact_percent=(
            100.0 * cds_frame_intact_total / (cds_frame_total - cds_frame_partial_total - cds_frame_pseudogene_total)
            if (cds_frame_total - cds_frame_partial_total - cds_frame_pseudogene_total)
            else float("nan")
        ),
    ), cds_rows


def aggregate_genome(summaries: List[RepliconSummary]) -> Dict[str, float]:
    total_length = sum(s.length_bp for s in summaries)

    if total_length:
        genome_gc = sum(
            s.gc_percent * s.length_bp
            for s in summaries
            if not math.isnan(s.gc_percent)
        ) / total_length
    else:
        genome_gc = float("nan")

    out = {
        "replicons_total": len(summaries),
        "chromosome_or_other_replicons": sum(1 for s in summaries if s.record_type == "chromosome_or_other"),
        "plasmids_total": sum(1 for s in summaries if s.record_type == "plasmid"),
        "length_bp_total": total_length,
        "gc_percent_weighted": genome_gc,

        "genes_total": sum(s.genes_total for s in summaries),
        "cds_total": sum(s.cds_total for s in summaries),
        "cds_with_translation_total": sum(s.cds_with_translation for s in summaries),
        "pseudogenes_total": sum(s.pseudogenes_total for s in summaries),

        "trna_total": sum(s.trna_total for s in summaries),
        "rrna_total": sum(s.rrna_total for s in summaries),
        "rrna_5s_total": sum(s.rrna_5s for s in summaries),
        "rrna_16s_total": sum(s.rrna_16s for s in summaries),
        "rrna_23s_total": sum(s.rrna_23s for s in summaries),
        "rrna_other_total": sum(s.rrna_other for s in summaries),
        "rrna_complete_total": sum(s.rrna_complete for s in summaries),
        "rrna_partial_total": sum(s.rrna_partial for s in summaries),

        "ncrna_total": sum(s.ncrna_total for s in summaries),
        "tmrna_total": sum(s.tmrna_total for s in summaries),
        "other_rna_total": sum(s.other_rna_total for s in summaries),

        "annotated_bases_bp_total": sum(s.annotated_bases_bp for s in summaries),
        "coding_bases_bp_total": sum(s.coding_bases_bp for s in summaries),

        "intergenic_spacers_n_total": sum(s.intergenic_spacers_n for s in summaries),
        "intergenic_spacers_bp_total": sum(s.intergenic_spacers_total_bp for s in summaries),

        "annotated_overlaps_n_total": sum(s.annotated_overlaps_n for s in summaries),
        "annotated_overlaps_bp_total": sum(s.annotated_overlaps_total_bp for s in summaries),

        "cds_frame_total": sum(s.cds_frame_total for s in summaries),
        "cds_frame_intact_total": sum(s.cds_frame_intact_total for s in summaries),
        "cds_frame_partial_total": sum(s.cds_frame_partial_total for s in summaries),
        "cds_frame_pseudogene_total": sum(s.cds_frame_pseudogene_total for s in summaries),
        "cds_frame_problematic_total": sum(s.cds_frame_problematic_total for s in summaries),
        "cds_frame_error_total": sum(s.cds_frame_error_total for s in summaries),
        "cds_frame_internal_stop_total": sum(s.cds_frame_internal_stop_total for s in summaries),
        "cds_frame_length_error_total": sum(s.cds_frame_length_error_total for s in summaries),
        "cds_frame_missing_start_total": sum(s.cds_frame_missing_start_total for s in summaries),
        "cds_frame_missing_stop_total": sum(s.cds_frame_missing_stop_total for s in summaries),
        "cds_frame_ambiguous_codon_total": sum(s.cds_frame_ambiguous_codon_total for s in summaries),
        "cds_frame_translation_mismatch_total": sum(s.cds_frame_translation_mismatch_total for s in summaries),
    }

    if total_length:
        out["coding_density_percent"] = 100.0 * out["coding_bases_bp_total"] / total_length
    else:
        out["coding_density_percent"] = float("nan")

    def weighted_avg(attr: str, weight_attr: str = "length_bp") -> float:
        nums = []
        dens = []
        for s in summaries:
            val = getattr(s, attr)
            w = getattr(s, weight_attr)
            if not math.isnan(val):
                nums.append(val * w)
                dens.append(w)
        return sum(nums) / sum(dens) if dens and sum(dens) else float("nan")

    out["avg_gene_length_bp_weighted"] = weighted_avg("avg_gene_length_bp")
    out["avg_cds_length_bp_weighted"] = weighted_avg("avg_cds_length_bp")

    out["intergenic_spacers_mean_bp"] = (
        out["intergenic_spacers_bp_total"] / out["intergenic_spacers_n_total"]
        if out["intergenic_spacers_n_total"] else float("nan")
    )
    out["annotated_overlaps_mean_bp"] = (
        out["annotated_overlaps_bp_total"] / out["annotated_overlaps_n_total"]
        if out["annotated_overlaps_n_total"] else float("nan")
    )

    out["genes_per_mbp"] = (out["genes_total"] / total_length * 1e6) if total_length else float("nan")
    out["cds_per_mbp"] = (out["cds_total"] / total_length * 1e6) if total_length else float("nan")
    out["trna_per_mbp"] = (out["trna_total"] / total_length * 1e6) if total_length else float("nan")
    out["rrna_operon_proxy_min"] = min(
        out["rrna_5s_total"], out["rrna_16s_total"], out["rrna_23s_total"]
    )

    eligible = out["cds_frame_total"] - out["cds_frame_partial_total"] - out["cds_frame_pseudogene_total"]
    out["cds_frame_non_partial_non_pseudogene_intact_percent"] = (
        100.0 * out["cds_frame_intact_total"] / eligible if eligible else float("nan")
    )

    return out


# ------------------------------- output -------------------------------- #

def write_replicon_tsv(summaries: List[RepliconSummary], path: str) -> None:
    if not summaries:
        raise ValueError("No summaries to write")

    fieldnames = list(summaries[0].__dataclass_fields__.keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for s in summaries:
            formatted = {k: format_value(v) for k, v in s.__dict__.items()}
            writer.writerow(formatted)


def write_cds_tsv(cds_rows: List[CDSFrameIntegrity], path: str) -> None:
    fieldnames = list(CDSFrameIntegrity.__dataclass_fields__.keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in cds_rows:
            writer.writerow({k: format_value(v) for k, v in row.__dict__.items()})


def write_genome_tsv(genome_summary: Dict[str, float], path: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in genome_summary.items():
            writer.writerow([key, format_value(value)])


def print_human_readable(genome_summary: Dict[str, float]) -> None:
    print("Genome annotation summary")
    print("--------------")
    for key, value in genome_summary.items():
        print(f"{key}: {format_value(value)}")


# ------------------------------- main ---------------------------------- #

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Evaluates annotated CDS reading-frame integrity in a bacterial genome.")
    p.add_argument(
        "-i", "--input",
        required=True,
        help="Input GenBank flatfile (.gb, .gbk, .gbff)"
    )
    p.add_argument(
        "-o", "--out-prefix",
        default=None,
        help="Output prefix (default: input filename without extension)"
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    in_path = args.input
    out_prefix = args.out_prefix
    if out_prefix is None:
        base = os.path.basename(in_path)
        out_prefix = re.sub(r"\.(gb|gbk|gbff|genbank)$", "", base, flags=re.IGNORECASE)

    records = list(SeqIO.parse(in_path, "genbank"))
    if not records:
        raise SystemExit(f"No GenBank records found in {in_path!r}")

    record_summaries = [summarize_record(r) for r in records]
    summaries = [summary for summary, _ in record_summaries]
    cds_rows = [row for _, rows in record_summaries for row in rows]
    genome_summary = aggregate_genome(summaries)

    replicon_tsv = f"{out_prefix}.replicons.tsv"
    genome_tsv = f"{out_prefix}.genome.tsv"
    cds_tsv = f"{out_prefix}.cds_frame_integrity.tsv"

    write_replicon_tsv(summaries, replicon_tsv)
    write_genome_tsv(genome_summary, genome_tsv)
    write_cds_tsv(cds_rows, cds_tsv)
    print_human_readable(genome_summary)

    print("\nWrote:")
    print(f"  {replicon_tsv}")
    print(f"  {genome_tsv}")
    print(f"  {cds_tsv}")


if __name__ == "__main__":
    main()
