#!/usr/bin/env python3

# Compare annotated genes between two GenBank files by identifier and genomic position.
#
# A gene is eligible for comparison if, and only if, it has an identical matching
# feature of type CDS, tRNA, or rRNA in the same genome.
#
# "Identical matching" means:
# - same record ID
# - same chosen identifier
# - same start
# - same end
# - same strand
#
# Identifiers are chosen with this priority:
# standard_name > locus_tag > gene > protein_id
#
# For standard_name, the last word is removed before comparison. This is useful
# when otherwise identical annotations end with different suffixes such as
# "gene" versus "CDS".
#
# Between the two genomes, two eligible genes are considered a match when:
# 1) they have the same identifier
# 2) the absolute difference between their start positions is within a
#    user-defined threshold (default: 500 bp)
#
# Genes that do not have a matching CDS, tRNA, or rRNA annotation are not compared.
# Instead, they are reported as warnings.
#
# The script:
# - parses gene, CDS, tRNA, and rRNA annotations from both GenBank files
# - keeps only gene features that have an identical CDS, tRNA, or rRNA annotation
# - compares eligible genes between the two genomes
# - writes a CSV table of matched and unmatched eligible genes
# - writes a CSV table of warnings for unsupported gene annotations

from Bio import SeqIO
import csv
import argparse
import re
from collections import defaultdict


IDENTIFIER_PRIORITY = ("standard_name", "locus_tag", "gene", "protein_id")
SUPPORTING_TYPES = {"CDS", "tRNA", "rRNA"}


def normalize_standard_name(name):
    """
    Normalize standard_name by removing the last word.

    Examples:
        "abc gene" -> "abc"
        "abc CDS"  -> "abc"

    Args:
        name (str): Raw standard_name value.

    Returns:
        str: Normalized standard_name.
    """
    parts = name.strip().split()
    if len(parts) > 1:
        return " ".join(parts[:-1])
    return name.strip()


def normalize_rna_label_suffix(value):
    """
    Remove a trailing standalone RNA label token used in some gene annotations.

    Examples:
        "thrT tRNA" -> "thrT"
        "rrsA rRNA" -> "rrsA"
    """
    return re.sub(r"\s+(?:trna|rrna)\s*$", "", value, flags=re.IGNORECASE).strip()


def get_feature_identifier(feature):
    """
    Return the best available identifier for a feature.

    Priority:
    standard_name > locus_tag > gene > protein_id

    For standard_name, the last word is removed before comparison.

    Args:
        feature: A Biopython SeqFeature object.

    Returns:
        tuple[str, str] | tuple[None, None]:
            (identifier_type, identifier_value) or (None, None)
    """
    qualifiers = feature.qualifiers

    for key in IDENTIFIER_PRIORITY:
        values = qualifiers.get(key, [])
        if values and values[0].strip():
            value = values[0].strip()

            if key == "standard_name":
                value = normalize_standard_name(value)

            # Gene annotations may include a trailing "tRNA"/"rRNA" token while
            # corresponding tRNA/rRNA features often omit it.
            value = normalize_rna_label_suffix(value)

            return key, value

    return None, None


def extract_features(file_path):
    """
    Parse a GenBank file and extract gene/CDS/tRNA/rRNA features with identifiers.

    Args:
        file_path (str): Path to a GenBank file.

    Returns:
        dict[str, list[dict]]: Features grouped by feature type.
    """
    feature_buckets = defaultdict(list)

    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type not in {"gene", "CDS", "tRNA", "rRNA"}:
                    continue

                identifier_type, identifier = get_feature_identifier(feature)
                if identifier is None:
                    continue

                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand

                feature_buckets[feature.type].append({
                    "record_id": record.id,
                    "feature_type": feature.type,
                    "identifier_type": identifier_type,
                    "identifier": identifier,
                    "start": start,
                    "end": end,
                    "strand": strand,
                })

    return feature_buckets


def validate_genes(features_by_type):
    """
    Keep only gene features that have an identical CDS, tRNA, or rRNA annotation.

    "Identical" means same:
    - record_id
    - identifier
    - start
    - end
    - strand

    Args:
        features_by_type (dict[str, list[dict]]): Parsed features by type.

    Returns:
        tuple[list[dict], list[dict]]:
            eligible_genes, warning_genes
    """
    supporting_index = {}
    for feature_type in SUPPORTING_TYPES:
        for feat in features_by_type.get(feature_type, []):
            key = (
                feat["record_id"],
                feat["identifier"],
                feat["start"],
                feat["end"],
                feat["strand"],
            )
            supporting_index.setdefault(key, set()).add(feature_type)

    eligible_genes = []
    warning_genes = []

    for gene in features_by_type.get("gene", []):
        key = (
            gene["record_id"],
            gene["identifier"],
            gene["start"],
            gene["end"],
            gene["strand"],
        )

        matched_support_types = sorted(supporting_index.get(key, set()))
        if matched_support_types:
            gene_copy = dict(gene)
            gene_copy["supporting_feature_types"] = ",".join(matched_support_types)
            eligible_genes.append(gene_copy)
        else:
            gene_copy = dict(gene)
            gene_copy["warning"] = (
                "No identical matching CDS, tRNA, or rRNA annotation found"
            )
            warning_genes.append(gene_copy)

    return eligible_genes, warning_genes


def group_by_identifier(genes):
    """
    Group gene records by identifier value.

    Args:
        genes (list[dict]): Eligible gene records.

    Returns:
        dict[str, list[dict]]
    """
    grouped = defaultdict(list)
    for gene in genes:
        grouped[gene["identifier"]].append(gene)
    return grouped


def match_genes_by_position(genes1, genes2, max_start_diff=500):
    """
    Match eligible genes between two genomes using:
    1. identical identifier value
    2. start-position difference within max_start_diff base pairs

    Greedy one-to-one matching is used per identifier.

    Args:
        genes1 (list[dict]): Eligible genes from genome 1.
        genes2 (list[dict]): Eligible genes from genome 2.
        max_start_diff (int): Maximum allowed absolute difference in start positions.

    Returns:
        tuple[list[dict], list[dict], list[dict]]:
            matched_rows, unmatched_genome1, unmatched_genome2
    """
    grouped1 = group_by_identifier(genes1)
    grouped2 = group_by_identifier(genes2)

    all_identifiers = sorted(set(grouped1) | set(grouped2))

    matched_rows = []
    unmatched_genome1 = []
    unmatched_genome2 = []

    for identifier in all_identifiers:
        entries1 = sorted(grouped1.get(identifier, []), key=lambda x: x["start"])
        entries2 = sorted(grouped2.get(identifier, []), key=lambda x: x["start"])

        used2 = set()

        for g1 in entries1:
            best_j = None
            best_diff = None

            for j, g2 in enumerate(entries2):
                if j in used2:
                    continue

                diff = abs(g1["start"] - g2["start"])
                if diff <= max_start_diff:
                    if best_diff is None or diff < best_diff:
                        best_diff = diff
                        best_j = j

            if best_j is not None:
                g2 = entries2[best_j]
                used2.add(best_j)

                matched_rows.append({
                    "identifier": identifier,
                    "identifier_type_genome1": g1["identifier_type"],
                    "identifier_type_genome2": g2["identifier_type"],
                    "record_id_genome1": g1["record_id"],
                    "record_id_genome2": g2["record_id"],
                    "supporting_types_genome1": g1.get("supporting_feature_types", ""),
                    "supporting_types_genome2": g2.get("supporting_feature_types", ""),
                    "start_genome1": g1["start"],
                    "start_genome2": g2["start"],
                    "end_genome1": g1["end"],
                    "end_genome2": g2["end"],
                    "strand_genome1": g1["strand"],
                    "strand_genome2": g2["strand"],
                    "start_difference_bp": abs(g1["start"] - g2["start"]),
                    "position_match_within_threshold": True,
                    "status": "matched",
                })
            else:
                unmatched_genome1.append(g1)

        for j, g2 in enumerate(entries2):
            if j not in used2:
                unmatched_genome2.append(g2)

    return matched_rows, unmatched_genome1, unmatched_genome2


def build_output_table(matched_rows, unmatched_genome1, unmatched_genome2, threshold):
    """
    Build a unified output table for CSV export.

    Args:
        matched_rows (list[dict]): Successfully matched gene pairs.
        unmatched_genome1 (list[dict]): Unmatched eligible genes from genome 1.
        unmatched_genome2 (list[dict]): Unmatched eligible genes from genome 2.
        threshold (int): Start-position threshold in base pairs.

    Returns:
        list[dict]: Table rows.
    """
    rows = []

    rows.extend(matched_rows)

    for g1 in unmatched_genome1:
        rows.append({
            "identifier": g1["identifier"],
            "identifier_type_genome1": g1["identifier_type"],
            "identifier_type_genome2": "",
            "record_id_genome1": g1["record_id"],
            "record_id_genome2": "",
            "supporting_types_genome1": g1.get("supporting_feature_types", ""),
            "supporting_types_genome2": "",
            "start_genome1": g1["start"],
            "start_genome2": "",
            "end_genome1": g1["end"],
            "end_genome2": "",
            "strand_genome1": g1["strand"],
            "strand_genome2": "",
            "start_difference_bp": "",
            "position_match_within_threshold": False,
            "status": f"unique_to_genome1_or_no_match_within_{threshold}bp",
        })

    for g2 in unmatched_genome2:
        rows.append({
            "identifier": g2["identifier"],
            "identifier_type_genome1": "",
            "identifier_type_genome2": g2["identifier_type"],
            "record_id_genome1": "",
            "record_id_genome2": g2["record_id"],
            "supporting_types_genome1": "",
            "supporting_types_genome2": g2.get("supporting_feature_types", ""),
            "start_genome1": "",
            "start_genome2": g2["start"],
            "end_genome1": "",
            "end_genome2": g2["end"],
            "strand_genome1": "",
            "strand_genome2": g2["strand"],
            "start_difference_bp": "",
            "position_match_within_threshold": False,
            "status": f"unique_to_genome2_or_no_match_within_{threshold}bp",
        })

    return rows


def write_csv(rows, output_csv):
    """
    Write comparison rows to a CSV file.

    Args:
        rows (list[dict]): Output rows.
        output_csv (str): Output file path.
    """
    fieldnames = [
        "identifier",
        "identifier_type_genome1",
        "identifier_type_genome2",
        "record_id_genome1",
        "record_id_genome2",
        "supporting_types_genome1",
        "supporting_types_genome2",
        "start_genome1",
        "start_genome2",
        "end_genome1",
        "end_genome2",
        "strand_genome1",
        "strand_genome2",
        "start_difference_bp",
        "position_match_within_threshold",
        "status",
    ]

    with open(output_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_warnings_csv(warning_rows, output_csv):
    """
    Write warning rows to a CSV file.

    Args:
        warning_rows (list[dict]): Gene features lacking supporting annotations.
        output_csv (str): Output file path.
    """
    fieldnames = [
        "record_id",
        "identifier",
        "identifier_type",
        "start",
        "end",
        "strand",
        "warning",
    ]

    with open(output_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for row in warning_rows:
            writer.writerow({
                "record_id": row["record_id"],
                "identifier": row["identifier"],
                "identifier_type": row["identifier_type"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "warning": row["warning"],
            })


def print_warnings(label, warning_rows, max_examples=10):
    """
    Print warnings for genes without matching CDS/tRNA/rRNA annotations.

    Args:
        label (str): Genome label.
        warning_rows (list[dict]): Warning gene rows.
        max_examples (int): Maximum number of examples to print.
    """
    print(
        f"\nWarnings for {label}: {len(warning_rows)} gene feature(s) lack an identical CDS, tRNA, or rRNA annotation."
    )

    for row in warning_rows[:max_examples]:
        print(
            f"  - {row['record_id']} | {row['identifier_type']}={row['identifier']} "
            f"| {row['start']}-{row['end']} | strand={row['strand']}"
        )

    if len(warning_rows) > max_examples:
        print(f"  ... {len(warning_rows) - max_examples} additional warning(s) not shown")


def print_summary(
    genes1,
    genes2,
    warnings1,
    warnings2,
    matched_rows,
    unmatched_genome1,
    unmatched_genome2,
    threshold,
):
    """
    Print a short comparison summary.
    """
    print("Comparison Summary")
    print("------------------")
    print(f"Eligible genes in Genome 1: {len(genes1)}")
    print(f"Eligible genes in Genome 2: {len(genes2)}")
    print(f"Warnings in Genome 1: {len(warnings1)}")
    print(f"Warnings in Genome 2: {len(warnings2)}")
    print(f"Matched genes (same identifier and within {threshold} bp): {len(matched_rows)}")
    print(f"Unmatched eligible genes from Genome 1: {len(unmatched_genome1)}")
    print(f"Unmatched eligible genes from Genome 2: {len(unmatched_genome2)}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare eligible genes from two GenBank genomes by identifier and "
            "start-position proximity, and warn about gene annotations lacking "
            "matching CDS/tRNA/rRNA annotations."
        )
    )
    parser.add_argument("genbank1", help="Path to genome 1 GenBank file")
    parser.add_argument("genbank2", help="Path to genome 2 GenBank file")
    parser.add_argument(
        "-x",
        "--max-start-diff",
        type=int,
        default=500,
        help="Maximum allowed difference in gene start position (default: 500 bp)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="gene_comparison_table.csv",
        help="Output CSV file for comparison results (default: gene_comparison_table.csv)",
    )
    parser.add_argument(
        "--warnings-output",
        default="gene_annotation_warnings.csv",
        help="Output CSV file for warning rows (default: gene_annotation_warnings.csv)",
    )

    args = parser.parse_args()

    features1 = extract_features(args.genbank1)
    features2 = extract_features(args.genbank2)

    genes1, warnings1 = validate_genes(features1)
    genes2, warnings2 = validate_genes(features2)

    matched_rows, unmatched_genome1, unmatched_genome2 = match_genes_by_position(
        genes1,
        genes2,
        max_start_diff=args.max_start_diff,
    )

    output_rows = build_output_table(
        matched_rows,
        unmatched_genome1,
        unmatched_genome2,
        args.max_start_diff,
    )

    all_warnings = warnings1 + warnings2

    write_csv(output_rows, args.output)
    write_warnings_csv(all_warnings, args.warnings_output)

    print_summary(
        genes1,
        genes2,
        warnings1,
        warnings2,
        matched_rows,
        unmatched_genome1,
        unmatched_genome2,
        args.max_start_diff,
    )
    print_warnings("Genome 1", warnings1)
    print_warnings("Genome 2", warnings2)

    print(f"\nDetailed comparison table written to: {args.output}")
    print(f"Warning table written to: {args.warnings_output}")


if __name__ == "__main__":
    main()
