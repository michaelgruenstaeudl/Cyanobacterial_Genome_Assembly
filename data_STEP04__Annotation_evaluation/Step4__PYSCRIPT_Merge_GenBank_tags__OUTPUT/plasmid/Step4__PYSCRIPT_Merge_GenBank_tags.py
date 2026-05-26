#!/usr/bin/env python3

import argparse
import copy
from Bio import SeqIO

TRANSL_TABLE = 11
NEARBY_BP = 10

SKIP_QUALIFIERS = {"standard_name"}

GENE_SKIP_QUALIFIERS = {
    "protein_id",
    "translation",
    "transl_table",
    "codon_start",
    "product",
}


def log(msg):
    print(msg, flush=True)


def normalize_gene_feature(feature):
    if feature.type == "Gene":
        feature.type = "gene"


def coords(feature):
    return (
        int(feature.location.start),
        int(feature.location.end),
        feature.location.strand,
    )


def exact_key(feature):
    start, end, strand = coords(feature)
    return feature.type, start, end, strand


def is_near_duplicate_same_type(pgap_feature, donor_feature):
    pgap_start, pgap_end, pgap_strand = coords(pgap_feature)
    donor_start, donor_end, donor_strand = coords(donor_feature)

    return (
        pgap_feature.type == donor_feature.type
        and pgap_strand == donor_strand
        and abs(pgap_start - donor_start) <= NEARBY_BP
        and abs(pgap_end - donor_end) <= NEARBY_BP
    )


def add_valid_cds_translation(feature, record):
    if feature.type != "CDS":
        return

    dna = feature.extract(record.seq)

    try:
        protein = dna.translate(table=TRANSL_TABLE, cds=True)
        log("      CDS translated with cds=True")
    except Exception as e:
        log(f"      cds=True failed: {e}")
        protein = dna.translate(table=TRANSL_TABLE, to_stop=True)
        log("      CDS translated with to_stop=True fallback")

    feature.qualifiers["transl_table"] = [str(TRANSL_TABLE)]
    feature.qualifiers["translation"] = [str(protein)]


def remove_skipped_qualifiers(feature):
    for qualifier in SKIP_QUALIFIERS:
        if qualifier in feature.qualifiers:
            del feature.qualifiers[qualifier]
            log(f"    Removed skipped qualifier /{qualifier}")

    if feature.type == "gene":
        for qualifier in GENE_SKIP_QUALIFIERS:
            if qualifier in feature.qualifiers:
                del feature.qualifiers[qualifier]
                log(f"    Removed invalid gene qualifier /{qualifier}")


def transfer_missing_qualifiers(pgap_feature, donor_feature):
    added = 0
    transferred_to_cds = False

    for qualifier, values in donor_feature.qualifiers.items():
        if qualifier in SKIP_QUALIFIERS:
            log(f"    Skipped /{qualifier}; never transferred")
            continue

        if pgap_feature.type == "gene" and qualifier in GENE_SKIP_QUALIFIERS:
            log(f"    Skipped /{qualifier}; not valid for gene feature")
            continue

        if qualifier not in pgap_feature.qualifiers:
            pgap_feature.qualifiers[qualifier] = list(values)
            added += 1
            log(f"    Added qualifier /{qualifier}")

            if pgap_feature.type == "CDS":
                transferred_to_cds = True
        else:
            log(f"    Already present /{qualifier}")

    return added, transferred_to_cds


def merge_and_add_annotations(pgap_record, donor_record):
    pgap_by_exact = {}

    for pgap_feature in pgap_record.features:
        normalize_gene_feature(pgap_feature)
        remove_skipped_qualifiers(pgap_feature)
        pgap_by_exact.setdefault(exact_key(pgap_feature), []).append(pgap_feature)

    added_qualifiers = 0
    added_features = 0
    skipped_nearby = 0

    for donor_feature in donor_record.features:
        normalize_gene_feature(donor_feature)
        remove_skipped_qualifiers(donor_feature)

        key = exact_key(donor_feature)

        if key in pgap_by_exact:
            pgap_feature = pgap_by_exact[key][0]
            log(f"Exact same-type match found for {donor_feature.type}")

            added, transferred_to_cds = transfer_missing_qualifiers(
                pgap_feature,
                donor_feature,
            )
            added_qualifiers += added

            if transferred_to_cds:
                add_valid_cds_translation(pgap_feature, pgap_record)

            continue

        nearby_same_type = any(
            is_near_duplicate_same_type(pgap_feature, donor_feature)
            for pgap_feature in pgap_record.features
        )

        if nearby_same_type:
            skipped_nearby += 1
            log(
                f"Skipped donor {donor_feature.type}; same-type PGAP feature "
                f"has both start and stop within {NEARBY_BP} bp"
            )
            continue

        new_feature = copy.deepcopy(donor_feature)
        normalize_gene_feature(new_feature)
        remove_skipped_qualifiers(new_feature)

        if new_feature.type == "CDS":
            add_valid_cds_translation(new_feature, pgap_record)

        pgap_record.features.append(new_feature)
        pgap_by_exact.setdefault(exact_key(new_feature), []).append(new_feature)

        added_features += 1
        log(f"Added new {new_feature.type} feature")

    return added_qualifiers, added_features, skipped_nearby


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge annotations from a donor GenBank file into a PGAP GenBank file. "
            "Exact same-type coordinate matches receive missing qualifiers. "
            "Donor features are added only when no PGAP feature of the same type "
            "has both start and stop coordinates within 10 bp. "
            "/standard_name is never transferred. Protein-specific qualifiers "
            "are not added to gene features."
        )
    )

    parser.add_argument("pgap_file")
    parser.add_argument("donor_file")
    parser.add_argument("output_file")

    args = parser.parse_args()

    log("Reading input files")

    pgap_records = list(SeqIO.parse(args.pgap_file, "genbank"))
    donor_records = list(SeqIO.parse(args.donor_file, "genbank"))

    if len(pgap_records) != len(donor_records):
        raise ValueError("Record count mismatch between files")

    total_qualifiers = 0
    total_features = 0
    total_skipped_nearby = 0

    for pgap_record, donor_record in zip(pgap_records, donor_records):
        q, f, s = merge_and_add_annotations(pgap_record, donor_record)
        total_qualifiers += q
        total_features += f
        total_skipped_nearby += s

    log("Writing output file")
    SeqIO.write(pgap_records, args.output_file, "genbank")

    log("Done.")
    log(f"Total qualifiers added: {total_qualifiers}")
    log(f"Total features added: {total_features}")
    log(
        "Total donor features skipped due to nearby same-type PGAP feature: "
        f"{total_skipped_nearby}"
    )
    log(f"Output written to: {args.output_file}")


if __name__ == "__main__":
    main()
