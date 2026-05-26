#!/usr/bin/env python3

import argparse
import copy
import difflib
from Bio import SeqIO

TRANSL_TABLE = 11
NEARBY_BP = 10
SAME_TAG_NEARBY_BP = 100
SIMILAR_PROTEIN_MAX_GAP_BP = 300
SIMILAR_PROTEIN_MIN_RATIO = 0.90
SIMILAR_PROTEIN_MAX_LEN_DIFF_FRAC = 0.20

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


def pair_key(feature):
    start, end, strand = coords(feature)
    return start, end, strand


def pair_key_no_strand(feature):
    start, end, _ = coords(feature)
    return start, end


def is_near_duplicate_same_type(pgap_feature, donor_feature):
    pgap_start, pgap_end, pgap_strand = coords(pgap_feature)
    donor_start, donor_end, donor_strand = coords(donor_feature)

    return (
        pgap_feature.type == donor_feature.type
        and pgap_strand == donor_strand
        and abs(pgap_start - donor_start) <= NEARBY_BP
        and abs(pgap_end - donor_end) <= NEARBY_BP
    )


def has_shared_qualifier_value(feature_a, feature_b, qualifier_name):
    values_a = feature_a.qualifiers.get(qualifier_name)
    values_b = feature_b.qualifiers.get(qualifier_name)

    if not values_a or not values_b:
        return False

    return bool(set(values_a) & set(values_b))


def is_invalid_gene_name(value):
    normalized = str(value).strip().lower()
    return normalized == "hypothetical protein cds"


def is_near_duplicate_same_tag(pgap_feature, donor_feature):
    pgap_start, pgap_end, _ = coords(pgap_feature)
    donor_start, donor_end, _ = coords(donor_feature)

    if not (
        abs(pgap_start - donor_start) <= SAME_TAG_NEARBY_BP
        and abs(pgap_end - donor_end) <= SAME_TAG_NEARBY_BP
    ):
        return False

    if has_shared_qualifier_value(pgap_feature, donor_feature, "gene"):
        return True
    if has_shared_qualifier_value(pgap_feature, donor_feature, "protein_id"):
        return True
    if has_shared_qualifier_value(pgap_feature, donor_feature, "locus_tag"):
        return True
    if has_shared_qualifier_value(pgap_feature, donor_feature, "standard_name"):
        return True

    return False


def interval_gap(start_a, end_a, start_b, end_b):
    if end_a <= start_b:
        return start_b - end_a
    if end_b <= start_a:
        return start_a - end_b
    return 0


def normalize_translation_string(value):
    return "".join(str(value).split()).rstrip("*")


def get_cds_translation(feature, record):
    translation_values = feature.qualifiers.get("translation")
    if translation_values:
        cleaned = normalize_translation_string(translation_values[0])
        if cleaned:
            return cleaned

    dna = feature.extract(record.seq)
    try:
        protein = dna.translate(table=TRANSL_TABLE, cds=True)
    except Exception:
        protein = dna.translate(table=TRANSL_TABLE, to_stop=True)

    return normalize_translation_string(str(protein))


def get_feature_translation(feature, record, cds_features):
    if feature.type == "CDS":
        return get_cds_translation(feature, record)

    if feature.type != "gene":
        return ""

    feature_start, feature_end, feature_strand = coords(feature)
    for cds_feature in cds_features:
        cds_start, cds_end, cds_strand = coords(cds_feature)
        if (
            feature_start == cds_start
            and feature_end == cds_end
            and feature_strand == cds_strand
        ):
            return get_cds_translation(cds_feature, record)

    return ""


def are_proteins_too_similar(seq_a, seq_b):
    if not seq_a or not seq_b:
        return False

    longer = max(len(seq_a), len(seq_b))
    shorter = min(len(seq_a), len(seq_b))
    if longer == 0:
        return False

    if (longer - shorter) / longer > SIMILAR_PROTEIN_MAX_LEN_DIFF_FRAC:
        return False

    ratio = difflib.SequenceMatcher(a=seq_a, b=seq_b, autojunk=False).ratio()
    return ratio >= SIMILAR_PROTEIN_MIN_RATIO


def is_similar_cds_duplicate(
    donor_feature,
    donor_record,
    donor_cds_features,
    pgap_cds_feature,
    pgap_record,
):
    donor_start, donor_end, donor_strand = coords(donor_feature)
    pgap_start, pgap_end, pgap_strand = coords(pgap_cds_feature)

    if donor_strand != pgap_strand:
        return False

    if (
        interval_gap(donor_start, donor_end, pgap_start, pgap_end)
        > SIMILAR_PROTEIN_MAX_GAP_BP
    ):
        return False

    donor_translation = get_feature_translation(
        donor_feature,
        donor_record,
        donor_cds_features,
    )
    pgap_translation = get_cds_translation(pgap_cds_feature, pgap_record)

    return are_proteins_too_similar(donor_translation, pgap_translation)


def get_skip_reason(
    donor_feature,
    pgap_record,
    donor_record,
    donor_cds_features,
    pgap_cds_features,
):
    nearby_same_type = any(
        is_near_duplicate_same_type(pgap_feature, donor_feature)
        for pgap_feature in pgap_record.features
    )
    if nearby_same_type:
        return (
            f"same-type PGAP feature has both start and stop within "
            f"{NEARBY_BP} bp"
        )

    nearby_same_tag = any(
        is_near_duplicate_same_tag(pgap_feature, donor_feature)
        for pgap_feature in pgap_record.features
    )
    if nearby_same_tag:
        return (
            "PGAP feature with same gene, protein_id, locus_tag, or "
            f"standard_name has both start and stop within {SAME_TAG_NEARBY_BP} bp"
        )

    if donor_feature.type in {"CDS", "gene"}:
        similar_cds_duplicate = any(
            is_similar_cds_duplicate(
                donor_feature,
                donor_record,
                donor_cds_features,
                pgap_cds_feature,
                pgap_record,
            )
            for pgap_cds_feature in pgap_cds_features
        )
        if similar_cds_duplicate:
            return (
                "translation is highly similar to nearby PGAP CDS within "
                f"{SIMILAR_PROTEIN_MAX_GAP_BP} bp"
            )

    return None


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
        if "gene" in feature.qualifiers:
            original_gene_values = feature.qualifiers["gene"]
            filtered_gene_values = [
                value
                for value in original_gene_values
                if not is_invalid_gene_name(value)
            ]

            if not filtered_gene_values:
                del feature.qualifiers["gene"]
                log("    Removed invalid gene qualifier value from /gene")
            elif len(filtered_gene_values) != len(original_gene_values):
                feature.qualifiers["gene"] = filtered_gene_values
                log("    Removed invalid gene qualifier value from /gene")

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

        if pgap_feature.type == "gene" and qualifier == "gene":
            valid_gene_values = [
                value for value in values if not is_invalid_gene_name(value)
            ]
            if not valid_gene_values:
                log("    Skipped /gene; invalid CDS-like gene name")
                continue
            values = valid_gene_values

        if qualifier not in pgap_feature.qualifiers:
            pgap_feature.qualifiers[qualifier] = list(values)
            added += 1
            log(f"    Added qualifier /{qualifier}")

            if pgap_feature.type == "CDS":
                transferred_to_cds = True
        else:
            log(f"    Already present /{qualifier}")

    return added, transferred_to_cds


def remove_unpaired_gene_cds_features(record):
    paired_by_coords = {}
    gene_keys = set()
    for feature in record.features:
        normalize_gene_feature(feature)
        if feature.type == "gene":
            gene_keys.add(pair_key_no_strand(feature))

        if feature.type not in {"gene", "CDS"}:
            continue

        key = pair_key(feature)
        paired_by_coords.setdefault(key, {"gene": 0, "CDS": 0})
        paired_by_coords[key][feature.type] += 1

    filtered_features = []
    removed_count = 0
    added_gene_for_rna_count = 0

    for feature in record.features:
        normalize_gene_feature(feature)
        feature_type_lower = feature.type.lower()
        if feature_type_lower in {"trna", "rrna"}:
            key = pair_key_no_strand(feature)
            if key not in gene_keys:
                new_gene = copy.deepcopy(feature)
                new_gene.type = "gene"
                normalize_gene_feature(new_gene)
                remove_skipped_qualifiers(new_gene)

                filtered_features.append(new_gene)
                gene_keys.add(key)
                added_gene_for_rna_count += 1
                log(
                    f"Added matching gene for unpaired {feature.type}; "
                    "used same coordinates and strand"
                )

            filtered_features.append(feature)
            continue

        if feature.type not in {"gene", "CDS"}:
            filtered_features.append(feature)
            continue

        key = pair_key(feature)
        counts = paired_by_coords.get(key, {"gene": 0, "CDS": 0})
        if counts["gene"] > 0 and counts["CDS"] > 0:
            filtered_features.append(feature)
            continue

        removed_count += 1
        log(
            f"Removed unpaired {feature.type}; no same-coordinate "
            f"gene/CDS counterpart"
        )

    if removed_count or added_gene_for_rna_count:
        record.features = filtered_features

    return removed_count, added_gene_for_rna_count


def enforce_rna_gene_pairs(record):
    gene_keys = {
        pair_key_no_strand(feature)
        for feature in record.features
        if feature.type == "gene"
    }

    added_count = 0
    features_to_add = []

    for feature in record.features:
        feature_type_lower = feature.type.lower()
        if feature_type_lower not in {"trna", "rrna"}:
            continue

        key = pair_key_no_strand(feature)
        if key in gene_keys:
            continue

        new_gene = copy.deepcopy(feature)
        new_gene.type = "gene"
        normalize_gene_feature(new_gene)
        remove_skipped_qualifiers(new_gene)

        features_to_add.append(new_gene)
        gene_keys.add(key)
        added_count += 1
        log(
            f"Added missing gene for {feature.type} in verification pass; "
            "used same start/stop coordinates"
        )

    if features_to_add:
        record.features.extend(features_to_add)

    return added_count


def merge_and_add_annotations(pgap_record, donor_record):
    pgap_by_exact = {}
    pgap_cds_features = []
    donor_cds_features = []
    donor_pairs = {}

    for pgap_feature in pgap_record.features:
        normalize_gene_feature(pgap_feature)
        pgap_by_exact.setdefault(exact_key(pgap_feature), []).append(pgap_feature)
        if pgap_feature.type == "CDS":
            pgap_cds_features.append(pgap_feature)

    for donor_feature in donor_record.features:
        normalize_gene_feature(donor_feature)
        if donor_feature.type == "CDS":
            donor_cds_features.append(donor_feature)

        if donor_feature.type in {"gene", "CDS"}:
            key = pair_key(donor_feature)
            donor_pairs.setdefault(key, {"gene": [], "CDS": []})
            donor_pairs[key][donor_feature.type].append(donor_feature)

    added_qualifiers = 0
    added_features = 0
    skipped_nearby = 0

    for donor_feature in donor_record.features:
        donor_pair_partner = None
        if donor_feature.type in {"gene", "CDS"}:
            key = pair_key(donor_feature)
            pair_bucket = donor_pairs.get(key, {})
            partner_type = "CDS" if donor_feature.type == "gene" else "gene"
            partner_candidates = pair_bucket.get(partner_type, [])
            if partner_candidates:
                donor_pair_partner = partner_candidates[0]

            if donor_pair_partner is None:
                skipped_nearby += 1
                log(
                    f"Skipped donor {donor_feature.type}; no same-length "
                    f"donor {partner_type} counterpart for paired transfer"
                )
                continue

            donor_exact_exists = exact_key(donor_feature) in pgap_by_exact
            partner_exact_exists = exact_key(donor_pair_partner) in pgap_by_exact

            donor_skip_reason = None if donor_exact_exists else get_skip_reason(
                donor_feature,
                pgap_record,
                donor_record,
                donor_cds_features,
                pgap_cds_features,
            )
            partner_skip_reason = None if partner_exact_exists else get_skip_reason(
                donor_pair_partner,
                pgap_record,
                donor_record,
                donor_cds_features,
                pgap_cds_features,
            )

            if donor_skip_reason or partner_skip_reason:
                skipped_nearby += 1
                reason = donor_skip_reason or partner_skip_reason
                log(
                    f"Skipped donor {donor_feature.type}; paired gene/CDS "
                    f"transfer failed because {reason}"
                )
                continue

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

        skip_reason = get_skip_reason(
            donor_feature,
            pgap_record,
            donor_record,
            donor_cds_features,
            pgap_cds_features,
        )
        if skip_reason:
            skipped_nearby += 1
            log(
                f"Skipped donor {donor_feature.type}; {skip_reason}"
            )
            continue

        new_feature = copy.deepcopy(donor_feature)
        normalize_gene_feature(new_feature)
        remove_skipped_qualifiers(new_feature)

        if new_feature.type == "CDS":
            add_valid_cds_translation(new_feature, pgap_record)
            pgap_cds_features.append(new_feature)

        pgap_record.features.append(new_feature)
        pgap_by_exact.setdefault(exact_key(new_feature), []).append(new_feature)

        added_features += 1
        log(f"Added new {new_feature.type} feature")

    removed_unpaired, added_rna_pair_genes = remove_unpaired_gene_cds_features(
        pgap_record
    )
    verification_added_rna_pair_genes = enforce_rna_gene_pairs(pgap_record)

    return (
        added_qualifiers,
        added_features,
        skipped_nearby,
        removed_unpaired,
        added_rna_pair_genes + verification_added_rna_pair_genes,
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge annotations from a donor GenBank file into a PGAP GenBank file. "
            "Exact same-type coordinate matches receive missing qualifiers. "
            "Donor features are added only when no PGAP feature of the same type "
            "has both start and stop coordinates within 10 bp, and when no PGAP "
            "feature shares the same /gene, /protein_id, /locus_tag, or "
            "/standard_name qualifier with both start and stop coordinates "
            "within 100 bp. Donor gene/CDS features are also skipped when their "
            "protein translation is highly similar to a same-strand nearby PGAP "
            "CDS feature. Donor gene/CDS features are transferred only as "
            "same-coordinate pairs (or not transferred at all). "
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
    total_removed_unpaired = 0
    total_added_rna_pair_genes = 0

    for pgap_record, donor_record in zip(pgap_records, donor_records):
        q, f, s, r, a = merge_and_add_annotations(pgap_record, donor_record)
        total_qualifiers += q
        total_features += f
        total_skipped_nearby += s
        total_removed_unpaired += r
        total_added_rna_pair_genes += a

    log("Writing output file")
    SeqIO.write(pgap_records, args.output_file, "genbank")

    log("Done.")
    log(f"Total qualifiers added: {total_qualifiers}")
    log(f"Total features added: {total_features}")
    log(
        "Total donor features skipped by duplication/pairing filters: "
        f"{total_skipped_nearby}"
    )
    log(
        "Total orphan gene/CDS features removed in final cleanup: "
        f"{total_removed_unpaired}"
    )
    log(
        "Total gene features added to pair unpaired tRNA/rRNA: "
        f"{total_added_rna_pair_genes}"
    )
    log(f"Output written to: {args.output_file}")


if __name__ == "__main__":
    main()
