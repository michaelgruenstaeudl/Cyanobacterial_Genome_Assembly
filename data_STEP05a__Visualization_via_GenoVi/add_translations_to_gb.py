#!/usr/bin/env python3
from __future__ import annotations

import argparse
from typing import Optional, Tuple

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


def _get_int_qual(feature: SeqFeature, key: str, default: int) -> int:
    vals = feature.qualifiers.get(key)
    if not vals:
        return default
    try:
        return int(vals[0])
    except Exception:
        return default


def _infer_translation_for_cds(
    feature: SeqFeature,
    record_seq,
    default_table: int = 11,
    keep_stop: bool = False,
) -> Optional[str]:
    # Extract nucleotide sequence for the feature (handles join/complement)
    try:
        nt = feature.extract(record_seq)
    except Exception:
        return None

    if nt is None or len(nt) == 0:
        return None

    codon_start = _get_int_qual(feature, "codon_start", 1)
    if codon_start < 1:
        codon_start = 1
    if codon_start > 3:
        codon_start = 1

    table = _get_int_qual(feature, "transl_table", default_table)

    # Apply codon_start (1-based)
    nt = nt[codon_start - 1 :]

    # Trim to multiple of 3
    if len(nt) < 3:
        return None
    trim_len = (len(nt) // 3) * 3
    nt = nt[:trim_len]
    if len(nt) < 3:
        return None

    # Translate. to_stop removes terminal stop; internal stops remain as '*'
    try:
        aa = nt.translate(table=table, to_stop=not keep_stop)
    except Exception:
        return None

    s = str(aa)
    if not s:
        return None
    return s


def add_translations(
    in_gb: str,
    out_gb: str,
    default_table: int = 11,
    overwrite: bool = False,
) -> Tuple[int, int, int]:
    total_cds = 0
    added = 0
    skipped = 0

    records = []
    for rec in SeqIO.parse(in_gb, "genbank"):
        if rec.seq is None or len(rec.seq) == 0:
            raise RuntimeError(
                f"Record {rec.id} has no sequence. Cannot infer translations without ORIGIN sequence."
            )

        for feat in rec.features:
            if feat.type != "CDS":
                continue
            total_cds += 1

            has_translation = "translation" in feat.qualifiers and len(feat.qualifiers["translation"]) > 0
            if has_translation and not overwrite:
                skipped += 1
                continue

            aa = _infer_translation_for_cds(
                feat,
                rec.seq,
                default_table=default_table,
                keep_stop=False,
            )
            if aa is None:
                continue

            feat.qualifiers["translation"] = [aa]
            if "transl_table" not in feat.qualifiers:
                feat.qualifiers["transl_table"] = [str(default_table)]
            if "codon_start" not in feat.qualifiers:
                feat.qualifiers["codon_start"] = ["1"]
            added += 1

        records.append(rec)

    SeqIO.write(records, out_gb, "genbank")
    return total_cds, added, skipped


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Add /translation qualifiers to CDS features in a GenBank file using Biopython."
    )
    ap.add_argument("-i", "--input", required=True, help="Input GenBank file")
    ap.add_argument("-o", "--output", required=True, help="Output GenBank file with translations added")
    ap.add_argument("--table", type=int, default=11, help="Default genetic code table (default: 11 for bacteria)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing /translation qualifiers")
    args = ap.parse_args()

    total_cds, added, skipped = add_translations(
        args.input, args.output, default_table=args.table, overwrite=args.overwrite
    )

    print(f"CDS features seen: {total_cds}")
    print(f"Translations added: {added}")
    print(f"CDS skipped (already had translation): {skipped}")
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
