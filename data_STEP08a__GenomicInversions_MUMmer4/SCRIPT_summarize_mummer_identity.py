#!/usr/bin/env python3

import csv
import sys


def fasta_length(path):
    length = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def summarize_coords(path, reference_length):
    total = 0
    weighted_identity = 0.0
    forward = 0
    reverse = 0

    with open(path, encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 7:
                continue

            len1 = int(row[4])
            pidy = float(row[6])
            s2 = int(row[2])
            e2 = int(row[3])

            total += len1
            weighted_identity += len1 * pidy

            if s2 < e2:
                forward += len1
            elif s2 > e2:
                reverse += len1

    identity = weighted_identity / total if total else 0.0
    divergence = 100.0 - identity
    aligned_fraction = 100.0 * total / reference_length if reference_length else 0.0
    forward_percent = 100.0 * forward / total if total else 0.0
    reverse_percent = 100.0 * reverse / total if total else 0.0

    return total, identity, divergence, aligned_fraction, forward, forward_percent, reverse, reverse_percent


if len(sys.argv) != 5:
    sys.exit(
        "Usage:\n"
        "  python SCRIPT_summarize_mummer_identity.py "
        "REF.fasta QRY.fasta original.coords.tsv inversion_corrected.coords.tsv"
    )


ref_fasta = sys.argv[1]
qry_fasta = sys.argv[2]
original_coords = sys.argv[3]
corrected_coords = sys.argv[4]

ref_len = fasta_length(ref_fasta)
qry_len = fasta_length(qry_fasta)
length_difference = qry_len - ref_len

print(f"reference_chromosome_fasta: {ref_fasta}")
print(f"input_chromosome_fasta: {qry_fasta}")
print("=" * 70)
print("Input chromosome length comparison")
print("=" * 70)
print(f"Reference chromosome length: {ref_len} bp")
print(f"Query chromosome length:     {qry_len} bp")
print(f"Length difference:           {length_difference} bp")

corrected_identity = None
corrected_divergence = None
corrected_aligned_fraction = None

for label, coords in [
    ("Original chromosome orientation", original_coords),
    ("Inversion-corrected chromosome orientation", corrected_coords),
]:
    total, identity, divergence, aligned_fraction, forward, forward_percent, reverse, reverse_percent = summarize_coords(
        coords,
        ref_len,
    )

    print()
    print("=" * 70)
    print(label)
    print("=" * 70)
    print(f"Total aligned reference bases: {total} bp")
    print(f"Forward-orientation aligned bases: {forward} bp ({forward_percent:.2f}%)")
    print(f"Reverse-orientation aligned bases: {reverse} bp ({reverse_percent:.2f}%)")

    if label == "Inversion-corrected chromosome orientation":
        corrected_identity = identity
        corrected_divergence = divergence
        corrected_aligned_fraction = aligned_fraction

print()
print("=" * 70)
print("Summary values after inversion correction")
print("=" * 70)
print(f"Alignment-length-weighted nucleotide identity: {corrected_identity:.4f}%")
print(f"Alignment-length-weighted nucleotide divergence: {corrected_divergence:.4f}%")
print(f"Reference-aligned fraction: {corrected_aligned_fraction:.4f}%")
