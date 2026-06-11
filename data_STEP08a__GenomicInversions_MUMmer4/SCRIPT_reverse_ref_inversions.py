#!/usr/bin/env python3

from pathlib import Path
import csv
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    records = {}
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    records[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            records[name] = "".join(seq)
    return records


def write_fasta(records, path, width=80):
    with open(path, "w") as out:
        for name, seq in records.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def reverse_inversions(fasta_in, csv_in, fasta_out):
    records = read_fasta(fasta_in)

    with open(csv_in) as f:
        reader = csv.DictReader(f)
        intervals = list(reader)

    for row in intervals:
        chrom = row["reference"]
        start = int(row["start"])
        end = int(row["end"])

        if chrom not in records:
            raise ValueError(f"Chromosome {chrom} not found in FASTA.")

        seq = records[chrom]

        # CSV coordinates are 1-based inclusive.
        left = seq[:start - 1]
        middle = seq[start - 1:end]
        right = seq[end:]

        records[chrom] = left + revcomp(middle) + right

    write_fasta(records, fasta_out)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(
            "Usage:\n"
            "  python SCRIPT_reverse_ref_inversions.py "
            "REF.fasta six_major_inversions_ref.csv REF_reoriented.fasta"
        )

    reverse_inversions(
        fasta_in=Path(sys.argv[1]),
        csv_in=Path(sys.argv[2]),
        fasta_out=Path(sys.argv[3]),
    )
