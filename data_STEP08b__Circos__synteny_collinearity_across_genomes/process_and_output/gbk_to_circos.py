#!/usr/bin/env python3
from pathlib import Path
from Bio import SeqIO

def safe_id(s: str) -> str:
    # Circos ids should be simple: letters, numbers, underscore
    out = []
    for ch in s:
        if ch.isalnum() or ch == "_":
            out.append(ch)
        elif ch in " .:-|/":
            out.append("_")
    return "".join(out).strip("_") or "seq"

def main(gbk_path: str, prefix: str):
    gbk = Path(gbk_path)
    records = list(SeqIO.parse(str(gbk), "genbank"))
    if not records:
        raise SystemExit(f"No records found in {gbk}")

    # Write FASTA (one entry per record/replicon)
    fasta_path = Path(f"{prefix}.fna")
    SeqIO.write(records, str(fasta_path), "fasta")

    # Circos karyotype + simple gene track
    karyo_path = Path(f"{prefix}.karyotype.txt")
    genes_path = Path(f"{prefix}.genes.txt")

    with karyo_path.open("w", encoding="utf-8") as kf, genes_path.open("w", encoding="utf-8") as gf:
        for i, rec in enumerate(records, start=1):
            rec_id_raw = rec.id or rec.name or f"{prefix}_{i}"
            rec_id = safe_id(f"{prefix}_{rec_id_raw}")
            length = len(rec.seq)

            # Karyotype line format:
            # chr - <id> <label> <start> <end> <color>
            label = f"{prefix}:{rec_id_raw}"
            color = "vdgrey" if prefix.lower().endswith("a") else "vlgrey"
            kf.write(f"chr - {rec_id} {label} 0 {length} {color}\n")

            # Genes (very simple): write CDS as intervals for a heatmap/histogram track
            # Format for a "text" or "highlight"-style file depends on Circos plot type,
            # but for many plots: <chr> <start> <end> <value>
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                loc = feat.location
                if loc is None:
                    continue
                start = int(loc.start)
                end = int(loc.end)
                if end <= start:
                    continue
                gf.write(f"{rec_id}\t{start}\t{end}\t1\n")

    print(f"Wrote: {fasta_path}")
    print(f"Wrote: {karyo_path}")
    print(f"Wrote: {genes_path}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        raise SystemExit("Usage: gbk_to_circos.py <genome.gbk> <prefix>\nExample: gbk_to_circos.py genomeA.gbk A")
    main(sys.argv[1], sys.argv[2])
