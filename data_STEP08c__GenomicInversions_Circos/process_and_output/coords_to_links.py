#!/usr/bin/env python3
from pathlib import Path

def main(coords_tsv: str, out_links: str, min_len: int = 5000, min_id: float = 95.0):
    inp = Path(coords_tsv)
    out = Path(out_links)

    with inp.open("r", encoding="utf-8") as f, out.open("w", encoding="utf-8") as g:
        for line in f:
            line = line.strip()
            if not line or line.startswith("[") or line.startswith("="):
                continue

            # show-coords -T output is tab-separated. Typical columns include:
            # 0:S1 1:E1 2:S2 3:E2 4:LEN1 5:LEN2 6:%IDY  ...  last two: REF_QRY sequence IDs
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            try:
                s1 = int(parts[0]); e1 = int(parts[1])
                s2 = int(parts[2]); e2 = int(parts[3])
                l1 = int(parts[4])
                pid = float(parts[6])
                ref = parts[-2].strip()
                qry = parts[-1].strip()
            except Exception:
                continue

            if l1 < min_len or pid < min_id:
                continue

            # Circos wants: <chr1> <start1> <end1> <chr2> <start2> <end2>
            # Keep start <= end for plotting
            a1, b1 = (s1, e1) if s1 <= e1 else (e1, s1)
            a2, b2 = (s2, e2) if s2 <= e2 else (e2, s2)

            # Detect inversion in second genome
            is_inversion = s2 > e2
            inv_flag = "1" if is_inversion else "0"

            # Determine inversion
            inv = 1 if (s2 > e2) else 0
            g.write(f"{ref}\t{a1}\t{b1}\t{qry}\t{a2}\t{b2}\tinv={inv}\n")



    print(f"Wrote links: {out}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        raise SystemExit("Usage: coords_to_links.py <A_vs_B.coords.tsv> <links.txt> [min_len] [min_id]\n"
                         "Example: coords_to_links.py A_vs_B.coords.tsv links.txt 1000 90")
    min_len = int(sys.argv[3]) if len(sys.argv) >= 4 else 5000
    min_id  = float(sys.argv[4]) if len(sys.argv) >= 5 else 95.0
    main(sys.argv[1], sys.argv[2], min_len=min_len, min_id=min_id)
