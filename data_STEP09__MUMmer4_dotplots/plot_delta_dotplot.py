#!/usr/bin/env python3
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt

COLOR_MAP = {
    "+": "blue",
    "-": "red",
}

LINE_WIDTH = 2.4   # triple original
POINT_SIZE = 6.0   # triple original
BP_TO_MBP = 1e6

def usage():
    sys.stderr.write(
        "Usage: plot_delta_dotplot.py <filtered.delta> <output_basename>\n"
    )
    sys.exit(1)

def read_delta(delta_path: str):
    rows = []
    ref_name = None
    qry_name = None

    with open(delta_path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    if not lines:
        raise ValueError("Delta file is empty")

    lines = lines[1:]

    header_re = re.compile(r"^>(\S+)\s+(\S+)\s+(\d+)\s+(\d+)")

    for line in lines:
        line = line.strip()
        if not line or line == "0":
            continue

        m = header_re.match(line)
        if m:
            ref_name = m.group(1)
            qry_name = m.group(2)
            continue

        parts = line.split()
        if len(parts) == 7 and all(p.lstrip("-").isdigit() for p in parts):
            rs, re_, qs, qe, err = map(int, parts[:5])
            rows.append(
                {
                    "rs": rs / BP_TO_MBP,
                    "re": re_ / BP_TO_MBP,
                    "qs": qs / BP_TO_MBP,
                    "qe": qe / BP_TO_MBP,
                    "strand": "+" if (qe - qs) > 0 else "-",
                }
            )

    if not rows:
        raise ValueError("No alignment records parsed")

    return pd.DataFrame(rows), ref_name, qry_name

def dotplot(df, ref_name, qry_name, out_base):
    fig, ax = plt.subplots(figsize=(10, 10))

    for strand, g in df.groupby("strand", sort=True):
        color = COLOR_MAP[strand]
        for _, r in g.iterrows():
            ax.plot(
                [r["rs"], r["re"]],
                [r["qs"], r["qe"]],
                linewidth=LINE_WIDTH,
                color=color,
            )
        ax.scatter(g["rs"], g["qs"], s=POINT_SIZE, alpha=0.5, color=color)
        ax.scatter(g["re"], g["qe"], s=POINT_SIZE, alpha=0.5, color=color)

    ax.set_xlabel(f"{ref_name} (Mbp)")
    ax.set_ylabel(f"{qry_name} (Mbp)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.3, alpha=0.4)

    fig.tight_layout()
    fig.savefig(f"{out_base}.png", dpi=300)
    fig.savefig(f"{out_base}.svg")
    plt.close(fig)

def main():
    if len(sys.argv) != 3:
        usage()

    delta_file = sys.argv[1]
    out_base = sys.argv[2]

    df, ref_name, qry_name = read_delta(delta_file)
    dotplot(df, ref_name, qry_name, out_base)

if __name__ == "__main__":
    main()
