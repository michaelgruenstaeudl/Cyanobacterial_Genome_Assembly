#!/usr/bin/env python3
import sys
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
        "Usage: plot_coords_dotplot.py <coords.tsv> <output_basename>\n"
    )
    sys.exit(1)

def read_coords_tsv(path: str):
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)

    if df.shape[1] < 6:
        raise ValueError("Unexpected show-coords format")

    df = df.rename(columns={
        0: "rs", 1: "re",
        2: "qs", 3: "qe",
        df.shape[1] - 2: "ref_name",
        df.shape[1] - 1: "qry_name",
    })

    for c in ["rs", "re", "qs", "qe"]:
        df[c] = pd.to_numeric(df[c], errors="raise") / BP_TO_MBP

    df["strand"] = df.apply(
        lambda r: "+" if (r["qe"] - r["qs"]) > 0 else "-", axis=1
    )

    ref_name = df["ref_name"].iloc[0]
    qry_name = df["qry_name"].iloc[0]

    return df, ref_name, qry_name

def dotplot_segments(df, ref_name, qry_name, out_base):
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

    coords_tsv = sys.argv[1]
    out_base = sys.argv[2]

    df, ref_name, qry_name = read_coords_tsv(coords_tsv)
    dotplot_segments(df, ref_name, qry_name, out_base)

if __name__ == "__main__":
    main()
