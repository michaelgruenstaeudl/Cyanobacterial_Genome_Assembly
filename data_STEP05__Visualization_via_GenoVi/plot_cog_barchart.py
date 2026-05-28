#!/usr/bin/env python3
"""Create a GenoVi-compatible COG bar chart from the COG classification CSV."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Force all fonts to Arial, Helvetica, or sans-serif
import matplotlib as mpl
mpl.rcParams['font.family'] = ['sans-serif']
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Liberation Sans', 'Helvetica', 'Arial', 'sans-serif']
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['svg.fonttype'] = 'none'  # Prefer text, not paths


COG_NAMES = {
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "Unclassified": "Unclassified",
}

COG_COLOR_MAP = {
    "D": "#637BB7",
    "M": "#2659A8",
    "N": "#6497B0",
    "O": "#2084AE",
    "T": "#42A9B3",
    "U": "#126974",
    "V": "#61C3A6",
    "W": "#259775",
    "Z": "#15793C",
    "J": "#A660A7",
    "K": "#9D3F97",
    "L": "#8574B5",
    "X": "#4F489E",
    "C": "#AAD382",
    "E": "#7DB040",
    "F": "#B2B36D",
    "G": "#8F8A2F",
    "H": "#EBBC86",
    "I": "#AF7E35",
    "P": "#DB8856",
    "Q": "#C46426",
    "R": "#626262",
    "S": "#909090",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a COG bar chart SVG from a GenoVi COG classification CSV."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="Input COG classification CSV path.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output SVG path (default: <input_basename>_COG_barplot.svg).",
    )
    return parser.parse_args()


def build_dataframe(csv_path: Path) -> pd.DataFrame:
    raw = pd.read_csv(csv_path)
    replicon_row = raw[raw.iloc[:, 0] == "Replicon"]
    total_row = raw[raw.iloc[:, 0] == "Total"]
    if replicon_row.empty or total_row.empty:
        raise ValueError("CSV must contain rows starting with 'Replicon' and 'Total'.")

    df = pd.DataFrame(
        {
            "COG": replicon_row.iloc[0, 1:].dropna().astype(str).str.strip().values,
            "Count": total_row.iloc[0, 1:].dropna().astype(int).values,
        }
    )
    df["Label"] = df["COG"].apply(lambda x: f"{COG_NAMES.get(x, x)} [{x}]")
    df["Color"] = df["COG"].map(COG_COLOR_MAP).fillna("#B0B0B0")
    return df


def main() -> None:
    args = parse_args()
    csv_path = args.input.resolve()
    if not csv_path.exists():
        raise FileNotFoundError(f"Input file not found: {csv_path}")

    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("default")

    df = build_dataframe(csv_path)
    fig, ax = plt.subplots(figsize=(14, 6))
    bars = ax.bar(df["Label"], df["Count"], color=df["Color"], width=0.8)

    # Force legend font to sans-serif
    legend = ax.legend(["COG Category"], loc="upper right", frameon=False, prop={"family": "sans-serif", "size": 14})
    for text in legend.get_texts():
        text.set_fontfamily("sans-serif")

    # Explicitly set fontfamily for all text (use 'sans-serif' for portability)
    ax.set_title("COG Category", fontsize=14, fontfamily="sans-serif")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(
        df["Label"], rotation=65, ha="right", va="top", rotation_mode="anchor", fontfamily="sans-serif"
    )
    ax.tick_params(axis="x", labelsize=8, pad=10)
    ax.margins(x=0.02)
    plt.subplots_adjust(bottom=0.36)
    y_max = max(float(df["Count"].max()), 1.0)
    ax.set_ylim(0, y_max * 1.15)

    for bar in bars:
        h = bar.get_height()
        ax.annotate(
            f"{int(h)}",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            fontfamily="sans-serif"
        )

    plt.tight_layout()
    out_svg = args.output or csv_path.with_name(f"{csv_path.stem}_COG_barplot.svg")
    fig.savefig(out_svg, format="svg", dpi=300, bbox_inches="tight")
    print(f"Saved SVG to: {out_svg}")


if __name__ == "__main__":
    main()
