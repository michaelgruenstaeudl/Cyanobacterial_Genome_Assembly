#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

COLOR_MAP = {
    "+": "blue",
    "-": "red",
}

LINE_WIDTH = 2.4   # triple original
POINT_SIZE = 6.0   # triple original
INVERSION_LABEL_LINE_WIDTH = 0.5
INVERSION_LABEL_LINE_ALPHA = 0.35
INVERSION_LABEL_TEXT_ALPHA = 0.9
INVERSION_LABEL_TEXT_SIZE = 9
INVERSION_LABEL_TEXT_DX = 2
INVERSION_LABEL_TEXT_DY = 3
INVERSION_LABEL_TEXT_DY_TIGHT = 10
INVERSION_LABEL_TIGHT_GAP_BP = 70000
INVERSION_NAME_TEXT_SIZE = 10
INVERSION_NAME_TEXT_DX_LOWER = 72
INVERSION_NAME_TEXT_DX_UPPER = 18
AXIS_LABEL_TEXT_SIZE = 13
TICK_LABEL_TEXT_SIZE = 11
BP_TO_MBP = 1e6
INVERSION_CSV = "six_major_inversions_ref.csv"

def usage():
    sys.stderr.write(
        "Usage: plot_coords_dotplot.py <coords.tsv> <output_basename>\n"
    )
    sys.exit(1)

def resolve_coords_path(path: str) -> Path:
    p = Path(path)
    if p.exists():
        return p

    here = Path.cwd()
    hints = sorted(x.name for x in here.glob("*.coords.tsv"))
    msg = [f"Input file not found: {path}"]
    if hints:
        msg.append("Available *.coords.tsv files in current directory:")
        msg.extend(f"  - {h}" for h in hints)
    raise FileNotFoundError("\n".join(msg))

def read_coords_tsv(path: str):
    resolved = resolve_coords_path(path)
    df = pd.read_csv(resolved, sep="\t", header=None, comment="#", dtype=str)

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

def read_inversion_intervals(path: Path):
    df = pd.read_csv(path)
    required = {"start", "end"}
    if not required.issubset(df.columns):
        raise ValueError(f"Inversion CSV must contain columns: {sorted(required)}")

    intervals = []
    for _, row in df.iterrows():
        start = float(row["start"]) / BP_TO_MBP
        end = float(row["end"]) / BP_TO_MBP
        intervals.append((start, end))

    return intervals

def choose_red_segment(df, start_mbp, end_mbp):
    red_df = df[df["strand"] == "-"]
    overlaps = red_df.apply(
        lambda r: max(0.0, min(float(r["re"]), end_mbp) - max(float(r["rs"]), start_mbp)),
        axis=1,
    )
    if overlaps.empty or overlaps.max() <= 0.0:
        return None
    return red_df.loc[overlaps.idxmax()]

def format_bp_label(value_mbp):
    return f"{int(round(value_mbp * BP_TO_MBP)):,} bp"

def label_vertical_offsets(start_mbp, end_mbp):
    gap_bp = abs(end_mbp - start_mbp) * BP_TO_MBP
    if gap_bp <= INVERSION_LABEL_TIGHT_GAP_BP:
        return INVERSION_LABEL_TEXT_DY_TIGHT, -INVERSION_LABEL_TEXT_DY_TIGHT
    return INVERSION_LABEL_TEXT_DY, -INVERSION_LABEL_TEXT_DY

def dotplot_segments(df, ref_name, qry_name, out_base, inversion_intervals):
    fig, ax = plt.subplots(figsize=(10, 10))
    max_x = float(df["re"].max())
    max_y = float(df["qe"].max())
    mid_x = max_x / 2.0
    lower_group_x = mid_x + (max_x * 0.005)
    upper_group_x = max_x * 0.41

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

    for inversion_index, (start_mbp, end_mbp) in enumerate(inversion_intervals):
        segment = choose_red_segment(df, start_mbp, end_mbp)
        if segment is None:
            continue

        is_lower_group = inversion_index < 3
        label_x_pos = lower_group_x if is_lower_group else upper_group_x
        start_text_dy, end_text_dy = label_vertical_offsets(start_mbp, end_mbp)
        inversion_name = f"Inversion {inversion_index + 1}"
        inversion_name_y = (start_mbp + end_mbp) / 2.0
        inversion_name_dx = INVERSION_NAME_TEXT_DX_LOWER if is_lower_group else -INVERSION_NAME_TEXT_DX_UPPER
        inversion_name_ha = "left" if is_lower_group else "right"

        ax.plot(
            [start_mbp, mid_x],
            [start_mbp, start_mbp],
            linewidth=INVERSION_LABEL_LINE_WIDTH,
            color=COLOR_MAP["-"],
            alpha=INVERSION_LABEL_LINE_ALPHA,
        )
        ax.annotate(
            format_bp_label(start_mbp),
            xy=(label_x_pos, start_mbp),
            xytext=(INVERSION_LABEL_TEXT_DX, start_text_dy),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=INVERSION_LABEL_TEXT_SIZE,
            color=COLOR_MAP["-"],
            alpha=INVERSION_LABEL_TEXT_ALPHA,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.2),
        )
        ax.plot(
            [end_mbp, mid_x],
            [end_mbp, end_mbp],
            linewidth=INVERSION_LABEL_LINE_WIDTH,
            color=COLOR_MAP["-"],
            alpha=INVERSION_LABEL_LINE_ALPHA,
        )
        ax.annotate(
            format_bp_label(end_mbp),
            xy=(label_x_pos, end_mbp),
            xytext=(INVERSION_LABEL_TEXT_DX, end_text_dy),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=INVERSION_LABEL_TEXT_SIZE,
            color=COLOR_MAP["-"],
            alpha=INVERSION_LABEL_TEXT_ALPHA,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.2),
        )
        ax.annotate(
            inversion_name,
            xy=(label_x_pos, inversion_name_y),
            xytext=(inversion_name_dx, 0),
            textcoords="offset points",
            ha=inversion_name_ha,
            va="center",
            fontsize=INVERSION_NAME_TEXT_SIZE,
            color=COLOR_MAP["-"],
            fontweight="bold",
            alpha=INVERSION_LABEL_TEXT_ALPHA,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.2),
        )

    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    ax.set_xlabel("Limnothrix sp. BL-A-16 (Mbp)", fontsize=AXIS_LABEL_TEXT_SIZE)
    ax.set_ylabel("Bacterial chromosome assembled here (Mbp)", fontsize=AXIS_LABEL_TEXT_SIZE)
    ax.set_xlim(0.0, max_x)
    ax.set_ylim(0.0, max_y)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, which="major", linewidth=0.6, color="0.55", alpha=0.8)
    ax.grid(True, which="minor", linewidth=0.3, color="0.75", alpha=0.5)
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_TEXT_SIZE)

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
    inversion_csv = Path(__file__).with_name(INVERSION_CSV)
    inversion_intervals = read_inversion_intervals(inversion_csv)
    dotplot_segments(df, ref_name, qry_name, out_base, inversion_intervals)

if __name__ == "__main__":
    main()
