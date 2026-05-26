### K-mer spectrum analysis for bacterial genome assembly

#### Preparation
```bash
module load Jellyfish
```

#### K-mer counting with Jellyfish on Illumina data

```bash

# Parameters
K=21
THREADS=8
HASH="1G"

# Output directory
ILM_OUT=~/data/Limnothrix/05_kmer_spectrum/KmerSpectrum_Illumina_k${K}
mkdir -p "$ILM_OUT"
cd "$ILM_OUT"

# Input reads (edit as needed; FASTQ or FASTQ.GZ both work with zcat -f)
ILLUMINA_READS=(
  ~/data/Limnothrix/02_processed_reads/Illumina_filt_R1_paired.fastq.gz
  ~/data/Limnothrix/02_processed_reads/Illumina_filt_R2_paired.fastq.gz
)

# Count canonical k-mers (-C) and write Jellyfish database
# zcat -f works for both .fastq and .fastq.gz
jellyfish count -C -m "$K" -s "$HASH" -t "$THREADS" \
  /dev/fd/0 -o "illumina_reads_k${K}.jf" < <(zcat -f "${ILLUMINA_READS[@]}")

# Make k-mer multiplicity histogram
jellyfish histo -t "$THREADS" "illumina_reads_k${K}.jf" > "illumina_reads_k${K}.histo"
```

#### K-mer counting with Jellyfish on ONT data

```bash

# Parameters
K=17
THREADS=8
HASH="1G"

# Output directory
ONT_OUT=~/data/Limnothrix/05_kmer_spectrum/KmerSpectrum_ONT_k${K}
mkdir -p "$ONT_OUT"
cd "$ONT_OUT"

# Input reads (edit as needed; FASTQ or FASTQ.GZ)
ONT_READS=(
  ~/data/Limnothrix/02_processed_reads/Nanopore_filtered.q75.fastq.gz
)

# Count canonical k-mers (-C) and write Jellyfish database
jellyfish count -C -m "$K" -s "$HASH" -t "$THREADS" \
  /dev/fd/0 -o "ont_reads_k${K}.jf" < <(zcat -f "${ONT_READS[@]}")

# Make k-mer multiplicity histogram
jellyfish histo -t "$THREADS" "ont_reads_k${K}.jf" > "ont_reads_k${K}.histo"
```


#### Visualization of Jellyfish via Python's Seaborn

```python
import os
import tempfile
import shutil

import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use("Agg")  # robust non-interactive backend
import matplotlib.pyplot as plt


def load_histo(path: str) -> pd.DataFrame:
    """
    Read a Jellyfish .histo file (two columns: multiplicity, n_kmers) into a DataFrame.
    """
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["multiplicity", "n_kmers"])
    df["multiplicity"] = pd.to_numeric(df["multiplicity"], errors="coerce")
    df["n_kmers"] = pd.to_numeric(df["n_kmers"], errors="coerce")
    df = df.dropna().sort_values("multiplicity")
    return df


def copy_to_tmp(src_path: str) -> str:
    """
    Copy a file to /tmp and return the /tmp path (often faster on mounted filesystems).
    """
    tmp_path = os.path.join(tempfile.gettempdir(), os.path.basename(src_path))
    shutil.copy2(src_path, tmp_path)
    return tmp_path


def plot_two_spectra(
    df_illumina: pd.DataFrame,
    df_ont: pd.DataFrame,
    out_svg: str,
    title: str = "Read k-mer spectra",
    x_max: int | None = 200,
    log_y: bool = True,
) -> None:
    """
    Plot Illumina and ONT k-mer spectra on the same axes and save as SVG.

    x_max:
      - set to an integer (e.g., 200) to zoom for readability
      - set to None to plot the entire x-range (can be very wide)
    """
    df_i = df_illumina.copy()
    df_o = df_ont.copy()

    if x_max is not None:
        df_i = df_i[df_i["multiplicity"] <= x_max]
        df_o = df_o[df_o["multiplicity"] <= x_max]

    df_i["platform"] = "Illumina"
    df_o["platform"] = "Oxford Nanopore"
    df_all = pd.concat([df_i, df_o], ignore_index=True)

    sns.set_theme(style="whitegrid")

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    sns.lineplot(
        data=df_all,
        x="multiplicity",
        y="n_kmers",
        hue="platform",
        linewidth=2,
        ax=ax,
    )

    ax.set_title(title)
    ax.set_xlabel("k-mer multiplicity")
    ax.set_ylabel("Number of distinct k-mers")

    if log_y:
        ax.set_yscale("log")

    tmp_svg = os.path.join(tempfile.gettempdir(), os.path.basename(out_svg))
    fig.savefig(tmp_svg, format="svg")
    plt.close(fig)

    os.makedirs(os.path.dirname(out_svg), exist_ok=True)
    shutil.copy2(tmp_svg, out_svg)


# ---------------- Parameters ----------------
K_ILLUMINA = 21
K_ONT = 17

# Set to None to plot full multiplicity range; keep as 200/500 for a classic zoomed spectrum view
X_MAX = 200

# Use absolute paths if you want to avoid dependence on current working directory:
# ILM_BASE = os.path.expanduser(f"~/data/Limnothrix/05_kmer_spectrum/KmerSpectrum_Illumina_k{K_ILLUMINA}")
# ONT_BASE = os.path.expanduser(f"~/data/Limnothrix/05_kmer_spectrum/KmerSpectrum_ONT_k{K_ONT}")

ILM_BASE = os.path.expanduser(f"05_kmer_spectrum/KmerSpectrum_Illumina_k{K_ILLUMINA}")
ONT_BASE = os.path.expanduser(f"05_kmer_spectrum/KmerSpectrum_ONT_k{K_ONT}")

illumina_histo = os.path.join(ILM_BASE, f"illumina_reads_k{K_ILLUMINA}.histo")
ont_histo = os.path.join(ONT_BASE, f"ont_reads_k{K_ONT}.histo")

# Output (choose a combined output location; here: Illumina directory)
combined_svg = os.path.join(
    ILM_BASE, f"kmer_spectrum_Illumina_k{K_ILLUMINA}_ONT_k{K_ONT}_combined.svg"
)

# ---------------- Load data (full histograms) ----------------
if not os.path.isfile(illumina_histo):
    raise FileNotFoundError(f"Histogram not found: {illumina_histo}")
if not os.path.isfile(ont_histo):
    raise FileNotFoundError(f"Histogram not found: {ont_histo}")

illumina_histo_tmp = copy_to_tmp(illumina_histo)
ont_histo_tmp = copy_to_tmp(ont_histo)

df_illumina = load_histo(illumina_histo_tmp)
df_ont = load_histo(ont_histo_tmp)

# ---------------- Plot combined ----------------
plot_two_spectra(
    df_illumina=df_illumina,
    df_ont=df_ont,
    out_svg=combined_svg,
    title=f"k-mer spectra (Illumina k={K_ILLUMINA}, ONT k={K_ONT})",
    x_max=X_MAX,
    log_y=True,
)

print("Wrote:")
print(" ", combined_svg)

```