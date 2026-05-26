### Collinearity between two genomes as mummer dotplots

#### Install MUMmer4, Run nucmer and produce plot-ready tables

##### Installing MUMmer4
```bash

# Installation of MUMmer4
# mamba create -y -n mummer4_env -c conda-forge -c bioconda mummer4 python

# Activating MUmer4 as conda package
conda activate mummer4_env

# Checking that mummer available
command -v nucmer  >/dev/null 2>&1 || { echo "ERROR: nucmer not found"; exit 1; }
command -v delta-filter >/dev/null 2>&1 || { echo "ERROR: delta-filter not found"; exit 1; }
command -v show-coords >/dev/null 2>&1 || { echo "ERROR: show-coords not found"; exit 1; }
```

##### Running nucmer and producing plot-ready tables
__Important__: Work in a directory with no spaces and use simple filenames (no special characters other than underscores)

```bash

REF="Limnothrix_sp_BL_A_16_CP166615.fasta"
QRY="FinalAssembly_Bactopia.fasta"
PFX="Limnothrix_sp_BLA16_vs_BactopiaAssembly"

## STEP 1. RUNNING NUCMER ALIGNMENT
# --maxmatch is common for whole-genome dotplots
# -p sets output prefix: creates ${PFX}.delta, etc.
nucmer --maxmatch -p "${PFX}" "${REF}" "${QRY}"

## STEP 2. FILTER THE DELTA TO EMPHASIZE COLLINEARITY AND INVERSIONS
# For bacterial genomes, start with min aligned block length 1000 or 5000
# -r -q keeps best reciprocal alignments (helps reduce clutter), as in the Hippocamplus post. :contentReference[oaicite:2]{index=2}
# Adjust -l if you want more/less detail.
delta-filter -r -q -l 1000 "${PFX}.delta" > "${PFX}.filt.delta"

echo "Tip: if the plot is too dense, raise the filter length, for example:"
echo "  delta-filter -r -q -l 5000 ${PFX}.delta > ${PFX}.filt.delta"

## STEP 3. EXPORTING COORDINATES TABLE FOR PLOTTING
# show-coords parses .delta and outputs coordinate summaries. :contentReference[oaicite:3]{index=3}
# Flags:
#  -H  no header
#  -T  tab-delimited
#  -r  sort by reference
#  -c  include percent coverage columns (handy, optional)
show-coords -H -T -r -c "${PFX}.filt.delta" > "${PFX}.coords.tsv"

```
