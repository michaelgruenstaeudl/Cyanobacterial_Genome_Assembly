### Collinearity between two genomes as mummer dotplots

#### Installing and activating MUMmer4
```bash
# Install mamba into the base conda environment
conda install -y -n base -c conda-forge mamba

# Create a dedicated MUMmer4 environment
mamba create -y -n mummer4_env -c conda-forge -c bioconda mummer4 bedtools python pandas matplotlib

# Activate the MUMmer4 environment
conda activate mummer4_env

# Check that the required MUMmer4 programs are available
command -v nucmer       >/dev/null 2>&1 || { echo "ERROR: nucmer not found"; exit 1; }
command -v delta-filter >/dev/null 2>&1 || { echo "ERROR: delta-filter not found"; exit 1; }
command -v show-coords  >/dev/null 2>&1 || { echo "ERROR: show-coords not found"; exit 1; }
command -v bedtools     >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }
```

---

#### Running nucmer and producing plot-ready tables
__Important__: Work in a directory with no spaces and use simple filenames (no special characters other than underscores)

```bash
cd data_STEP08a__GenomicInversions_MUMmer4

REF="./input/Limnothrix_sp_BL_A_16_CP166615.fasta"
QRY="./input/FinalAssembly_bactGenome_corrected.fasta"
PFX="Limnothrix_sp_BLA16_vs_BacterialChr__MUMmer4"

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

---

#### Numeric inference of inversion locations

##### Identifying inversion locations
```bash
show-coords -rclTH "${PFX}.filt.delta" \
| awk -F'\t' '$3 > $4 && $5 >= 1000 && $7 >= 95 {
    s=$1; e=$2; if (s>e) {t=s; s=e; e=t}
    # BED-like: chrom start end (tab-delimited)
    print $12"\t"s"\t"e
}' \
| sort -k1,1 -k2,2n \
| bedtools merge -i - -d 10000 \
| awk 'BEGIN{OFS="\t"} {len=$3-$2+1; print $1,$2,$3,len}' \
| sort -k4,4nr \
| head -n 6 \
| sort -k1,1 -k2,2n \
| awk 'BEGIN{OFS=","; print "reference,start,end,length"} {print $1,$2,$3,$4}' \
> six_major_inversions_ref.csv
```

---

#### Visualization of dotplots via matplotlib

##### Plot from `show-coords` TSV
```bash
python SCRIPT_MUMmer4_plot_coords_dotplot.py "${PFX}.coords.tsv" "${PFX}.dotplot"
```