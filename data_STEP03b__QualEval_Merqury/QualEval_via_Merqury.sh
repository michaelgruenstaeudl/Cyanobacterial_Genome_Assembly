#!/usr/bin/env bash
set -euo pipefail

eval "$(/homes/mgruenstaeudl/miniconda3/bin/conda shell.bash hook)"
conda activate merqury

# ----------------------------------------
# Input and output files
# ----------------------------------------
#GENOME="bactGenome"
GENOME="plasmid"

R1=~/HorseThiefReservoir_Fall2024/02b_backmapping/Limnothrix_Unicycler_R1.fastq.gz
R2=~/HorseThiefReservoir_Fall2024/02b_backmapping/Limnothrix_Unicycler_R2.fastq.gz
ASM=~/HorseThiefReservoir_Fall2024/02b_backmapping/FinalAssembly_${GENOME}_corrected.fasta

OUT=${GENOME}_Illumina_MERQURY_report
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"
MEMORY=8

# ----------------------------------------
# Run analysis
# ----------------------------------------
mkdir -p "$OUT"
cd "$OUT"

mkdir -p logs
PREFIX="$(basename "$OUT")_$(date +%Y-%m-%d)"

# Required by several Merqury installations
export MERQURY="${MERQURY:-$CONDA_PREFIX/share/merqury}"

# STEP 1. Build Illumina read k-mer database
meryl count \
  k="$K" \
  threads="$THREADS" \
  memory="$MEMORY" \
  "$R1" "$R2" \
  output "${PREFIX}.reads.meryl" \
  > logs/01_meryl_reads.log 2>&1

# STEP 2. Run full Merqury workflow
# Produces QV, error rate, completeness, spectra-cn plots, spectra-asm plots,
# histograms, and assembly/read k-mer comparison outputs.
"$MERQURY/merqury.sh" \
  "${PREFIX}.reads.meryl" \
  "$ASM" \
  "$PREFIX" \
  > logs/02_merqury_full.log 2>&1

# STEP 3. Also run qv.sh explicitly for a compact QV table
"$MERQURY/eval/qv.sh" \
  "${PREFIX}.reads.meryl" \
  "$ASM" \
  "${PREFIX}.qv_only" \
  > logs/03_qv_only.log 2>&1

# STEP 4. Collect the most useful tabular outputs
{
  echo "Merqury report: $PREFIX"
  echo "Genome: $GENOME"
  echo "Assembly: $ASM"
  echo "Read k-mer size: $K"
  echo

  echo "=== Consensus QV / error rate ==="
  cat "${PREFIX}.qv" 2>/dev/null || true
  echo

  echo "=== Explicit qv.sh output ==="
  cat "${PREFIX}.qv_only.qv" 2>/dev/null || true
  echo

  echo "=== K-mer completeness ==="
  cat "${PREFIX}.completeness.stats" 2>/dev/null || true
  echo

  echo "=== Main output files ==="
  ls -lh
} > "${PREFIX}.summary.txt"
