#!/usr/bin/env bash
set -euo pipefail

eval "$(/homes/mgruenstaeudl/miniconda3/bin/conda shell.bash hook)"
conda activate merqury

module load BWA
module load SAMtools

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

mkdir -p logs mapped_reads
PREFIX="$(basename "$OUT")_$(date +%Y-%m-%d)"

# Required by several Merqury installations
export MERQURY="${MERQURY:-$CONDA_PREFIX/share/merqury}"

# ----------------------------------------
# STEP 1. Extract reads that map to the selected genome
# ----------------------------------------

# Index assembly for BWA if needed
if [[ ! -f "${ASM}.bwt" ]]; then
  bwa index "$ASM" > logs/00_bwa_index.log 2>&1
fi

# Map Illumina reads to the selected assembly
bwa mem -t "$THREADS" \
  "$ASM" \
  "$R1" \
  "$R2" \
  2> logs/01_bwa_mem.log | \
samtools view -@ "$THREADS" -b -F 4 -o mapped_reads/${GENOME}.mapped.bam -

# Sort and index mapped BAM
samtools sort -@ "$THREADS" \
  -o mapped_reads/${GENOME}.mapped.sorted.bam \
  mapped_reads/${GENOME}.mapped.bam \
  > logs/02_samtools_sort.log 2>&1

samtools index \
  mapped_reads/${GENOME}.mapped.sorted.bam \
  > logs/03_samtools_index.log 2>&1

# Convert mapped alignments back to paired FASTQ files
samtools fastq -@ "$THREADS" \
  -1 mapped_reads/${GENOME}.mapped.R1.fastq.gz \
  -2 mapped_reads/${GENOME}.mapped.R2.fastq.gz \
  -0 /dev/null \
  -s /dev/null \
  -n \
  mapped_reads/${GENOME}.mapped.sorted.bam \
  > logs/04_samtools_fastq.log 2>&1

MAPPED_R1=mapped_reads/${GENOME}.mapped.R1.fastq.gz
MAPPED_R2=mapped_reads/${GENOME}.mapped.R2.fastq.gz

# ----------------------------------------
# STEP 2. Build Illumina read k-mer database
# ----------------------------------------
meryl count \
  k="$K" \
  threads="$THREADS" \
  memory="$MEMORY" \
  "$MAPPED_R1" "$MAPPED_R2" \
  output "${PREFIX}.reads.meryl" \
  > logs/05_meryl_reads.log 2>&1

# ----------------------------------------
# STEP 3. Run full Merqury workflow
# ----------------------------------------
"$MERQURY/merqury.sh" \
  "${PREFIX}.reads.meryl" \
  "$ASM" \
  "$PREFIX" \
  > logs/06_merqury_full.log 2>&1

# ----------------------------------------
# STEP 4. Also run qv.sh explicitly for a compact QV table
# ----------------------------------------
"$MERQURY/eval/qv.sh" \
  "${PREFIX}.reads.meryl" \
  "$ASM" \
  "${PREFIX}.qv_only" \
  > logs/07_qv_only.log 2>&1

# ----------------------------------------
# STEP 5. Collect the most useful tabular outputs
# ----------------------------------------
{
  echo "Merqury report: $PREFIX"
  echo "Genome: $GENOME"
  echo "Assembly: $ASM"
  echo "Original R1: $R1"
  echo "Original R2: $R2"
  echo "Mapped R1: $MAPPED_R1"
  echo "Mapped R2: $MAPPED_R2"
  echo "Read k-mer size: $K"
  echo

  echo "=== Mapping summary ==="
  samtools flagstat mapped_reads/${GENOME}.mapped.sorted.bam
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
