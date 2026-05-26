#!/bin/bash
#SBATCH --job-name=merqury_compare_two
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10

set -euo pipefail

eval "$(/homes/mgruenstaeudl/miniconda3/bin/conda shell.bash hook)"
conda activate merqury

#--- INPUT -------------------------------------------------------------
R1=/homes/mgruenstaeudl/data/Limnothrix/02_processed_reads/Illumina_filt_R1_paired.fastq.gz
R2=/homes/mgruenstaeudl/data/Limnothrix/02_processed_reads/Illumina_filt_R2_paired.fastq.gz

ASM1=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Mercury/input/BacterialGenome_Bactopia.fasta
ASM2=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Mercury/input/PlasmidGenome_Bactopia.fasta

OUT=merqury_compare
PREFIX="$(basename "$OUT")_$(date +%Y-%m-%d)"
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"

#--- RUN ---------------------------------------------------------------
mkdir -p "$OUT"
cd "$OUT"

# Merqury writes relative logs/, so keep it inside $OUT
mkdir -p logs
mkdir -p figs

# Build the meryl database from paired-end reads
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1" "$R2" \
  output "${PREFIX}.meryl"

# Run Merqury comparison (two assemblies)
merqury.sh "${PREFIX}.meryl" "$ASM1" "$ASM2" "$PREFIX"

#--- VISUALIZE ---------------------------------------------------------

# Plot only BacterialGenome
Rscript "$MERQURY/plot/plot_spectra_cn.R" \
  -f "${PREFIX}.BacterialGenome_Bactopia.spectra-cn.hist" \
  -o "figs/${PREFIX}.BacterialGenome_Bactopia.spectra-cn" \
  -x 10 -y 6 -m 750 -t line -p

# Plot only PlasmidGenome
Rscript "$MERQURY/plot/plot_spectra_cn.R" \
  -f "${PREFIX}.PlasmidGenome_Bactopia.spectra-cn.hist" \
  -o "figs/${PREFIX}.PlasmidGenome_Bactopia.spectra-cn" \
  -x 10 -y 6 -m 6000 -t line -p

# Plot both genomes simultaneously (custom R script)
Rscript "$MERQURY/plot/plot_spectra_cn_cuberootlogscale.R" \
  -f "${PREFIX}.spectra-cn.hist" \
  -o "figs/${PREFIX}.BothGenomes.spectra-cn" \
  -x 10 -y 6 -m 2500 -t fill -p
