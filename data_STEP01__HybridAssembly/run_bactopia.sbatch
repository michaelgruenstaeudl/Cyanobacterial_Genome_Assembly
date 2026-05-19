#!/bin/bash
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --job-name=bactopia_HorseThief
#SBATCH --mail-type=all
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00

# Load software
source /homes/mgruenstaeudl/.bashrc
conda activate bactopia

# Run Bactopia
INPATH=/homes/mgruenstaeudl/hybrid_assembly/

bactopia \
  --sample Limnothrix_Unicycler \
  --r1 ${INPATH}/Illumina_filt_R1_paired.fastq.gz \
  --r2 ${INPATH}/Illumina_filt_R2_paired.fastq.gz \
  --ont ${INPATH}/Nanopore_filtered.q75.fastq.gz \
  --hybrid \
  --genome_size 4600000 \
  --coverage 0 \
  --use_bbmap \
  --phix /homes/mgruenstaeudl/references/phix174.fasta \
  --phix_k 31 \
  --maq 20 \
  --minlength 75 \
  --trimq 15 \
  --min_contig_len 500 \
  --min_contig_cov 2 \
  --unicycler_mode conservative \
  --nfconfig /homes/mgruenstaeudl/hybrid_assembly/Limnothrix_denovo_assembly.conf
