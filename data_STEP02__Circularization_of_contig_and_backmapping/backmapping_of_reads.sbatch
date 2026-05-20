#!/bin/bash
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --job-name=backmapping_Illumina
#SBATCH --mail-type=all
#SBATCH --mem=12G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00

bash backmapping_of_reads.sh
