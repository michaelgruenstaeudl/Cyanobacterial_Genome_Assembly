#!/bin/bash
#SBATCH --job-name=merqury_plasm
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --mail-type=all
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10

bash QualEval_via_Merqury.sh
