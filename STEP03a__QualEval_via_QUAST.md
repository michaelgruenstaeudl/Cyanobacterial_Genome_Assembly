#### Installation
```bash
## Requires Python3 !

# Clone QUAST if not already present
if [ ! -d "quast" ]; then
    git clone https://github.com/ablab/quast.git
fi

# Navigate to QUAST directory
cd quast
```

#### Run QUAST on both assemblies
```bash
# ----------------------------------------
# Input and output files
# ----------------------------------------
GENOME="bactGenome"
#GENOME="plasmid"

ILLUM_R1=~/HorseThiefReservoir_Fall2024/02b_backmapping/Limnothrix_Unicycler_R1.fastq.gz
ILLUM_R2=~/HorseThiefReservoir_Fall2024/02b_backmapping/Limnothrix_Unicycler_R2.fastq.gz
ONT_READS=~/HorseThiefReservoir_Fall2024/02b_backmapping/Limnothrix_Unicycler_ONT.fastq.gz
ASM=~/HorseThiefReservoir_Fall2024/02b_backmapping/FinalAssembly_${GENOME}_corrected.fasta
ILLUM_BACKM=~/HorseThiefReservoir_Fall2024/02b_backmapping/${GENOME}.illumina.sorted.bam
ONT_BACKM=~/HorseThiefReservoir_Fall2024/02b_backmapping/${GENOME}.ont.sorted.bam

python3 quast/quast.py $ASM \
  -o ${GENOME}_Illumina_QUAST_report \
  -t 8 \
  --k-mer-stats \
  --report-all-metrics \
  --pe1 $ILLUM_R1 \
  --pe2 $ILLUM_R2 \
  --bam $ILLUM_BACKM

python3 quast/quast.py $ASM \
  -o ${GENOME}_Nanopore_QUAST_report \
  -t 8 \
  --k-mer-stats \
  --report-all-metrics \
  --nanopore $ONT_READS \
  --bam $ONT_BACKM
```