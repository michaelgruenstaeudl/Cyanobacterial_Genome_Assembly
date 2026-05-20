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
GENOME="bactGenome"
#GENOME="plasmid"

quast.py FinalAssembly_${GENOME}_corrected.fasta \
  -o ${GENOME}_Illumina_QUAST_report \
  -t 16 \
  --k-mer-stats \
  --report-all-metrics \
  --pe1 Limnothrix_Unicycler_R1.fastq.gz \
  --pe2 Limnothrix_Unicycler_R2.fastq.gz \
  --bam ${GENOME}.illumina.sorted.bam

quast.py FinalAssembly_${GENOME}_corrected.fasta \
  -o ${GENOME}_Nanopore_QUAST_report \
  -t 16 \
  --k-mer-stats \
  --report-all-metrics \
  --nanopore Limnothrix_Unicycler_ONT.fastq.gz \
  --bam ${GENOME}.ont.sorted.bam
```