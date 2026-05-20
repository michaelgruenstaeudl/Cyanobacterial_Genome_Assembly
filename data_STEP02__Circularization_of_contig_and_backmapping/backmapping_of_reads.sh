module load BWA
module load SAMtools

INDIR="~/HorseThiefReservoir_Fall2024/02b_backmapping"

#GENOME=bactGenome
GENOME=plasmid
# ----------------------------------------
# Input files
# ----------------------------------------
ASSEMBLY="${INDIR}/FinalAssembly_${GENOME}_corrected.fasta"
ILLUMINA_R1="${INDIR}/Limnothrix_Unicycler_R1.fastq.gz"
ILLUMINA_R2="${INDIR}/Limnothrix_Unicycler_R2.fastq.gz"
ONT_READS="${INDIR}/Limnothrix_Unicycler_ONT.fastq.gz"

# ----------------------------------------
# Index assembly
# ----------------------------------------
bwa index ${ASSEMBLY}

# ----------------------------------------
# Illumina backmapping
# ----------------------------------------
bwa mem -t 8 \
    ${ASSEMBLY} \
    ${ILLUMINA_R1} \
    ${ILLUMINA_R2} | \
samtools sort -@ 8 -m 4G -o ${GENOME}.illumina.sorted.bam
samtools index ${GENOME}.illumina.sorted.bam

# ----------------------------------------
# ONT backmapping
# ----------------------------------------
bwa mem -x ont2d -t 8 \
    ${ASSEMBLY} \
    ${ONT_READS} | \
samtools sort -@ 8 -m 4G -o ${GENOME}.ont.sorted.bam
samtools index ${GENOME}.ont.sorted.bam
