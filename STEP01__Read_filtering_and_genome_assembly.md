### Cyanobacterial Genome Assembly (Hybrid Assembly)

#### Installation of the necessary tools
```
# Installation of filtlong
cd ~
mkdir -p git
mkdir -p bin
cd git
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
cp bin/filtlong ~/bin/filtlong
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Installation of bactopia
mamba create -n bactopia -c conda-forge -c bioconda bactopia
```

#### Filtering the sequence reads
```
# Combining Illumina reads
LOC_ILLUMINA_DAT1_R1=~/data/HorseThief/HorseThief_Illumina_Jan2025_run1_R1.fastq.gz  # Change!
LOC_ILLUMINA_DAT1_R2=~/data/HorseThief/HorseThief_Illumina_Jan2025_run1_R2.fastq.gz  # Change!
LOC_ILLUMINA_DAT2_R1=~/data/HorseThief/HorseThief_Illumina_Jan2025_run2_R1.fastq.gz  # Change!
LOC_ILLUMINA_DAT2_R2=~/data/HorseThief/HorseThief_Illumina_Jan2025_run2_R2.fastq.gz  # Change!

cat $LOC_ILLUMINA_DAT1_R1 $LOC_ILLUMINA_DAT2_R1 > Illumina_all_R1.fastq.gz
cat $LOC_ILLUMINA_DAT1_R2 $LOC_ILLUMINA_DAT2_R2 > Illumina_all_R2.fastq.gz

# Filtering the Illumina reads using \texttt{Trimmomatic
module load Trimmomatic

# Note: This will take 20-30 min!
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred33 \
  Illumina_all_R1.fastq.gz Illumina_all_R2.fastq.gz \
  Illumina_filt_R1_paired.fastq.gz Illumina_filt_R1_unpaired.fastq.gz \
  Illumina_filt_R2_paired.fastq.gz Illumina_filt_R2_unpaired.fastq.gz \
  ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  
rm Illumina_all_R1.fastq.gz Illumina_all_R2.fastq.gz
```

```
# Combining Nanopore reads
LOC_NANOPORE_DAT1=~/data/HorseThief/Nanopore_2025_06_12__fastq_pass/  # Change!
LOC_NANOPORE_DAT2=~/data/HorseThief/Nanopore_2025_06_24__fastq_pass/  # Change!

zcat $LOC_NANOPORE_DAT1/*.fastq.gz >  Nanopore_all.fastq
zcat $LOC_NANOPORE_DAT2/*.fastq.gz >> Nanopore_all.fastq

# Filtering the Nanopore reads using filtlong
filtlong --min_length 1 --keep_percent 75 Nanopore_all.fastq | 
	gzip > Nanopore_filtered.q75.fastq.gz

rm Nanopore_all.fastq
```

#### Hybrid DeNovo Assembly using Bactopia

##### Assembly primarily through Illumina reads
```
conda activate bactopia

# Under this option, Bactopia creates a hybrid assembly using Unicycler to assemble the short reads first, then bridging gaps with long reads.
# Good for low-coverage noisy long-reads.

bactopia \
   --sample HorseThief_Unicycler \
   --r1 Illumina_filt_R1_paired.fastq.gz \
   --r2 Illumina_filt_R2_paired.fastq.gz \
   --ont Nanopore_filtered.q75.fastq.gz \
   --hybrid
   --phix /homes/mgruenstaeudl/references/phix174.fasta   # Removal of any PhiX spike-ins
```

##### Assembly primarily through Oxford Nanopore reads
```
conda activate bactopia

# Under this option, Bactopia creates a hybrid assembly using Dragonflye to assemble the long-reads first, then polishing the assembly with the short-reads.
# Good for high-coverage high-quality long-reads.

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
```
