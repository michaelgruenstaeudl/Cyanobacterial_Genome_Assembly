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
LOC_ILLUMINA_DAT1_R1=~/data/Limnothrix/Limnothrix_Illumina_Jan2025_run1_R1.fastq.gz  # Change!
LOC_ILLUMINA_DAT1_R2=~/data/Limnothrix/Limnothrix_Illumina_Jan2025_run1_R2.fastq.gz  # Change!
LOC_ILLUMINA_DAT2_R1=~/data/Limnothrix/Limnothrix_Illumina_Jan2025_run2_R1.fastq.gz  # Change!
LOC_ILLUMINA_DAT2_R2=~/data/Limnothrix/Limnothrix_Illumina_Jan2025_run2_R2.fastq.gz  # Change!

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
LOC_NANOPORE_DAT1=~/data/Limnothrix/Nanopore_2025_06_12__fastq_pass/  # Change!
LOC_NANOPORE_DAT2=~/data/Limnothrix/Nanopore_2025_06_24__fastq_pass/  # Change!

zcat $LOC_NANOPORE_DAT1/*.fastq.gz >  Nanopore_all.fastq
zcat $LOC_NANOPORE_DAT2/*.fastq.gz >> Nanopore_all.fastq

# Filtering the Nanopore reads using filtlong
filtlong --min_length 1 --keep_percent 75 Nanopore_all.fastq | 
	gzip > Nanopore_filtered.q75.fastq.gz

rm Nanopore_all.fastq
```

#### De Novo Assembly using Bactopia

##### Option that calls Unicycler to generate the assembly
```
conda activate bactopia

# Under this option, Bactopia creates a hybrid assembly using Unicycler to assemble the short reads first, then bridging gaps with long reads.
# Good for low-coverage noisy long-reads.

# Will take 30-60 min!
bactopia \
   --sample Limnothrix_Unicycler \
   --r1 Illumina_filt_R1_paired.fastq.gz \
   --r2 Illumina_filt_R2_paired.fastq.gz \
   --ont Nanopore_filtered.q75.fastq.gz \
   --hybrid
```

##### Option that calls Dragonflye to generate the assembly
```
conda activate bactopia

# Under this option, Bactopia creates a hybrid assembly using Dragonflye to assemble the long-reads first, then polishing the assembly with the short-reads.
# Good for high-coverage high-quality long-reads.

# Will take 30-60 min!
bactopia \
   --sample Limnothrix_Dragonflye \
   --r1 Illumina_filt_R1_paired.fastq.gz \
   --r2 Illumina_filt_R2_paired.fastq.gz \
   --ont Nanopore_filtered.q75.fastq.gz \
   --short_polish
```
