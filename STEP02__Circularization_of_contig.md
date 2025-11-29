### Converting assemblies to complete bacterial genomes

#### Using circlator to evaluate overlap of linear contig
```
# Create environment (Python 3.7 + Circlator 1.5.5 + older OpenSSL for compatibility)
conda create -n circlator_env -c bioconda -c conda-forge python=3.7 circlator=1.5.5 openssl=1.0

# Activate environment
conda activate circlator_env

# Run Circlator
circlator fixstart Limnothrix_polished_genome.fa fixstart_out

# Check the output
head fixstart_out.fasta

# Confirm overlap of contig ends (500 bp)
tail -c 1000 fixstart_out.fasta | head -c 500 > end.fa && head -c 500 fixstart_out.fasta > start.fa && diff start.fa end.fa && echo "Checked 500 bp overlap at contig ends"

```
