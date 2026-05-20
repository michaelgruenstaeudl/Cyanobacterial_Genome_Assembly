### Inference of QV for genome assembly via CheckM2

#### Installation of CheckM2 on Beocat
```bash
conda create -n checkm2 -c bioconda -c conda-forge checkm2 -y
conda activate checkm2
```

#### Run CheckM2 on both assemblies
```bash
conda activate checkm2
checkm2 database --download

# ----------------------------------------
# Input and output files
# ----------------------------------------

GENOME="bactGenome"
#GENOME="plasmid"

ASM=FinalAssembly_${GENOME}_corrected.fasta

# ----------------------------------------
# Run analyses
# ----------------------------------------
checkm2 predict \
    --threads 8 \
    --input $ASM \
    --output-directory checkm2_${GENOME}_output
```