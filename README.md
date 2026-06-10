# Cyanobacterial Genome Assembly

Scripts and workflows for the **assembly, polishing, quality assessment, and visualization** of complete cyanobacterial genomes.

---

### Overview

This repository contains step-by-step protocols and scripts for generating high-quality cyanobacterial genome assemblies from sequencing data. The workflow covers hybrid assembly, contig circularization, quality evaluation, sequence annotation, and multiple visualization approaches commonly used in comparative and structural genomics.

---

### 01. Genome Assembly
- [Hybrid genome assembly](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP01__Read_filtering_and_genome_assembly.md)

---

### 02. Contig circularization and Backmapping
- [Contig circularization](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP02__Circularization_of_contig_and_backmapping.md)

---

### 03. Assembly Quality Assessment
- [Quality evaluation via QUAST](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP03a__QualEval_via_QUAST.md)

- [Quality evaluation via Merqury](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP03b__QualEval_via_Merqury.md)

- [Quality evaluation via CheckM2](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP03c__QualEval_via_CheckM2.md)

---

### 04. Adding and evaluating genome sequence annotations

#### Generate sequence annotations via PGAP
```bash
$ head Limnothrix_sp_HT2024_plasmid_GENOME.yaml
fasta:
    class: File
    location: Limnothrix_sp_HT2024_plasmid_GENOME.fasta
submol:
    class: File
    location: Limnothrix_sp_HT2024_plasmid_SUBMOL_PRE.yaml

$ head Limnothrix_sp_HT2024_plasmid_SUBMOL_PRE.yaml
topology: 'circular'
location: plasmid
organism:
    genus_species: 'Limnothrix sp.'
    strain: 'HorseThiefRes_Aug2024'
contact_info:
    last_name: 'Michael'
    first_name: 'Gruenstaeudl'
    email: 'm_gruenstaeudl@fhsu.edu'
    organization: 'Fort Hays State University'

$ head Limnothrix_sp_HT2024_plasmid_GENOME.fasta
>Limnothrix sp. HorseThiefRes_Aug2024 plasmid, complete sequence [plasmid-name=unnamed]
ATGCCTGTGCAGGACACAGACGAGGCTATACCTGTCAGTGTACCTGAATCCTCGAATTCC
TCAATCCCATTAGCTGAGGTCTCACCCCGAGACAAGCCGTGGGATAAGCATCGCGCTAAT
...

$ pgap.py -r Limnothrix_sp_HT2024_plasmid_GENOME.yaml
Output will be placed in: /mnt/c/Users/m_gruenstaeudl/OneDrive - Fort Hays State University/Desktop/temp/plasmid/output
PGAP version 2025-05-06.build7983 is up to date.
PGAP completed successfully.
```

#### Merging sequence annotations of RAST and PGAP

[Merging PAGP-generated and own annotations](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/data_STEP04__Annotation_evaluation/Step4__PYSCRIPT_Merge_GenBank_tags.py)
```bash
python Step4__Merge_GenBank_tags.py \
    Limnothrix_sp_HT2024_plasmid_PGAP.gbk \
    Limnothrix_sp_HT2024_plasmid_GENOME.gb \
    Limnothrix_sp_HT2024_plasmid_MERGED.gb
```


#### Evaluating quality and correcting sequence annotations

- [Evaluate if the reading frames of the genes of the genome are intact](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/data_STEP04__Annotation_evaluation/Step1__PYSCRIPT_Evaluate_reading_frames_of_genes.py)
```bash
# Bacterial genome
python Step1__PYSCRIPT_Evaluate_reading_frames_of_genes.py -i Limnothrix_sp_HT2024_Bactopia.gb -o Limnothrix_sp_HT2024_Bactopia_ANNOTATION-INFO
# Bacterial plasmid
python Step1__PYSCRIPT_Evaluate_reading_frames_of_genes.py -i Limnothrix_sp_HT2024_plasmid.gb -o Limnothrix_sp_HT2024_plasmid_ANNOTATION-INFO
```

- [Compare gene set of two input genomes by gene name and start-position proximity](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/data_STEP04__Annotation_evaluation/Step2__PYSCRIPT_Compare_genes_by_name_and_position.py)
```bash
# Bacterial genome
python Step2__PYSCRIPT_Compare_genes_by_name_and_position.py Limnothrix_sp_HT2024_Bactopia.gb Limnothrix_sp_HT2024_bacass.gb --max-start-diff 500
```

- [Standardize the annotations of a bacterial genome](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/data_STEP04__Annotation_evaluation/Step3__PYSCRIPT_Standardize_annotations_of_bacterial_genome.py)
This script ensures that every `CDS` and every `gene` annotation contain at least a `gene`-tag as well as a `product`-tag. The `gene`-tag contains the four-letter gene abbreviation. The full behaviour of the script is as follows:
```bash
# Bacterial genome
python Step3__PYSCRIPT_Standardize_annotations_of_bacterial_genome.py input.gb output.gb
grep -v '/gene="unknown_gene"' output.gb | grep -v "/locus_tag=" > output_final.gb
# Bacterial plasmid
python Step3__PYSCRIPT_Standardize_annotations_of_bacterial_genome.py Limnothrix_sp_HT2024_plasmid.gb Limnothrix_sp_HT2024_plasmid_TMP.gb
grep -v '/gene="unknown_gene"' Limnothrix_sp_HT2024_plasmid_TMP.gb | grep -v "/locus_tag=" > Limnothrix_sp_HT2024_plasmid_FINAL.gb
```

##### Behavior Table
| Situation                                                                                         | Action                                                                                                                   | Log style                                                             |
| ------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------- |
| CDS already has valid `gene` and `product`                                                        | Copy missing values to the paired `gene` feature; CDS remains authoritative                                              | White if nothing changes, yellow if synchronization changes something |
| Only the `gene` feature has valid `gene` and `product`                                            | Copy missing values to the paired `CDS` feature                                                                          | Yellow if this changes the CDS, otherwise white                       |
| `gene` and `CDS` disagree                                                                         | Prefer CDS values and record the conflict in the report                                                                  | Yellow if resolved successfully, red if still unresolved              |
| Valid `product` exists but `gene` is missing                                                      | Try local mapping first; otherwise query UniProt **cyanobacteria** by product to infer the most common gene abbreviation | Yellow on successful resolution, red on failure                       |
| Valid `gene` exists but `product` is missing                                                      | Try local mapping first; otherwise query UniProt **cyanobacteria** by gene to infer the most common product description  | Yellow on successful resolution, red on failure                       |
| `standard_name` is present                                                                        | Use it as supporting information during local resolution and as the displayed name in logs                               | Shown as `standard_name` instead of `locus_tag`                       |
| `standard_name` is `hypothetical protein CDS` or `hypothetical protein gene`                      | Do **not** query UniProt and do **not** log the annotation; still standardize `/product` to `hypothetical protein`       | No log output                                                         |
| Product is already `hypothetical protein` for one of those hypothetical-standard-name annotations | Leave it as `hypothetical protein` or rewrite it to the same standard form                                               | No log output                                                         |
| Nothing reliable can be inferred                                                                  | Fall back to unresolved values such as `unknown_gene` or remaining missing information, and record it in the report      | Red                                                                   |
| Annotation does not change                                                                        | Keep existing values as they are                                                                                         | White                                                                 |

##### Priority Order

| Priority | Rule                                                                                                                                                                                              |
| -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1        | Prefer existing CDS qualifiers over gene qualifiers                                                                                                                                               |
| 2        | Copy missing values across paired `gene` and `CDS` features                                                                                                                                       |
| 3        | Apply the local mapping table                                                                                                                                                                     |
| 4        | Query UniProt cyanobacteria by product if `gene` is missing                                                                                                                                       |
| 5        | Query UniProt cyanobacteria by gene if `product` is missing                                                                                                                                       |
| 6        | Skip logging and UniProt lookup for annotations whose `standard_name` is `hypothetical protein CDS` or `hypothetical protein gene`, but still standardize their product to `hypothetical protein` |
| 7        | Record conflicts and unresolved cases in the report                                                                                                                                               |

---

### 05. Gene-Level Visualization
- [Visualization of gene location via GenoVi](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP05__Visualization_via_GenoVi.md)

<div style="display:flex; justify-content:center; gap:2%;">
    <img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP05__Visualization_via_GenoVi/Limnothrix_sp_HT2024_chromosome__combined_vertical.png" style="width:48%;">
    <img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP05__Visualization_via_GenoVi/Limnothrix_sp_HT2024_plasmid__combined_vertical.png" style="width:48%;">
</div>


---

### 06. Coverage Visualization
- [Visualization of sequencing coverage via Circleator](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP06__Coverage_Viz_Circleator.md)

<div style="display:flex; justify-content:center; gap:2%;">
    <img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP06__Coverage_Viz_Circleator/output/chromosome.circleator.final.svg" style="width:48%;">
    <img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP06__Coverage_Viz_Circleator/output/plasmid.circleator.final.svg" style="width:4.8%;">
</div>

---

### 07. k-mer Spectrum Analysis
- [Visualization of k-mer spectra via Jellyfish](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP07a__Kmer_spectrum_Jellyfish.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP07a__Kmer_spectrum_Jellyfish/kmer_spectrum_Illumina_k21_ONT_k17_combined.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

###### **TO DO:** Extend the x-axis in this figure.

- [Visualization of k-mer spectra via Merqury](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP07b__Kmer_spectrum_Merqury.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP07b__Kmer_spectrum_Merqury/output/merqury_compare_2026-01-30.BothGenomes.spectra-cn.fl.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

###### **TO DO:** Create the same figure based on the Nanopore data; maybe create a facet plot for both.

---

### 08. Whole-Genome Alignment
- Show dotplots produced by MUMmer4

- [Visualization of MUMmer4 results as dotplots](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP08a__MUMmer4_dotplots.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP08a__GenomicInversions_MUMmer4/Limnothrix_sp_BLA16_vs_BacterialChr__MUMmer4.dotplot.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">


- Show inversions within the assembly using Circos

- [Show synteny and collinearity between *Limnothrix* B-16 and the assembly using Circos](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP08c__Circos_collinearity_visualization.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP08c__GenomicInversions_Circos/process_and_output/circos/circos.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

---

### 10. Metagenomic analysis of 16S rRNA amplicon data using QIIME2
- [Metagenomic analysis](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP10__QIIME2_analysis.md)

| Metric             | HorseThief_lake_water_2 | HorseThief_in_cultivation |
| ------------------ | ----------------------- | ------------------------- |
| Raw reads          | 54,376                  | 91,377                    |
| ASV method         | DADA2                   | DADA2                     |
| Number of ASVs     | 201                     | 272                       |
| Number of genera   | 53                      | 67                        |
| Number of families | 40                      | 46                        |
| Shannon diversity  | 6.570                   | 6.341                     |
| Pielou evenness    | 0.859                   | 0.784                     |
| Simpson diversity  | 0.983                   | 0.972                     |


---

### 11. Preparing submission to GenBank

- [Submission preparation to GenBank via Geneious submission plugin](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP11__Submission_via_Geneious_plugin.md)

---

## Notes

- Image files are stored in the corresponding `data_STEPXX__*` directories.
- All workflows assume Linux environments and standard bioinformatics toolchains.
- Individual steps can be adapted to other bacterial genomes with minimal modification.
