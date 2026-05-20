### Inference of QV for genome assembly via Mercury

#### Installation of Mercury on Beocat
```bash
# Load compiler
module load GCC
# Installation of Mercury
conda config --set solver classic
conda config --set channel_priority strict
env -u LD_LIBRARY_PATH conda create -n merqury -c conda-forge -c bioconda merqury --solver=classic
# Activation of Mercury
conda activate merqury
```

#### Running Merqury

##### Run script `run_QualEval_via_Merqury.sh` on Beocat
```bash
#!/bin/bash
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --job-name=mercury_plasmid
#SBATCH --mail-type=all
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10

bash run_QualEval_via_Merqury.sh
```

##### File Hygiene
Note: This script runs from OUTSIDE $OUT (do not cd).

```bash
OUT=plasmid_Illumina_MERQURY_report

shopt -s nullglob

# Compress all meryl databases created by the QV workflow:
#   *.reads.meryl  *.asm.meryl  *.asm.only.meryl
for mdb in "$OUT"/*.meryl; do
  if [[ -d "$mdb" ]]; then
    tarball="${mdb}.tar.gz"

    echo "  -> Compressing $(basename "$mdb")"
    tar -czf "$tarball" -C "$OUT" "$(basename "$mdb")"

    # Verify archive exists and is non-empty
    if [[ -s "$tarball" ]]; then
      rm -rf "$mdb"
      echo "     Removed original $(basename "$mdb")"
    else
      echo "     ERROR: Failed to create $(basename "$tarball"); keeping original" >&2
    fi
  fi
done

shopt -u nullglob
```