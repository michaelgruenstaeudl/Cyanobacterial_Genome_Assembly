### Visualize synteny and collinearity across assemblies using Circos

#### Installation
```bash
# Installation of Circos and MUMmer
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

conda create -n circos_compare -c conda-forge -c bioconda \
  python=3.11 biopython mummer4 circos

conda activate circos_compare

conda install -y -c conda-forge mummer
```


#### Set up a folder with the two input genomes
```bash
mkdir -p work
cd work
GENOME1=Corrected_Limnothrix_sp_BL-A-16_CP166615.gbk
GENOME2=FinalAssembly_Bactopia__corrected_withTranslations.gbk
```

#### STEP 1. Converting GenBank file to set of input files expected by Circos
A script named `gbk_to_circos.py` is used to convert GenBank file to the following files expected by Circos:
- FASTA file 
- Circos karyotype information 
- simple gene track

##### gbk_to_circos.py
```python

#!/usr/bin/env python3
from pathlib import Path
from Bio import SeqIO

def safe_id(s: str) -> str:
    # Circos ids should be simple: letters, numbers, underscore
    out = []
    for ch in s:
        if ch.isalnum() or ch == "_":
            out.append(ch)
        elif ch in " .:-|/":
            out.append("_")
    return "".join(out).strip("_") or "seq"

def main(gbk_path: str, prefix: str):
    gbk = Path(gbk_path)
    records = list(SeqIO.parse(str(gbk), "genbank"))
    if not records:
        raise SystemExit(f"No records found in {gbk}")

    # Write FASTA (one entry per record/replicon)
    fasta_path = Path(f"{prefix}.fna")
    SeqIO.write(records, str(fasta_path), "fasta")

    # Circos karyotype + simple gene track
    karyo_path = Path(f"{prefix}.karyotype.txt")
    genes_path = Path(f"{prefix}.genes.txt")

    with karyo_path.open("w", encoding="utf-8") as kf, genes_path.open("w", encoding="utf-8") as gf:
        for i, rec in enumerate(records, start=1):
            rec_id_raw = rec.id or rec.name or f"{prefix}_{i}"
            rec_id = safe_id(f"{prefix}_{rec_id_raw}")
            length = len(rec.seq)

            # Karyotype line format:
            # chr - <id> <label> <start> <end> <color>
            label = f"{prefix}:{rec_id_raw}"
            color = "vdgrey" if prefix.lower().endswith("a") else "vlgrey"
            kf.write(f"chr - {rec_id} {label} 0 {length} {color}\n")

            # Genes (very simple): write CDS as intervals for a heatmap/histogram track
            # Format for a "text" or "highlight"-style file depends on Circos plot type,
            # but for many plots: <chr> <start> <end> <value>
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                loc = feat.location
                if loc is None:
                    continue
                start = int(loc.start)
                end = int(loc.end)
                if end <= start:
                    continue
                gf.write(f"{rec_id}\t{start}\t{end}\t1\n")

    print(f"Wrote: {fasta_path}")
    print(f"Wrote: {karyo_path}")
    print(f"Wrote: {genes_path}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        raise SystemExit("Usage: gbk_to_circos.py <genome.gbk> <prefix>\nExample: gbk_to_circos.py genomeA.gbk A")
    main(sys.argv[1], sys.argv[2])
```

##### Running `gbk_to_circos.py`
```bash
python gbk_to_circos.py $GENOME1 genome1
python gbk_to_circos.py $GENOME2 genome2
```

#### STEP 2. Align genomes (MUMmer) and convert to Circos links

```bash
nucmer --maxmatch -p genome1_vs_genome2 genome1.fna genome2.fna
show-coords -rclT genome1_vs_genome2.delta > genome1_vs_genome2.coords.tsv
```

#### STEP 3. Convert coordinates to links
Conversion is conducted by Python script `coords_to_links.py`:

##### Script `coords_to_links.py`
```python

#!/usr/bin/env python3
from pathlib import Path

def main(coords_tsv: str, out_links: str, min_len: int = 5000, min_id: float = 95.0):
    inp = Path(coords_tsv)
    out = Path(out_links)

    with inp.open("r", encoding="utf-8") as f, out.open("w", encoding="utf-8") as g:
        for line in f:
            line = line.strip()
            if not line or line.startswith("[") or line.startswith("="):
                continue

            # show-coords -T output is tab-separated. Typical columns include:
            # 0:S1 1:E1 2:S2 3:E2 4:LEN1 5:LEN2 6:%IDY  ...  last two: REF_QRY sequence IDs
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            try:
                s1 = int(parts[0]); e1 = int(parts[1])
                s2 = int(parts[2]); e2 = int(parts[3])
                l1 = int(parts[4])
                pid = float(parts[6])
                ref = parts[-2].strip()
                qry = parts[-1].strip()
            except Exception:
                continue

            if l1 < min_len or pid < min_id:
                continue

            # Circos wants: <chr1> <start1> <end1> <chr2> <start2> <end2>
            # Keep start <= end for plotting
            a1, b1 = (s1, e1) if s1 <= e1 else (e1, s1)
            a2, b2 = (s2, e2) if s2 <= e2 else (e2, s2)

            # Detect inversion in second genome
            is_inversion = s2 > e2
            inv_flag = "1" if is_inversion else "0"

            # Determine inversion
            inv = 1 if (s2 > e2) else 0
            g.write(f"{ref}\t{a1}\t{b1}\t{qry}\t{a2}\t{b2}\tinv={inv}\n")



    print(f"Wrote links: {out}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        raise SystemExit("Usage: coords_to_links.py <A_vs_B.coords.tsv> <links.txt> [min_len] [min_id]\n"
                         "Example: coords_to_links.py A_vs_B.coords.tsv links.txt 1000 90")
    min_len = int(sys.argv[3]) if len(sys.argv) >= 4 else 5000
    min_id  = float(sys.argv[4]) if len(sys.argv) >= 5 else 95.0
    main(sys.argv[1], sys.argv[2], min_len=min_len, min_id=min_id)
```

##### Running `coords_to_links.py`
```bash
python coords_to_links.py genome1_vs_genome2.coords.tsv links.txt 5000 95
```

Important: the `ref` and `qry` IDs in `links.txt` must match the IDs in `karyotype.txt`.
If your IDs do not match, the easiest fix is to ensure the FASTA headers written by `gbk_to_circos.py` match the karyotype IDs. The script prefixes records as `genome1_<recordid>` and `genome2_<recordid>`, and those same ids should appear in the FASTA headers.


#### STEP 4. Correcting the IDs in the karyotype files
The correction is conducted by Python script `fix_karyotype_ids_from_links.py`:

##### Script `fix_karyotype_ids_from_links.py`
```python

#!/usr/bin/env python3

from pathlib import Path
import sys

def get_link_ids(links_path: Path):
    with links_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 4:
                return fields[0], fields[3]
    return None, None


def fix_karyotype(karyo_in: Path, karyo_out: Path, new_id: str):
    with karyo_in.open() as fin, karyo_out.open("w") as fout:
        for line in fin:
            if line.startswith("chr"):
                parts = line.rstrip("\n").split()
                # Circos karyotype format:
                # chr - <id> <label> <start> <end> <color>
                parts[2] = new_id
                fout.write(" ".join(parts) + "\n")
            else:
                fout.write(line)


def main():
    links = Path("links.txt")
    if not links.exists():
        sys.exit("ERROR: links.txt not found")

    id1, id2 = get_link_ids(links)
    if id1 is None:
        print("links.txt is empty, no changes made")
        return

    print(f"Using IDs from links.txt: {id1}, {id2}")

    fix_karyotype(
        Path("genome1.karyotype.txt"),
        Path("genome1.karyotype.fixed.txt"),
        id1
    )

    fix_karyotype(
        Path("genome2.karyotype.txt"),
        Path("genome2.karyotype.fixed.txt"),
        id2
    )

    print("Wrote:")
    print("  genome1.karyotype.fixed.txt")
    print("  genome2.karyotype.fixed.txt")


if __name__ == "__main__":
    main()

```

##### Run script `fix_karyotype_ids_from_links.py`
```bash
python fix_karyotype_ids_from_links.py
```


#### STEP 5. Prepare inputs, configure Circos, run Circos

##### Assemble input files for Circos in specific folder `circos`
```bash
mkdir -p circos
cat genome1.karyotype.fixed.txt genome2.karyotype.fixed.txt > karyotype.txt
cp karyotype.txt links.txt genome1.genes.txt genome2.genes.txt circos/
cd circos
```

##### Manually create a minimal Circos configuration file `circos.conf`
```conf
karyotype = karyotype.txt

<image>
dir   = .
file  = circos
png   = yes
svg   = yes
radius = 1500p
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>

chromosomes_units = 100000
chromosomes_display_default = yes

# Put genome1 and genome2 on opposite halves
# If you have multiple replicons, list them explicitly in the order you want.
# Example:
# chromosomes = genome1_NC_000000_1;genome2_NC_000000_1
# chromosomes_order = genome1_NC_000000_1,genome2_NC_000000_1
# For a quick start, you can omit and let Circos place everything.

<ideogram>
<spacing>
default = 0.02r
</spacing>
radius    = 0.85r
thickness = 20p
fill      = yes
stroke_color = dgrey
stroke_thickness = 1p
show_label = yes
label_font = default
label_radius = 1.02r
label_size = 18p
</ideogram>

# Optional: plot gene density-like bars (here just CDS intervals as value=1)
<plots>
<plot>
type = histogram
file = genome1.genes.txt
r1 = 0.84r
r0 = 0.72r
min = 0
max = 1
thickness = 1
extend_bin = no
</plot>

<plot>
type = histogram
file = genome2.genes.txt
r1 = 0.28r
r0 = 0.16r
min = 0
max = 1
thickness = 1
extend_bin = no
</plot>
</plots>


<links>
  <link>
    file   = links.txt
    radius = 0.90r

    color = blue

    <rules>
      <rule>
        condition = var(inv) == 1
        color     = red
      </rule>
    </rules>

  </link>
</links>
```

##### Run Circos
```bash
# Activate circos
conda activate circos_compare
# Run circos
circos -conf circos.conf
```

The final image is saved as ./circos/circos.png
