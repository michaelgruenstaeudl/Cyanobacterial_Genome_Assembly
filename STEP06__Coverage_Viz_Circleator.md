### Visualization of sequencing coverage via Circleator

##### Installation of Circleator and dependecies

###### On Debian
```bash
# Install dependencies via conda
conda activate base
conda install -n base --solver=libmamba \
-c conda-forge -c bioconda \
perl \
perl-app-cpanminus \
perl-bioperl-core \
perl-bioperl \
perl-bio-featureio \
perl-clone \
perl-json \
perl-log-log4perl \
perl-svg \
perl-text-csv \
perl-module-build \
perl-cache-cache \
perl-digest-sha1 \
vcftools \
samtools
# Install remaining BioPerl modules if necessary
cpanm --force Bio::Perl
# Install Circleator
curl -L -o Circleator-1.0.2.tar.gz \
  https://github.com/jonathancrabtree/Circleator/archive/refs/tags/1.0.2.tar.gz
tar xzvf Circleator-1.0.2.tar.gz
cd Circleator-1.0.2
which perl
perl Build.PL
./Build
./Build test
./Build install
```

---

##### Correction of LOCUS name

Before running Circleator, ensure the GenBank LOCUS name and BAM sequence name are short and identical.
If `samtools` prints `libtinfow`/`libncursesw` "no version information available", that warning is usually harmless.

```bash
# Build Circleator-ready GenBank and BAM files (skip any missing inputs)
IN_DIR="backmapping_results"
GB_DIR="input"

for GENOME in chromosome plasmid; do
  SHORT_ID="Limnothrix"

  IN_GB="${GB_DIR}/Limnothrix_sp_HT2024_${GENOME}.gb"
  if [ ! -f "$IN_GB" ]; then
    echo "SKIP: $IN_GB not found"
    continue
  fi

  cp "$IN_GB" "${GENOME}.circleator.gb"
  sed -i -E "s/^(LOCUS[[:space:]]+)[^[:space:]]+/\1${SHORT_ID}/" "${GENOME}.circleator.gb"

  for PLATFORM in illumina ont; do
    IN_BAM="${IN_DIR}/${GENOME}.${PLATFORM}.sorted.bam"
    if [ ! -f "$IN_BAM" ]; then
      echo "SKIP: $IN_BAM not found"
      continue
    fi

    samtools view -H "$IN_BAM" \
    | awk -v sid="$SHORT_ID" 'BEGIN{done=0} /^@SQ/ && !done {sub(/SN:[^ \t]+/,"SN:" sid); done=1} {print}' \
    > "${GENOME}.${PLATFORM}.header.sam"

    samtools reheader "${GENOME}.${PLATFORM}.header.sam" "$IN_BAM" \
    > "${GENOME}.${PLATFORM}.circleator.bam"

    samtools index "${GENOME}.${PLATFORM}.circleator.bam"
  done
done

rm -f chromosome.illumina.header.sam chromosome.ont.header.sam
rm -f plasmid.illumina.header.sam plasmid.ont.header.sam
```

---

##### Configuration file for Circleator
Circleator draws figures based on the information in its configuration file.
The following script creates one config file per genome with both Illumina and ONT coverage tracks.

```bash
for GENOME in chromosome plasmid; do
  SHORT_ID="Limnothrix"

  [ -f "${GENOME}.illumina.circleator.bam" ] || continue
  [ -f "${GENOME}.ont.circleator.bam" ] || continue

  cat > "${GENOME}.circleator.conf" <<EOF
## Circleator config with both coverage tracks
small-cgap
## Inner coverage ring (Illumina, red)
new cov_illumina graph 0.14 graph-function=BAMCoverage,bam-file=${GENOME}.illumina.circleator.bam,bam-seqid=${SHORT_ID},graph-min=0,graph-max=data_max,window-size=5000,heightf=0.25,opacity=0.85,color1=#e41a1c

small-cgap
## Inner coverage ring (ONT, blue)
new cov_ont graph 0.14 graph-function=BAMCoverage,bam-file=${GENOME}.ont.circleator.bam,bam-seqid=${SHORT_ID},graph-min=0,graph-max=data_max,window-size=10000,heightf=0.25,opacity=0.85,color1=#3aa0d5
EOF
done
```

##### Running Circleator
```bash
for GENOME in chromosome plasmid; do
  SHORT_ID="Limnothrix"
  [ -f "${GENOME}.circleator.gb" ] || continue
  [ -f "${GENOME}.circleator.conf" ] || continue

  N1=$(samtools view -c -F 4 "${GENOME}.illumina.circleator.bam" "$SHORT_ID" 2>/dev/null || echo 0)
  N2=$(samtools view -c -F 4 "${GENOME}.ont.circleator.bam" "$SHORT_ID" 2>/dev/null || echo 0)

  if [ "$N1" -gt 0 ] || [ "$N2" -gt 0 ]; then
    circleator --data="${GENOME}.circleator.gb" --config="${GENOME}.circleator.conf" > "${GENOME}.circleator.svg"
  else
    echo "SKIP: ${GENOME} has no mapped reads on $SHORT_ID"
  fi
done
```

---

##### Correcting the raw SVG
Circleator produces an SVG in the old SVG 1.0 format that is not rendered correctly in today's SVG editors. Hence, the output of Circleator must be rasterized first using [Apache Batik](https://xmlgraphics.apache.org/batik/download.html).

For that, download the [Batik binary](https://www.apache.org/dyn/closer.cgi?filename=/xmlgraphics/batik/binaries/batik-bin-1.19.tar.gz&action=download) into the same directopry as the output of Circleator and unzip it.

```bash
# Set BATIK_HOME and confirm you are in the right place
export BATIK_HOME="$PWD/batik-1.19"
ls "$BATIK_HOME"/lib | head

# Convert all Circleator SVG outputs to PDF
for svg in *circleator.svg; do
  [ -e "$svg" ] || continue
  java -cp "$BATIK_HOME/lib/*:$BATIK_HOME/extensions/*:$BATIK_HOME/batik-rasterizer-1.19.jar" \
    org.apache.batik.apps.rasterizer.Main \
    -m application/pdf \
    -scriptSecurityOff \
    -d "${svg%.svg}.pdf" \
    "$svg"
done
```

---

##### Generate a color legend using Python
Circleator does not produce any legends for the figures it produces. Hence, the legends must be generated separately by the user.

```bash
pip install svgwrite
python3 legend/Circleator_legend_maker.py
```

---

##### Combine figure with legend

```bash
# Convert Batik PDF outputs back to SVG (requires pdftocairo from poppler-utils)
# Debian/Ubuntu: sudo apt install -y poppler-utils
for GENOME in chromosome plasmid; do
  [ -f "${GENOME}.circleator.pdf" ] || continue
  pdftocairo -svg "${GENOME}.circleator.pdf" "${GENOME}.circleator.frompdf.svg"
done

# Combine each figure SVG with the legend SVG into one final SVG
pip install svgutils
python3 - <<'PY'
from pathlib import Path
import re
from svgutils.transform import fromfile, SVGFigure

def px(value: str) -> float:
  m = re.match(r"\s*([0-9.]+)", value or "0")
  return float(m.group(1)) if m else 0.0

legend_candidates = [Path("circleator_legend.svg"), Path("legend/circleator_legend.svg")]
legend_path = next((p for p in legend_candidates if p.exists()), None)
if legend_path is None:
  raise SystemExit("Legend SVG not found (expected circleator_legend.svg or legend/circleator_legend.svg)")

for genome in ("chromosome", "plasmid"):
  fig_path = Path(f"{genome}.circleator.frompdf.svg")
  if not fig_path.exists():
    print(f"SKIP: {fig_path} not found")
    continue

  fig = fromfile(str(fig_path))
  leg = fromfile(str(legend_path))
  w1, h1 = map(px, fig.get_size())
  w2, h2 = map(px, leg.get_size())

  gap = 20
  out_w = w1 + gap + w2
  out_h = max(h1, h2)

  out = SVGFigure(f"{out_w}px", f"{out_h}px")
  fig_root = fig.getroot()
  leg_root = leg.getroot()
  fig_root.moveto(0, (out_h - h1) / 2)
  leg_root.moveto(w1 + gap, (out_h - h2) / 2)
  out.append([fig_root, leg_root])

  out_path = f"{genome}.circleator.final.svg"
  out.save(out_path)
  fig_path.unlink()
  print(f"Wrote {out_path}")
PY
```

---

##### Hygiene cleanup

Remove intermediate files produced during this workflow, while keeping SVG outputs and Circleator config files.

```bash
for GENOME in chromosome plasmid; do
  rm -f "${GENOME}.illumina.header.sam" "${GENOME}.ont.header.sam"
  rm -f "${GENOME}.illumina.circleator.bam" "${GENOME}.ont.circleator.bam"
  rm -f "${GENOME}.illumina.circleator.bam.bai" "${GENOME}.ont.circleator.bam.bai"
  rm -f "${GENOME}.circleator.gb" "${GENOME}.circleator.pdf" "${GENOME}.circleator.frompdf.svg"
done

# Optional: remove downloaded Batik bundle if no longer needed
# rm -rf batik-1.19 Circleator-1.0.2 Circleator-1.0.2.tar.gz
```

