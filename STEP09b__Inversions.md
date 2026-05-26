

#### Visualization of dotplots via matplotlib

##### Plot from `show-coords` TSV

```bash
python plot_coords_dotplot.py Limnothrix_sp_BLA16_vs_BactopiaAssembly.coords.tsv Limnothrix_sp_BLA16_vs_BactopiaAssembly.dotplot
```

##### Plot directly from filtered `.delta` file

```bash
python plot_delta_dotplot.py Limnothrix_sp_BLA16_vs_BactopiaAssembly.filt.delta Limnothrix_sp_BLA16_vs_BactopiaAssembly.delta_dotplot
```

#### Numeric inference of inversion locations

##### Identifying inversion locations
```bash
conda install -c bioconda bedtools
```

```bash
show-coords -rclTH Limnothrix_sp_BLA16_vs_BactopiaAssembly.filt.delta \
| awk -F'\t' '$3 > $4 && $5 >= 1000 && $7 >= 95 {
    s=$1; e=$2; if (s>e) {t=s; s=e; e=t}
    # BED-like: chrom start end (tab-delimited)
    print $12"\t"s"\t"e
}' \
| sort -k1,1 -k2,2n \
| bedtools merge -i - -d 10000 \
| awk 'BEGIN{OFS="\t"} {len=$3-$2+1; print $1,$2,$3,len}' \
| sort -k4,4nr \
| head -n 6 \
| sort -k1,1 -k2,2n \
| awk 'BEGIN{OFS=","; print "reference,start,end,length"} {print $1,$2,$3,$4}' \
> six_major_inversions_ref.csv
```

##### RESULTS
```text
Limnothrix_sp_BL_A_16_CP166615	82994	113322	30329
Limnothrix_sp_BL_A_16_CP166615	592628	736045	143418
Limnothrix_sp_BL_A_16_CP166615	1756173	1817785	61613
Limnothrix_sp_BL_A_16_CP166615	2674925	2733309	58385
Limnothrix_sp_BL_A_16_CP166615	3085082	3131342	46261
Limnothrix_sp_BL_A_16_CP166615	4038055	4191007	152953
```