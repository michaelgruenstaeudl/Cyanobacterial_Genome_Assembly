### Barchart-style gene tally compatible with GenoVi visualization

```python
import os, pandas as pd, matplotlib.pyplot as plt, tkinter as tk
from tkinter import filedialog

# Minimal seaborn-style white grid
plt.style.use("seaborn-v0_8-whitegrid")

# ----------------------------
# Select input CSV interactively
# ----------------------------
root = tk.Tk(); root.withdraw()
csv_path = filedialog.askopenfilename(
    title="Select COG classification CSV file",
    filetypes=[("CSV files","*.csv"),("All files","*.*")]
)
if not csv_path:
    raise FileNotFoundError("No input file selected.")

# ----------------------------
# Read CSV and extract Replicon / Total rows
# ----------------------------
raw = pd.read_csv(csv_path)
replicon_row = raw[raw.iloc[:,0] == "Replicon"]
total_row     = raw[raw.iloc[:,0] == "Total"]
if replicon_row.empty or total_row.empty:
    raise ValueError("CSV must contain rows starting with 'Replicon' and 'Total'.")

# Build tidy table: one COG per row
df = pd.DataFrame({
    "COG":   replicon_row.iloc[0,1:].dropna().astype(str).str.strip().values,
    "Count": total_row.iloc[0,1:].dropna().astype(int).values
})

# ----------------------------
# COG functional descriptions
# ----------------------------
cog_names = {
"D":"Cell cycle control, cell division, chromosome partitioning",
"M":"Cell wall/membrane/envelope biogenesis","N":"Cell motility",
"O":"Posttranslational modification, protein turnover, chaperones",
"T":"Signal transduction mechanisms",
"U":"Intracellular trafficking, secretion, and vesicular transport",
"V":"Defense mechanisms","W":"Extracellular structures","Y":"Nuclear structure",
"Z":"Cytoskeleton","A":"RNA processing and modification",
"B":"Chromatin structure and dynamics",
"J":"Translation, ribosomal structure and biogenesis","K":"Transcription",
"L":"Replication, recombination and repair","X":"Mobilome: prophages, transposons",
"C":"Energy production and conversion","E":"Amino acid transport and metabolism",
"F":"Nucleotide transport and metabolism","G":"Carbohydrate transport and metabolism",
"H":"Coenzyme transport and metabolism","I":"Lipid transport and metabolism",
"P":"Inorganic ion transport and metabolism",
"Q":"Secondary metabolites biosynthesis, transport and catabolism",
"R":"General function prediction only","S":"Function unknown",
"Unclassified":"Unclassified"
}

# Build x-axis labels as "description [COG]"
df["Label"] = df["COG"].apply(
    lambda x: f"{cog_names.get(x, x)} [{x}]"
)

# ----------------------------
# GenoVi "strong" palette (resolved from Circos config)
# ----------------------------
cog_color_map = {
"D":"#637BB7","M":"#2659A8","N":"#6497B0","O":"#2084AE","T":"#42A9B3","U":"#126974",
"V":"#61C3A6","W":"#259775","Z":"#15793C","J":"#A660A7","K":"#9D3F97",
"L":"#8574B5","X":"#4F489E","C":"#AAD382","E":"#7DB040","F":"#B2B36D",
"G":"#8F8A2F","H":"#EBBC86","I":"#AF7E35","P":"#DB8856","Q":"#C46426",
"R":"#626262","S":"#909090"
}
df["Color"] = df["COG"].map(cog_color_map).fillna("#B0B0B0")

# ----------------------------
# Bar plot
# ----------------------------
fig, ax = plt.subplots(figsize=(14,6))
bars = ax.bar(df["Label"], df["Count"], color=df["Color"], width=0.8)

ax.set_title("COG Category", fontsize=14)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(
    df["Label"], rotation=65,
    ha="right", va="top", rotation_mode="anchor"
)
ax.tick_params(axis="x", labelsize=8, pad=10)
ax.margins(x=0.02); plt.subplots_adjust(bottom=0.36)

# Annotate counts
for b in bars:
    h = b.get_height()
    ax.text(
        b.get_x()+b.get_width()/2, h+5, f"{int(h)}",
        ha="center", va="bottom", fontsize=8
    )

# ----------------------------
# Save SVG next to input CSV
# ----------------------------
plt.tight_layout()
out_svg = os.path.join(
    os.path.dirname(csv_path),
    f"{os.path.splitext(os.path.basename(csv_path))[0]}_COG_barplot.svg"
)
fig.savefig(out_svg, format="svg", dpi=300, bbox_inches="tight")
print(f"Saved SVG to: {out_svg}")

plt.show()
```