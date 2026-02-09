# motif-anchor-analysis

Investigates whether SNPs within CTCF motifs explain discrepancies between empirical (ChIP-seq) and predicted (AlphaGenome) allele-specific CTCF binding at loop anchors in GM12878 and HepG2 cell lines.

## Pipeline

The master script `SCR_sample_vs_motif-strength.sh` runs every step in order:

```bash
./SCR_sample_vs_motif-strength.sh <output_folder>
```

### Step 1 — Sample anchors (`SCR_from-colab.py`)

Reads full ChIP-seq-vs-AlphaGenome comparison tables for GM12878 and HepG2. Filters to anchors outside a low-signal exclusion zone, then randomly samples N "good" anchors (predicted and observed agree) and N "bad" anchors (along the x-axis only, i.e. observed differential binding but no predicted difference). Outputs a combined TSV of selected anchors plus BED files for each cell line.

### Step 2 — Intersect with CTCF motifs (`SCR_intersect_motifs_with_anchors.sh`)

Uses `bedtools intersect` to overlap the selected anchor BED files with diploid CTCF motif calls (MA0139.1), producing per-cell-line TSVs of anchors annotated with overlapping motif SNP details (ref/alt scores, delta scores, sequences).

### Step 3 — Analyze SNP effects on motifs (`SCR_analyze-snp-motif-anchors.py`)

Parses the intersect output for each genome (GM or HG2). Classifies every motif-overlapping SNP by absolute delta score into three tiers (< 1, < 5, >= 5) and summarizes each anchor by its strongest overlapping SNP effect. Outputs a per-anchor summary and a detailed per-SNP table.

### Step 4 — Plot results (`SCR_plot_examples.py`)

For each cell line, generates a 5-panel figure:

1. Scatter of empirical vs predicted log2 fold-change, colored by variant
2. Box plot of predicted log2 diff grouped by SNP effect tier
3. Scatter of SNP effect vs predicted log2 diff
4. Confusion heatmap of SNP effect vs prediction quality (all rows)
5. Confusion heatmap of SNP effect vs prediction quality (unique anchors only)

Also produces a combined heatmap across both cell lines.

### Step 5 — Combine details (`FUNC_make-file-with-details.py`)

Joins the per-SNP detail tables from both cell lines with the selected-anchor metadata into a single combined TSV.

## Utility script

- **`FUNC_get_variants.py`** — Extracts variant IDs (chrom_pos_ref_alt) from the motif intersection file for downstream lookups.

## Data files

- `GM12878_5-0.1_outside.tsv` / `HepG2_5-0.1_outside.tsv` — Pre-filtered anchors outside the exclusion zone with large empirical log2 differences, used as the "bad prediction" pool when `have_bad=yes`.

## Requirements

- Python 3 with `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`
- `bedtools`
