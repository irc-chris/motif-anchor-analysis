# motif-anchor-analysis

Investigates whether SNPs within CTCF motifs explain discrepancies between empirical (ChIP-seq) and predicted (AlphaGenome) allele-specific CTCF binding at loop anchors in GM12878 and HepG2 cell lines.

## Pipeline

The master script `SCR_motif-strength_vs_prediction.sh` runs every step in order:

```bash
./SCR_motif-strength_vs_prediction.sh <output_folder_prefix>
```

The script accepts configurable thresholds for filtering at the top:
- `motif_threshold=8` — Minimum motif effect delta score to include in final filtered outputs and last heatmap
- `chip_diff_threshold=3` — ChIP-seq log2 fold-change threshold for analyzing big changes
- `empirical_threshold=15` — Empirical either homolog threshold for anchor selection
- `bad_ag_threshold=0.1` — AlphaGenome log2 maximum for "bad" predictions to be from 0
- `good_ag_threshold=1.0` — AlphaGenome log2 maximum for "good" predictions to be within from empirical

Output folder is automatically named with these thresholds appended (e.g., `prefix-c15-cd3-pb0.1-pg1.0-m8`).

### Step 1 — Sample anchors (`SCR_from-claud.py`)

Reads full ChIP-seq-vs-AlphaGenome comparison tables for GM12878 and HepG2. Filters anchors by the configured thresholds, then randomly samples anchors into "good" (predicted and observed agree) and "bad" (observed differential binding but no predicted difference) categories. Outputs a combined TSV of selected anchors plus BED files for each cell line.

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

### Step 5 — Filter and process output files

Filters combined anchor-SNP data files by `motif_threshold` using `awk`, producing two additional TSVs containing only high-effect motifs. Additionally:

- Reorders columns in all TSVs using `../prep_scripts/reorder-cols.py` for consistency
- Sorts all TSVs by genomic position using `../prep_scripts/sorting.py`
- Generates reference-allele read-bias plots with `SCR_plot_ref_reads.py`

## Utility scripts

- **`FUNC_get_variants.py`** — Extracts variant IDs (chrom_pos_ref_alt) from the motif intersection file for downstream lookups.
- **`SCR_from-claud.py`** — Alternative anchor sampling script (currently used; replaces `SCR_from-colab.py`)

## Data files

- `GM12878_5-0.1_outside.tsv` / `HepG2_5-0.1_outside.tsv` — Pre-filtered anchors outside the exclusion zone with large empirical log2 differences, used as the "bad prediction" pool when `have_bad=yes`.

## Requirements

- Python 3 with `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`
- `bedtools`
