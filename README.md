# motif-anchor-analysis

Analyzes the relationship between SNP-disrupted CTCF motifs and loop anchor prediction accuracy. The pipeline selects a sample of loop anchors (both well-predicted and poorly-predicted), intersects them with CTCF motif variants (AlphaGenome output), classifies the strength of each SNP's effect on the motif, and produces summary tables and plots.

---

## File naming conventions

| Prefix | Meaning |
|--------|---------|
| `SCR_` | **Script** — runnable entry point, takes command-line arguments |
| `FUNC_` | **Function/utility** — helper module called by scripts, or a one-off utility with hardcoded paths |

---

## Pipeline overview

The full pipeline is driven by a single shell script:

```
SCR_sample_vs_motif-strength.sh <output_folder>
```

Internally it runs the following steps in order:

```
SCR_sample_vs_motif-strength.sh <folder>
 ├── SCR_from-colab.py                        # step 1 – sample anchors
 ├── SCR_intersect_motifs_with_anchors.sh      # step 2 – bedtools intersect
 ├── SCR_analyze-snp-motif-anchors.py GM       # step 3 – SNP effect scoring (GM12878)
 ├── SCR_analyze-snp-motif-anchors.py HG2      # step 3 – SNP effect scoring (HepG2)
 ├── SCR_plot_examples.py                      # step 4 – plots
 ├── FUNC_make-file-with-details.py            # step 5 – combine outputs
 └── ../prep_scripts/sorting.py  (external)    # step 6 – sort all TSVs
```

---

## Scripts (`SCR_`)

### `SCR_sample_vs_motif-strength.sh`

**Entry point — runs the full pipeline.**

```bash
bash SCR_sample_vs_motif-strength.sh <output_folder>
```

| Positional arg | Description |
|---|---|
| `$1` | Output folder to create and write all results into |

Creates the output folder, then calls every downstream script in sequence. All other scripts write into this folder.

---

### `SCR_from-colab.py`

**Selects a random sample of anchors** from pre-computed ChIP-seq vs AlphaGenome comparison data, split into "bad" (poorly predicted) and "good" (well predicted) groups.

```bash
python SCR_from-colab.py <folder> <have_bad>
```

| Positional arg | Description |
|---|---|
| `sys.argv[1]` | Output folder |
| `sys.argv[2]` | `yes` to load pre-computed bad anchors from the repo's `.tsv` data files; `no` to sample them live from the full dataset |

**Key parameters (hardcoded at top of script):**

| Variable | Default | Meaning |
|---|---|---|
| `exclude` | `1` | Half-width of the "near-zero" exclusion square (anchors within ±1 on both axes are filtered out before sampling) |
| `shift` | `1` | Half-width of the diagonal band that defines "between the lines" (good prediction) |
| `N` | `5` | Number of anchors sampled per group per genome |

**Reads:**
- `/mnt/altnas/work/ishawnia/empvsag/full-CTCF-HG2-15.tsv`
- `/mnt/altnas/work/ishawnia/empvsag/full-CTCF-GM12878-15.tsv`
- `GM12878_5-0.1_outside.tsv` *(only when `have_bad=yes`)*
- `HepG2_5-0.1_outside.tsv` *(only when `have_bad=yes`)*

**Outputs (written to `<folder>`):**

| File | Description |
|---|---|
| `motif_anchor_analysis_selected_anchors.tsv` | All sampled anchors (good + bad, both genomes) with full metadata |
| `gm12878_selected_anchors_variants.bed` | BED file of selected GM12878 anchors |
| `hepg2_selected_anchors_variants.bed` | BED file of selected HepG2 anchors |

---

### `SCR_intersect_motifs_with_anchors.sh`

**Intersects selected anchor BED files with CTCF motif variant files** using `bedtools intersect`.

```bash
bash SCR_intersect_motifs_with_anchors.sh <gm_bed> <hg2_bed> <output_folder>
```

| Positional arg | Description |
|---|---|
| `$1` | GM12878 anchor BED file |
| `$2` | HepG2 anchor BED file |
| `$3` | Output folder |

**Motif files used (hardcoded, on GPFS):**
- GM12878: `FULL_LCL_ALPHAGENOME_RUNS_10.31.25/GM12878/combined/CTCF_MA0139.1.bed`
- HepG2: `HepG2/diploidGenomeOutput/combined/CTCF_MA0139.1.bed`

**Outputs (written to `<folder>`):**

| File | Description |
|---|---|
| `GM_anchors_with_motifs.tsv` | Intersection result for GM12878 |
| `HG2_anchors_with_motifs.tsv` | Intersection result for HepG2 |

Both files share the same columns:
`CHR, POS1, POS2, motif_chr, motif_start, motif_end, motif_id, snp_pos, ref, alt, genotype, variant_type, ref_score, alt_score, delta_score, ref_seq, alt_seq`

---

### `SCR_analyze-snp-motif-anchors.py`

**Classifies the effect of each SNP on its overlapping CTCF motif** and summarises per anchor.

```bash
python SCR_analyze-snp-motif-anchors.py <TYPE> <folder>
```

| Positional arg | Description |
|---|---|
| `sys.argv[1]` | Genome identifier: `GM` or `HG2` |
| `sys.argv[2]` | Folder containing the intersect TSV and to write output into |

**Motif effect classification:**

| `motif_effect` value | Condition |
|---|---|
| `0` | `abs(delta_score)` < 1 — negligible effect |
| `1` | 1 ≤ `abs(delta_score)` < 5 — moderate effect |
| `2` | `abs(delta_score)` ≥ 5 — strong effect |

Per-anchor summary uses the **strongest** effect seen across all overlapping SNPs.

**Outputs (written to `<folder>`):**

| File | Description |
|---|---|
| `{TYPE}_JUST-snp-analysis.tsv` | One row per anchor: `CHR, POS1, POS2, n_motif_snps, max_abs_delta, strongest_effect, anchor_summary, GENOME` |
| `{TYPE}_snp_motif_details.tsv` | One row per anchor–SNP–motif combination, with all scores and the per-anchor summary columns merged in |

---

### `SCR_plot_examples.py`

**Generates diagnostic plots** comparing SNP motif effect to anchor prediction quality.

```bash
python SCR_plot_examples.py <folder>
```

| Positional arg | Description |
|---|---|
| `sys.argv[1]` | Folder containing all upstream outputs |

Runs independently for each genome (`GM`, `HG2`), producing a 5-panel figure:

| Panel | Content |
|---|---|
| 0 | Scatter: empirical vs predicted log₂ fold-change, colored by variant |
| 1 | Boxplot: predicted log₂ diff grouped by strongest SNP effect category |
| 2 | Scatter: predicted log₂ diff vs strongest SNP effect |
| 3 | Confusion heatmap: SNP effect vs `btwn_lines` (all rows) |
| 4 | Confusion heatmap: SNP effect vs `btwn_lines` (unique anchors only) |

Also produces a combined heatmap across both genomes.

**Outputs (written to `<folder>`):**

| File | Description |
|---|---|
| `GM_motif_anchor_snp_effects_summary.png` | 5-panel figure for GM12878 |
| `HG2_motif_anchor_snp_effects_summary.png` | 5-panel figure for HepG2 |
| `btwn_lines_vs_strongest_effect.png` | Combined heatmap (both genomes) |

---

## Helper/utility files (`FUNC_`)

`FUNC_` files are either reusable utility modules called by the pipeline or standalone one-off helpers. They are **not** meant to be run as top-level entry points except where noted.

---

### `FUNC_make-file-with-details.py`

**Joins the two per-genome SNP motif detail files with the selected anchor metadata** into a single combined TSV. Called directly by `SCR_sample_vs_motif-strength.sh`.

```bash
python FUNC_make-file-with-details.py <hg2_details> <gm_details> <anchor_details> <output_file>
```

| Positional arg | Description |
|---|---|
| `sys.argv[1]` | `HG2_snp_motif_details.tsv` |
| `sys.argv[2]` | `GM_snp_motif_details.tsv` |
| `sys.argv[3]` | `motif_anchor_analysis_selected_anchors.tsv` |
| `sys.argv[4]` | Output path (default: `combined2.tsv`) |

Both detail files are inner-joined to the anchor metadata on `CHR, POS1, POS2`, then concatenated.

**Output:**

| File | Description |
|---|---|
| `combined_details.tsv` | Combined HepG2 + GM12878 SNP motif details with anchor metadata |

---

### `FUNC_get_variants.py`

**One-off utility** — extracts variant identifiers from a `GM_anchors_with_motifs.tsv` file and writes a flat variant table. Paths are **hardcoded** and must be edited before use.

**Outputs:**
- `GM_variants.tsv` — columns: `variant_id, CHROM, POS, REF, ALT`

---

## Included data files

| File | Description |
|---|---|
| `GM12878_5-0.1_outside.tsv` | Pre-computed set of "bad" (poorly-predicted, outside square) GM12878 anchors, used when `have_bad=yes` |
| `HepG2_5-0.1_outside.tsv` | Same for HepG2 |

These files let the pipeline reproduce a fixed set of "bad" examples without re-running the full upstream comparison.

---

## Output folder structure

After a complete run, `<folder>/` will contain:

```
<folder>/
├── motif_anchor_analysis_selected_anchors.tsv   # selected anchors (all groups)
├── gm12878_selected_anchors_variants.bed         # GM12878 anchor BED
├── hepg2_selected_anchors_variants.bed           # HepG2 anchor BED
├── GM_anchors_with_motifs.tsv                    # bedtools intersect output – GM
├── HG2_anchors_with_motifs.tsv                   # bedtools intersect output – HepG2
├── GM_JUST-snp-analysis.tsv                      # per-anchor SNP summary – GM
├── HG2_JUST-snp-analysis.tsv                     # per-anchor SNP summary – HepG2
├── GM_snp_motif_details.tsv                      # per-SNP detail rows – GM
├── HG2_snp_motif_details.tsv                     # per-SNP detail rows – HepG2
├── combined_details.tsv                          # merged across both genomes
├── GM_motif_anchor_snp_effects_summary.png       # 5-panel plot – GM
├── HG2_motif_anchor_snp_effects_summary.png      # 5-panel plot – HepG2
└── btwn_lines_vs_strongest_effect.png            # combined heatmap
```

---

## Dependencies

- Python ≥ 3.8
- `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`
- `bedtools` (must be on `PATH`)
- Access to GPFS motif files and `/mnt/altnas` data paths (cluster environment)
