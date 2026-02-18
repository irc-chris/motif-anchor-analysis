#!/usr/bin/env python3
import pandas as pd
import sys

# ---------------------------- 
# Input / output 
# ---------------------------- 
TYPE = sys.argv[1]  # "GM" or "HG2"
folder = sys.argv[2]

INPUT = f"{folder}/2-{TYPE}_anchors_with_motifs.tsv"
OUTPUT = f"{folder}/3-{TYPE}_JUST-snp-analysis.tsv"
DETAIL_OUTPUT = f"{folder}/3-{TYPE}_snp_motif_details.tsv"

# Column names for input file
COLUMNS = [
    "CHR", "POS1", "POS2",
    "motif_chr", "motif_start", "motif_end", "motif_id",
    "snp_pos", "ref", "alt", "phase", "better_strand",
    "delta_score", "h1_score", "h2_score", "ref_seq", "alt_seq"
]

# Columns to round to 3 decimals
SCORE_COLUMNS = ['delta_score', 'h1_score', 'h2_score', 'abs_delta']    

def classify_motif_effect(abs_delta):
    """Classify motif effect based on absolute delta score."""
    if pd.isna(abs_delta):
        return "No SNP in motif"
    elif abs_delta < 1:
        return 0
    elif abs_delta < 5:
        return 1
    else:
        return 2


def anchor_sentence(row):
    """Generate description for anchor based on SNP effects."""
    if row["n_motif_snps"] == 0:
        return "Does not overlap a CTCF motif SNP."
    elif row["strongest_effect"] == 2:
        return "SNP that strongly alters CTCF motif strength."
    elif row["strongest_effect"] == 1:
        return "SNP that moderately affects CTCF motif strength."
    else:
        return "SNP with no meaningful effect on CTCF motif strength."


# ---------------------------- 
# Load and prepare data
# ---------------------------- 
df = pd.read_csv(INPUT, sep="\t", header=0)
df.columns = COLUMNS

# Convert scores to float and calculate absolute delta
df['delta_score'] = df['delta_score'].replace('.', 0).astype(float)
df['h1_score'] = df['h1_score'].replace('.', 0).astype(float)
df['h2_score'] = df['h2_score'].replace('.', 0).astype(float)
df["abs_delta"] = df["delta_score"].abs()

# Classify motif effects
df["motif_effect"] = df["abs_delta"].apply(classify_motif_effect)

# ---------------------------- 
# Per-anchor summary
# ---------------------------- 
anchor_summary = (
    df.groupby(["CHR", "POS1", "POS2"])
    .agg(
        n_motif_snps=("delta_score", "count"),
        max_abs_delta=("abs_delta", "max"),
        strongest_effect=("motif_effect", lambda x: (
            2 if 2 in x.values else 1 if 1 in x.values else 0
        ))
    )
    .reset_index()
)

# Add summary descriptions
anchor_summary["anchor_summary"] = anchor_summary.apply(anchor_sentence, axis=1)
anchor_summary["GENOME"] = TYPE

# ---------------------------- 
# Merge summary back to details
# ---------------------------- 
df = df.merge(
    anchor_summary[["CHR", "POS1", "POS2", "n_motif_snps", "strongest_effect", "anchor_summary"]],
    on=["CHR", "POS1", "POS2"],
    how="left"
)
df["GENOME"] = TYPE

# Round score columns to 3 decimals
for col in SCORE_COLUMNS:
    if col in df.columns:
        df[col] = df[col].replace('.', 0).astype(float).round(3)

# ---------------------------- 
# Sort data
# ---------------------------- 
anchor_summary = anchor_summary.sort_values(by=['CHR', 'POS1', 'POS2'])
df = df.sort_values(by=['CHR', 'POS1', 'POS2', 'motif_chr', 'motif_start', 'motif_end', 'snp_pos'])

# ---------------------------- 
# Save outputs
# ---------------------------- 
anchor_summary.to_csv(OUTPUT, sep="\t", index=False)
df.to_csv(DETAIL_OUTPUT, sep="\t", index=False)

print(f"Saved anchor summary to {OUTPUT}")
print(f"Saved SNP motif details to {DETAIL_OUTPUT}")