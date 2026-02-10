#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

# ----------------------------
# Input / output
# ----------------------------
TYPE = sys.argv[1]  # "GM" or "HG2"
folder = sys.argv[2]  # folder to read / write files
INPUT = f"{folder}/2-{TYPE}_anchors_with_motifs.tsv"
OUTPUT = f"{folder}/3-{TYPE}_JUST-snp-analysis.tsv"

# ----------------------------
# Load intersect output
# ----------------------------
df = pd.read_csv(INPUT, sep="\t", header=0)

df.columns = [
    "CHR", "POS1", "POS2",
    "motif_chr", "motif_start", "motif_end",
    "motif_id",
    "snp_pos",
    "ref", "alt", "genotype", "variant_type",
    "ref_score", "alt_score", "delta_score",
    "ref_seq", "alt_seq"
]

# ----------------------------
# Absolute delta score
# ----------------------------
df['delta_score'] = (
    df['delta_score']
    .replace('.', 0)
    .astype(float)
)

df["abs_delta"] = df["delta_score"].abs()

# ----------------------------
# Motif effect classification (STEP 3 + 4)
# ----------------------------
def classify_motif_effect(row):
    if pd.isna(row["delta_score"]):
        return "No SNP in motif"
    elif row["abs_delta"] < 1:
        return 0
    elif row["abs_delta"] < 5:
        return 1
    else:
        return 2

df["motif_effect"] = df.apply(classify_motif_effect, axis=1)

# ----------------------------
# Per-anchor summary
# ----------------------------
anchor_summary = (
    df
    .groupby(["CHR", "POS1", "POS2"])
    .agg(
        n_motif_snps=("delta_score", "count"),
        max_abs_delta=("abs_delta", "max"),
        strongest_effect=("motif_effect", lambda x: (
            2 if 2 in x.values
            else 1 if 1 in x.values
            else 0
        ))
    )
    .reset_index()
)

# ----------------------------
# One-sentence description per anchor
# ----------------------------
def anchor_sentence(row):
    if row["n_motif_snps"] == 0:
        return "Does not overlap a CTCF motif SNP."
    elif row["strongest_effect"] == 2:
        return "SNP that strongly alters CTCF motif strength."
    elif row["strongest_effect"] == 1:
        return "SNP that moderately affects CTCF motif strength."
    else:
        return "SNP with no meaningful effect on CTCF motif strength."

anchor_summary["anchor_summary"] = anchor_summary.apply(anchor_sentence, axis=1)

print(df.shape, 'before')
print(df)
print(anchor_summary)

# Merge n_motif_snps back to df
df = df.merge(
    anchor_summary[["CHR", "POS1", "POS2", "n_motif_snps", "strongest_effect", "anchor_summary"]], 
    on=["CHR", "POS1", "POS2"], 
    how="left"
)

print(df.shape, 'after')
print(df)

df["anchor_summary"] = df.apply(anchor_sentence, axis=1)
df["GENOME"] = TYPE
df['alt_score'] = df['alt_score'].replace('.', 0).astype(float).round(3)
df['ref_score'] = df['ref_score'].replace('.', 0).astype(float).round(3)
df['delta_score'] = df['delta_score'].replace('.', 0).astype(float).round(3)
df['abs_delta'] = df['abs_delta'].replace('.', 0).astype(float).round(3)

### sort by column {CHR, POS1, POS2}
anchor_summary = anchor_summary.sort_values(by=['CHR', 'POS1', 'POS2'])
df = df.sort_values(by=['CHR', 'POS1', 'POS2', 'motif_chr', 'motif_start', 'motif_end', 'snp_pos'])

anchor_summary["GENOME"] = TYPE

# ----------------------------
# Save output
# ----------------------------
anchor_summary.to_csv(OUTPUT, sep="\t", index=False, header=True)
df.to_csv(f"{folder}/3-{TYPE}_snp_motif_details.tsv", sep="\t", index=False, header=True)

print(f"Saved SNP motif analysis to {OUTPUT}")
