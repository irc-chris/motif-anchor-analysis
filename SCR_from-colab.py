import pandas as pd
import sklearn.cluster as sc
import matplotlib.pyplot as plt
import numpy as np
import sys

### Seeing area between lines.
exclude = 1
shift = 1
N = 9
folder = sys.argv[1]  # folder to read / write files
have_bad = sys.argv[2]
sample_bad = sys.argv[3]
sample_good = sys.argv[4]
GM_EXP = sys.argv[5] if len(sys.argv) > 5 else None
print(f"Arguments: folder={folder}, have_bad={have_bad}, sample_bad={sample_bad}, sample_good={sample_good}, GM_EXP={GM_EXP}")

def get_okay(df):
  mask = (
      (df["DIFF_LOG2_data2"] >= df["DIFF_LOG2_data1"] - shift) &
      (df["DIFF_LOG2_data2"] <= df["DIFF_LOG2_data1"] + shift)
  )

  df.loc[:, "btwn_lines"] = mask.astype(int)
  count = df[df['btwn_lines']==1].count()['btwn_lines']
  return df, count

def get_outside_square(df):
  mask = (
      (df["DIFF_LOG2_data1"] >= -exclude) & (df["DIFF_LOG2_data1"] <= exclude) &
      (df["DIFF_LOG2_data2"] >= -exclude) & (df["DIFF_LOG2_data2"] <= exclude)
  )
  return df[~mask]

def analysis(df, i):
  df1, count1 = get_okay(df)
  df2 = get_outside_square(df)
  df2, count2 = get_okay(df2)

  percentage = 100 * count2 / df2.shape[0]
  axes[i].scatter(df2['DIFF_LOG2_data1'], df2['DIFF_LOG2_data2'], c=df2['btwn_lines'])
  axes[i].set_title(f'CTCF Experiment {names[i]} ChipSeq Log2Diff vs AG\n n={count2}/{df2.shape[0]},  {percentage:.3f}%')
  return df2

dfHg2 = pd.read_csv('/mnt/altnas/work/ishawnia/empvsag/full-CTCF-HG2-15.tsv', sep='\t')
dfGM = pd.read_csv('/mnt/altnas/work/ishawnia/empvsag/full-CTCF-GM12878-15.tsv', sep='\t')

# Add 'along-x' field to dfHg2
dfHg2.loc[:, 'along-x'] = ((dfHg2['DIFF_LOG2_data2'] >= -1) & (dfHg2['DIFF_LOG2_data2'] <= 1)).astype(int)

# Add 'along-x' field to dfGM
dfGM.loc[:, 'along-x'] = ((dfGM['DIFF_LOG2_data2'] >= -1) & (dfGM['DIFF_LOG2_data2'] <= 1)).astype(int)

dfHg2_outside = get_outside_square(dfHg2)
dfHg2_outside['GENOME'] = 'HepG2'
dfGM_outside = get_outside_square(dfGM)
dfGM_outside['GENOME'] = 'GM12878'

####### TAKE ONLY ONE GM Experiment

if have_bad == 'yes':
    dfGM_bad = pd.read_csv('/mnt/altnas/work/ishawnia/empvsag/motif-anchor-analysis/GM12878_5-0.1_outside.tsv', sep='\t')
    dfHG2_bad = pd.read_csv('/mnt/altnas/work/ishawnia/empvsag/motif-anchor-analysis/HepG2_5-0.1_outside.tsv', sep='\t')
    if GM_EXP:
      dfGM_bad = dfGM_bad[dfGM_bad['transcription_factor_data1'].str.endswith(GM_EXP)]
    if sample_bad == 'yes':
      hepg2_random_bad = dfHG2_bad.sample(n=N, replace=False)
      gm_random_bad = dfGM_bad.sample(n=N, replace=False)
    else:
      hepg2_random_bad = dfHG2_bad
      gm_random_bad = dfGM_bad
    hepg2_random_bad['GENOME'] = 'HepG2'
    gm_random_bad['GENOME'] = 'GM12878'
    hepg2_random_bad['along-x'] = 1
    gm_random_bad['along-x'] = 1
else:
  h_outside_bad = dfHg2_outside[(dfHg2_outside['along-x'] == 1) & (abs(dfHg2_outside['DIFF_LOG2_data1']) > 3)]
  g_outside_bad = dfGM_outside[(dfGM_outside['along-x'] == 1) & (abs(dfGM_outside['DIFF_LOG2_data1']) > 3)]
  if GM_EXP:
      g_outside_bad = g_outside_bad[g_outside_bad['transcription_factor_data1'].str.endswith(GM_EXP)]
  if sample_bad == 'yes':
    hepg2_random_bad = h_outside_bad.sample(n=N, replace=False)
    gm_random_bad = g_outside_bad.sample(n=N, replace=False)
  else:
    hepg2_random_bad = h_outside_bad
    gm_random_bad = g_outside_bad

h_outside_good = dfHg2_outside[(dfHg2_outside['btwn_lines'] == 1) & (abs(dfHg2_outside['DIFF_LOG2_data1']) > 3)]
g_outside_good = dfGM_outside[(dfGM_outside['btwn_lines'] == 1) & (abs(dfGM_outside['DIFF_LOG2_data1']) > 3)]
if GM_EXP:
      g_outside_good = g_outside_good[g_outside_good['transcription_factor_data1'].str.endswith(GM_EXP)]
  
if sample_good == 'yes':
  hepg2_random_good = h_outside_good.sample(n=N, replace=False)
  gm_random_good = g_outside_good.sample(n=N, replace=False)
else:
  hepg2_random_good = h_outside_good
  gm_random_good = g_outside_good

print("\n#####")
print(hepg2_random_bad)
print("\n#####")
print(hepg2_random_good)
print("\n#####")
print(gm_random_bad)
print("\n#####")
print(gm_random_good)

combined_anchors = pd.concat([hepg2_random_bad, hepg2_random_good, gm_random_bad, gm_random_good])

just_anchor_variants = combined_anchors[['CHR', 'POS1', 'POS2', 'GENOME']]
hepg2_variants = just_anchor_variants[just_anchor_variants['GENOME'] == 'HepG2'].drop(columns=['GENOME']).drop_duplicates()
hepg2_variants.to_csv(f'{folder}/1-hepg2_selected_anchors_variants.bed', sep='\t', index=False, header=['#chrom', 'start', 'end'])
gm_variants = just_anchor_variants[just_anchor_variants['GENOME'] == 'GM12878'].drop(columns=['GENOME']).drop_duplicates()
gm_variants.to_csv(f'{folder}/1-gm12878_selected_anchors_variants.bed', sep='\t', index=False, header=['#chrom', 'start', 'end'])

combined_anchors.to_csv(f'{folder}/1-motif_anchor_analysis_selected_anchors.tsv', sep='\t', index=False)
