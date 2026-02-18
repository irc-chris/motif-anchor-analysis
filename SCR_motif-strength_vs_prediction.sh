folder="$1"
motif_threshold=8
chip_diff_threshold=3
empirical_threshold=15
bad_ag_threshold=0.1
good_ag_threshold=1.0

folder="${folder}-c${empirical_threshold}-cd${chip_diff_threshold}-pb${bad_ag_threshold}-pg${good_ag_threshold}-m${motif_threshold}"
mkdir -p "$folder"

Experiment= # DKV or AKB or DRZ

# python SCR_from-colab.py $1 yes no no $Experiment $chip_diff_threshold
python SCR_from-claud.py $folder $chip_diff_threshold $bad_ag_threshold $good_ag_threshold
./SCR_intersect_motifs_with_anchors.sh  $folder/1-gm12878_selected_anchors_variants.bed $folder/1-hepg2_selected_anchors_variants.bed $folder
python SCR_analyze-snp-motif-anchors.py GM $folder
python SCR_analyze-snp-motif-anchors.py HG2 $folder

python SCR_plot_examples.py $folder $chip_diff_threshold $empirical_threshold $motif_threshold $bad_ag_threshold $good_ag_threshold $Experiment 
# python FUNC_make-file-with-details.py $folder/3-HG2_snp_motif_details.tsv $folder/3-GM_snp_motif_details.tsv $folder/1-motif_anchor_analysis_selected_anchors.tsv $folder/4-combined_details.tsv
awk -v threshold=$motif_threshold '$29 > threshold || $30 > threshold' $folder/6-GM_combined_anchor_snp_data.tsv > $folder/7-GM-combined_details_motifs$motif_threshold.tsv
awk -v threshold=$motif_threshold '$29 > threshold || $30 > threshold' $folder/6-HG2_combined_anchor_snp_data.tsv > $folder/7-HG2-combined_details_motifs$motif_threshold.tsv

python SCR_plot_ref_reads.py $folder/6-combined_anchor_snp_data.tsv

python ../prep_scripts/reorder-cols.py "$folder/6-GM_combined_anchor_snp_data.tsv"
python ../prep_scripts/reorder-cols.py "$folder/6-HG2_combined_anchor_snp_data.tsv"
python ../prep_scripts/reorder-cols.py "$folder/6-combined_anchor_snp_data.tsv"

for file in $folder/*.tsv; do
    python ../prep_scripts/sorting.py "$file" "$file"
done