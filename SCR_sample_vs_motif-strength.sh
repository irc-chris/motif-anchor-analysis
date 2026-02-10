folder="$1"
mkdir -p "$folder"
motif_threshold=5

python SCR_from-colab.py $1 yes no no
./SCR_intersect_motifs_with_anchors.sh  $folder/1-gm12878_selected_anchors_variants.bed $folder/1-hepg2_selected_anchors_variants.bed $folder
python SCR_analyze-snp-motif-anchors.py GM $folder
python SCR_analyze-snp-motif-anchors.py HG2 $folder
python SCR_plot_examples.py $folder
python FUNC_make-file-with-details.py $folder/3-HG2_snp_motif_details.tsv $folder/3-GM_snp_motif_details.tsv $folder/1-motif_anchor_analysis_selected_anchors.tsv $folder/4-combined_details.tsv
awk -v threshold=$motif_threshold '$29 > threshold || $30 > threshold' $folder/4-combined_details_GM.tsv > $folder/5-combined_details_GM-motifs$motif_threshold.tsv
awk -v threshold=$motif_threshold '$29 > threshold || $30 > threshold' $folder/4-combined_details_HG2.tsv > $folder/5-combined_details_HG2-motifs$motif_threshold.tsv
    
for file in $folder/*.tsv; do
    python ../prep_scripts/sorting.py "$file" "$file"
done