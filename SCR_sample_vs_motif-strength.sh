folder="$1"
mkdir -p "$folder"

python SCR_from-colab.py $1 yes
./SCR_intersect_motifs_with_anchors.sh  $folder/gm12878_selected_anchors_variants.bed $folder/hepg2_selected_anchors_variants.bed $folder
python SCR_analyze-snp-motif-anchors.py GM $folder
python SCR_analyze-snp-motif-anchors.py HG2 $folder
python SCR_plot_examples.py $folder
python FUNC_make-file-with-details.py $folder/HG2_snp_motif_details.tsv $folder/GM_snp_motif_details.tsv $folder/motif_anchor_analysis_selected_anchors.tsv $folder/combined_details.tsv

for file in $folder/*.tsv; do
    python ../prep_scripts/sorting.py "$file" "$file"
done