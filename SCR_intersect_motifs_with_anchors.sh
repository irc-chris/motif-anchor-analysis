folder=$3

echo -e "CHR\tPOS1\tPOS2\tmotif_chr\tmotif_start\tmotif_end\tmotif_id\tsnp_pos\tref\talt\tgenotype\tvariant_type\tref_score\talt_score\tdelta_score\tref_seq\talt_seq" > $folder/2-GM_anchors_with_motifs.tsv
echo -e "CHR\tPOS1\tPOS2\tmotif_chr\tmotif_start\tmotif_end\tmotif_id\tsnp_pos\tref\talt\tgenotype\tvariant_type\tref_score\talt_score\tdelta_score\tref_seq\talt_seq" > $folder/2-HG2_anchors_with_motifs.tsv

bedtools intersect \
  -a $1 \
  -b /gpfs0/work/suhas/temp_processing/p-e_paper_2022/figure2/diploid_motifs/FULL_LCL_ALPHAGENOME_RUNS_10.31.25/GM12878/combined/CTCF_MA0139.1.bed \
  -wa -wb \
  >> $folder/2-GM_anchors_with_motifs.tsv


bedtools intersect \
  -a $2 \
  -b /gpfs0/work/suhas/temp_processing/p-e_paper_2022/figure2/diploid_motifs/HepG2/diploidGenomeOutput/combined/CTCF_MA0139.1.bed \
  -wa -wb \
  >> $folder/2-HG2_anchors_with_motifs.tsv