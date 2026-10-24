input=~/Treeshrew/FileSave/bed_by_celltype_for_traning/bed_after_sort/C8.1.sort.bed
C8_1_pr1_bed=~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C8Epithelial_cells/data/cut_site_bigwig/C8.1.pr1.bed
C8_1_pr2_bed=~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C8Epithelial_cells/data/cut_site_bigwig/C8.1.pr2.bed
C8_1_pr1_sort_bed=~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C8Epithelial_cells/data/cut_site_bigwig/C8.1.pr1.sort.bed
C8_1_pr2_sort_bed=~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C8Epithelial_cells/data/cut_site_bigwig/C8.1.pr2.sort.bed

awk -v pr1="$C8_1_pr1_bed" \
    -v pr2="$C8_1_pr2_bed" \
    'BEGIN{srand(123)}
     {
        if (rand() < 0.5) print > pr1;
        else print > pr2;
     }' "$input"

sort -k1,1 -k2,2n $C8_1_pr1_bed > $C8_1_pr1_sort_bed
sort -k1,1 -k2,2n $C8_1_pr2_bed > $C8_1_pr2_sort_bed