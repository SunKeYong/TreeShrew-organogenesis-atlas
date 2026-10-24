dir="~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/data/region_for_predict/vista/heart/hs_mm_ht/filter_common_name"
region_bed="${dir}/heart_tissue_VISTA_HsToTS_05.filter.bed"
chrom_sizes="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.chrom.sizes"
output_bed="${dir}/heart_tissue_VISTA_HsToTS_05.modified.bed"


#-i srcwinnum将原始区间的名称（第 4 列）与该窗口在原始区间内的顺序编号结合起来!!!!!!!!!!!!!!!!!!
bedtools makewindows -b $region_bed -w 1000 -s 750 -i srcwinnum | \
awk -v OFS="\t" '{print $1,$2,$3,$4,".",".",".",".",".",int(($3-$2)/2)}' > $output_bed
