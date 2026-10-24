TS_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/heart_tissue_VISTA_HsToTS_05.bed"
TS_modified_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/heart_tissue_VISTA_HsToTS_05.modified.bed"
name_list="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/heart_tissue_VISTA_common_name.txt"
filter_bed="~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/data/region_for_predict/vista/heart/hs_mm_ht/filter_common_name/heart_tissue_VISTA_HsToTS_05.filter.bed"

awk 'BEGIN {
    map["CM065820.1"] = "chr1"
    map["CM065821.1"] = "chr2"
    map["CM065822.1"] = "chr3"
    map["CM065823.1"] = "chr4"
    map["CM065824.1"] = "chr5"
    map["CM065825.1"] = "chr6"
    map["CM065826.1"] = "chr7"
    map["CM065827.1"] = "chr8"
    map["CM065828.1"] = "chr9"
    map["CM065829.1"] = "chr10"
    map["CM065830.1"] = "chr11"
    map["CM065831.1"] = "chr12"
    map["CM065832.1"] = "chr13"
    map["CM065833.1"] = "chr14"
    map["CM065834.1"] = "chr15"
    map["CM065835.1"] = "chr16"
    map["CM065836.1"] = "chr17"
    map["CM065837.1"] = "chr18"
    map["CM065838.1"] = "chr19"
    map["CM065839.1"] = "chr20"
    map["CM065840.1"] = "chr21"
    map["CM065841.1"] = "chr22"
    map["CM065842.1"] = "chr23"
    map["CM065843.1"] = "chr24"
    map["CM065844.1"] = "chr25"
    map["CM065845.1"] = "chr26"
    map["CM065846.1"] = "chr27"
    map["CM065847.1"] = "chr28"
    map["CM065848.1"] = "chr29"
    map["CM065849.1"] = "chr30"
    map["CM065850.1"] = "chrX"
    OFS="\t"  # 设置输出字段分隔符为制表符
}
{
    if($1 in map) {
        $1 = map[$1]
        print
    } else {
        # 如果染色体不在映射中，输出到单独文件以便检查
        print $0 > "/dev/stderr"
    }
}' $TS_bed > $TS_modified_bed 2> ~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/unmapped_chromosomes.txt


#根据第四列（名称）与给定列表的匹配来过滤BED文件
awk 'NR==FNR {a[$0]; next} $4 in a' $name_list $TS_modified_bed > $filter_bed

