input_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/heart_tissue_VISTA_HsToTS_05.bed"
chain_file="~/software/chain_files/treeshrew/GCA_033439345.1ToHg38.over.chain.gz"
output_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/heart_tissue_VISTA_TSbackHs_05.bed"
unmapped_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/heart_tissue_VISTA_TSbackHs_05.unmapped.bed"


#minMatch默认为0.95.参考文献改为0.5捕获更多的同源序列
~/software/liftOver_files/liftOver $input_bed $chain_file $output_bed $unmapped_bed -minMatch=0.5


# 提取第四列的所有值
awk '{print $4}' ~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/heart_tissue_VISTA_TSbackHs_05.bed | sort | uniq > ~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/liftover/heart_tissue_VISTA_common_name.txt
