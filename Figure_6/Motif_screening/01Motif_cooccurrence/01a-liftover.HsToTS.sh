input_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/heart_tissue_VISTA.bed"
chain_file="~/software/chain_files/treeshrew/hg38ToGCA_033439345.1.over.chain.gz"
output_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/heart_tissue_VISTA_HsToTS.bed"
unmapped_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/data/region_for_predict/vista/heart/heart_tissue_VISTA_HsToTS.unmapped.bed"


#minMatch默认为0.95.参考文献改为0.5捕获更多的同源序列
~/software/liftOver_files/liftOver $input_bed $chain_file $output_bed $unmapped_bed -minMatch=0.5
