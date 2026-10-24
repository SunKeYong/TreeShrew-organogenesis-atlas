#
setwd("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step3_peak_plot")
#
library(ArchR)
#
load("../step2_ArchR/human_ArchR/Hs_markerPeaks.RData")
load("../step2_ArchR/Treeshrew_ArchR/TS_markerPeaks.RData")
#

heatmapPeaks <- markerHeatmap(
  seMarker = Hs_markerPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#



#
Hs_markerPeaks@assays@data$Mean
Hs_markerPeaks@NAMES
#
Hs_markerList <- getMarkers(Hs_markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
TS_markerList <- getMarkers(TS_markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
#
Hs_markerList$B_cells
#
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(rownames(as.data.frame(Hs_markerList$B_cells))[1:1000]),],
                   color = paletteContinuous(),
                   scale = "row",
                   cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[Hs_markerList$B_cells$idx,],
                   color = paletteContinuous(),
                   scale = "row",
                   cluster_rows = F,cluster_cols = F)
#
#                                                     
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(rownames(as.data.frame(Hs_markerList$B_cells))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$CLP))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$Erythroid_cells))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$HSPCs))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$HSPCs_MPO))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$Kupffer_cells))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$Macrophages))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$Mast_cells))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$Megakaryocytes))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$Microglia))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$Monocytes))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$NK_cells))[1:100],
                                                     rownames(as.data.frame(Hs_markerList$T_cells))[1:100]),],
                   color = paletteContinuous(),
                   scale = "row",
                   cluster_rows = F,cluster_cols = F)
#
#
#
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(rownames(as.data.frame(TS_markerList$B_cells))[1:100],
                                                     rownames(as.data.frame(TS_markerList$CLP))[1:100],
                                                     rownames(as.data.frame(TS_markerList$Erythroid_cells))[1:100],
                                                     rownames(as.data.frame(TS_markerList$HSPCs))[1:100],
                                                     rownames(as.data.frame(TS_markerList$HSPCs_MPO))[1:100],
                                                     rownames(as.data.frame(TS_markerList$Kupffer_cells))[1:100],
                                                     rownames(as.data.frame(TS_markerList$Macrophages))[1:100],
                                                     rownames(as.data.frame(TS_markerList$Mast_cells))[1:100],
                                                     rownames(as.data.frame(TS_markerList$Megakaryocytes))[1:100],
                                                     rownames(as.data.frame(TS_markerList$Microglia))[1:100],
                                                     rownames(as.data.frame(TS_markerList$Monocytes))[1:100],
                                                     rownames(as.data.frame(TS_markerList$NK_cells))[1:100],
                                                     rownames(as.data.frame(TS_markerList$T_cells))[1:100]),],
                   color = paletteContinuous(),
                   scale = "row",
                   cluster_rows = F,cluster_cols = F)
#
#
#
#
#
######### ######### ######### #
######### ######### ######### #
######### ######### ######### #
#
TS_markerList$B_cells
#
TS_markerList$B_cells$idx <- rownames(TS_markerList$B_cells)
TS_markerList$CLP$idx <- rownames(TS_markerList$CLP)
TS_markerList$Erythroid_cells$idx <- rownames(TS_markerList$Erythroid_cells)
TS_markerList$HSPCs$idx <- rownames(TS_markerList$HSPCs)
TS_markerList$HSPCs_MPO$idx <- rownames(TS_markerList$HSPCs_MPO)
TS_markerList$Kupffer_cells$idx <- rownames(TS_markerList$Kupffer_cells)
TS_markerList$Macrophages$idx <- rownames(TS_markerList$Macrophages)
TS_markerList$Mast_cells$idx <- rownames(TS_markerList$Mast_cells)
TS_markerList$Megakaryocytes$idx <- rownames(TS_markerList$Megakaryocytes)
TS_markerList$Microglia$idx <- rownames(TS_markerList$Microglia)
TS_markerList$Monocytes$idx <- rownames(TS_markerList$Monocytes)
TS_markerList$NK_cells$idx <- rownames(TS_markerList$NK_cells)
TS_markerList$T_cells$idx <- rownames(TS_markerList$T_cells)
#
#
Hs_markerList$B_cells$idx <- rownames(Hs_markerList$B_cells)
Hs_markerList$CLP$idx <- rownames(Hs_markerList$CLP)
Hs_markerList$Erythroid_cells$idx <- rownames(Hs_markerList$Erythroid_cells)
Hs_markerList$HSPCs$idx <- rownames(Hs_markerList$HSPCs)
Hs_markerList$HSPCs_MPO$idx <- rownames(Hs_markerList$HSPCs_MPO)
Hs_markerList$Kupffer_cells$idx <- rownames(Hs_markerList$Kupffer_cells)
Hs_markerList$Macrophages$idx <- rownames(Hs_markerList$Macrophages)
Hs_markerList$Mast_cells$idx <- rownames(Hs_markerList$Mast_cells)
Hs_markerList$Megakaryocytes$idx <- rownames(Hs_markerList$Megakaryocytes)
Hs_markerList$Microglia$idx <- rownames(Hs_markerList$Microglia)
Hs_markerList$Monocytes$idx <- rownames(Hs_markerList$Monocytes)
Hs_markerList$NK_cells$idx <- rownames(Hs_markerList$NK_cells)
Hs_markerList$T_cells$idx <- rownames(Hs_markerList$T_cells)
#
#
#
#
#
#
#
#
#
hs_orth_peak_df <- read.table("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/orth.human.raw.sort.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(hs_orth_peak_df) <- c("chr", "start", "end", "peak_id")
#
ts_orth_peak_df <- read.table("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/orth.ts.raw.sort.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(ts_orth_peak_df) <- c("chr", "start", "end", "peak_id")
#
head(hs_orth_peak_df)
head(ts_orth_peak_df)
#
dim(hs_orth_peak_df)
dim(ts_orth_peak_df)
#
hs_orth_peak_df$hs_orth_idx <- as.character(1:108744)
ts_orth_peak_df$ts_orth_idx <- as.character(1:108744)
#
head(hs_orth_peak_df)
head(ts_orth_peak_df)
#
# Hs_markerPeaks
# Hs_markerList
# TS_markerPeakssx
# TS_markerList
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$B_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$B_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);
ts_sp_tmp <- setdiff(ts_tmp,co_tmp);
hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_B_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_B_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_B_cells <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_B_cells <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_B_cells <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_B_cells <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_B_cells <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_B_cells <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[hs_co_tmp_num_B_cells,],scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画人的数据
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[ts_co_tmp_num_B_cells,],scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$CLP)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$CLP)) )$peak_id
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_CLP.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_CLP.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
ts_co_tmp_num_CLP <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_CLP <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_CLP <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_CLP <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_CLP <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_CLP <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Erythroid_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Erythroid_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_Erythroid_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_Erythroid_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_Erythroid_cells <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_Erythroid_cells <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_Erythroid_cells <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_Erythroid_cells <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_Erythroid_cells <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_Erythroid_cells <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$HSPCs)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$HSPCs)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_HSPCs.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_HSPCs.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_HSPCs <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_HSPCs <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_HSPCs <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_HSPCs <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_HSPCs <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_HSPCs <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$HSPCs_MPO)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$HSPCs_MPO)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_HSPCs_MPO.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_HSPCs_MPO.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_HSPCs_MPO <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_HSPCs_MPO <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_HSPCs_MPO <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_HSPCs_MPO <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_HSPCs_MPO <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_HSPCs_MPO <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Kupffer_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Kupffer_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_Kupffer_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_Kupffer_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_Kupffer_cells <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_Kupffer_cells <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_Kupffer_cells <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_Kupffer_cells <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_Kupffer_cells <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_Kupffer_cells <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Macrophages)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Macrophages)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_Macrophages.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_Macrophages.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_Macrophages <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_Macrophages <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_Macrophages <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_Macrophages <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_Macrophages <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_Macrophages <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Mast_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Mast_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_Mast_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_Mast_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_Mast_cells <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_Mast_cells <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_Mast_cells <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_Mast_cells <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_Mast_cells <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_Mast_cells <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Megakaryocytes)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Megakaryocytes)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_Megakaryocytes.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_Megakaryocytes.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_Megakaryocytes <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_Megakaryocytes <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_Megakaryocytes <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_Megakaryocytes <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_Megakaryocytes <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_Megakaryocytes <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Microglia)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Microglia)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_Microglia.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_Microglia.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_Microglia <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_Microglia <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_Microglia <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_Microglia <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_Microglia <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_Microglia <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Monocytes)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Monocytes)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_Monocytes.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_Monocytes.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_Monocytes <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_Monocytes <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_Monocytes <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_Monocytes <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_Monocytes <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_Monocytes <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$NK_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$NK_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_NK_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_NK_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_NK_cells <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_NK_cells <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_NK_cells <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_NK_cells <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_NK_cells <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_NK_cells <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$T_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$T_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
write.table(subset(ts_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/ts_orth_co_T_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(subset(hs_orth_peak_df, peak_id %in% co_tmp ),file = "./Marker_heatmap/hs_orth_co_T_cells.peak.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ts_co_tmp_num_T_cells <- subset(ts_orth_peak_df, peak_id %in% co_tmp )$ts_orth_idx
hs_co_tmp_num_T_cells <- subset(hs_orth_peak_df, peak_id %in% co_tmp )$hs_orth_idx
ts_sp_tmp_num_T_cells <- subset(ts_orth_peak_df, peak_id %in% ts_sp_tmp )$ts_orth_idx
hs_sp_tmp_num_T_cells <- subset(hs_orth_peak_df, peak_id %in% hs_sp_tmp )$hs_orth_idx
ts_hs_tmp_num_T_cells <- subset(ts_orth_peak_df, peak_id %in% hs_sp_tmp )$ts_orth_idx
hs_ts_tmp_num_T_cells <- subset(hs_orth_peak_df, peak_id %in% ts_sp_tmp )$hs_orth_idx
#
#
co_num <- c(length(ts_co_tmp_num_B_cells),
            length(ts_co_tmp_num_CLP),
            length(ts_co_tmp_num_Erythroid_cells),
            length(ts_co_tmp_num_HSPCs),
            length(ts_co_tmp_num_HSPCs_MPO),
            length(ts_co_tmp_num_Kupffer_cells),
            length(ts_co_tmp_num_Macrophages),
            length(ts_co_tmp_num_Mast_cells),
            length(ts_co_tmp_num_Megakaryocytes),
            length(ts_co_tmp_num_Microglia),
            length(ts_co_tmp_num_Monocytes),
            length(ts_co_tmp_num_NK_cells),
            length(ts_co_tmp_num_T_cells))
#
hs_sp_num <- c(length(hs_sp_tmp_num_B_cells),
               length(hs_sp_tmp_num_CLP),
               length(hs_sp_tmp_num_Erythroid_cells),
               length(hs_sp_tmp_num_HSPCs),
               length(hs_sp_tmp_num_HSPCs_MPO),
               length(hs_sp_tmp_num_Kupffer_cells),
               length(hs_sp_tmp_num_Macrophages),
               length(hs_sp_tmp_num_Mast_cells),
               length(hs_sp_tmp_num_Megakaryocytes),
               length(hs_sp_tmp_num_Microglia),
               length(hs_sp_tmp_num_Monocytes),
               length(hs_sp_tmp_num_NK_cells),
               length(hs_sp_tmp_num_T_cells))
#
#
ts_sp_num <- c(length(ts_sp_tmp_num_B_cells),
               length(ts_sp_tmp_num_CLP),
               length(ts_sp_tmp_num_Erythroid_cells),
               length(ts_sp_tmp_num_HSPCs),
               length(ts_sp_tmp_num_HSPCs_MPO),
               length(ts_sp_tmp_num_Kupffer_cells),
               length(ts_sp_tmp_num_Macrophages),
               length(ts_sp_tmp_num_Mast_cells),
               length(ts_sp_tmp_num_Megakaryocytes),
               length(ts_sp_tmp_num_Microglia),
               length(ts_sp_tmp_num_Monocytes),
               length(ts_sp_tmp_num_NK_cells),
               length(ts_sp_tmp_num_T_cells))
#
#
barplot(co_num/c(hs_sp_num+co_num))
barplot(co_num/c(ts_sp_num+co_num))
#
barplot(co_num)
barplot(hs_sp_num)
barplot(ts_sp_num)
#
sum(co_num)
sum(hs_sp_num)
sum(ts_sp_num)
##
#
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(hs_co_tmp_num_B_cells,
                                                     hs_co_tmp_num_CLP,
                                                     hs_co_tmp_num_Erythroid_cells,
                                                     hs_co_tmp_num_HSPCs,
                                                     hs_co_tmp_num_HSPCs_MPO,
                                                     hs_co_tmp_num_Kupffer_cells,
                                                     hs_co_tmp_num_Macrophages,
                                                     hs_co_tmp_num_Mast_cells,
                                                     hs_co_tmp_num_Megakaryocytes,
                                                     hs_co_tmp_num_Microglia,
                                                     hs_co_tmp_num_Monocytes,
                                                     hs_co_tmp_num_NK_cells,
                                                     hs_co_tmp_num_T_cells),],
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画人的数据
#
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(ts_co_tmp_num_B_cells,
                                                     ts_co_tmp_num_CLP,
                                                     ts_co_tmp_num_Erythroid_cells,
                                                     ts_co_tmp_num_HSPCs,
                                                     ts_co_tmp_num_HSPCs_MPO,
                                                     ts_co_tmp_num_Kupffer_cells,
                                                     ts_co_tmp_num_Macrophages,
                                                     ts_co_tmp_num_Mast_cells,
                                                     ts_co_tmp_num_Megakaryocytes,
                                                     ts_co_tmp_num_Microglia,
                                                     ts_co_tmp_num_Monocytes,
                                                     ts_co_tmp_num_NK_cells,
                                                     ts_co_tmp_num_T_cells),],
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据
#
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(hs_sp_tmp_num_B_cells,
                                                     hs_sp_tmp_num_CLP,
                                                     hs_sp_tmp_num_Erythroid_cells,
                                                     hs_sp_tmp_num_HSPCs,
                                                     hs_sp_tmp_num_HSPCs_MPO,
                                                     hs_sp_tmp_num_Kupffer_cells,
                                                     hs_sp_tmp_num_Macrophages,
                                                     hs_sp_tmp_num_Mast_cells,
                                                     hs_sp_tmp_num_Megakaryocytes,
                                                     hs_sp_tmp_num_Microglia,
                                                     hs_sp_tmp_num_Monocytes,
                                                     hs_sp_tmp_num_NK_cells,
                                                     hs_sp_tmp_num_T_cells),],
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-开放，同源树鼩不开放，画人的数据
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(ts_hs_tmp_num_B_cells,
                                                     ts_hs_tmp_num_CLP,
                                                     ts_hs_tmp_num_Erythroid_cells,
                                                     ts_hs_tmp_num_HSPCs,
                                                     ts_hs_tmp_num_HSPCs_MPO,
                                                     ts_hs_tmp_num_Kupffer_cells,
                                                     ts_hs_tmp_num_Macrophages,
                                                     ts_hs_tmp_num_Mast_cells,
                                                     ts_hs_tmp_num_Megakaryocytes,
                                                     ts_hs_tmp_num_Microglia,
                                                     ts_hs_tmp_num_Monocytes,
                                                     ts_hs_tmp_num_NK_cells,
                                                     ts_hs_tmp_num_T_cells),],
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-开放，同源树鼩不开放，画树鼩的数据
#
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(hs_ts_tmp_num_B_cells,
                                                     hs_ts_tmp_num_CLP,
                                                     hs_ts_tmp_num_Erythroid_cells,
                                                     hs_ts_tmp_num_HSPCs,
                                                     hs_ts_tmp_num_HSPCs_MPO,
                                                     hs_ts_tmp_num_Kupffer_cells,
                                                     hs_ts_tmp_num_Macrophages,
                                                     hs_ts_tmp_num_Mast_cells,
                                                     hs_ts_tmp_num_Megakaryocytes,
                                                     hs_ts_tmp_num_Microglia,
                                                     hs_ts_tmp_num_Monocytes,
                                                     hs_ts_tmp_num_NK_cells,
                                                     hs_ts_tmp_num_T_cells),],
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-不开放，同源树鼩开放，画人的数据
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(ts_sp_tmp_num_B_cells,
                                                     ts_sp_tmp_num_CLP,
                                                     ts_sp_tmp_num_Erythroid_cells,
                                                     ts_sp_tmp_num_HSPCs,
                                                     ts_sp_tmp_num_HSPCs_MPO,
                                                     ts_sp_tmp_num_Kupffer_cells,
                                                     ts_sp_tmp_num_Macrophages,
                                                     ts_sp_tmp_num_Mast_cells,
                                                     ts_sp_tmp_num_Megakaryocytes,
                                                     ts_sp_tmp_num_Microglia,
                                                     ts_sp_tmp_num_Monocytes,
                                                     ts_sp_tmp_num_NK_cells,
                                                     ts_sp_tmp_num_T_cells),],
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-开放，同源树鼩不开放，画树鼩的数据
#
#
#
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(hs_sp_tmp_num_B_cells[1:100],
                                                     hs_sp_tmp_num_CLP[1:100],
                                                     hs_sp_tmp_num_Erythroid_cells[1:100],
                                                     hs_sp_tmp_num_HSPCs[1:100],
                                                     hs_sp_tmp_num_HSPCs_MPO[1:100],
                                                     hs_sp_tmp_num_Kupffer_cells[1:100],
                                                     hs_sp_tmp_num_Macrophages[1:100],
                                                     hs_sp_tmp_num_Mast_cells[1:100],
                                                     hs_sp_tmp_num_Megakaryocytes[1:100],
                                                     hs_sp_tmp_num_Microglia[1:100],
                                                     hs_sp_tmp_num_Monocytes[1:100],
                                                     hs_sp_tmp_num_NK_cells[1:100],
                                                     hs_sp_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-开放，同源树鼩不开放，画人的数据
#
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(ts_hs_tmp_num_B_cells[1:100],
                                                     ts_hs_tmp_num_CLP[1:100],
                                                     ts_hs_tmp_num_Erythroid_cells[1:100],
                                                     ts_hs_tmp_num_HSPCs[1:100],
                                                     ts_hs_tmp_num_HSPCs_MPO[1:100],
                                                     ts_hs_tmp_num_Kupffer_cells[1:100],
                                                     ts_hs_tmp_num_Macrophages[1:100],
                                                     ts_hs_tmp_num_Mast_cells[1:100],
                                                     ts_hs_tmp_num_Megakaryocytes[1:100],
                                                     ts_hs_tmp_num_Microglia[1:100],
                                                     ts_hs_tmp_num_Monocytes[1:100],
                                                     ts_hs_tmp_num_NK_cells[1:100],
                                                     ts_hs_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 树鼩-开放，同源人不开放，画人的数据
########
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(hs_co_tmp_num_B_cells[1:50],
                                                     hs_co_tmp_num_CLP[1:50],
                                                     hs_co_tmp_num_Erythroid_cells[1:50],
                                                     hs_co_tmp_num_HSPCs[1:50],
                                                     hs_co_tmp_num_HSPCs_MPO[1:50],
                                                     hs_co_tmp_num_Kupffer_cells[1:50],
                                                     hs_co_tmp_num_Macrophages[1:50],
                                                     hs_co_tmp_num_Mast_cells[1:50],
                                                     hs_co_tmp_num_Megakaryocytes[1:50],
                                                     hs_co_tmp_num_Microglia[1:50],
                                                     hs_co_tmp_num_Monocytes[1:32],
                                                     hs_co_tmp_num_NK_cells[1:50],
                                                     hs_co_tmp_num_T_cells[1:50]),],
                   color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画人的数据
#
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(ts_co_tmp_num_B_cells[1:100],
                                                     ts_co_tmp_num_CLP[1:100],
                                                     ts_co_tmp_num_Erythroid_cells[1:100],
                                                     ts_co_tmp_num_HSPCs[1:100],
                                                     ts_co_tmp_num_HSPCs_MPO[1:100],
                                                     ts_co_tmp_num_Kupffer_cells[1:100],
                                                     ts_co_tmp_num_Macrophages[1:100],
                                                     ts_co_tmp_num_Mast_cells[1:100],
                                                     ts_co_tmp_num_Megakaryocytes[1:100],
                                                     ts_co_tmp_num_Microglia[1:100],
                                                     ts_co_tmp_num_Monocytes[1:32],
                                                     ts_co_tmp_num_NK_cells[1:100],
                                                     ts_co_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据
#

#

#

###### ###### #
###### ###### #
###### ###### #
#

pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(hs_sp_tmp_num_B_cells[1:100],
                                                     hs_sp_tmp_num_CLP[1:100],
                                                     hs_sp_tmp_num_Erythroid_cells[1:100],
                                                     hs_sp_tmp_num_HSPCs[1:100],
                                                     hs_sp_tmp_num_HSPCs_MPO[1:100],
                                                     hs_sp_tmp_num_Kupffer_cells[1:100],
                                                     hs_sp_tmp_num_Macrophages[1:100],
                                                     hs_sp_tmp_num_Mast_cells[1:100],
                                                     hs_sp_tmp_num_Megakaryocytes[1:100],
                                                     hs_sp_tmp_num_Microglia[1:100],
                                                     hs_sp_tmp_num_Monocytes[1:100],
                                                     hs_sp_tmp_num_NK_cells[1:100],
                                                     hs_sp_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous("solarExtra")[5:220],raster = TRUE, raster_device = "png",
                   show_rownames = F,show_colnames = F,
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-开放，同源树鼩不开放，画人的数据 # hs_open_ts_closed_plot_hs
#
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(ts_hs_tmp_num_B_cells[1:100],
                                                     ts_hs_tmp_num_CLP[1:100],
                                                     ts_hs_tmp_num_Erythroid_cells[1:100],
                                                     ts_hs_tmp_num_HSPCs[1:100],
                                                     ts_hs_tmp_num_HSPCs_MPO[1:100],
                                                     ts_hs_tmp_num_Kupffer_cells[1:100],
                                                     ts_hs_tmp_num_Macrophages[1:100],
                                                     ts_hs_tmp_num_Mast_cells[1:100],
                                                     ts_hs_tmp_num_Megakaryocytes[1:100],
                                                     ts_hs_tmp_num_Microglia[1:100],
                                                     ts_hs_tmp_num_Monocytes[1:100],
                                                     ts_hs_tmp_num_NK_cells[1:100],
                                                     ts_hs_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous("solarExtra")[5:220],raster = TRUE, raster_device = "png",
                   show_rownames = F,show_colnames = F,
                   scale = "row",cluster_rows = F,cluster_cols = F) # 树鼩-开放，同源人不开放，画人的数据 # ts_open_hs_closed_plot_hs

#
#
#
#
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(ts_sp_tmp_num_B_cells[1:100],
                                                     ts_sp_tmp_num_CLP[1:100],
                                                     ts_sp_tmp_num_Erythroid_cells[1:100],
                                                     ts_sp_tmp_num_HSPCs[1:100],
                                                     ts_sp_tmp_num_HSPCs_MPO[1:100],
                                                     ts_sp_tmp_num_Kupffer_cells[1:100],
                                                     ts_sp_tmp_num_Macrophages[1:100],
                                                     ts_sp_tmp_num_Mast_cells[1:100],
                                                     ts_sp_tmp_num_Megakaryocytes[1:100],
                                                     ts_sp_tmp_num_Microglia[1:100],
                                                     ts_sp_tmp_num_Monocytes[1:100],
                                                     ts_sp_tmp_num_NK_cells[1:100],
                                                     ts_sp_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous("solarExtra")[5:220],raster = TRUE, raster_device = "png",
                   show_rownames = F,show_colnames = F,
                   scale = "row",cluster_rows = F,cluster_cols = F) # 树鼩-开放，同源人不开放，画树鼩的数据  # ts_open_hs_closed_plot_ts
#
#
#
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(hs_ts_tmp_num_B_cells[1:100],
                                                     hs_ts_tmp_num_CLP[1:100],
                                                     hs_ts_tmp_num_Erythroid_cells[1:100],
                                                     hs_ts_tmp_num_HSPCs[1:100],
                                                     hs_ts_tmp_num_HSPCs_MPO[1:100],
                                                     hs_ts_tmp_num_Kupffer_cells[1:100],
                                                     hs_ts_tmp_num_Macrophages[1:100],
                                                     hs_ts_tmp_num_Mast_cells[1:100],
                                                     hs_ts_tmp_num_Megakaryocytes[1:100],
                                                     hs_ts_tmp_num_Microglia[1:100],
                                                     hs_ts_tmp_num_Monocytes[1:100],
                                                     hs_ts_tmp_num_NK_cells[1:100],
                                                     hs_ts_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous("solarExtra")[5:220],raster = TRUE, raster_device = "png",
                   show_rownames = F,show_colnames = F,
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-开放，同源树鼩不开放，画树鼩的数据  # hs_open_ts_closed_plot_ts
#






########
pheatmap::pheatmap(Hs_markerPeaks@assays@data$Mean[c(hs_co_tmp_num_B_cells[1:100],
                                                     hs_co_tmp_num_CLP[1:100],
                                                     hs_co_tmp_num_Erythroid_cells[1:100],
                                                     hs_co_tmp_num_HSPCs[1:100],
                                                     hs_co_tmp_num_HSPCs_MPO[1:100],
                                                     hs_co_tmp_num_Kupffer_cells[1:100],
                                                     hs_co_tmp_num_Macrophages[1:100],
                                                     hs_co_tmp_num_Mast_cells[1:100],
                                                     hs_co_tmp_num_Megakaryocytes[1:100],
                                                     hs_co_tmp_num_Microglia[1:100],
                                                     hs_co_tmp_num_Monocytes[1:43],
                                                     hs_co_tmp_num_NK_cells[1:100],
                                                     hs_co_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous("solarExtra")[5:220],raster = TRUE, raster_device = "png",
                   show_rownames = F,show_colnames = F,
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画人的数据  # ts_hs_co_plot_hs
#
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(ts_co_tmp_num_B_cells[1:100],
                                                     ts_co_tmp_num_CLP[1:100],
                                                     ts_co_tmp_num_Erythroid_cells[1:100],
                                                     ts_co_tmp_num_HSPCs[1:100],
                                                     ts_co_tmp_num_HSPCs_MPO[1:100],
                                                     ts_co_tmp_num_Kupffer_cells[1:100],
                                                     ts_co_tmp_num_Macrophages[1:100],
                                                     ts_co_tmp_num_Mast_cells[1:100],
                                                     ts_co_tmp_num_Megakaryocytes[1:100],
                                                     ts_co_tmp_num_Microglia[1:100],
                                                     ts_co_tmp_num_Monocytes[1:43],
                                                     ts_co_tmp_num_NK_cells[1:100],
                                                     ts_co_tmp_num_T_cells[1:100]),],
                   color = paletteContinuous("solarExtra")[5:220],raster = TRUE, raster_device = "png",
                   show_rownames = F,show_colnames = F,
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据 # ts_hs_co_plot_ts

#
save.image("CS_immune_plot.save.RData")
#

#
#
co_num <- c(ts_co_tmp_num_B_cells,
            ts_co_tmp_num_CLP,
            ts_co_tmp_num_Erythroid_cells,
            ts_co_tmp_num_HSPCs,
            ts_co_tmp_num_HSPCs_MPO,
            ts_co_tmp_num_Kupffer_cells,
            ts_co_tmp_num_Macrophages,
            ts_co_tmp_num_Mast_cells,
            ts_co_tmp_num_Megakaryocytes,
            ts_co_tmp_num_Microglia,
            ts_co_tmp_num_Monocytes,
            ts_co_tmp_num_NK_cells,
            ts_co_tmp_num_T_cells)
#
hs_co_num <- c(hs_co_tmp_num_B_cells,
            hs_co_tmp_num_CLP,
            hs_co_tmp_num_Erythroid_cells,
            hs_co_tmp_num_HSPCs,
            hs_co_tmp_num_HSPCs_MPO,
            hs_co_tmp_num_Kupffer_cells,
            hs_co_tmp_num_Macrophages,
            hs_co_tmp_num_Mast_cells,
            hs_co_tmp_num_Megakaryocytes,
            hs_co_tmp_num_Microglia,
            hs_co_tmp_num_Monocytes,
            hs_co_tmp_num_NK_cells,
            hs_co_tmp_num_T_cells)
#
hs_sp_num <- c(hs_sp_tmp_num_B_cells,
               hs_sp_tmp_num_CLP,
               hs_sp_tmp_num_Erythroid_cells,
               hs_sp_tmp_num_HSPCs,
               hs_sp_tmp_num_HSPCs_MPO,
               hs_sp_tmp_num_Kupffer_cells,
               hs_sp_tmp_num_Macrophages,
               hs_sp_tmp_num_Mast_cells,
               hs_sp_tmp_num_Megakaryocytes,
               hs_sp_tmp_num_Microglia,
               hs_sp_tmp_num_Monocytes,
               hs_sp_tmp_num_NK_cells,
               hs_sp_tmp_num_T_cells)
#
#
ts_sp_num <- c(ts_sp_tmp_num_B_cells,
               ts_sp_tmp_num_CLP,
               ts_sp_tmp_num_Erythroid_cells,
               ts_sp_tmp_num_HSPCs,
               ts_sp_tmp_num_HSPCs_MPO,
               ts_sp_tmp_num_Kupffer_cells,
               ts_sp_tmp_num_Macrophages,
               ts_sp_tmp_num_Mast_cells,
               ts_sp_tmp_num_Megakaryocytes,
               ts_sp_tmp_num_Microglia,
               ts_sp_tmp_num_Monocytes,
               ts_sp_tmp_num_NK_cells,
               ts_sp_tmp_num_T_cells)
#
length(co_num) # 9843
length(hs_co_num) # 9843
length(hs_sp_num) # 25985
length(ts_sp_num) # 59867
#
length(unique(co_num)) # 8104
length(unique(hs_co_num)) # 8104
length(unique(hs_sp_num)) # 20244
length(unique(ts_sp_num)) # 43772
#
#
############# # ############## # ############## # ############## # ############## # ############## # ############## # ############## # ############## # #
############# # ############## # ############## # ############## # ############## # ############## # ############## # ############## # ############## # #
############# # ############## # ############## # ############## # ############## # ############## # ############## # ############## # ############## # #
############# # ############## # ############## # ############## # ############## # ############## # ############## # ############## # ############## # #
############# # ############## # ############## # ############## # ############## # ############## # ############## # ############## # ############## # #
############# # ############## # ############## # ############## # ############## # ############## # ############## # ############## # ############## # #
# 
ts_orth_co.all.list <- read.csv("Marker_heatmap/ts_orth_co.cat.bed",header = F,sep = "\t")
head(ts_orth_co.all.list)
#
ts_orth_co.all.list$PEAKS <- paste(paste(ts_orth_co.all.list$V1,ts_orth_co.all.list$V2,sep=":"),ts_orth_co.all.list$V3,sep = "-")
#
head(HE_bedX)
#
TS_bedX.co <- subset(TS_bedX, PEAKS %in% ts_orth_co.all.list$PEAKS)
table(TS_bedX.co$ArchR_peakType)
table(TS_bed$ArchR_peakType)
table(TS_bedX$ArchR_peakType)
#
table(TS_bedX.co$ArchR_peakType)/sum(table(TS_bedX.co$ArchR_peakType))
table(TS_bed$ArchR_peakType)/sum(table(TS_bed$ArchR_peakType))
table(TS_bedX$ArchR_peakType)/sum(table(TS_bedX$ArchR_peakType))
#
TS_marker_long <- TS_markerList@listData %>% map_dfr(as.data.frame, .id = "celltype")
dim(TS_marker_long)
TS_marker_long <- TS_marker_long[which(TS_marker_long$Log2FC >= 0.5 & TS_marker_long$FDR <= 0.1),]
dim(TS_marker_long)
TS_marker_long$PEAKS <- paste(paste(TS_marker_long$seqnames,TS_marker_long$start,sep=":"),TS_marker_long$end,sep = "-")
length(unique(TS_marker_long$PEAKS))
#
TS_marker_long2 <- subset(TS_bed, PEAKS %in% TS_marker_long$PEAKS)
dim(TS_marker_long2)
#

#
#
dim(TS_bed) # 所有peaks # 1238436
dim(TS_bedX) # 免疫细胞 所有peaks # 270659
dim(TS_bed_keep2) # 108744 conserved
dim(TS_marker_long2) # marker peak # 48191
dim(TS_bedX.co)  # 8104 # conserved celltype specific的
#
#
barplot(c(1238436,270659,108744,48191,8104))
#
TS_ratio_counts <- cbind(table(TS_bed$ArchR_peakType)/sum(table(TS_bed$ArchR_peakType)),
table(TS_bedX$ArchR_peakType)/sum(table(TS_bedX$ArchR_peakType)),
table(TS_bed_keep2$ArchR_peakType)/sum(table(TS_bed_keep2$ArchR_peakType)),
table(TS_marker_long2$ArchR_peakType)/sum(table(TS_marker_long2$ArchR_peakType)),
table(TS_bedX.co$ArchR_peakType)/sum(table(TS_bedX.co$ArchR_peakType)))
#
colnames(TS_ratio_counts) <- c("All","immune","conserved","marker","specific")
#
#
# 加载必要的包
library(ggplot2)
library(tidyr)
library(tidyr)
library(tibble)
## 转为长格式
TS_long <- TS_ratio_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "Region") %>%
  pivot_longer(cols = -Region, names_to = "Category", values_to = "Ratio")

# 设置顺序
TS_long$Category <- factor(TS_long$Category, 
                           levels = c("All","immune","conserved","marker","specific"))

# 绘图
ggplot(TS_long, aes(x = Category, y = Ratio, fill = Region)) +
  geom_col(position = "stack") +
  labs(title = "TS Ratio Composition", x = "Category", y = "Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
################ ####。#### #
################ ####。#### #
################ ####。#### #
############# # ############## # ############## # ############## # ############## # ############## # ############## # ############## # ############## # #
############# # ############## # ############## # ############## # ############## # ############## # ############## # ############## # ############## # #
# 
hs_orth_co.all.list <- read.csv("Marker_heatmap/hs_orth_co.cat.bed",header = F,sep = "\t")
head(hs_orth_co.all.list)
dim(hs_orth_co.all.list)
#
hs_orth_co.all.list$PEAKS <- paste(paste(hs_orth_co.all.list$V1,hs_orth_co.all.list$V2,sep=":"),hs_orth_co.all.list$V3,sep = "-")
length(hs_orth_co.all.list$PEAKS)
length(unique(hs_orth_co.all.list$PEAKS))
#
head(HE_bed)
#
HE_bedX.co <- subset(HE_bed, PEAKS %in% hs_orth_co.all.list$PEAKS)
dim(HE_bedX.co)
head(HE_bedX.co)
#
HE_marker_long <- Hs_markerList@listData %>% map_dfr(as.data.frame, .id = "celltype")
dim(HE_marker_long)
HE_marker_long <- HE_marker_long[which(HE_marker_long$Log2FC >= 0.5 & HE_marker_long$FDR <= 0.1),]
dim(HE_marker_long)
HE_marker_long$PEAKS <- paste(paste(HE_marker_long$seqnames,HE_marker_long$start,sep=":"),HE_marker_long$end,sep = "-")
length(unique(HE_marker_long$PEAKS))
#
HE_marker_long2 <- subset(HE_bed, PEAKS %in% HE_marker_long$PEAKS)
dim(HE_marker_long2)
#
#
#
dim(HE_bed) # 所有peaks # 1016382
dim(HE_bedX) # 免疫细胞 所有peaks # 244562
dim(HE_bed_keep2) # 108744 conserved
dim(HE_marker_long2) # marker peak # 25529
dim(HE_bedX.co)  # 8104 # conserved celltype specific的
#
#
barplot(c(1016382,244562,108744,25529,8104))
#
HE_ratio_counts <- cbind(table(HE_bed$ArchR_peakType)/sum(table(HE_bed$ArchR_peakType)),
                         table(HE_bedX$ArchR_peakType)/sum(table(HE_bedX$ArchR_peakType)),
                         table(HE_bed_keep2$ArchR_peakType)/sum(table(HE_bed_keep2$ArchR_peakType)),
                         table(HE_marker_long2$ArchR_peakType)/sum(table(HE_marker_long2$ArchR_peakType)),
                         table(HE_bedX.co$ArchR_peakType)/sum(table(HE_bedX.co$ArchR_peakType)))
#
colnames(HE_ratio_counts) <- c("All","immune","conserved","marker","specific")
#
#
# 加载必要的包
library(ggplot2)
library(tidyr)
library(tidyr)
library(tibble)
## 转为长格式
HE_long <- HE_ratio_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "Region") %>%
  pivot_longer(cols = -Region, names_to = "Category", values_to = "Ratio")

# 设置顺序
HE_long$Category <- factor(HE_long$Category, 
                           levels = c("All","immune","conserved","marker","specific"))

# 绘图
ggplot(HE_long, aes(x = Category, y = Ratio, fill = Region)) +
  geom_col(position = "stack") +
  labs(title = "TS Ratio Composition", x = "Category", y = "Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
#


#


#





#



#


#





#



#


#





#



#


#





#



#





#



#


















































































































































































































###### ###### #
###### ###### #
###### ###### #
#
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[c(intersect(TS_markerList$B_cells$idx,ts_co_tmp_num_B_cells[1:150]),
                                                     intersect(TS_markerList$CLP$idx,ts_co_tmp_num_CLP[1:150]),
                                                     intersect(TS_markerList$Erythroid_cells$idx,ts_co_tmp_num_Erythroid_cells[1:150]),
                                                     intersect(TS_markerList$HSPCs$idx,ts_co_tmp_num_HSPCs[1:150]),
                                                     intersect(TS_markerList$HSPCs_MPO$idx,ts_co_tmp_num_HSPCs_MPO[1:150]),
                                                     intersect(TS_markerList$Kupffer_cells$idx,ts_co_tmp_num_Kupffer_cells[1:150]),
                                                     intersect(TS_markerList$Macrophages$idx,ts_co_tmp_num_Macrophages[1:150]),
                                                     intersect(TS_markerList$Mast_cells$idx,ts_co_tmp_num_Mast_cells[1:150]),
                                                     intersect(TS_markerList$Megakaryocytes$idx,ts_co_tmp_num_Megakaryocytes[1:150]),
                                                     intersect(TS_markerList$Microglia$idx,ts_co_tmp_num_Microglia[1:150]),
                                                     intersect(TS_markerList$Monocytes$idx,ts_co_tmp_num_Monocytes[1:150]),
                                                     intersect(TS_markerList$NK_cells$idx,ts_co_tmp_num_NK_cells[1:150]),
                                                     intersect(TS_markerList$T_cells$idx,ts_co_tmp_num_T_cells[1:150])),],color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据
#
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[intersect(TS_markerList$B_cells$idx,ts_co_tmp_num_B_cells),],color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[ts_co_tmp_num_B_cells,],color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[TS_markerList$B_cells$idx,],color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据
pheatmap::pheatmap(TS_markerPeaks@assays@data$Mean[TS_markerList$B_cells$idx[1:10],],color = paletteContinuous(),
                   scale = "row",cluster_rows = F,cluster_cols = F) # 人-树鼩都开放，画树鼩的数据
#
#
unique(c(unique(ts_co_tmp_num_B_cells),
         unique(ts_co_tmp_num_CLP),
         unique(ts_co_tmp_num_Erythroid_cells),
         unique(ts_co_tmp_num_HSPCs),
         unique(ts_co_tmp_num_HSPCs_MPO),
         unique(ts_co_tmp_num_Kupffer_cells),
         unique(ts_co_tmp_num_Macrophages),
         unique(ts_co_tmp_num_Mast_cells),
         unique(ts_co_tmp_num_Megakaryocytes),
         unique(ts_co_tmp_num_Microglia),
         unique(ts_co_tmp_num_Monocytes),
         unique(ts_co_tmp_num_NK_cells),
         unique(ts_co_tmp_num_T_cells)))
#
#
#
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$B_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$B_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp)
#
TS_markerList_unlist <- unlist(TS_markerList)
TS_markerList_unlist$celltype <- unlist(strsplit(rownames(TS_markerList_unlist),".",fixed=T))[seq(1,2*length(rownames(TS_markerList_unlist))[1],2)]
TS_markerList_unlist$rownames <- unlist(strsplit(rownames(TS_markerList_unlist),".",fixed=T))[seq(2,2*length(rownames(TS_markerList_unlist))[1],2)]
head(TS_markerList_unlist)
Hs_markerList_unlist <- unlist(Hs_markerList)
Hs_markerList_unlist$celltype <- unlist(strsplit(rownames(Hs_markerList_unlist),".",fixed=T))[seq(1,2*length(rownames(Hs_markerList_unlist))[1],2)]
Hs_markerList_unlist$rownames <- unlist(strsplit(rownames(Hs_markerList_unlist),".",fixed=T))[seq(2,2*length(rownames(Hs_markerList_unlist))[1],2)]
head(Hs_markerList_unlist)
#
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% TS_markerList_unlist$rownames )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% Hs_markerList_unlist$rownames )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp)
length(co_tmp)
length(unique(co_tmp))
#
#
#
#
#
#


#
TS_co_matrix <- TS_markerPeaks@assays@data$Mean[c(ts_co_tmp_num_B_cells,
                                                  ts_co_tmp_num_CLP,
                                                  ts_co_tmp_num_Erythroid_cells,
                                                  ts_co_tmp_num_HSPCs,
                                                  ts_co_tmp_num_HSPCs_MPO,
                                                  ts_co_tmp_num_Kupffer_cells,
                                                  ts_co_tmp_num_Macrophages,
                                                  ts_co_tmp_num_Mast_cells,
                                                  ts_co_tmp_num_Megakaryocytes,
                                                  ts_co_tmp_num_Microglia,
                                                  ts_co_tmp_num_Monocytes,
                                                  ts_co_tmp_num_NK_cells,
                                                  ts_co_tmp_num_T_cells),]
#
HS_co_matrix <- Hs_markerPeaks@assays@data$Mean[c(hs_co_tmp_num_B_cells,
                                                  hs_co_tmp_num_CLP,
                                                  hs_co_tmp_num_Erythroid_cells,
                                                  hs_co_tmp_num_HSPCs,
                                                  hs_co_tmp_num_HSPCs_MPO,
                                                  hs_co_tmp_num_Kupffer_cells,
                                                  hs_co_tmp_num_Macrophages,
                                                  hs_co_tmp_num_Mast_cells,
                                                  hs_co_tmp_num_Megakaryocytes,
                                                  hs_co_tmp_num_Microglia,
                                                  hs_co_tmp_num_Monocytes,
                                                  hs_co_tmp_num_NK_cells,
                                                  hs_co_tmp_num_T_cells),]
#
colnames(TS_co_matrix) <- paste("TS_",c("B_cells","CLP","Erythroid_cells","HSPCs","HSPCs_MPO","Kupffer_cells","Macrophages","Mast_cells",
                                        "Megakaryocytes","Microglia","Monocytes","NK_cells","T_cells"),sep = "_")
colnames(HS_co_matrix) <- paste("HS_",c("B_cells","CLP","Erythroid_cells","HSPCs","HSPCs_MPO","Kupffer_cells","Macrophages","Mast_cells",
                                        "Megakaryocytes","Microglia","Monocytes","NK_cells","T_cells"),sep = "_")
#
head(HS_co_matrix)
head(TS_co_matrix)
#
#
pheatmap::pheatmap(t(scale(t(TS_co_matrix))),cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(t(scale(t(HS_co_matrix))),cluster_rows = F,cluster_cols = F)
#
cbind(TS_co_matrix,HS_co_matrix)
pheatmap::pheatmap(cor(cbind(TS_co_matrix,HS_co_matrix)))
#
TS_co_matrix2 <- t(scale(t(TS_co_matrix)))
HS_co_matrix2 <- t(scale(t(HS_co_matrix)))
#
pheatmap::pheatmap(cor(cbind(TS_co_matrix2,HS_co_matrix2)),cluster_rows = F,cluster_cols = F)
#
cor_values <- cor(cbind(TS_co_matrix2,HS_co_matrix2))
pheatmap::pheatmap(cor_values[1:13,14:26],cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(cor_values[1:13,14:26],cluster_rows = T,cluster_cols = T)
pheatmap::pheatmap(cor_values[1:13,14:26],cluster_rows = T,cluster_cols = T)
#
pheatmap::pheatmap(cor_values[1:13,paste("HS_",c("HSPCs_MPO",
                                                 "Monocytes",
                                                 "Microglia",
                                                 "Kupffer_cells",
                                                 "Macrophages",
                                                 "Mast_cells",
                                                 "Erythroid_cells",
                                                 "Megakaryocytes",
                                                 "B_cells",
                                                 "NK_cells",
                                                 "T_cells",
                                                 "CLP",
                                                 "HSPCs"),sep = "_")],cluster_rows = T,cluster_cols = F)
#
#
rank_id_TS <- paste("TS_",c("CLP","HSPCs","NK_cells","T_cells","B_cells","pDCs","Mast_cells","Erythroid_cells","Megakaryocytes","Osteoclasts",
                            "Monocytes","Microglia","Kupffer_cells","Macrophages","cDCs","HSPCs_MPO"),sep = "_")
rank_id_HS <- paste("HS_",c("CLP","HSPCs","NK_cells","T_cells","B_cells","pDCs","Mast_cells","Erythroid_cells","Megakaryocytes","Osteoclasts",
                            "Monocytes","Microglia","Kupffer_cells","Macrophages","cDCs","HSPCs_MPO"),sep = "_")
#
pheatmap::pheatmap(cor_values[rank_id_TS,rank_id_HS],cluster_rows = F,cluster_cols = F)
#
cor_values2 <- cor(cbind(TS_co_matrix2,HS_co_matrix2),method = c("spearman"))
pheatmap::pheatmap(cor_values2[rank_id_TS,rank_id_HS],cluster_rows = F,cluster_cols = F)
#
pheatmap::pheatmap(cor_values2[c(1,3,4,5,6,7,8,9,10,11,13,16),c(1,3,4,5,6,7,8,9,10,11,13,16)+16],cluster_rows = T,cluster_cols = T)
#
#
#
#
#
#
#
library(ChIPseeker)
# 从 BED 文件读取
peak_ann <- "/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis/setp1_peak/orth.human.sort.final.bed"
peak_ann <- readPeakFile(peak_ann)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(peak_ann, 
                         tssRegion = c(-2000, 200),   # TSS上下游区域定义为promoter
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")      # 添加基因ID转换（可选）
# 转换为数据框
anno_df <- as.data.frame(peakAnno)
head(anno_df)
# 查看注释统计（每个peak被分配到哪个区域）
table(anno_df$annotation)
#
table(substring(anno_df$annotation,1,5))
table(substring(anno_df$annotation,1,5))/sum(table(substring(anno_df$annotation,1,5)))
#
#
# 3' UT  5' UT  Dista  Downs  Exon   Intro  Promo 
# 12241   2205  67410    359  18290 147861  47779 
#
anno_df_tmp <- subset(anno_df, V4 %in% rownames(TS_co_matrix2))
table(substring(anno_df_tmp$annotation,1,5))
table(substring(anno_df_tmp$annotation,1,5))/sum(table(substring(anno_df_tmp$annotation,1,5)))
#
# 3' UT 5' UT Dista Downs Exon  Intro Promo 
# 283    65  1523    16   430  3755  2296 
#
47779/296145
2296/8368
#