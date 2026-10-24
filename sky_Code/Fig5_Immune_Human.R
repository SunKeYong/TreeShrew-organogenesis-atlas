setwd("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step2_ArchR/human_ArchR")
#
#
.libPaths("/home/sunkeyong/R/x86_64-pc-linux-gnu-library/4.4")
library(ArchR) 
library(Seurat) 
library(dplyr)
library(GenomeInfoDb)
library(ensembldb)
library(ggplot2)
library(patchwork)
set.seed(1234)
addArchRThreads(threads = 8)
library(BSgenome) 
#
addArchRGenome("hg38")
addArchRThreads(threads = 16) 
#
ArrowFiles <- c("/home/sunkeyong/Human_em/Cross_species/Immune_Human/ArchR/arrow_files/Human_ery.arrow",
                "/home/sunkeyong/Human_em/Cross_species/Immune_Human/ArchR/arrow_files/Human_immune.arrow")
#
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "CS_Human_Immune_ArchRoutput",
  copyArrows = TRUE
)
#
#
projHeme1$cells <- unlist(strsplit(projHeme1$cellNames,"#"))[seq(2,2*dim(projHeme1@cellColData)[1],2)]
head(projHeme1$cells)
#
load("/home/sunkeyong/Human_em/Cross_species/Immune_Human/H_WBC_Ery.seurat.type.final.RData")
DimPlot(H_WBC_Ery.seurat.type.final,reduction = "UMAP_ArchR2",label=T)+NoLegend()
#
H_WBC_Ery.seurat.type.final.meta <- as.data.frame(H_WBC_Ery.seurat.type.final@meta.data)
head(H_WBC_Ery.seurat.type.final.meta)
H_WBC_Ery.seurat.type.final.meta$ID_defined <- as.character(H_WBC_Ery.seurat.type.final@active.ident)
#
Before_axis <- as.data.frame(H_WBC_Ery.seurat.type.final@reductions$UMAP_ArchR2@cell.embeddings)
head(Before_axis)
#
H_WBC_Ery.seurat.type.final.meta <- H_WBC_Ery.seurat.type.final.meta[projHeme1$cells,]
Before_axis <- Before_axis[projHeme1$cells,]
#
#
projHeme1$ID_defined <- as.character(H_WBC_Ery.seurat.type.final.meta$ID_defined)
projHeme1$embryo <- as.character(H_WBC_Ery.seurat.type.final.meta$embryo)
projHeme1$time <- as.character(H_WBC_Ery.seurat.type.final.meta$time)
projHeme1$batch <- as.character(H_WBC_Ery.seurat.type.final.meta$batch)
projHeme1$Before_axis_UMAP1 <- as.character(Before_axis$UMAPArchR2_1)
projHeme1$Before_axis_UMAP2 <- as.character(Before_axis$UMAPArchR2_2)
#
p <- plotBrowserTrack(ArchRProj = projHeme1, groupBy = "ID_defined", geneSymbol = "HBB", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),)
grid::grid.newpage()
grid::grid.draw(p$HBB)
##
#
#
#
#
save.image("CS_Human_Immune.RData")
#
#















#
HT_keep_celltypes <- c("B cells",
                       "Border-associated macrophages",
                       "Common lymphoid progenitor",
                       "Erythrocyte",
                       "Haematopoietic stem and progenitor cells",
                       "Haematopoietic stem and progenitor cells (Mpo+)",
                       "Kupffer cells",
                       "Mast cells",
                       "Megakaryocyte",
                       "Microglia",
                       "Monocytes",
                       "Natural killer cells",
                       "T cells")
#
H_WBC_Ery.seurat.type.final$ID_defined <- as.character(H_WBC_Ery.seurat.type.final@active.ident)
DimPlot(H_WBC_Ery.seurat.type.final,reduction = "UMAP_ArchR2",label=T)+NoLegend()
#
H_WBC_Ery.seurat.type.final_meta <- as.data.frame(H_WBC_Ery.seurat.type.final@meta.data)
head(H_WBC_Ery.seurat.type.final_meta)
H_WBC_Ery.seurat.type.final_meta_keep <- subset(H_WBC_Ery.seurat.type.final_meta, ID_defined %in% HT_keep_celltypes)
head(H_WBC_Ery.seurat.type.final_meta_keep)
dim(H_WBC_Ery.seurat.type.final_meta)
dim(H_WBC_Ery.seurat.type.final_meta_keep)
#
#
CS_Human_Immune.seurat.type_keep2 <- subset(H_WBC_Ery.seurat.type.final,cells=rownames(H_WBC_Ery.seurat.type.final_meta_keep)) 
CS_Human_Immune.seurat.type_keep2@active.ident <- as.factor(CS_Human_Immune.seurat.type_keep2$ID_defined)
table(CS_Human_Immune.seurat.type_keep2@active.ident)
#
DimPlot(CS_Human_Immune.seurat.type_keep2,reduction = "UMAP_ArchR2",label=T)+NoLegend()
#
#
#
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"T cells"="T_cells")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Border-associated macrophages"="Macrophages")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Common lymphoid progenitor"="CLP")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Erythrocyte"="Erythroid_cells")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Haematopoietic stem and progenitor cells"="HSPCs")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Kupffer cells"="Kupffer_cells")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Mast cells"="Mast_cells")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Megakaryocyte"="Megakaryocytes")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Microglia"="Microglia")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Monocytes"="Monocytes")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Haematopoietic stem and progenitor cells (Mpo+)"="HSPCs_MPO")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"Natural killer cells"="NK_cells")
CS_Human_Immune.seurat.type_keep2 <- RenameIdents(CS_Human_Immune.seurat.type_keep2,"B cells"="B_cells")
#
CS_Human_Immune.seurat.type_keep2$ID_defined_short <- as.character(CS_Human_Immune.seurat.type_keep2@active.ident)
#
DimPlot(CS_Human_Immune.seurat.type_keep2,reduction = "UMAP_ArchR2",group.by = "ID_defined",
        label = T,repel = F,raster = F,raster.dpi=c(1512, 1512))+NoLegend()
DimPlot(CS_Human_Immune.seurat.type_keep2,reduction = "UMAP_ArchR2",group.by = "ID_defined_short",
        label = T,repel = F,raster = F,raster.dpi=c(1512, 1512))+NoLegend()
#
#
table(CS_Human_Immune.seurat.type_keep2$ID_defined_short)
#
################
#
#
peak_df <- read.table("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/orth.human.raw.sort.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(peak_df) <- c("chr", "start", "end", "peak_id")
# 转换为 GRanges
hs_orth_gr <- GRanges(
  seqnames = peak_df$chr,
  ranges = IRanges(start = peak_df$start, end = peak_df$end),
  peak_id = peak_df$peak_id
)
#
hs_orth_gr
addArchRThreads(threads = 2) 
#
projHeme1_meta <- as.data.frame(projHeme1@cellColData)
projHeme1_meta$cellNames <- rownames(projHeme1_meta)
projHeme1_meta$barcodes <- unlist(strsplit(projHeme1_meta$cellNames,"#"))[seq(2,2*length(projHeme1_meta$cellNames),2)]
rownames(projHeme1_meta) <- projHeme1_meta$barcodes
#
projHeme1_meta_keep <- projHeme1_meta[Cells(CS_Human_Immune.seurat.type_keep2),]
#
projHeme1$barcodes <- unlist(strsplit(projHeme1$cellNames,"#"))[seq(2,2*length(projHeme1$cellNames),2)]
projHeme1_keep <- projHeme1[projHeme1_meta_keep$cellNames,]
projHeme1_keep
CS_Human_Immune.seurat.type_keep2
#
CS_Human_Immune.seurat.type_keep2_meta <- as.data.frame(CS_Human_Immune.seurat.type_keep2@meta.data)
CS_Human_Immune.seurat.type_keep2_meta2 <- CS_Human_Immune.seurat.type_keep2_meta[projHeme1_keep$barcodes,]
#
projHeme1_keep$ID_defined2 <- as.character(CS_Human_Immune.seurat.type_keep2_meta2$ID_defined_short)
#
table(projHeme1_keep$ID_defined2)
#
addArchRThreads(threads = 8)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "HBB", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,  
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = getPeakSet(projHeme1_keep))
grid::grid.newpage()
grid::grid.draw(p$HBB)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "MPO", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = getPeakSet(projHeme1_keep))
grid::grid.newpage()
grid::grid.draw(p$MPO)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CD79A", 
                      upstream = 60000,downstream = 60000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = getPeakSet(projHeme1_keep))
grid::grid.newpage()
grid::grid.draw(p$CD79A)
#
#
#
#
#
#
addArchRThreads(threads = 5)
#
#
projHeme1_keep <- addPeakSet(projHeme1_keep,peakSet = hs_orth_gr,force=TRUE)
#
save.image("CS_Human_Immune.RData")
#
projHeme1_keep <- addPeakMatrix(projHeme1_keep,force = T)
getAvailableMatrices(projHeme1_keep)
#
markerPeaks <- getMarkerFeatures(
  ArchRProj = projHeme1_keep, 
  useMatrix = "PeakMatrix", 
  groupBy = "ID_defined2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#
save.image("CS_Human_Immune.RData")
#
#
plotEmbedding(ArchRProj = projHeme1_keep, colorBy = "cellColData", name = "ID_defined2", embedding = "UMAP2")
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "HBB", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),)
grid::grid.newpage()
grid::grid.draw(p$HBB)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CSF1R", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),)
grid::grid.newpage()
grid::grid.draw(p$CSF1R)
##
#
#
#
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList$B_cells
markerList$Macrophages
markerList$pDCs
markerList$Osteoclasts
markerList$Monocytes
markerList$Microglia
#
markerList_GR <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList_GR
markerList_GR$Erythroid
#
#
heatmapPeaks <- markerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
markerPeaks@elementMetadata
#
pheatmap::pheatmap(markerPeaks@assays@data$Mean[c(rownames(as.data.frame(markerList$B_cells))[1:100],
                                                  rownames(as.data.frame(markerList$cDCs))[1:100],
                                                  rownames(as.data.frame(markerList$CLP))[1:100],
                                                  rownames(as.data.frame(markerList$Erythroid_cells))[1:100],
                                                  rownames(as.data.frame(markerList$HSPCs))[1:100],
                                                  rownames(as.data.frame(markerList$HSPCs_MPO))[1:100],
                                                  rownames(as.data.frame(markerList$Kupffer_cells))[1:100],
                                                  rownames(as.data.frame(markerList$Macrophages))[1:100],
                                                  rownames(as.data.frame(markerList$Mast_cells))[1:100],
                                                  rownames(as.data.frame(markerList$Megakaryocytes))[1:100],
                                                  rownames(as.data.frame(markerList$Microglia))[1:100],
                                                  rownames(as.data.frame(markerList$Monocytes))[1:100],
                                                  rownames(as.data.frame(markerList$NK_cells))[1:100],
                                                  rownames(as.data.frame(markerList$Osteoclasts))[1:100],
                                                  rownames(as.data.frame(markerList$pDCs))[1:100],
                                                  rownames(as.data.frame(markerList$T_cells))[1:100]),],
                   scale = "row",
                   cluster_rows = F,cluster_cols = F)
#
#
#
#
#
#
save.image("CS_Human_Immune.RData")
#
Hs_markerPeaks <- markerPeaks
Hs_markerList <- markerList
#
save(Hs_markerPeaks,file = "./Hs_markerPeaks.RData")
#save(Hs_markerList,file = "./Marker_heatmap/Hs_markerList.RData")
#
hs_orth_gr
#
#
#



#
#
features =  getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid_cells"]
#
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "HBB", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,  
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = GRanges(hs_orth_peak_df_co_celltype_sorted))
grid::grid.newpage()
grid::grid.draw(p$HBB)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "ARMC2", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,  
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = GRanges(hs_orth_peak_df_co_celltype_sorted))
grid::grid.newpage()
grid::grid.draw(p$ARMC2)
#
min(table(projHeme1_keep$ID_defined2))
#
#
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "HBB", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,  
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = GRanges(subset(hs_orth_peak_df, peak_id %in% co_tmp )))
grid::grid.newpage()
grid::grid.draw(p$HBB)
#
############ 
### B_cells
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$B_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$B_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CCND3", 
                      useGroups = c("B_cells","Kupffer_cells","Macrophages","Erythroid_cells","Megakaryocytes"), 
                      upstream = 140000,downstream = 20000,scCellsMax = 660,tileSize = 50,
                      minCells = 660,maxCells = 3500,ylim = c(0.05,0.9999999),features = GRanges(subset(hs_orth_peak_df, peak_id %in% co_tmp )))
grid::grid.newpage()
grid::grid.draw(p$CCND3)
#
#
############ 
### Erythroid_cells
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Erythroid_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Erythroid_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "HBB",, 
                      upstream = 40000,downstream = 90000,scCellsMax = 660,tileSize = 50,
                      minCells = 660,maxCells = 3500,ylim = c(0.05,0.9999999),features = GRanges(subset(hs_orth_peak_df, peak_id %in% co_tmp )))
grid::grid.newpage()
grid::grid.draw(p$HBB)
#
#
#
addArchRThreads(threads = 8)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CD19",
                      pal = c("#FDCDAC","#CDE8B9","#DE77AE","#00BFFF","#C2A5CF","#32CD32","#3375EF"),
                      useGroups = c("B_cells","T_cells","Kupffer_cells","Macrophages","Microglia","Erythroid_cells","Megakaryocytes"), 
                      upstream = 16000,downstream = 86000,scCellsMax = 6600,tileSize = 50,
                      minCells = 780,maxCells = 3500,ylim = c(0.05,0.9999999))
grid::grid.newpage()
grid::grid.draw(p$CD19)
#
#