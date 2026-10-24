setwd("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step2_ArchR/Treeshrew_ArchR")
#
#
##
.libPaths("/home/sunkeyong/R/x86_64-pc-linux-gnu-library/4.4")
#
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
.libPaths("/home/sunkeyong/R/x86_64-pc-linux-gnu-library/4.4")
##
load("/home/sunkeyong/softeware/ref_tree_shrew/ArchR_annotation_Storage/TreeShrew.genomeAnnotation.RData")
load("/home/sunkeyong/softeware/ref_tree_shrew/ArchR_annotation_Storage/TreeShrew.geneAnnotation.RData")
#
library(BSgenome.tupBel1.UCSC.TreeShrew)
TreeShrew.genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.tupBel1.UCSC.TreeShrew)
#
addArchRThreads(threads = 2) 
#
load("/home/sunkeyong/Human_em/Cross_species/Immune_TreeShrew/TS_immune.seurat.type.final.RData")
TS_immune.seurat.type.final
TS_immune.seurat.type.final.meta <- as.data.frame(TS_immune.seurat.type.final@meta.data)
#
ArrowFiles <- c("/home/sunkeyong/Human_em/Cross_species/Immune_TreeShrew/ArchR/TS_Immune.arrow")
addArchRThreads(threads = 8)
#
#
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "TreeShrew",
  geneAnnotation = TreeShrew.geneAnnotation,
  genomeAnnotation = TreeShrew.genomeAnnotation,
  copyArrows = T
)
#
#
save.image("CS_TreeShrew_Immune.RData")
#
#
projHeme1
TS_immune.seurat.type.final
table(TS_immune.seurat.type.final$Defined_ID)
#
#
HT_keep_celltypes2 <- c("Pre-B cells",
                        "Border-associated macrophages",
                        "Common lymphoid progenitor",
                        "Erythroid progenitor",
                        "HSPCs",
                        "Monocytic myeloid-derived suppressor cells",
                        "Kupffer cells",
                        "Mast cells",
                        "Megakaryocytes",
                        "Microglia",
                        "Monocytes",
                        "NK cells",
                        "T cells",
                        "Pre-T cells",
                        "Activated T cells")
##
#
table(TS_immune.seurat.type.final@active.ident)
#
#
TS_immune.seurat.type.final_keep <- subset(TS_immune.seurat.type.final,idents=HT_keep_celltypes2)
DimPlot(TS_immune.seurat.type.final,reduction = "UMAP_ArchRB",label = T)+NoLegend()
DimPlot(TS_immune.seurat.type.final_keep,reduction = "UMAP_ArchRB",label = T)+NoLegend()
#
#
#
CS_TShrew_Immune.seurat.type_keep2 <- TS_immune.seurat.type.final_keep
#
#
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"T cells"="T_cells")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Pre-T cells"="T_cells")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Activated T cells"="T_cells")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Border-associated macrophages"="Macrophages")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Common lymphoid progenitor"="CLP")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Erythroid progenitor"="Erythroid_cells")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"HSPCs"="HSPCs")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Kupffer cells"="Kupffer_cells")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Mast cells"="Mast_cells")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Megakaryocytes"="Megakaryocytes")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Microglia"="Microglia")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Monocytes"="Monocytes")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Monocytic myeloid-derived suppressor cells"="HSPCs_MPO")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"NK cells"="NK_cells")
CS_TShrew_Immune.seurat.type_keep2 <- RenameIdents(CS_TShrew_Immune.seurat.type_keep2,"Pre-B cells"="B_cells")
#
DimPlot(CS_TShrew_Immune.seurat.type_keep2,reduction = "UMAP_ArchR",group.by = "Defined_ID",
        label = T,repel = F,raster = F,raster.dpi=c(1512, 1512))+NoLegend()
#
CS_TShrew_Immune.seurat.type_keep2$ID_defined2 <- as.character(CS_TShrew_Immune.seurat.type_keep2@active.ident)
#
DimPlot(CS_TShrew_Immune.seurat.type_keep2,reduction = "UMAP_ArchR",group.by = "ID_defined2",
        label = T,repel = F,raster = F,raster.dpi=c(1512, 1512))+NoLegend()
#
#
#
#











#
peak_df <- read.table("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/orth.ts.raw.sort.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(peak_df) <- c("chr", "start", "end", "peak_id")
# 转换为 GRanges
ts_orth_gr <- GRanges(
  seqnames = peak_df$chr,
  ranges = IRanges(start = peak_df$start, end = peak_df$end),
  peak_id = peak_df$peak_id
)
#
addArchRThreads(threads = 2) 
#
#
#
projHeme1$barcodes <- unlist(strsplit(projHeme1$cellNames,"#"))[seq(2,2*length(projHeme1$cellNames),2)]
projHeme1.meta <- as.data.frame(projHeme1@cellColData)
projHeme1.meta$cellNames <- projHeme1$cellNames
rownames(projHeme1.meta) <- projHeme1.meta$barcodes
projHeme1.meta_tmp <- projHeme1.meta[Cells(CS_TShrew_Immune.seurat.type_keep2),]
projHeme1_keep <- projHeme1[projHeme1.meta_tmp$cellNames,]
#
#
CS_TShrew_Immune.seurat.type_keep2_meta <- as.data.frame(CS_TShrew_Immune.seurat.type_keep2@meta.data)
CS_TShrew_Immune.seurat.type_keep2_meta <- CS_TShrew_Immune.seurat.type_keep2_meta[projHeme1_keep$barcodes,]
projHeme1_keep$ID_defined2 <- as.character(CS_TShrew_Immune.seurat.type_keep2_meta$ID_defined2)
#
addArchRThreads(threads = 12) 
#
projHeme1_keep <- addPeakSet(projHeme1_keep,peakSet = ts_orth_gr,force=TRUE)
projHeme1_keep <- addPeakMatrix(projHeme1_keep,force = T)
getAvailableMatrices(projHeme1_keep)
#
tail(ts_orth_gr)
#
save.image("CS_TreeShrew_Immune.RData")
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "HBE1", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = getPeakSet(projHeme1_keep))
grid::grid.newpage()
grid::grid.draw(p$HBE1)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CD79A", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = getPeakSet(projHeme1_keep))
grid::grid.newpage()
grid::grid.draw(p$CD79A)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CD24", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = getPeakSet(projHeme1_keep))
grid::grid.newpage()
grid::grid.draw(p$CD24)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CSF1R", 
                      upstream = 160000,downstream = 160000,scCellsMax = 500,
                      minCells = 200,maxCells = 3500,ylim = c(0.01,0.9999999),features = getPeakSet(projHeme1_keep))
grid::grid.newpage()
grid::grid.draw(p$CSF1R)
##
#
markerPeaks <- getMarkerFeatures(
  ArchRProj = projHeme1_keep, 
  useMatrix = "PeakMatrix", 
  groupBy = "ID_defined2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#
save.image("CS_TreeShrew_Immune.RData")
#
#
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList$B_cells
markerList$Macrophages
markerList$CLP
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
heatmapPeaks <- markerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
markerPeaks@elementMetadata
#
pheatmap::pheatmap(markerPeaks@assays@data$Mean[c(rownames(as.data.frame(markerList$B_cells))[1:100],
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
                                                  rownames(as.data.frame(markerList$T_cells))[1:100]),],
                   scale = "row",
                   cluster_rows = F,cluster_cols = F)
#

#

#
#
#
TS_markerPeaks <- markerPeaks
#TS_markerList <- markerList
#
save(TS_markerPeaks,file = "./TS_markerPeaks.RData")
#save(TS_markerList,file = "./Marker_heatmap/TS_markerList.RData")
#
save.image("CS_TreeShrew_Immune.RData")
#
#
#
#




##############################
##############################
table(projHeme1_keep$ID_defined2)
min(table(projHeme1_keep$ID_defined2))
#
### B_cells
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$B_cells)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$B_cells)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CCND3", 
                      useGroups = c("B_cells","Kupffer_cells","Macrophages","Erythroid_cells","Megakaryocytes"), 
                      upstream = 20000,downstream = 140000,scCellsMax = 900,tileSize = 50,
                      minCells = 900,maxCells = 3500,ylim = c(0.05,0.9999999),features = GRanges(subset(ts_orth_peak_df, peak_id %in% co_tmp )))
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
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CCND3", 
                      useGroups = c("B_cells","Kupffer_cells","Macrophages","Erythroid_cells","Megakaryocytes"), 
                      upstream = 20000,downstream = 140000,scCellsMax = 900,tileSize = 50,
                      minCells = 900,maxCells = 3500,ylim = c(0.05,0.9999999),features = GRanges(subset(ts_orth_peak_df, peak_id %in% co_tmp )))
grid::grid.newpage()
grid::grid.draw(p$CCND3)
#
############ 
### Megakaryocytes
ts_tmp <- subset(ts_orth_peak_df, ts_orth_idx %in% rownames(as.data.frame(TS_markerList$Megakaryocytes)) )$peak_id
hs_tmp <- subset(hs_orth_peak_df, hs_orth_idx %in% rownames(as.data.frame(Hs_markerList$Megakaryocytes)) )$peak_id
co_tmp <- intersect(ts_tmp,hs_tmp);ts_sp_tmp <- setdiff(ts_tmp,co_tmp);hs_sp_tmp <- setdiff(hs_tmp,co_tmp)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CCND3", 
                      useGroups = c("B_cells","Kupffer_cells","Macrophages","Erythroid_cells","Megakaryocytes"), 
                      upstream = 20000,downstream = 140000,scCellsMax = 900,tileSize = 50,
                      minCells = 900,maxCells = 3500,ylim = c(0.05,0.9999999),features = GRanges(subset(ts_orth_peak_df, peak_id %in% co_tmp )))
grid::grid.newpage()
grid::grid.draw(p$CCND3)
#
#
#
addArchRThreads(threads = 8)
#
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CD19", 
                      upstream = 20000,downstream = 140000,scCellsMax = 9000,tileSize = 50,
                      minCells = 900,maxCells = 3500,ylim = c(0.05,0.9999999))
grid::grid.newpage()
grid::grid.draw(p$CD19)
#
p <- plotBrowserTrack(ArchRProj = projHeme1_keep, groupBy = "ID_defined2", geneSymbol = "CD19", 
                      pal = c("#FDCDAC","#CDE8B9","#DE77AE","#00BFFF","#C2A5CF","#32CD32","#3375EF"),
                      useGroups = c("B_cells","T_cells","Kupffer_cells","Macrophages","Microglia","Erythroid_cells","Megakaryocytes"), 
                      upstream = 16000,downstream = 66000,scCellsMax = 9000,tileSize = 50,
                      minCells = 900,maxCells = 3500,ylim = c(0.05,0.9999999))
grid::grid.newpage()
grid::grid.draw(p$CD19)
#
#



########################
########################
########################
########################
########################
########################
########################
########################
########################
########################
load("/home/sunkeyong/Human_em/Cross_species/Immune_TreeShrew/TS_immune.seurat.type.final.RData")
#
TS_immune.seurat.type.final
#
#
table(TS_immune.seurat.type.final@active.ident)
#
#
DimPlot(TS_immune.seurat.type.final,
        reduction = "UMAP_ArchRB",label = T,raster = F,pt.size = 0.00001,alpha = 0.5)+NoLegend()
#
TS_immune.seurat.type.final <- RenameIdents(TS_immune.seurat.type.final,"Activated T cells"="T cells")
TS_immune.seurat.type.final <- RenameIdents(TS_immune.seurat.type.final,"Pre-T cells"="T cells")
#
#
TS_immune.seurat.type.final$Defined_ID2 <- as.character(TS_immune.seurat.type.final@active.ident)
TS_immune.seurat.type.final@active.ident <- as.factor(TS_immune.seurat.type.final$Defined_ID2)
table(TS_immune.seurat.type.final@active.ident )
#
color_use <- c("#DDAD3B","#00BFFF","#9970AB","#B15928",## 1-4
               "#87CEFA","#3288BD","#32CD32","#48D1CC",## 5-8
               "#FF33CC","#DE77AE","#B3CDE3","#F4CAE4",## 9-12
               "#3375EF","#C2A5CF","#EE82EE","#FF7F00",## 13-16
               "#9E0142","#6495ED","#00FFFF","#FDCDAC",## 17-20
               "#CDE8B9")
#
DimPlot(TS_immune.seurat.type.final,cols = color_use,
        reduction = "UMAP_ArchRB",label = T,raster = F,pt.size = 0.00001,alpha = 0.5)+NoLegend()
#
#
p <- DimPlot(TS_immune.seurat.type.final,cols = color_use,label = F,
             repel = T,pt.size = 0.00000000002,
             reduction = "UMAP_ArchRB",
             shuffle=T,raster=F,alpha = 0.3,
             raster.dpi = c(3312, 3312))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
#
ggsave("../Fig_save/TS.immune.UMAPP.stage.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
#
#
#
load("/home/sunkeyong/Human_em/Cross_species/Immune_Human/H_WBC_Ery.seurat.type.final.RData")
#
H_WBC_Ery.seurat.type.final
#
DimPlot(H_WBC_Ery.seurat.type.final,reduction = "UMAP_ArchR2",label=T)+NoLegend()
#
table(H_WBC_Ery.seurat.type.final@active.ident)
H_WBC_Ery.seurat.type.final$ID_use <- as.character(H_WBC_Ery.seurat.type.final@active.ident)
#
H_WBC_Ery.seurat.type.final@active.ident <- as.factor(H_WBC_Ery.seurat.type.final$ID_use)
table(H_WBC_Ery.seurat.type.final@active.ident)
#
# H_WBC_Ery.seurat.type.final <- RenameIdents(H_WBC_Ery.seurat.type.final,"Terminal megakaryocyte"="Megakaryocyte")
#
#
color_use_human <- c("#FDCDAC","#00BFFF","#9970AB","#87CEFA",
                     "#32CD32","#FF33CC","#B15928","#DE77AE",
                     "#B3CDE3","#3375EF","#C2A5CF","#EE82EE",
                     "#9E0142","#6495ED","#00FFFF","#CDE8B9","#3375EF")
#
DimPlot(H_WBC_Ery.seurat.type.final,cols = color_use_human,
        reduction = "UMAP_ArchR2",label = T,raster = F,pt.size = 0.00001,alpha = 1)+NoLegend()
#
#
p <- DimPlot(H_WBC_Ery.seurat.type.final,cols = color_use_human,label = F,
             repel = T,pt.size = 0.00000000002,
             reduction = "UMAP_ArchR2",
             shuffle=T,raster=F,alpha = 0.35,
             raster.dpi = c(3312, 3312))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
#
ggsave("../Fig_save/HE.immune.UMAPP.stage.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
#
#
#
HEi_color <- names(table(H_WBC_Ery.seurat.type.final@active.ident))
HEi_color
#
HEi_color <- as.data.frame(cbind(HEi_color,color_use_human))
head(HEi_color)
colnames(HEi_color) <- c("cell_type","color")
#
rownames(HEi_color) <- HEi_color$cell_type
#
HEi_color
#
#



#