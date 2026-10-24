# cell type define 
# taking immune white blood cell as example
# using the data from Qiu Nature 2019 mouse organogenesis datasets for quick cell type label transfer
# and using the marker genes from panglaoDB cellmarker etc datasets
#
setwd("CellType_define/Immune_cells")
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
#
inputFiles <- c("White_blood_cells.bed.gz")
names(inputFiles) <- c("White_blood_cells")
#
#
#addArchRChrPrefix(chrPrefix = FALSE)
#
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  geneAnnotation = TreeShrew.geneAnnotation,
  genomeAnnotation = TreeShrew.genomeAnnotation,
  filterTSS = 1,
  filterFrags = 10,
  addTileMat = T,
  offsetPlus = 0,
  offsetMinus = 0,
  addGeneScoreMat = T,
  force = T
)
#
#
#
#

#
addArchRThreads(threads = 8)
#
#
proj_WBC <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "TreeShrew",
  geneAnnotation = TreeShrew.geneAnnotation,
  genomeAnnotation = TreeShrew.genomeAnnotation,
  copyArrows = T
)
#
proj_WBC <- addIterativeLSI(
  ArchRProj = proj_WBC,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 3, 
  sampleCellsPre = 50000,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(2,3), 
    sampleCells = 50000, 
    n.start = 10,
    maxClusters = 100), 
  corCutOff = 0.6,
  UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE,
                    fast_sgd = TRUE),
  nPlot = 50000,
  varFeatures = 150000, 
  dimsToUse = 1:50
)
#
#
#
#
proj_WBC <- addUMAP(
  ArchRProj = proj_WBC, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
#
#
#
proj_WBC <- addClusters(
  input = proj_WBC,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_1",
  resolution = 1,
  maxClusters = 100
)
#
#
proj_WBC <- addClusters(
  input = proj_WBC,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_2",
  resolution = 2,
  maxClusters = 100
)
#
#
proj_WBC <- addClusters(
  input = proj_WBC,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_3",
  resolution = 3,
  maxClusters = 100
)
#
plotEmbedding(ArchRProj = proj_WBC, colorBy = "cellColData", name = "Clusters_1", baseSize = 3, size = 0.01, embedding = "UMAP")+NoLegend()
#
proj_WBC$Barcodes <- unlist(strsplit(proj_WBC$cellNames,"#"))[seq(2,2*length(proj_WBC$cellNames),2)]
proj_WBC.meta <- as.data.frame(proj_WBC@cellColData)
proj_WBC.meta$cellNames <- rownames(proj_WBC.meta)
head(proj_WBC.meta)
#
#
x.num <- c(rep(1,25000000),rep(0,25000000),rep(2,25000000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(proj_WBC.meta$Barcodes)),
                      sample(x.num,length(proj_WBC.meta$Barcodes)),
                      sample(x.num,length(proj_WBC.meta$Barcodes)),
                      sample(x.num,length(proj_WBC.meta$Barcodes)),
                      sample(x.num,length(proj_WBC.meta$Barcodes)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- proj_WBC.meta$Barcodes
#
TS_immune.seurat.type <- CreateSeuratObject(counts = t(x.num.sample), assay = "SPATAC", project = "E17.5", min.cells = 0, min.features = 0)
TS_immune.seurat.type
TS_immune.seurat.type <- NormalizeData(TS_immune.seurat.type, normalization.method = "LogNormalize", scale.factor = 10000)
TS_immune.seurat.type <- FindVariableFeatures(TS_immune.seurat.type, selection.method = "vst", nfeatures = 4)
TS_immune.seurat.type <- ScaleData(TS_immune.seurat.type, features = rownames(TS_immune.seurat.type))
variablegene <- VariableFeatures(object = TS_immune.seurat.type)
TS_immune.seurat.type <- RunPCA(TS_immune.seurat.type, features = variablegene,npcs =2,reduction.name = "UMAP_SnapATAC2",reduction.key = "UMAP_SnapATAC2")
TS_immune.seurat.type <- RunPCA(TS_immune.seurat.type, features = variablegene,npcs =2,reduction.name = "UMAP_ArchR",reduction.key = "UMAP_ArchR")
TS_immune.seurat.type
#
#
s.UMAP <- proj_WBC@embeddings$UMAP$df
rownames(s.UMAP ) <- proj_WBC$Barcodes
colnames(s.UMAP ) <- c("UMAPArchR_1","UMAPArchR_2")
head(s.UMAP)
TS_immune.seurat.type@reductions$UMAP_ArchR@cell.embeddings <- as.matrix(s.UMAP)
#
DimPlot(TS_immune.seurat.type,label = T,reduction = "UMAP_ArchR")+NoLegend()
#
TS_immune.seurat.type$Clusters_1 <- as.character(proj_WBC$Clusters_1)
TS_immune.seurat.type$Clusters_2 <- as.character(proj_WBC$Clusters_2)
TS_immune.seurat.type$Clusters_3 <- as.character(proj_WBC$Clusters_3)
#
DimPlot(TS_immune.seurat.type,group.by = "Clusters_1",label = T,reduction = "UMAP_ArchR")+NoLegend()
#
#
#
#
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "ACTB", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$ACTB)
#
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "CD3E", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$CD3E)
#
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "CD19", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$CD19)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "CD79A", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$CD79A)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "Kit", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$KIT)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "PF4", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$PF4)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "CCR2", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$CCR2)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "MMP9", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$MMP9)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "BCL11B", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$BCL11B)
#
addArchRThreads(threads = 16)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "CD3G", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$CD3G)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "CD8A", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$CD8A)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_1", 
                      geneSymbol = "GAPDH", upstream = 350000,downstream = 450000)
grid::grid.newpage()
grid::grid.draw(p$GAPDH)
#
#
#
#
#
#
addArchRThreads(threads = 6)
markersGS <- getMarkerFeatures(
  ArchRProj = proj_WBC, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters_1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C26
markerList$C27
markerList$C20
markerList$C23
markerList$C7[1:10,]
markerList$C8[1:10,]
#
#
#
load("/home/sunkeyong/Multispecies_analysis/M_TreeShrew/Analysis_T0808/Five_part/SDpups_E14E15E16.DS.RData")
SDpups_E14E15E16.DS
#
SDpups_E14E15E16.DS@active.ident <- as.factor(SDpups_E14E15E16.DS$major_tra)
table(SDpups_E14E15E16.DS@active.ident)
SDpups_E14E15E16.DS.sub <- subset(SDpups_E14E15E16.DS,idents=c("B_cells","Definitive_erythroid","Megakaryocytes",
                                                               "Mast_cells","T_cells","Primitive_erythroid",
                                                               "White_blood_cells"))
SDpups_E14E15E16.DS.sub
#
rm(SDpups_E14E15E16.DS)
gc()
#
#
#
addArchRThreads(threads = 2)
#
save.image("MultiS_TS.RData")
#
proj_WBC <- addGeneIntegrationMatrix(
  ArchRProj = proj_WBC, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_E14E15E16",
  reducedDims = "IterativeLSI",
  seRNA = SDpups_E14E15E16.DS.sub,
  addToArrow = F,
  force= TRUE,
  sampleCellsATAC = 40000,
  sampleCellsRNA = 50477,
  dimsToUse = 1:50,
  scaleDims = NULL,
  corCutOff = 0.75,
  groupRNA = "cell_type_major_tra",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore")
#
#
save.image("MultiS_TS.RData")
#
proj_WBC$p.Group <- unlist(strsplit(proj_WBC$predictedGroup,"__"))[seq(1,2*length(proj_WBC$predictedGroup),2)]
proj_WBC$p.Lineage <- unlist(strsplit(proj_WBC$predictedGroup,"__"))[seq(2,2*length(proj_WBC$predictedGroup),2)]
#
#
plotEmbedding(ArchRProj = proj_WBC, colorBy = "cellColData", name = "predictedScore", baseSize = 3, size = 0.01, embedding = "UMAP")+NoLegend()
plotEmbedding(ArchRProj = proj_WBC, colorBy = "cellColData", name = "p.Group", baseSize = 3, size = 0.01, embedding = "UMAP")
plotEmbedding(ArchRProj = proj_WBC, colorBy = "cellColData", name = "p.Lineage", baseSize = 3, size = 0.01, embedding = "UMAP")
#
#
#
#
DimPlot(TS_immune.seurat.type,group.by = "Clusters_1",label = T,reduction = "UMAP_ArchR")+NoLegend()
TS_immune.seurat.type@active.ident <- as.factor(TS_immune.seurat.type$Clusters_1)
Doublet_cells <- Cells(subset(TS_immune.seurat.type,idents=c("C18","C20","C21","C11","C28")))
#
TS_immune.seurat.type$cell_kepp <- "cell_kepp"
TS_immune.seurat.type@active.ident <- as.factor(TS_immune.seurat.type$cell_kepp)
table(TS_immune.seurat.type@active.ident)
Idents(TS_immune.seurat.type,cells=Doublet_cells) <- "Doublet_cells"
#
TS_immune.seurat.type$cell_kepp <- as.character(TS_immune.seurat.type@active.ident)
DimPlot(TS_immune.seurat.type,group.by = "cell_kepp",label = T,reduction = "UMAP_ArchR")+NoLegend()
#
TS_immune.seurat.type.meta <- as.data.frame(TS_immune.seurat.type@meta.data)
cell_keeps <- rownames(TS_immune.seurat.type.meta[which(TS_immune.seurat.type.meta$cell_kepp=="cell_kepp"),])
length(cell_keeps)
table(TS_immune.seurat.type$cell_kepp)
#
#
#
proj_WBC.metaB <- as.data.frame(proj_WBC@cellColData)
proj_WBC.metaB$cellNames <- proj_WBC$cellNames
head(proj_WBC.metaB)
#
rownames(proj_WBC.metaB) <- proj_WBC.metaB$Barcodes
proj_WBC.metaB <- proj_WBC.metaB[cell_keeps,]
cell_keeps.ArchR <- proj_WBC.metaB$cellNames
#
proj_WBC.keep <- proj_WBC[cell_keeps.ArchR,]
proj_WBC.keep
#
#
#
#
save.image("MultiS_TS.RData")
#
#
proj_WBC.keep <- addIterativeLSI(
  ArchRProj = proj_WBC.keep,
  useMatrix = "TileMatrix", 
  name = "IterativeLSIB", 
  iterations = 3, 
  sampleCellsPre = 50000,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(2,3), 
    sampleCells = 50000, 
    n.start = 10,
    maxClusters = 100), 
  corCutOff = 0.6,
  UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE,
                    fast_sgd = TRUE),
  nPlot = 50000,
  varFeatures = 150000, 
  dimsToUse = 1:50)
#
save.image("MultiS_TS.RData")
#
#
#
proj_WBC.keep <- addUMAP(
  ArchRProj = proj_WBC.keep, 
  reducedDims = "IterativeLSIB", 
  name = "UMAPB", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
#
save.image("MultiS_TS.RData")
#
#

#
proj_WBC.keep <- addClusters(
  input = proj_WBC.keep,
  reducedDims = "IterativeLSIB",
  method = "Seurat",
  name = "Clusters_1B",
  resolution = 1,
  maxClusters = 100
)
#
save.image("MultiS_TS.RData")
#
#
#
proj_WBC.keep
SDpups_E14E15E16.DS.sub
#
proj_WBC.keep <- addGeneIntegrationMatrix(
  ArchRProj = proj_WBC.keep, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_E14E15E16B",
  reducedDims = "IterativeLSIB",
  seRNA = SDpups_E14E15E16.DS.sub,
  addToArrow = F,
  force= TRUE,
  sampleCellsATAC = 70000,
  sampleCellsRNA = 50477,
  nGenes = 3000,
  dimsToUse = 1:50,
  scaleDims = NULL,
  corCutOff = 0.75,
  groupRNA = "cell_type_major_tra",
  nameCell = "predictedCellB",
  nameGroup = "predictedGroupB",
  nameScore = "predictedScoreB")
#
#
save.image("MultiS_TS.RData")
#
plotEmbedding(ArchRProj = proj_WBC.keep, colorBy = "cellColData", name = "Clusters_1B", baseSize = 3, size = 0.01, embedding = "UMAPB")
plotEmbedding(ArchRProj = proj_WBC.keep, colorBy = "cellColData", name = "Clusters_2B", baseSize = 3, size = 0.01, embedding = "UMAPB")+NoLegend()
#
plotEmbedding(ArchRProj = proj_WBC.keep, colorBy = "cellColData", name = "predictedScore", baseSize = 3, size = 0.01, embedding = "UMAPB")+NoLegend()
plotEmbedding(ArchRProj = proj_WBC.keep, colorBy = "cellColData", name = "p.Group", baseSize = 3, size = 0.01, embedding = "UMAPB")
plotEmbedding(ArchRProj = proj_WBC.keep, colorBy = "cellColData", name = "p.Lineage", baseSize = 3, size = 0.01, embedding = "UMAPB")
#
proj_WBC.keep$p.GroupB <- unlist(strsplit(proj_WBC.keep$predictedGroupB,"__"))[seq(1,2*length(proj_WBC.keep$predictedGroupB),2)]
proj_WBC.keep$p.LineageB <- unlist(strsplit(proj_WBC.keep$predictedGroupB,"__"))[seq(2,2*length(proj_WBC.keep$predictedGroupB),2)]
#
plotEmbedding(ArchRProj = proj_WBC.keep, colorBy = "cellColData", name = "p.GroupB", baseSize = 3, size = 0.01, embedding = "UMAPB")
plotEmbedding(ArchRProj = proj_WBC.keep, colorBy = "cellColData", name = "p.LineageB", baseSize = 3, size = 0.01, embedding = "UMAPB")
#
#
#
VlnPlot(SDpups_E14E15E16.DS.sub,features = "MPO",group.by = "cell_type",pt.size = 0)+NoLegend()
VlnPlot(SDpups_E14E15E16.DS.sub,features = "EPX",group.by = "cell_type",pt.size = 0)+NoLegend()
#
p <- plotBrowserTrack(ArchRProj = proj_WBC.keep, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_2B", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "MZB1", upstream = 100000,downstream = 100000)
grid::grid.newpage()
grid::grid.draw(p$MZB1)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC.keep, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_2B", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "CD3E", upstream = 100000,downstream = 100000)
grid::grid.newpage()
grid::grid.draw(p$CD3E)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC.keep, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_2B", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "CD8B", upstream = 100000,downstream = 100000)
grid::grid.newpage()
grid::grid.draw(p$CD8B)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC.keep, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_2B", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "IL2RB", upstream = 100000,downstream = 100000)
grid::grid.newpage()
grid::grid.draw(p$IL2RB)
#
#
p <- plotBrowserTrack(ArchRProj = proj_WBC.keep, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_2B", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "S100B", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$S100B)
#
#
TS_immune.seurat.type
DimPlot(TS_immune.seurat.type,group.by = "cell_kepp",label = T,reduction = "UMAP_ArchR")+NoLegend()
#
TS_immune.seurat.type.keep <- subset(TS_immune.seurat.type,cells=proj_WBC.keep$Barcodes)
TS_immune.seurat.type.keep$Clusters_2B <- proj_WBC.keep$Clusters_2B
TS_immune.seurat.type.keep$p.GroupB <- proj_WBC.keep$p.GroupB
#
TS_immune.seurat.type.keep <- RunPCA(TS_immune.seurat.type.keep, features = variablegene,npcs =2,reduction.name = "UMAP_ArchRB",reduction.key = "UMAP_ArchRB")
TS_immune.seurat.type.keep
#
#
s.UMAP <- proj_WBC.keep@embeddings$UMAPB$df
rownames(s.UMAP ) <- proj_WBC.keep$Barcodes
colnames(s.UMAP ) <- c("UMAPArchRB_1","UMAPArchRB_2")
head(s.UMAP)
TS_immune.seurat.type.keep@reductions$UMAP_ArchRB@cell.embeddings <- as.matrix(s.UMAP)
#
#
TS_immune.seurat.type.keep
#
DimPlot(TS_immune.seurat.type.keep,group.by = "Clusters_2B",label = T,reduction = "UMAP_ArchRB")+NoLegend()
DimPlot(TS_immune.seurat.type.keep,group.by = "Clusters_3",label = T,reduction = "UMAP_ArchRB")+NoLegend()
DimPlot(TS_immune.seurat.type.keep,group.by = "p.GroupB",label = T,reduction = "UMAP_ArchRB")+NoLegend()
#
#
View(as.data.frame(table(TS_immune.seurat.type.keep$Clusters_2B,TS_immune.seurat.type.keep$p.GroupB)))
table(TS_immune.seurat.type.keep$p.GroupB)
#
#
#
#
TS_immune.seurat.type.keep_backup <- TS_immune.seurat.type.keep
#
defined_ID.list <- read.csv("./Anno.txt",sep = "\t")
head(defined_ID.list)
# Clusters_2B	sub_ID
# C1	T cells
# C2	T cells
# C3	Activated T cells
# C4	Activated T cells
# C5	B cell progenitors
# C6	B cell progenitors
# C7	B cell progenitors
# C8	B cell progenitors
# C9	B cells
# C10	B cells
# C11	B cells
# C12	Megakaryocytes
# C13	Megakaryocytes
# C14	Megakaryocyte-erythroid progenitor
# C15	Erythroid progenitor
# C16	Osteoclasts
# C17	Kupffer cells
# C18	Kupffer cells
# C19	Kupffer cells
# C20	Microglia
# C21	Border-associated macrophages
# C22	Border-associated macrophages
# C23	Border-associated macrophages
# C24	Border-associated macrophages
# C25	Border-associated macrophages
# C26	Plasmacytoid dendritic cells
# C27	Monocytes
# C28	Conventional dendritic cells
# C29	NK cells
# C30	Monocytic myeloid-derived suppressor cells
# C31	Eosinophils
# C32	Granulocytes
# C33	Mast cells
# C34	Common myeloid progenitors
# C35	HSPCs
# C36	Common lymphoid progenitor
# C37	C37
# C38	Pre-T cells
# C39	NK cells
#
TS.snap.seurat.type.meta <- as.data.frame(TS_immune.seurat.type.keep@meta.data)
TS.snap.seurat.type.meta$Barcodes <- rownames(TS.snap.seurat.type.meta)
head(TS.snap.seurat.type.meta)
#
TS.snap.seurat.type.meta <- merge(TS.snap.seurat.type.meta,defined_ID.list,by="Clusters_2B")
head(TS.snap.seurat.type.meta)
rownames(TS.snap.seurat.type.meta) <- TS.snap.seurat.type.meta$Barcodes
TS.snap.seurat.type.meta <- TS.snap.seurat.type.meta[Cells(TS_immune.seurat.type.keep),]
#
TS_immune.seurat.type.keep$Defined_ID <- TS.snap.seurat.type.meta$sub_ID
#
DimPlot(TS_immune.seurat.type.keep,group.by = "Defined_ID",label = T,reduction = "UMAP_ArchRB")+NoLegend()
#
#
TS_immune.seurat.type.keep@active.ident <- as.factor(TS_immune.seurat.type.keep$Defined_ID)
TS_immune.seurat.type.keep2 <- subset(TS_immune.seurat.type.keep,idents="C37",invert=T)
TS_immune.seurat.type.keep2 <- RenameIdents(TS_immune.seurat.type.keep2,"B cell progenitors"="Pre-B cells")
#
TS_immune.seurat.type.keep2$Defined_ID <- as.character(TS_immune.seurat.type.keep2@active.ident)
DimPlot(TS_immune.seurat.type.keep2,group.by = "Defined_ID",label = T,reduction = "UMAP_ArchRB")+NoLegend()
#
#
table(TS_immune.seurat.type.keep2$Defined_ID)
#
#
#
DimPlot(TS_immune.seurat.type.keep2,group.by = "Defined_ID",
        cols = TS_immune.color,
        label = T,reduction = "UMAP_ArchRB")+NoLegend()
#
#
TS_immune.seurat.type.final <- TS_immune.seurat.type.keep2
#
levels(TS_immune.seurat.type.final) <- c("Activated T cells",                         
                                    "B cells",                                   
                                    "Border-associated macrophages",
                                    "Common lymphoid progenitor",                
                                    "Common myeloid progenitors",                
                                    "Conventional dendritic cells",              
                                    "Eosinophils",                               
                                    "Erythroid progenitor",                      
                                    "Granulocytes",                              
                                    "HSPCs",                                     
                                    "Kupffer cells",                            
                                    "Mast cells",                                
                                    "Megakaryocyte-erythroid progenitor",        
                                    "Megakaryocytes",                            
                                    "Microglia",                                 
                                    "Monocytes",                                 
                                    "Monocytic myeloid-derived suppressor cells",
                                    "NK cells",                                  
                                    "Osteoclasts",                               
                                    "Plasmacytoid dendritic cells",              
                                    "Pre-B cells",                               
                                    "Pre-T cells",                               
                                    "T cells")
#
#
#
TS_immune.color <- c("#0099FF", "#6A3D9A", "#B3E2CD", "#E6AB02", "#48D1CC", ##1-5
                     "#00FFFF", "#EE82EE", "#FDCDAC", "#99CCCC", "#ABD9E9", ##6-10
                     "#1A9850", "#9933CC", "#B3EE3A", "#87CEFA", "#E08214", ##11-15
                     "#9970AB", "#DE77AE", "#F1E2CC", "#0000FF", "#9E0142", ##16-20
                     "#66C2A5", "#FF00CC", "#276419")
#
#
plotEmbedding(ArchRProj = proj_WBC.keep, colorBy = "cellColData", name = "Clusters_2B", baseSize = 3, size = 0.01, embedding = "UMAPB")+NoLegend()
proj_WBC.keep
#
idxPass1 <- which(proj_WBC.keep$Clusters_2B == "C2")
idxPass2 <- which(proj_WBC.keep$Clusters_2B == "C4")
#
cellsPass <- proj_WBC.keep$cellNames[c(idxPass2,idxPass1)]
proj_WBC.keepx <- proj_WBC.keep[cellsPass, ]
#
plotEmbedding(ArchRProj = proj_WBC.keepx, colorBy = "cellColData", name = "Clusters_2B", baseSize = 3, size = 0.01, embedding = "UMAPB")+NoLegend()
#
addArchRThreads(threads = 6)
markersGS <- getMarkerFeatures(
  ArchRProj = proj_WBC.keepx, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters_2B",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C2[1:10,]
markerList$C4[1:10,]
#
p <- plotBrowserTrack(ArchRProj = proj_WBC.keepx, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_2B", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "GZMAli1", upstream = 100000,downstream = 100000)
grid::grid.newpage()
grid::grid.draw(p$GZMAli1)
#
p <- plotBrowserTrack(ArchRProj = proj_WBC, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Clusters_2B", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "MPO", upstream = 100000,downstream = 100000)
grid::grid.newpage()
grid::grid.draw(p$MPO)
#








####################################################################################################
####################################################################################################
####################################################################################################
#
#
# plot figure
# 
DimPlot(TS_immune.seurat.type.final,group.by = "p.GroupB",label = T,reduction = "UMAP_ArchRB")+NoLegend()
DimPlot(TS_immune.seurat.type.final,group.by = "Defined_ID",label = T,reduction = "UMAP_ArchRB")+NoLegend()
#
levels(TS_immune.seurat.type.final) <- c("Activated T cells",                         
                                         "B cells",                                   
                                         "Border-associated macrophages",
                                         "Common lymphoid progenitor",                
                                         "Common myeloid progenitors",                
                                         "Conventional dendritic cells",              
                                         "Eosinophils",                               
                                         "Erythroid progenitor",                      
                                         "Granulocytes",                              
                                         "HSPCs",                                     
                                         "Kupffer cells",                            
                                         "Mast cells",                                
                                         "Megakaryocyte-erythroid progenitor",        
                                         "Megakaryocytes",                            
                                         "Microglia",                                 
                                         "Monocytes",                                 
                                         "Monocytic myeloid-derived suppressor cells",
                                         "NK cells",                                  
                                         "Osteoclasts",                               
                                         "Plasmacytoid dendritic cells",              
                                         "Pre-B cells",                               
                                         "Pre-T cells",                               
                                         "T cells")
#
TS_immune.color <- c("#0099FF", "#6A3D9A", "#B3E2CD", "#E6AB02", "#48D1CC", ##1-5
                     "#00FFFF", "#EE82EE", "#FDCDAC", "#99CCCC", "#ABD9E9", ##6-10
                     "#1A9850", "#9933CC", "#B3EE3A", "#87CEFA", "#E08214", ##11-15
                     "#9970AB", "#DE77AE", "#F1E2CC", "#0000FF", "#9E0142", ##16-20
                     "#66C2A5", "#FF00CC", "#276419")
#
#
TS_immune.seurat.type.final@active.ident <- as.factor(TS_immune.seurat.type.final$Defined_ID)
p <- DimPlot(TS_immune.seurat.type.final,
             cols = TS_immune.color,label = F,repel = T,pt.size = 0.00003,reduction = "UMAP_ArchRB",
             shuffle=T,raster=F,raster.dpi = c(1512, 1512))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./Figures/TS_immune.anno.umap1.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
ggsave("./Figures/TS_immune.anno.umap2.png",plot=p,device="png",dpi=400,units = "cm",width = 15,height = 15)
ggsave("./Figures/TS_immune.anno.umap3.png",plot=p,device="png",dpi=400,units = "cm",width = 18,height = 18)
ggsave("./Figures/TS_immune.anno.umap4.png",plot=p,device="png",dpi=400,units = "cm",width = 22,height = 22)
#
#
#
#
#
TS_immune.seurat.type.final.meta <- as.data.frame(TS_immune.seurat.type.final@meta.data)
proj_WBC.keep.meta <- as.data.frame(proj_WBC.keep@cellColData)
proj_WBC.keep.meta$cellNames <- proj_WBC.keep$cellNames
head(proj_WBC.keep.meta)
rownames(proj_WBC.keep.meta) <- proj_WBC.keep.meta$Barcodes
#
projArchR_TS_Final <- proj_WBC.keep[proj_WBC.keep.meta[rownames(TS_immune.seurat.type.final.meta),]$cellNames,]
projArchR_TS_Final
#
TS_immune.seurat.type.final.meta2 <- TS_immune.seurat.type.final.meta[projArchR_TS_Final$Barcodes,]
projArchR_TS_Final$Defined_ID <- as.character(TS_immune.seurat.type.final.meta2$Defined_ID)
#
plotEmbedding(ArchRProj = projArchR_TS_Final, colorBy = "cellColData", name = "Defined_ID", baseSize = 3, size = 0.01, embedding = "UMAPB")
#
#
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "MPO", upstream = 350000,downstream = 350000)
grid::grid.newpage()
grid::grid.draw(p$MPO)
#
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "ACTB", upstream = 100000,downstream = 100000)
grid::grid.newpage()
grid::grid.draw(p$ACTB)
#
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "MS4A1", upstream = 160000,downstream = 160000)
grid::grid.newpage()
grid::grid.draw(p$MS4A1)
#
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "HBE1", upstream = 100000,downstream = 20000)
grid::grid.newpage()
grid::grid.draw(p$HBE1)
#
#
save.image("MultiS_TS.RData")
save(TS_immune.seurat.type.final,file = "TS_immune.seurat.type.final.RData")
#
#
#
#
#
addArchRThreads(threads = 8)
#
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "MPO", upstream = 350000,downstream = 350000)
grid::grid.newpage()
grid::grid.draw(p$MPO)
##
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "PF4", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$PF4)
#
##
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "CD14", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$CD14)
#
##
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "CD14", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$CD14)
#
p <- plotBrowserTrack(ArchRProj = projArchR_TS_Final, maxCells = 1500,ylim = c(0.1,0.99999),groupBy = "Defined_ID", sizes = c(10, 1.5, 3, 1.5),
                      geneSymbol = "FCGR3A", upstream = 150000,downstream = 150000)
grid::grid.newpage()
grid::grid.draw(p$FCGR3A)
#
load("/home/sunkeyong/Multispecies_analysis/M_TreeShrew/Analysis_T0808/Subset_T250830/subset/Immune_cells/TS_immune.seurat.type.final.RData")
TS_immune.seurat.type.final
#
write.table(Cells(TS_immune.seurat.type.final),file = "/home/sunkeyong/Human_em/Cross_species/Immune_TreeShrew/bed/TS_WBC.Barcode.csv",quote = F, sep = "\t",col.names = F,row.names = F )
#
