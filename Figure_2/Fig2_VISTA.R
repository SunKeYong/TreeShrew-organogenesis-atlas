# 
setwd("VISTA_Analysis")
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
#

addArchRThreads(threads = 10) 
#

names(inputFiles) <- c("C1.1",
"C10.1",
"C11.1",
"C11.2",
"C11.3",
"C11.4",
"C12.1",
"C12.10",
"C12.11",
"C12.12",
"C12.2",
"C12.3",
"C12.4",
"C12.5",
"C12.6",
"C12.7",
"C12.8",
"C12.9",
"C13.1",
"C14.1",
"C14.10",
"C14.11",
"C14.12",
"C14.13",
"C14.14",
"C14.15",
"C14.16",
"C14.17",
"C14.18",
"C14.19",
"C14.2",
"C14.20",
"C14.21",
"C14.22",
"C14.23",
"C14.24",
"C14.3",
"C14.4",
"C14.5",
"C14.6",
"C14.7",
"C14.8",
"C14.9",
"C15.1",
"C15.2",
"C16.1",
"C16.2",
"C16.3",
"C16.4",
"C16.5",
"C16.6",
"C16.7",
"C17.1",
"C17.2",
"C17.3",
"C17.4",
"C18.1",
"C18.10",
"C18.11",
"C18.12",
"C18.13",
"C18.14",
"C18.15",
"C18.16",
"C18.17",
"C18.18",
"C18.19",
"C18.2",
"C18.20",
"C18.3",
"C18.4",
"C18.5",
"C18.6",
"C18.7",
"C18.8",
"C18.9",
"C19.1",
"C2.1",
"C2.2",
"C2.3",
"C2.4",
"C2.5",
"C2.6",
"C2.7",
"C20.1",
"C21.1",
"C22.1",
"C23.1",
"C23.10",
"C23.11",
"C23.12",
"C23.13",
"C23.2",
"C23.3",
"C23.4",
"C23.5",
"C23.6",
"C23.7",
"C23.8",
"C23.9",
"C24.1",
"C24.10",
"C24.11",
"C24.12",
"C24.13",
"C24.14",
"C24.15",
"C24.16",
"C24.17",
"C24.18",
"C24.19",
"C24.2",
"C24.20",
"C24.21",
"C24.22",
"C24.23",
"C24.3",
"C24.4",
"C24.5",
"C24.6",
"C24.7",
"C24.8",
"C24.9",
"C3.1",
"C3.2",
"C4.1",
"C5.1",
"C5.10",
"C5.11",
"C5.12",
"C5.13",
"C5.14",
"C5.15",
"C5.16",
"C5.17",
"C5.18",
"C5.19",
"C5.2",
"C5.20",
"C5.21",
"C5.22",
"C5.23",
"C5.24",
"C5.25",
"C5.26",
"C5.27",
"C5.3",
"C5.4",
"C5.5",
"C5.6",
"C5.7",
"C5.8",
"C5.9",
"C6.1",
"C6.2",
"C7.1",
"C7.10",
"C7.2",
"C7.3",
"C7.4",
"C7.5",
"C7.6",
"C7.7",
"C7.8",
"C7.9",
"C8.1",
"C8.10",
"C8.11",
"C8.12",
"C8.13",
"C8.14",
"C8.15",
"C8.16",
"C8.17",
"C8.18",
"C8.19",
"C8.2",
"C8.20",
"C8.21",
"C8.22",
"C8.23",
"C8.24",
"C8.25",
"C8.26",
"C8.27",
"C8.28",
"C8.3",
"C8.4",
"C8.5",
"C8.6",
"C8.7",
"C8.8",
"C8.9",
"C9.1",
"C9.2",
"C9.3",
"C9.4")
inputFiles = c("../Downsample_bed_files/C1.1.downsample.bed.gz",
"../Downsample_bed_files/C10.1.downsample.bed.gz",
"../Downsample_bed_files/C11.1.downsample.bed.gz",
"../Downsample_bed_files/C11.2.downsample.bed.gz",
"../Downsample_bed_files/C11.3.downsample.bed.gz",
"../Downsample_bed_files/C11.4.downsample.bed.gz",
"../Downsample_bed_files/C12.1.downsample.bed.gz",
"../Downsample_bed_files/C12.10.downsample.bed.gz",
"../Downsample_bed_files/C12.11.downsample.bed.gz",
"../Downsample_bed_files/C12.12.downsample.bed.gz",
"../Downsample_bed_files/C12.2.downsample.bed.gz",
"../Downsample_bed_files/C12.3.downsample.bed.gz",
"../Downsample_bed_files/C12.4.downsample.bed.gz",
"../Downsample_bed_files/C12.5.downsample.bed.gz",
"../Downsample_bed_files/C12.6.downsample.bed.gz",
"../Downsample_bed_files/C12.7.downsample.bed.gz",
"../Downsample_bed_files/C12.8.downsample.bed.gz",
"../Downsample_bed_files/C12.9.downsample.bed.gz",
"../Downsample_bed_files/C13.1.downsample.bed.gz",
"../Downsample_bed_files/C14.1.downsample.bed.gz",
"../Downsample_bed_files/C14.10.downsample.bed.gz",
"../Downsample_bed_files/C14.11.downsample.bed.gz",
"../Downsample_bed_files/C14.12.downsample.bed.gz",
"../Downsample_bed_files/C14.13.downsample.bed.gz",
"../Downsample_bed_files/C14.14.downsample.bed.gz",
"../Downsample_bed_files/C14.15.downsample.bed.gz",
"../Downsample_bed_files/C14.16.downsample.bed.gz",
"../Downsample_bed_files/C14.17.downsample.bed.gz",
"../Downsample_bed_files/C14.18.downsample.bed.gz",
"../Downsample_bed_files/C14.19.downsample.bed.gz",
"../Downsample_bed_files/C14.2.downsample.bed.gz",
"../Downsample_bed_files/C14.20.downsample.bed.gz",
"../Downsample_bed_files/C14.21.downsample.bed.gz",
"../Downsample_bed_files/C14.22.downsample.bed.gz",
"../Downsample_bed_files/C14.23.downsample.bed.gz",
"../Downsample_bed_files/C14.24.downsample.bed.gz",
"../Downsample_bed_files/C14.3.downsample.bed.gz",
"../Downsample_bed_files/C14.4.downsample.bed.gz",
"../Downsample_bed_files/C14.5.downsample.bed.gz",
"../Downsample_bed_files/C14.6.downsample.bed.gz",
"../Downsample_bed_files/C14.7.downsample.bed.gz",
"../Downsample_bed_files/C14.8.downsample.bed.gz",
"../Downsample_bed_files/C14.9.downsample.bed.gz",
"../Downsample_bed_files/C15.1.downsample.bed.gz",
"../Downsample_bed_files/C15.2.downsample.bed.gz",
"../Downsample_bed_files/C16.1.downsample.bed.gz",
"../Downsample_bed_files/C16.2.downsample.bed.gz",
"../Downsample_bed_files/C16.3.downsample.bed.gz",
"../Downsample_bed_files/C16.4.downsample.bed.gz",
"../Downsample_bed_files/C16.5.downsample.bed.gz",
"../Downsample_bed_files/C16.6.downsample.bed.gz",
"../Downsample_bed_files/C16.7.downsample.bed.gz",
"../Downsample_bed_files/C17.1.downsample.bed.gz",
"../Downsample_bed_files/C17.2.downsample.bed.gz",
"../Downsample_bed_files/C17.3.downsample.bed.gz",
"../Downsample_bed_files/C17.4.downsample.bed.gz",
"../Downsample_bed_files/C18.1.downsample.bed.gz",
"../Downsample_bed_files/C18.10.downsample.bed.gz",
"../Downsample_bed_files/C18.11.downsample.bed.gz",
"../Downsample_bed_files/C18.12.downsample.bed.gz",
"../Downsample_bed_files/C18.13.downsample.bed.gz",
"../Downsample_bed_files/C18.14.downsample.bed.gz",
"../Downsample_bed_files/C18.15.downsample.bed.gz",
"../Downsample_bed_files/C18.16.downsample.bed.gz",
"../Downsample_bed_files/C18.17.downsample.bed.gz",
"../Downsample_bed_files/C18.18.downsample.bed.gz",
"../Downsample_bed_files/C18.19.downsample.bed.gz",
"../Downsample_bed_files/C18.2.downsample.bed.gz",
"../Downsample_bed_files/C18.20.downsample.bed.gz",
"../Downsample_bed_files/C18.3.downsample.bed.gz",
"../Downsample_bed_files/C18.4.downsample.bed.gz",
"../Downsample_bed_files/C18.5.downsample.bed.gz",
"../Downsample_bed_files/C18.6.downsample.bed.gz",
"../Downsample_bed_files/C18.7.downsample.bed.gz",
"../Downsample_bed_files/C18.8.downsample.bed.gz",
"../Downsample_bed_files/C18.9.downsample.bed.gz",
"../Downsample_bed_files/C19.1.downsample.bed.gz",
"../Downsample_bed_files/C2.1.downsample.bed.gz",
"../Downsample_bed_files/C2.2.downsample.bed.gz",
"../Downsample_bed_files/C2.3.downsample.bed.gz",
"../Downsample_bed_files/C2.4.downsample.bed.gz",
"../Downsample_bed_files/C2.5.downsample.bed.gz",
"../Downsample_bed_files/C2.6.downsample.bed.gz",
"../Downsample_bed_files/C2.7.downsample.bed.gz",
"../Downsample_bed_files/C20.1.downsample.bed.gz",
"../Downsample_bed_files/C21.1.downsample.bed.gz",
"../Downsample_bed_files/C22.1.downsample.bed.gz",
"../Downsample_bed_files/C23.1.downsample.bed.gz",
"../Downsample_bed_files/C23.10.downsample.bed.gz",
"../Downsample_bed_files/C23.11.downsample.bed.gz",
"../Downsample_bed_files/C23.12.downsample.bed.gz",
"../Downsample_bed_files/C23.13.downsample.bed.gz",
"../Downsample_bed_files/C23.2.downsample.bed.gz",
"../Downsample_bed_files/C23.3.downsample.bed.gz",
"../Downsample_bed_files/C23.4.downsample.bed.gz",
"../Downsample_bed_files/C23.5.downsample.bed.gz",
"../Downsample_bed_files/C23.6.downsample.bed.gz",
"../Downsample_bed_files/C23.7.downsample.bed.gz",
"../Downsample_bed_files/C23.8.downsample.bed.gz",
"../Downsample_bed_files/C23.9.downsample.bed.gz",
"../Downsample_bed_files/C24.1.downsample.bed.gz",
"../Downsample_bed_files/C24.10.downsample.bed.gz",
"../Downsample_bed_files/C24.11.downsample.bed.gz",
"../Downsample_bed_files/C24.12.downsample.bed.gz",
"../Downsample_bed_files/C24.13.downsample.bed.gz",
"../Downsample_bed_files/C24.14.downsample.bed.gz",
"../Downsample_bed_files/C24.15.downsample.bed.gz",
"../Downsample_bed_files/C24.16.downsample.bed.gz",
"../Downsample_bed_files/C24.17.downsample.bed.gz",
"../Downsample_bed_files/C24.18.downsample.bed.gz",
"../Downsample_bed_files/C24.19.downsample.bed.gz",
"../Downsample_bed_files/C24.2.downsample.bed.gz",
"../Downsample_bed_files/C24.20.downsample.bed.gz",
"../Downsample_bed_files/C24.21.downsample.bed.gz",
"../Downsample_bed_files/C24.22.downsample.bed.gz",
"../Downsample_bed_files/C24.23.downsample.bed.gz",
"../Downsample_bed_files/C24.3.downsample.bed.gz",
"../Downsample_bed_files/C24.4.downsample.bed.gz",
"../Downsample_bed_files/C24.5.downsample.bed.gz",
"../Downsample_bed_files/C24.6.downsample.bed.gz",
"../Downsample_bed_files/C24.7.downsample.bed.gz",
"../Downsample_bed_files/C24.8.downsample.bed.gz",
"../Downsample_bed_files/C24.9.downsample.bed.gz",
"../Downsample_bed_files/C3.1.downsample.bed.gz",
"../Downsample_bed_files/C3.2.downsample.bed.gz",
"../Downsample_bed_files/C4.1.downsample.bed.gz",
"../Downsample_bed_files/C5.1.downsample.bed.gz",
"../Downsample_bed_files/C5.10.downsample.bed.gz",
"../Downsample_bed_files/C5.11.downsample.bed.gz",
"../Downsample_bed_files/C5.12.downsample.bed.gz",
"../Downsample_bed_files/C5.13.downsample.bed.gz",
"../Downsample_bed_files/C5.14.downsample.bed.gz",
"../Downsample_bed_files/C5.15.downsample.bed.gz",
"../Downsample_bed_files/C5.16.downsample.bed.gz",
"../Downsample_bed_files/C5.17.downsample.bed.gz",
"../Downsample_bed_files/C5.18.downsample.bed.gz",
"../Downsample_bed_files/C5.19.downsample.bed.gz",
"../Downsample_bed_files/C5.2.downsample.bed.gz",
"../Downsample_bed_files/C5.20.downsample.bed.gz",
"../Downsample_bed_files/C5.21.downsample.bed.gz",
"../Downsample_bed_files/C5.22.downsample.bed.gz",
"../Downsample_bed_files/C5.23.downsample.bed.gz",
"../Downsample_bed_files/C5.24.downsample.bed.gz",
"../Downsample_bed_files/C5.25.downsample.bed.gz",
"../Downsample_bed_files/C5.26.downsample.bed.gz",
"../Downsample_bed_files/C5.27.downsample.bed.gz",
"../Downsample_bed_files/C5.3.downsample.bed.gz",
"../Downsample_bed_files/C5.4.downsample.bed.gz",
"../Downsample_bed_files/C5.5.downsample.bed.gz",
"../Downsample_bed_files/C5.6.downsample.bed.gz",
"../Downsample_bed_files/C5.7.downsample.bed.gz",
"../Downsample_bed_files/C5.8.downsample.bed.gz",
"../Downsample_bed_files/C5.9.downsample.bed.gz",
"../Downsample_bed_files/C6.1.downsample.bed.gz",
"../Downsample_bed_files/C6.2.downsample.bed.gz",
"../Downsample_bed_files/C7.1.downsample.bed.gz",
"../Downsample_bed_files/C7.10.downsample.bed.gz",
"../Downsample_bed_files/C7.2.downsample.bed.gz",
"../Downsample_bed_files/C7.3.downsample.bed.gz",
"../Downsample_bed_files/C7.4.downsample.bed.gz",
"../Downsample_bed_files/C7.5.downsample.bed.gz",
"../Downsample_bed_files/C7.6.downsample.bed.gz",
"../Downsample_bed_files/C7.7.downsample.bed.gz",
"../Downsample_bed_files/C7.8.downsample.bed.gz",
"../Downsample_bed_files/C7.9.downsample.bed.gz",
"../Downsample_bed_files/C8.1.downsample.bed.gz",
"../Downsample_bed_files/C8.10.downsample.bed.gz",
"../Downsample_bed_files/C8.11.downsample.bed.gz",
"../Downsample_bed_files/C8.12.downsample.bed.gz",
"../Downsample_bed_files/C8.13.downsample.bed.gz",
"../Downsample_bed_files/C8.14.downsample.bed.gz",
"../Downsample_bed_files/C8.15.downsample.bed.gz",
"../Downsample_bed_files/C8.16.downsample.bed.gz",
"../Downsample_bed_files/C8.17.downsample.bed.gz",
"../Downsample_bed_files/C8.18.downsample.bed.gz",
"../Downsample_bed_files/C8.19.downsample.bed.gz",
"../Downsample_bed_files/C8.2.downsample.bed.gz",
"../Downsample_bed_files/C8.20.downsample.bed.gz",
"../Downsample_bed_files/C8.21.downsample.bed.gz",
"../Downsample_bed_files/C8.22.downsample.bed.gz",
"../Downsample_bed_files/C8.23.downsample.bed.gz",
"../Downsample_bed_files/C8.24.downsample.bed.gz",
"../Downsample_bed_files/C8.25.downsample.bed.gz",
"../Downsample_bed_files/C8.26.downsample.bed.gz",
"../Downsample_bed_files/C8.27.downsample.bed.gz",
"../Downsample_bed_files/C8.28.downsample.bed.gz",
"../Downsample_bed_files/C8.3.downsample.bed.gz",
"../Downsample_bed_files/C8.4.downsample.bed.gz",
"../Downsample_bed_files/C8.5.downsample.bed.gz",
"../Downsample_bed_files/C8.6.downsample.bed.gz",
"../Downsample_bed_files/C8.7.downsample.bed.gz",
"../Downsample_bed_files/C8.8.downsample.bed.gz",
"../Downsample_bed_files/C8.9.downsample.bed.gz",
"../Downsample_bed_files/C9.1.downsample.bed.gz",
"../Downsample_bed_files/C9.2.downsample.bed.gz",
"../Downsample_bed_files/C9.3.downsample.bed.gz",
"../Downsample_bed_files/C9.4.downsample.bed.gz")
#
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
  addTileMat = F,
  offsetPlus = 0,
  offsetMinus = 0,
  nChunk = 2,
  addGeneScoreMat = T,
  force = T
)
#
#
#
#
proj_TS_Downsample <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "TreeShrew_ArchR_output",
  geneAnnotation = TreeShrew.geneAnnotation,
  genomeAnnotation = TreeShrew.genomeAnnotation,
  copyArrows = T
)
#
#
load("TS_final_annotated_peaks.RData")
#
TS_final_annotated_peaks
table(TS_final_annotated_peaks@seqnames)
#
proj_TS_Downsample <- addPeakSet(proj_TS_Downsample,peakSet = TS_final_annotated_peaks,force = TRUE)
proj_TS_Downsample <- addPeakMatrix(proj_TS_Downsample,  ceiling = 4,force = TRUE)
#
#
#
table(proj_TS_Downsample$sub_ID)
getAvailableMatrices(proj_TS_Downsample)
#
TS_counts <- getGroupSE(ArchRProj = proj_TS_Downsample, useMatrix = "PeakMatrix",
                        groupBy = "Sample",divideN = TRUE,scaleTo = 1e6,verbose = TRUE)
#
#
########################################################################
########################################################################
########################################################################
# 开始创建树鼩坐标的SummarizedExperiment
# 不做GWAS分析
TS_counts2 <- TS_counts@assays@data$PeakMatrix
head(TS_counts2)
dim(TS_counts2)
###
#
#options(scipen = 9999999999999999)  # 设置一个很大的数来禁用科学计数法
#options(scipen = 0)
TS_peak_tmp <- as.data.frame(TS_final_annotated_peaks@ranges)
# 查看结果
head(TS_peak_tmp)
dim(TS_peak_tmp)
#
#
#
head(TS_peak_tmp)
TS_peak_tmp2 <- cbind(as.character(TS_final_annotated_peaks@seqnames),TS_peak_tmp$start,TS_peak_tmp$end)
head(TS_peak_tmp2)
TS_peak_tmp2 <- as.data.frame(TS_peak_tmp2)
head(TS_peak_tmp2)
colnames(TS_peak_tmp2) <- c("chr","start","end")
head(TS_peak_tmp2)
#
grep("e", TS_peak_tmp2$start)
grep("e", TS_peak_tmp2$end)
#
#
peaks <- makeGRangesFromDataFrame(TS_peak_tmp2, seqnames = "chr", start.field = "start", end.field = "end")
#
#创建SummarizedExperiment对象
head(TS_counts2)
class(TS_counts2)
TS_counts2 <- as.matrix(TS_counts2)
dim(TS_counts2)
head(TS_counts2)
#
#
ID_rank <- read.csv("GWAS_celltype_rank.csv",sep = "\t")
rownames(ID_rank) <- ID_rank$Subtype_ID_Num
ID_rank <- ID_rank[colnames(TS_counts2),]
head(ID_rank)
colnames(TS_counts2) <- ID_rank$Subtype_ID_Num_type
#
SE <- SummarizedExperiment(assays = list(counts = TS_counts2),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(TS_counts2)))
#
SE
#
######################################################################################################################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
######################################################################################################################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# 做一下 VISTA 的分析，嘿嘿
# 
#
SE
#添加GC偏差
library(chromVAR)
library(gchromVAR)
#library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
#library(BSgenome.Hsapiens.UCSC.hg38)
#
library(chromVAR)
library(gchromVAR)
#library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
#library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(123)
register(MulticoreParam(2))
#
#
#
#
library(BSgenome.tupBel1.UCSC.TreeShrew)
SE <- addGCBias(SE, genome = BSgenome.tupBel1.UCSC.TreeShrew)
#
SEX_bg <- getBackgroundPeaks(object = SE, niterations = 200)
########################################################################
########################################################################
########################################################################
files <- list.files('/home/sunkeyong/Multispecies_analysis/M_TreeShrew/Analysis_T0808/Peak_analysis_T251220/VISTA/gChromVAR/')
files <- paste0("/home/sunkeyong/Multispecies_analysis/M_TreeShrew/Analysis_T0808/Peak_analysis_T251220/VISTA/gChromVAR/",files)
head(read.table(files[1]))
#
files
ukbb <- importBedScore(rowRanges(SE), files, colidx = 5)
str(ukbb)
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb, 
                                       background_peaks = SEX_bg)
zdf <- reshape2::melt(t(assays(ukbb_wDEV)[["z"]]))
head(zdf)
colnames(zdf) <- c("ct", "tr", "Zscore")
zdf$gchromVAR_pvalue <- pnorm(zdf$Zscore, lower.tail = FALSE)
head(zdf)
#
head(zdf)
VISTA_base_TS_Genome <- zdf
#
#
#
VISTA_base_TS_Genome
#
table(VISTA_base_TS_Genome$tr)
#
View(VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$tr=="human_bv"),])
View(VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$tr=="human_fb"),])
View(VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$tr=="human_eye"),])
View(VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$tr=="human_lb"),])
View(VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$tr=="mouse_bv"),])
#
View(VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$ct=="C3.2:Ventricular cardiomyocytes"),])
View(VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$ct=="C3.2:Ventricular cardiomyocytes"),])
#
#
ggplot(VISTA_base_TS_Genome, aes(x = ct, y = tr, fill = Zscore)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Z-score") +
  labs(x = "Cell Type (ct)", y = "Trait (tr)", 
       title = "Z-score Heatmap") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )
#
#
library(dplyr)
#
##
#
#
#
#
#
#
#
###################################################### # 
###################################################### # 
###################################################### # 
# 先做人的 VISTA 的 组织
# 
target_listX <- c("human_ba","human_bv","human_drg","human_ear","human_eye",
                  "human_fm","human_hb","human_ht","human_lb","human_fb",
                  "human_mb","human_nose","human_nt", 
                  "human_pan","human_som","human_tail")
#
target_listX <- c("human_ht", "human_som",
                  "human_lb","human_ba","human_bv","human_fm","human_nose",
                  "human_fb","human_mb","human_hb","human_nt")
#
Keep_cell_type <- c()
for (i in 1:length(target_listX)) {
  tmp <- VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$tr==target_listX[i]),]
  tmp <- tmp[order(tmp[, 3], decreasing = TRUE), ]
  tmp <- as.character(tmp$ct)[1:5]
  Keep_cell_type <- c(Keep_cell_type,tmp)
}
#
Keep_cell_type <- unique(Keep_cell_type)
#
VISTA_X <- assays(ukbb_wDEV)[["z"]]
head(VISTA_X)
VISTA_X_tmp <- VISTA_X[,Keep_cell_type]
head(VISTA_X_tmp)
#
target_listX2 <- c("human_ht", "human_som",
                   "human_lb","human_ba","human_bv","human_fm","human_nose",
                   "human_fb","human_mb","human_hb","human_nt",
                   "human_negtiveR50")
#
VISTA_X_tmp <- VISTA_X_tmp[target_listX2,]
VISTA_X_tmp
#View(VISTA_X_tmp)
#
#
library(pheatmap)
library(ComplexHeatmap)
library(colorspace)
library(RColorBrewer)
library(circlize)
#
min(VISTA_X_tmp)
max(VISTA_X_tmp)
col_fun <- colorRamp2(
  c(-5, 0, 26), 
  c("#01665e", "white", "red")
)
#
col_fun <- colorRamp2(
  c(-5, 1, 4, 12, 25),                    # 5个数值断点
  c("#01665e", "#5ab4ac", "white", "#ff9999", "red")  # 5个颜色
)
#pheatmap(VISTA_X_tmp, scale = "row",color=col_fun)
pheatmap(VISTA_X_tmp,color=col_fun,cluster_rows = F,cluster_cols = F)
#pheatmap(VISTA_X_tmp)
#
#
#
Human_Celltypes_Keep <- c("C3.2:Ventricular cardiomyocytes",
                          "C3.1:Atrial cardiomyocytes",     
                          "C14.17:Lung Pericytes (HEYL and TBX5+)",  
                          "C15.2:Myotubes",
                          "C14.3:Chondrocytes (ZIC+)",    
                          "C14.15:Limb mesenchyme progenitors",               
                          "C12.5:Nephron progenitors",      
                          "C14.5:Early chondrocytes",       
                          "C14.4:Dermatome",                
                          "C14.22:Pre-osteoblasts (SP7+)",  
                          "C14.20:Mesodermal progenitors",  
                          "C6.1:Definitive erythroblasts",  
                          "C24.13:Megakaryocyte-erythroid progenitor",        
                          "C8.2:Branchial arch epithelium", 
                          "C8.20:Olfactory related cells",  
                          "C8.23:Parathyroid chief cells",         
                          "C5.10:GABAergic cortical interneurons (LHX6+)",    
                          "C18.10:Hypothalamus (DLX2 and NKX2-1+)",           
                          "C11.2:Intermediate neuronal progenitors",          
                          "C18.20:Telencephalon",     
                          "C5.13:GABAergic neurons",        
                          "C5.17:Glutamatergic neurons (PITX2+)",             
                          "C18.12:Midbrain-hindbrain boundary",  
                          "C5.20:Neural progenitor cells (ROR1+)",            
                          "C5.23:Spinal cord dorsal and ventral progenitors", 
                          "C5.22:Purkinje cells",           
                          "C5.6:Cholinergic neurons",                         
                          "C2.3:Parasympathetic neurons_GFRAL+")  
#
VISTA_X_tmp_keep <- VISTA_X_tmp[,Human_Celltypes_Keep]
#
min(VISTA_X_tmp_keep)
max(VISTA_X_tmp_keep)
col_fun <- colorRamp2(
  c(-5, 1, 4, 12, 25),                    # 5个数值断点
  c("#01665e", "#5ab4ac", "white", "#ff9999", "red")  # 5个颜色
)
#
pheatmap(VISTA_X_tmp_keep,color=col_fun,cluster_rows = F,cluster_cols = F)
#
#

##
#
#
#
#
#
#
#
###################################################### # 
###################################################### # 
###################################################### # 
# 做小鼠的 VISTA 的 组织
#
#
target_listX <- c("mouse_ba","mouse_bv",        
                  "mouse_fb","mouse_fm","mouse_hb","mouse_ht","mouse_lb",        
                  "mouse_mb","mouse_negtive","mouse_nose","mouse_nt")
#
target_listX <- c("mouse_ht","mouse_ba","mouse_lb","mouse_fm",            
                  "mouse_fb","mouse_hb",
                  "mouse_mb","mouse_nt")
#
Keep_cell_type <- c()
for (i in 1:length(target_listX)) {
  tmp <- VISTA_base_TS_Genome[which(VISTA_base_TS_Genome$tr==target_listX[i]),]
  tmp <- tmp[order(tmp[, 3], decreasing = TRUE), ]
  tmp <- as.character(tmp$ct)[1:5]
  Keep_cell_type <- c(Keep_cell_type,tmp)
}
#
Keep_cell_type <- unique(Keep_cell_type)
#
VISTA_X <- assays(ukbb_wDEV)[["z"]]
head(VISTA_X)
VISTA_X_tmp <- VISTA_X[,Keep_cell_type]
head(VISTA_X_tmp)
#
target_listX2 <- c("mouse_ht","mouse_ba","mouse_lb","mouse_fm",            
                   "mouse_fb","mouse_hb",
                   "mouse_mb","mouse_nt",
                   "mouse_negtiveR20")
#
VISTA_X_tmp <- VISTA_X_tmp[target_listX2,]
VISTA_X_tmp
#View(VISTA_X_tmp)
#
#
library(pheatmap)
library(ComplexHeatmap)
library(colorspace)
library(RColorBrewer)
library(circlize)
#
min(VISTA_X_tmp)
max(VISTA_X_tmp)
col_fun <- colorRamp2(
  c(-3, 0, 12), 
  c("#01665e", "white", "red")
)
#
col_fun <- colorRamp2(
  c(-4, 0.5, 4, 8, 12),                    # 5个数值断点
  c("#01665e", "#5ab4ac", "white", "#ff9999", "red")  # 5个颜色
)
#pheatmap(VISTA_X_tmp, scale = "row",color=col_fun)
pheatmap(VISTA_X_tmp,color=col_fun,cluster_rows = F,cluster_cols = F)
#pheatmap(VISTA_X_tmp)
#
#
#
Mouse_Celltypes_Keep <- c("C3.1:Atrial cardiomyocytes",
                          "C3.2:Ventricular cardiomyocytes",
                          "C7.3:Endocardial cells",
                          "C14.9:Hepatic mesenchyme",
                          "C14.24:Smooth muscle cells",
                          "C14.4:Dermatome",
                          "C14.22:Pre-osteoblasts (SP7+)",
                          "C14.6:Facial mesenchyme",
                          "C14.10:Intermediate mesoderm (KBTBD13+)",
                          "C14.15:Limb mesenchyme progenitors",
                          "C14.20:Mesodermal progenitors",
                          "C14.23:Sclerotome",
                          "C14.2:Chondrocytes (POU3F3+)",
                          "C14.7:Fibroblasts",
                          "C5.10:GABAergic cortical interneurons (LHX6+)",
                          "C5.12:GABAergic cortical interneurons (ZNF503+)",
                          "C18.20:Telencephalon",
                          "C11.2:Intermediate neuronal progenitors",
                          "C5.8:GABAergic cortical interneurons (DLX6+)",
                          "C5.19:Motor neurons",
                          "C5.14:GABAergic neurons (GATA3+)",
                          "C5.17:Glutamatergic neurons (PITX2+)",
                          "C5.4:Cajal-Retzius cells (TBR1+)",
                          "C5.13:GABAergic neurons",
                          "C16.2:Neural crest (PNS glia)",
                          "C5.15:GABAergic neurons (PAX7+)",
                          "C5.1:Cajal-Retzius cells (BARHL2+)",
                          "C16.4:Satellite glial cells",
                          "C14.5:Early chondrocytes",
                          "C5.27:Thalamic neuronal precursors")  
#
VISTA_X_tmp_keep <- VISTA_X_tmp[,Mouse_Celltypes_Keep]
#
min(VISTA_X_tmp_keep)
max(VISTA_X_tmp_keep)
col_fun <- colorRamp2(
  c(-4, 0.5, 4, 8, 12),                    # 5个数值断点
  c("#01665e", "#5ab4ac", "white", "#ff9999", "red")  # 5个颜色
)
#
pheatmap(VISTA_X_tmp_keep,color=col_fun,cluster_rows = F,cluster_cols = F)
#
save.image("TS.VISTA.RData")
#