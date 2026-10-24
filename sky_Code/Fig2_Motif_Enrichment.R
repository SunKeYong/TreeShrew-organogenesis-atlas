#  topic cCREs TF motif enrichment using ArchR
#  chromVARmotifs CIS-BP data (v2)
# 
setwd("path/Peak_Motif")
# 
# 
## 
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
#
##
########### ############################################################### #
# 
addArchRThreads(threads = 8)
#
library(BSgenome.tupBel1.UCSC.TreeShrew)
TreeShrew.genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.tupBel1.UCSC.TreeShrew)
#
proj_TS_Downsample <- addGroupCoverages(ArchRProj = proj_TS_Downsample, 
                               minCells = 1600,
                               maxCells = 10500,
                               maxFragments = 35 * 10^10,
                               minReplicates = 2,
                               maxReplicates = 5,
                               sampleRatio = 0.8,
                               groupBy = "Sample",
                               useLabels=F,
                               force = TRUE)
#
#
library(TFBSTools)
library(JASPAR2024)
#
library(JASPAR2024) 
library(TFBSTools)
jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
PFM_list <- TFBSTools::getMatrixSet(sq24, list(species = "Homo sapiens", collection = "CORE"))
print(PFM_list)
#
PWM_list <- toPWM(PFM_list)# 检查结果print(PWM_list)
#
library(chromVARmotifs)
library(TFBSTools)
# 加载人类 CIS-BP 数据 (v2)
data("human_pwms_v2") 
#
#
#
proj_TS_Downsample <- addMotifAnnotations(ArchRProj = proj_TS_Downsample,motifPWMs=human_pwms_v2, name = "Motif",force = T)
proj_TS_Downsample <- addBgdPeaks(proj_TS_Downsample)
#
#
proj_TS_Downsample <- addDeviationsMatrix(
  ArchRProj = proj_TS_Downsample, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
plotVarDev <- getVarDeviations(proj_TS_Downsample, name = "MotifMatrix", plot = TRUE)
plotVarDev
#
#
#
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}
#
#
#
#
#
peaks_tmp <- paste(TS_final_annotated_peaks@seqnames,TS_final_annotated_peaks@ranges,sep = ":")
#
Target_topic <- read.csv("../Topics_otsu/topic_62.bed",sep = "\t",header = F)
Target_topic_ArchR <- cbind(Target_topic,paste(paste(Target_topic$V1,as.numeric(Target_topic$V2)+1,sep = ":"),as.numeric(Target_topic$V3-1),sep="-"))
colnames(Target_topic_ArchR) <- c("chr","start","end","new")
Target_topic_ArchR <- as.data.frame(Target_topic_ArchR)
co_peak <- intersect(Target_topic_ArchR$new,peaks_tmp)
Target_topic_ArchR <- subset(Target_topic_ArchR, new %in% co_peak)
Target_topic_ArchR2 <- makeGRangesFromDataFrame(Target_topic_ArchR, seqnames = "chr", start.field = "start", end.field = "end")
Target_topic_motif <- customEnrichment(
  ranges = Target_topic_ArchR2, 
  matches = getMatches(proj_TS_Downsample), 
  bgdPeaks = getBgdPeaks(proj_TS_Downsample, force=TRUE)
)
#
######### #
head(Target_topic_motif)
#
pwm <- getPeakAnnotation(proj_TS_Downsample, "Motif")$motifs[[Target_topic_motif$feature[1]]]
#pwm
#
ppm <- PWMatrixToProbMatrix(pwm)
#ppm
library(ggseqlogo)
plot <- ggseqlogo(ppm, method = "bits")
plot
#

#
#
save.image("TS.motif.enrichment.RData")
#
#
#
#
