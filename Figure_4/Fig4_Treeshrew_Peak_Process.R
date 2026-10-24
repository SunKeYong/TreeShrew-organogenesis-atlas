library(ChIPseeker)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ArchR)
#
setwd("/home/sunkeyong/Human_em/Cross_species/3_species")
#
TS_peak_raw <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/TS_peaks_by_celltypes.csv",sep = ",",header = T)
head(TS_peak_raw)
tail(TS_peak_raw$Peaks)
dim(TS_peak_raw)
#
TS_peak_raw$X <- paste("TreeShrew_Peak_rank",1:dim(TS_peak_raw)[1],sep = "")
head(TS_peak_raw)
#
TS_peak_list <- cbind(TS_peak_raw$Peaks,TS_peak_raw$X)
TS_peak_list <- as.data.frame(TS_peak_list)
TS_peak_list$X <- unlist(strsplit(TS_peak_list$V1,":"))[seq(2,2*dim(TS_peak_list)[1],2)]
TS_peak_list$chr <- unlist(strsplit(TS_peak_list$V1,":"))[seq(1,2*dim(TS_peak_list)[1],2)]
TS_peak_list$start <- unlist(strsplit(TS_peak_list$X,"-"))[seq(1,2*dim(TS_peak_list)[1],2)]
TS_peak_list$end <- unlist(strsplit(TS_peak_list$X,"-"))[seq(2,2*dim(TS_peak_list)[1],2)]
head(TS_peak_list)
#
TS_peak_list2 <- cbind(TS_peak_list$chr,TS_peak_list$start,TS_peak_list$end,TS_peak_list$V2)
TS_peak_list2 <- as.data.frame(TS_peak_list2)
head(TS_peak_list2)
table(TS_peak_list2$V1)
#
TS_peak_list2 <- subset(TS_peak_list2, !grepl("chrM", V1))
TS_peak_list2 <- subset(TS_peak_list2, !grepl("chrY", V1))
TS_peak_list2 <- subset(TS_peak_list2, !grepl("conti", V1))
#
head(TS_peak_list2)
table(TS_peak_list2$V1)
#
write.table(TS_peak_list2,file = "./TreeShrew_/TS_peak_list2.raw.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#
#
TreeShrew_keep <- TS_peak_list2
colnames(TreeShrew_keep) <- c("chr","snapATAC2_start","end","TreeShrew_Peak_rank")
TreeShrew_keep$ArchR_start <- sprintf("%.0f", as.numeric(TreeShrew_keep$snapATAC2_start) + 1)
head(TreeShrew_keep)
#
# 注释peaks
library(ArchR)
library(BSgenome.tupBel1.UCSC.TreeShrew)
#
#
load("/home/sunkeyong/softeware/ref_tree_shrew/ArchR_annotation_Storage/TreeShrew.genomeAnnotation.RData")
load("/home/sunkeyong/softeware/ref_tree_shrew/ArchR_annotation_Storage/TreeShrew.geneAnnotation.RData")
#
#
#
TreeShrew_keep_GRange <- makeGRangesFromDataFrame(TreeShrew_keep, seqnames = "chr", start.field = "ArchR_start", end.field = "end")
TreeShrew_keep_GRange$rank <- TreeShrew_keep$TreeShrew_Peak_rank
TreeShrew_keep_GRange
#
# 注释
annotated_peaks <- annotatePeakList(
  peaks = TreeShrew_keep_GRange,
  genomeAnnotation = TreeShrew.genomeAnnotation,
  geneAnnotation = TreeShrew.geneAnnotation)
#
annotated_peaks
#
table(annotated_peaks$peakType)
max(annotated_peaks$N)
hist(annotated_peaks$GC)
which(annotated_peaks$GC > 0.8)
which(annotated_peaks$GC < 0.2)
#
annotated_peaks
#
annotated_peaks.meta <- as.data.frame(annotated_peaks)
head(annotated_peaks.meta)
#
table(annotated_peaks.meta$rank==TreeShrew_keep$TreeShrew_Peak_rank)
#
TreeShrew_keep$ArchR_distToGeneStart <- annotated_peaks.meta$distToGeneStart
TreeShrew_keep$ArchR_nearestGene <- annotated_peaks.meta$nearestGene
TreeShrew_keep$ArchR_peakType <- annotated_peaks.meta$peakType
TreeShrew_keep$ArchR_distToTSS <- annotated_peaks.meta$distToTSS
TreeShrew_keep$ArchR_GC <- annotated_peaks.meta$GC
#
head(TreeShrew_keep)
#
table(TreeShrew_keep$ArchR_peakType)
#
#
library(ChIPseeker)
library(GenomicFeatures)
library(rtracklayer)
#
txdb <- makeTxDbFromGFF("/home/sunkeyong/softeware/ref_tree_shrew/TreeShrew_scATAC_ref/genes/genes.gtf", format = "gtf")
txdb
#
peakAnno <- annotatePeak(
  TreeShrew_keep_GRange,
  TxDb = txdb,
  tssRegion = c(-2000, 100),
  level = "gene"
)
#
head(as.data.frame(peakAnno))
df <- as.data.frame(peakAnno)
head(df)
table(df$annotation)
plotAnnoPie(peakAnno)
#
df$anno4 <- case_when(
  grepl("^Promoter", df$annotation) ~ "Promoter",
  grepl("^Exon", df$annotation) | grepl("UTR", df$annotation) ~ "Exonic",
  grepl("^Intron", df$annotation) ~ "Intronic",
  TRUE ~ "Distal"
)
#
table(df$anno4)
prop.table(table(df$anno4))
#
#
TreeShrew_keep$ChIPseeker_RawAnnotation <- df$annotation
TreeShrew_keep$ChIPseeker_anno <- df$anno4
#
head(TreeShrew_keep)
#
table(TreeShrew_keep$ArchR_peakType)/dim(TreeShrew_keep)[1]
table(TreeShrew_keep$ChIPseeker_anno)/dim(TreeShrew_keep)[1]
#
hist(TreeShrew_keep$ArchR_GC)
#
TreeShrew_keep_Final <- TreeShrew_keep[which(TreeShrew_keep$ArchR_GC <= 0.75),]
TreeShrew_keep_Final <- TreeShrew_keep_Final[which(TreeShrew_keep_Final$ArchR_GC >= 0.25),]
dim(TreeShrew_keep)
dim(TreeShrew_keep_Final)
#
#
hist(TreeShrew_keep_Final$ArchR_GC)
#
head(TreeShrew_keep_Final)
#
#
TreeShrew_keep_Final_bed_tmp <- cbind(TreeShrew_keep_Final$chr,TreeShrew_keep_Final$snapATAC2_start,TreeShrew_keep_Final$end,TreeShrew_keep_Final$TreeShrew_Peak_rank)
head(TreeShrew_keep_Final_bed_tmp)
#
write.table(TreeShrew_keep_Final_bed_tmp,file = "./TreeShrew_//TreeShrew_keep_Final_bed_tmp.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#
TOPA_TE_anno <- read.csv("./TreeShrew_/TreeShrew_keep_Final_bed_tmp.dfam.TEonly_overlap_final.txt",sep = "\t",header = F)
head(TOPA_TE_anno)
table(TOPA_TE_anno$V8)
table(TOPA_TE_anno$V11)
table(TOPA_TE_anno$V12)
dim(TOPA_TE_anno) # 1081613
TOPA_TE_anno$ratio <- TOPA_TE_anno$V16/TOPA_TE_anno$V17
#TOPA_TE_anno2 <- TOPA_TE_anno[which(TOPA_TE_anno$ratio > 0.5),]
TOPA_TE_anno2 <- TOPA_TE_anno[which(TOPA_TE_anno$V16 > 250),]
head(TOPA_TE_anno2)
dim(TOPA_TE_anno2) # 131846
#
TOPA_TE_anno2 <- cbind(TOPA_TE_anno2$V4,TOPA_TE_anno2$V11,TOPA_TE_anno2$V12,TOPA_TE_anno2$V8)
head(TOPA_TE_anno2)
colnames(TOPA_TE_anno2) <- c("TreeShrew_Peak_rank","repClass","repFamily","repName")
TOPA_TE_anno2 <- as.data.frame(TOPA_TE_anno2)
head(TOPA_TE_anno2)
#
length(unique(TOPA_TE_anno2$rank))
#
head(TreeShrew_keep_Final)
#
TreeShrew_keep_Final_TE <- merge(TreeShrew_keep_Final,TOPA_TE_anno2,by="TreeShrew_Peak_rank", all = TRUE)
dim(TreeShrew_keep_Final_TE)
head(TreeShrew_keep_Final_TE)
#
dim(TreeShrew_keep_Final_TE)
dim(TOPA_TE_anno2)
dim(TreeShrew_keep_Final)
length(unique(TOPA_TE_anno2$TreeShrew_Peak_rank))
length(unique(TreeShrew_keep_Final_TE$TreeShrew_Peak_rank))
#
#
#
#
TreeShrew_keep_Final_TE2 <- TreeShrew_keep_Final_TE
TreeShrew_keep_Final_TE2$repClass[is.na(TreeShrew_keep_Final_TE2$repClass)] <- "other"
TreeShrew_keep_Final_TE2$tmp_rep <- paste(TreeShrew_keep_Final_TE2$TreeShrew_Peak_rank,TreeShrew_keep_Final_TE2$repClass,sep = "__")
TreeShrew_keep_Final_TE3 <- TreeShrew_keep_Final_TE2[!duplicated(TreeShrew_keep_Final_TE2$tmp_rep), ]
head(TreeShrew_keep_Final_TE2)
dim(TreeShrew_keep_Final_TE3)
dim(TreeShrew_keep_Final_TE2)
dim(TreeShrew_keep_Final_TE)
#
table(TreeShrew_keep_Final_TE2$repClass,TreeShrew_keep_Final_TE2$ChIPseeker_anno)
table(TreeShrew_keep_Final_TE2$ChIPseeker_anno)
colSums(table(TreeShrew_keep_Final_TE2$repClass,TreeShrew_keep_Final_TE2$ChIPseeker_anno)[-4,])/table(TreeShrew_keep_Final_TE2$ChIPseeker_anno)
#
#
colSums(table(TreeShrew_keep_Final_TE3$repClass,TreeShrew_keep_Final_TE3$ChIPseeker_anno)[-4,])/table(TreeShrew_keep_Final_TE$ChIPseeker_anno)
#   Distal    Exonic   Intronic  Promoter 
# 0.13052682 0.03396170 0.10445029 0.05732334 
#
colSums(table(TreeShrew_keep_Final_TE2$repClass,TreeShrew_keep_Final_TE2$ChIPseeker_anno)[-4,])/table(TreeShrew_keep_Final_TE$ChIPseeker_anno)
#
# Distal    Exonic  Intronic  Promoter 
# 0.13055496 0.03396960 0.10447095 0.05732334
#
colSums(table(TreeShrew_keep_Final_TE3$repClass,TreeShrew_keep_Final_TE3$ArchR_peakType)[-4,])/table(TreeShrew_keep_Final_TE$ChIPseeker_anno)
#   Distal    Exonic   Intronic  Promoter 
# 0.14218401 0.01646292 0.09876537 0.04891080  
#
colSums(table(TreeShrew_keep_Final_TE2$repClass,TreeShrew_keep_Final_TE2$ArchR_peakType)[-4,])/table(TreeShrew_keep_Final_TE$ChIPseeker_anno)
#
# Distal    Exonic  Intronic  Promoter 
# 0.14221416 0.01647083 0.09878431 0.04891080 
#
#
#
#
#
head(TS_peak_raw)
dim(TS_peak_raw)
#
#
colnames(TS_peak_raw)
#
celltype_list <- colnames(TS_peak_raw)[3:199]
#
i <-130
peaks_tmpX <-  TS_peak_raw[c(1,i+2)]
colnames(peaks_tmpX) <- c("TreeShrew_Peak_rank","celltype")
head(peaks_tmpX)
peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
dim(peaks_tmpX)
peaks_tmpX2 <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% peaks_tmpX$TreeShrew_Peak_rank )
dim(peaks_tmpX2)
head(peaks_tmpX2)
table(peaks_tmpX2$repClass)/dim(peaks_tmpX)[1]
t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ChIPseeker_anno))/as.numeric(table(peaks_tmpX2$ChIPseeker_anno)))
#
peaks_tmpX3 <- peaks_tmpX2[which(peaks_tmpX2$repClass!="other"),]
length(peaks_tmpX3$TreeShrew_Peak_rank)/dim(peaks_tmpX)[1]
length(unique(peaks_tmpX3$TreeShrew_Peak_rank))/dim(peaks_tmpX)[1]
#
TE_ratio_ArchR <- c()
TE_ratio_ChipSeeker <- c()
TE_CT_ratio <- c()
#
for (i in 1:197) {
  peaks_tmpX <-  TS_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("TreeShrew_Peak_rank","celltype")
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  peaks_tmpX2 <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% peaks_tmpX$TreeShrew_Peak_rank )
  TE_ratio_ArchR <- cbind(TE_ratio_ArchR,t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ArchR_peakType))/as.numeric(table(peaks_tmpX2$ArchR_peakType))))
  TE_ratio_ChipSeeker <- cbind(TE_ratio_ChipSeeker,t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ChIPseeker_anno))/as.numeric(table(peaks_tmpX2$ChIPseeker_anno))))
  peaks_tmpX3 <- peaks_tmpX2[which(peaks_tmpX2$repClass!="other"),]
  TE_CT_ratio <- cbind(TE_CT_ratio,length(unique(peaks_tmpX3$TreeShrew_Peak_rank))/dim(peaks_tmpX)[1])
}
#
#
#
ID_rank <- read.csv("/home/sunkeyong/Multispecies_analysis/M_TreeShrew/Analysis_T0808/pseudo_bulk_analysis/bed_dedu/ArchR_showing/GWAS_celltype_rank.csv",sep = "\t")
head(ID_rank)
#
#
TE_ratio <- TE_ratio_ArchR
#
colnames(TE_ratio) <- paste(rep(celltype_list,each=4),colnames(TE_ratio) ,sep="_")
TE_ratio_Distal <- TE_ratio[,paste(ID_rank$Subtype_ID_Num,"Distal" ,sep="_")]
TE_ratio_Exonic <- TE_ratio[,paste(ID_rank$Subtype_ID_Num,"Exonic" ,sep="_")]
TE_ratio_Intronic <- TE_ratio[,paste(ID_rank$Subtype_ID_Num,"Intronic" ,sep="_")]
TE_ratio_Promoter <- TE_ratio[,paste(ID_rank$Subtype_ID_Num,"Promoter" ,sep="_")]
#
barplot(1-TE_ratio_Distal["other",])
barplot(1-TE_ratio_Exonic["other",])
barplot(1-TE_ratio_Intronic["other",])
barplot(1-TE_ratio_Promoter["other",])
#
#View(t(TE_ratio_Distal))
#
boxplot(1-TE_ratio_Distal["other",],1-TE_ratio_Intronic["other",],
        1-TE_ratio_Exonic["other",],
        1-TE_ratio_Promoter["other",])
#
wilcox.test(1-TE_ratio_Exonic["other",],
            1-TE_ratio_Promoter["other",])
#
# save.image("TreeShrew_TOPA_TE.RData")
#
median(1-TE_ratio_Distal["other",]) # 0.06385472 
median(1-TE_ratio_Intronic["other",]) # 0.04563519
median(1-TE_ratio_Exonic["other",]) # 0.008488523
median(1-TE_ratio_Promoter["other",]) # 0.01185436
#
head(TE_CT_ratio)
TE_CT_ratio <- as.data.frame(t(TE_CT_ratio))
rownames(TE_CT_ratio) <- celltype_list
TE_CT_ratio$celltype <- celltype_list
head(TE_CT_ratio)
#
#
TE_CT_ratio2 <-TE_CT_ratio
TE_CT_ratio2 <- TE_CT_ratio2[ID_rank$Subtype_ID_Num,]
head(TE_CT_ratio2)
#
barplot(TE_CT_ratio2$V1)
boxplot(TE_CT_ratio2$V1)
median(TE_CT_ratio2$V1) # 0.04355954
#
#
#
save.image("TreeShrew_TOPA_TE.RData")
##
#
#
#
head(TreeShrew_keep_Final_TE)
#
#
#
#
#
############
head(TreeShrew_keep_Final_TE)
#
TreeShrew_keep_Final_TE3 <- TreeShrew_keep_Final_TE
TreeShrew_keep_Final_TE3$repClass[is.na(TreeShrew_keep_Final_TE3$repClass)] <- "None"
TreeShrew_keep_Final_TE3$repFamily[is.na(TreeShrew_keep_Final_TE3$repFamily)] <- "None"
TreeShrew_keep_Final_TE3$repName[is.na(TreeShrew_keep_Final_TE3$repName)] <- "None"
head(TreeShrew_keep_Final_TE3)
#
length(TreeShrew_keep_Final_TE3$TreeShrew_Peak_rank)
length(unique(TreeShrew_keep_Final_TE3$TreeShrew_Peak_rank))
#
#
TOPA_lift_mouse <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_Mouse.classified.bed",sep = "\t",header = F)
TOPA_lift_Hg38 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_human.classified.bed",sep = "\t",header = F)
#
head(TOPA_lift_mouse)
head(TOPA_lift_Hg38)
#
dim(TOPA_lift_mouse)
dim(TOPA_lift_Hg38)
#
table(TOPA_lift_mouse$V4==TOPA_lift_Hg38$V4)
#
TOPA_lift <- TOPA_lift_Hg38
colnames(TOPA_lift) <- c("chr","snapATAC2_start","end","TreeShrew_Peak_rank","TreeShrew_to_Hg38_State")
TOPA_lift$TreeShrew_to_Mm_State <- TOPA_lift_mouse$V5
head(TOPA_lift)
#
table(TOPA_lift$TreeShrew_to_Hg38_State)
table(TOPA_lift$TreeShrew_to_Mm_State)
#
table(TOPA_lift$TreeShrew_to_Hg38_State)/sum(table(TOPA_lift$TreeShrew_to_Hg38_State))
# CA_conserved   CA_divergent Human_specific 
#   0.4323962    0.3304168    0.2371871 
labels <- paste("TS_to_human", c("CA_conserved","CA_divergent","TS_specific"),sep="-")
labels
pie(table(TOPA_lift$TreeShrew_to_Hg38_State)/sum(table(TOPA_lift$TreeShrew_to_Hg38_State)), 
    labels = labels, col = c("#EAA0A0","#46ADD1","#2F6EAF"))
#
table(TOPA_lift$TreeShrew_to_Mm_State)/sum(table(TOPA_lift$TreeShrew_to_Mm_State))
# CA_conserved   CA_divergent  Human_specific 
#    0.1737207    0.3073893    0.5188900 
labels <- paste("TS_to_Mm", c("CA_conserved","CA_divergent","TS_specific"),sep="-")
labels
pie(table(TOPA_lift$TreeShrew_to_Mm_State)/sum(table(TOPA_lift$TreeShrew_to_Mm_State)), 
    labels = labels, col = c("#EAA0A0","#46ADD1","#2F6EAF"))
#
#
### 
# Full peak set
table(TreeShrew_keep_Final_TE3$repClass)/sum(table(TreeShrew_keep_Final_TE3$repClass))
#      DNA       LINE        LTR       None       SINE 
# 0.006840732 0.054372190 0.034995769 0.893540708 0.010250601 
#
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Hg38_State=="TS_specific" & TOPA_lift$TreeShrew_to_Mm_State=="TS_specific"),]
targetTE <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE$repClass)/sum(table(targetTE$repClass))
#.    DNA       LINE        LTR       None       SINE 
# 0.00663802 0.12260705 0.06381316 0.77867091 0.02827085 
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Hg38_State=="TS_specific"),]
targetTE <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE$repClass)/sum(table(targetTE$repClass))
#.    DNA       LINE        LTR       None       SINE 
# 0.006777595 0.117033517 0.061866409 0.787756755 0.026565723
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Hg38_State=="CA_conserved"),]
targetTE <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE$repClass)/sum(table(targetTE$repClass))
#.     DNA       LINE        LTR       None       SINE 
# 0.005916022 0.023462404 0.020597765 0.947250675 0.002773135 
#
#
########
head(TS_peak_raw)
head(TOPA_lift)
#
i <-197
peaks_tmpX <-  TS_peak_raw[c(1,i+2)]
colnames(peaks_tmpX) <- c("TreeShrew_Peak_rank","celltype")
head(peaks_tmpX)
peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
dim(peaks_tmpX)
peaks_tmpX2 <- subset(TOPA_lift, TreeShrew_Peak_rank %in% peaks_tmpX$TreeShrew_Peak_rank )
table(peaks_tmpX2$TreeShrew_to_Hg38_State)/dim(peaks_tmpX2)[1]
table(peaks_tmpX2$TreeShrew_to_Mm_State)/dim(peaks_tmpX2)[1]
#
CA_to_Hg38_ratio <- c()
CA_to_Mm_ratio <- c()
#
for (i in 1:197) {
  peaks_tmpX <-  TS_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("TreeShrew_Peak_rank","celltype")
  head(peaks_tmpX)
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  dim(peaks_tmpX)
  peaks_tmpX2 <- subset(TOPA_lift, TreeShrew_Peak_rank %in% peaks_tmpX$TreeShrew_Peak_rank )
  CA_to_Hg38_ratio <- cbind(CA_to_Hg38_ratio,table(peaks_tmpX2$TreeShrew_to_Hg38_State)/dim(peaks_tmpX2)[1])
  CA_to_Mm_ratio <- cbind(CA_to_Mm_ratio,table(peaks_tmpX2$TreeShrew_to_Mm_State)/dim(peaks_tmpX2)[1])
  #
}
#
#
head(CA_to_Hg38_ratio)
#
colnames(CA_to_Hg38_ratio) <- celltype_list
CA_to_Hg38_ratio <- CA_to_Hg38_ratio[,ID_rank$Subtype_ID_Num]
pheatmap::pheatmap(CA_to_Hg38_ratio,cluster_rows = F,cluster_cols = F)
###
head(CA_to_Mm_ratio)
colnames(CA_to_Mm_ratio) <- celltype_list
CA_to_Mm_ratio <- CA_to_Mm_ratio[,ID_rank$Subtype_ID_Num]
pheatmap::pheatmap(CA_to_Mm_ratio,cluster_rows = F,cluster_cols = F)
#
#
save.image("TreeShrew_TOPA_TE.RData")
#
##################### #
#
head(TreeShrew_keep_Final_TE3)
table(TreeShrew_keep_Final_TE3$repClass)
table(TreeShrew_keep_Final_TE3$repFamily)
#

#
#














########
############## ####################### ####################### #
###################### ####################### ####################### #
##
head(TreeShrew_keep_Final)
head(TOPA_lift)
dim(TOPA_lift)
dim(TreeShrew_keep_Final)
#
#
table(TOPA_lift$TreeShrew_to_Hg38_State)
table(TOPA_lift$TreeShrew_to_Mm_State)
#
table(TreeShrew_keep_Final$ChIPseeker_anno)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "CA_conserved")$TreeShrew_Peak_rank)$ChIPseeker_anno)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "CA_divergent")$TreeShrew_Peak_rank)$ChIPseeker_anno)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "TS_specific")$TreeShrew_Peak_rank)$ChIPseeker_anno)
#
#
table(TreeShrew_keep_Final$ArchR_peakType)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "CA_conserved")$TreeShrew_Peak_rank)$ArchR_peakType)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "CA_divergent")$TreeShrew_Peak_rank)$ArchR_peakType)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "TS_specific")$TreeShrew_Peak_rank)$ArchR_peakType)
#
table(TreeShrew_keep_Final$ChIPseeker_anno)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Mm_State %in% "CA_conserved")$TreeShrew_Peak_rank)$ChIPseeker_anno)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Mm_State %in% "CA_divergent")$TreeShrew_Peak_rank)$ChIPseeker_anno)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Mm_State %in% "TS_specific")$TreeShrew_Peak_rank)$ChIPseeker_anno)
#
table(TreeShrew_keep_Final$ArchR_peakType)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Mm_State %in% "CA_conserved")$TreeShrew_Peak_rank)$ArchR_peakType)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Mm_State %in% "CA_divergent")$TreeShrew_Peak_rank)$ArchR_peakType)
table(subset(TreeShrew_keep_Final,  TreeShrew_Peak_rank %in% subset(TOPA_lift,TreeShrew_to_Mm_State %in% "TS_specific")$TreeShrew_Peak_rank)$ArchR_peakType)
#
#
###
####
i <- 1
peaks_tmpX <-  TS_peak_raw[c(1,i+2)]
colnames(peaks_tmpX) <- c("TreeShrew_Peak_rank","celltype")
peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
head(peaks_tmpX)
peaks_tmpX_A <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "CA_conserved")$TreeShrew_Peak_rank))
peaks_tmpX_B <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "CA_divergent")$TreeShrew_Peak_rank))
peaks_tmpX_C <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "TS_specific")$TreeShrew_Peak_rank))
#
peaks_tmpX_A2 <- table(peaks_tmpX_A$repClass)/dim(peaks_tmpX_A)[1]
peaks_tmpX_B2 <- table(peaks_tmpX_B$repClass)/dim(peaks_tmpX_B)[1]
peaks_tmpX_C2 <- table(peaks_tmpX_C$repClass)/dim(peaks_tmpX_C)[1]
TE_ratio_CA_conserved <- cbind(TE_ratio_CA_conserved,1-peaks_tmpX_A2[4])
TE_ratio_CA_divergent <- cbind(TE_ratio_CA_divergent,1-peaks_tmpX_B2[4])
TE_ratio_TS_specific <- cbind(TE_ratio_TS_specific,1-peaks_tmpX_C2[4])
#
#
TE_ratio_CA_conserved_hg38 <- c()
TE_ratio_CA_divergent_hg38 <- c()
TE_ratio_TS_specific_hg38 <- c()
#
TE_ratio_CA_conserved_mm10 <- c()
TE_ratio_CA_divergent_mm10 <- c()
TE_ratio_TS_specific_mm10 <- c()
#
#
for (i in 1:197) {
  peaks_tmpX <-  TS_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("TreeShrew_Peak_rank","celltype")
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  head(peaks_tmpX)
  peaks_tmpX_A <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "CA_conserved")$TreeShrew_Peak_rank))
  peaks_tmpX_B <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "CA_divergent")$TreeShrew_Peak_rank))
  peaks_tmpX_C <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Hg38_State %in% "TS_specific")$TreeShrew_Peak_rank))
  #
  peaks_tmpX_A2 <- table(peaks_tmpX_A$repClass)/dim(peaks_tmpX_A)[1]
  peaks_tmpX_B2 <- table(peaks_tmpX_B$repClass)/dim(peaks_tmpX_B)[1]
  peaks_tmpX_C2 <- table(peaks_tmpX_C$repClass)/dim(peaks_tmpX_C)[1]
  TE_ratio_CA_conserved_hg38 <- cbind(TE_ratio_CA_conserved_hg38,1-peaks_tmpX_A2[4])
  TE_ratio_CA_divergent_hg38 <- cbind(TE_ratio_CA_divergent_hg38,1-peaks_tmpX_B2[4])
  TE_ratio_TS_specific_hg38 <- cbind(TE_ratio_TS_specific_hg38,1-peaks_tmpX_C2[4])
  #
  peaks_tmpX_A <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Mm_State %in% "CA_conserved")$TreeShrew_Peak_rank))
  peaks_tmpX_B <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Mm_State %in% "CA_divergent")$TreeShrew_Peak_rank))
  peaks_tmpX_C <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% intersect(peaks_tmpX$TreeShrew_Peak_rank,subset(TOPA_lift,TreeShrew_to_Mm_State %in% "TS_specific")$TreeShrew_Peak_rank))
  #
  peaks_tmpX_A2 <- table(peaks_tmpX_A$repClass)/dim(peaks_tmpX_A)[1]
  peaks_tmpX_B2 <- table(peaks_tmpX_B$repClass)/dim(peaks_tmpX_B)[1]
  peaks_tmpX_C2 <- table(peaks_tmpX_C$repClass)/dim(peaks_tmpX_C)[1]
  TE_ratio_CA_conserved_mm10 <- cbind(TE_ratio_CA_conserved_mm10,1-peaks_tmpX_A2[4])
  TE_ratio_CA_divergent_mm10 <- cbind(TE_ratio_CA_divergent_mm10,1-peaks_tmpX_B2[4])
  TE_ratio_TS_specific_mm10 <- cbind(TE_ratio_TS_specific_mm10,1-peaks_tmpX_C2[4])
}
#
#
boxplot(as.numeric(TE_ratio_CA_conserved_hg38),as.numeric(TE_ratio_CA_divergent_hg38),as.numeric(TE_ratio_TS_specific_hg38),
        as.numeric(TE_ratio_CA_conserved_mm10),as.numeric(TE_ratio_CA_divergent_mm10),as.numeric(TE_ratio_TS_specific_mm10),
        TE_CT_ratio2$V1)
#
TS_TE_ratio_merge <- cbind(as.numeric(TE_ratio_CA_conserved_hg38),as.numeric(TE_ratio_CA_divergent_hg38),as.numeric(TE_ratio_TS_specific_hg38),
                           as.numeric(TE_ratio_CA_conserved_mm10),as.numeric(TE_ratio_CA_divergent_mm10),as.numeric(TE_ratio_TS_specific_mm10),
                           TE_CT_ratio2$V1)
TS_TE_ratio_merge <- as.data.frame(TS_TE_ratio_merge)
colnames(TS_TE_ratio_merge) <- c("CA_conserved_hg38","CA_divergent_hg38","TS_specific_hg38","CA_conserved_mm10","CA_divergent_mm10","TS_specific_mm10","TE_all_ratio")
head(TS_TE_ratio_merge)
#
save(TS_TE_ratio_merge,file = "./data_save_For_plot2/TS_TE_ratio_merge.RData")
#
plot(CA_to_Hg38_ratio[3,],TS_TE_ratio_merge$TE_all_ratio)
plot(CA_to_Mm_ratio[3,],TS_TE_ratio_merge$TE_all_ratio)
#
TS_CA_to_Hg38_ratio <- CA_to_Hg38_ratio
TS_CA_to_Mm_ratio <- CA_to_Mm_ratio
#
save(TS_CA_to_Hg38_ratio,file = "./data_save_For_plot2/TS_CA_to_Hg38_ratio.RData")
save(TS_CA_to_Mm_ratio,file = "./data_save_For_plot2/TS_CA_to_Mm_ratio.RData")
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























#
############ ###### #################################### #
############ ###### #################################### #
############ ###### #################################### #
# 计算 树鼩的Fig2B
TS_to_M_0.10 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_Mouse.Min0.10.bed",header = F,sep = "\t")
TS_to_M_0.25 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_Mouse.Min0.25.bed",header = F,sep = "\t")
TS_to_M_0.50 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_Mouse.Min0.50.bed",header = F,sep = "\t")
TS_to_M_0.75 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_Mouse.Min0.75.bed",header = F,sep = "\t")
#
TS_to_M_0.10$TreeShrew_Peak_rank <- TS_to_M_0.10$V4
TS_to_M_0.25$TreeShrew_Peak_rank <- TS_to_M_0.25$V4
TS_to_M_0.50$TreeShrew_Peak_rank <- TS_to_M_0.50$V4
TS_to_M_0.75$TreeShrew_Peak_rank <- TS_to_M_0.75$V4
#
TS_to_M_0.10 <- merge(TS_to_M_0.10,TOPA_lift,by="TreeShrew_Peak_rank", all =T)
TS_to_M_0.25 <- merge(TS_to_M_0.25,TOPA_lift,by="TreeShrew_Peak_rank", all =T)
TS_to_M_0.50 <- merge(TS_to_M_0.50,TOPA_lift,by="TreeShrew_Peak_rank", all =T)
TS_to_M_0.75 <- merge(TS_to_M_0.75,TOPA_lift,by="TreeShrew_Peak_rank", all =T)
#
rownames(TS_to_M_0.10) <- TS_to_M_0.10$TreeShrew_Peak_rank
rownames(TS_to_M_0.25) <- TS_to_M_0.25$TreeShrew_Peak_rank
rownames(TS_to_M_0.50) <- TS_to_M_0.50$TreeShrew_Peak_rank
rownames(TS_to_M_0.75) <- TS_to_M_0.75$TreeShrew_Peak_rank
#
TS_to_M_0.10 <- TS_to_M_0.10[TOPA_lift$TreeShrew_Peak_rank,]
TS_to_M_0.25 <- TS_to_M_0.25[TOPA_lift$TreeShrew_Peak_rank,]
TS_to_M_0.50 <- TS_to_M_0.50[TOPA_lift$TreeShrew_Peak_rank,]
TS_to_M_0.75 <- TS_to_M_0.75[TOPA_lift$TreeShrew_Peak_rank,]
#
head(TS_to_M_0.10)
head(TS_to_M_0.25)
head(TS_to_M_0.50)
head(TS_to_M_0.75)
#
TS_to_H_0.10 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_human.Min0.10.bed",header = F,sep = "\t")
TS_to_H_0.25 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_human.Min0.25.bed",header = F,sep = "\t")
TS_to_H_0.50 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_human.Min0.50.bed",header = F,sep = "\t")
TS_to_H_0.75 <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/TS_to_human.Min0.75.bed",header = F,sep = "\t")
#
TS_to_H_0.10$TreeShrew_Peak_rank <- TS_to_H_0.10$V4
TS_to_H_0.25$TreeShrew_Peak_rank <- TS_to_H_0.25$V4
TS_to_H_0.50$TreeShrew_Peak_rank <- TS_to_H_0.50$V4
TS_to_H_0.75$TreeShrew_Peak_rank <- TS_to_H_0.75$V4
#
TS_to_H_0.10 <- merge(TS_to_H_0.10,TOPA_lift,by="TreeShrew_Peak_rank", all =T)
TS_to_H_0.25 <- merge(TS_to_H_0.25,TOPA_lift,by="TreeShrew_Peak_rank", all =T)
TS_to_H_0.50 <- merge(TS_to_H_0.50,TOPA_lift,by="TreeShrew_Peak_rank", all =T)
TS_to_H_0.75 <- merge(TS_to_H_0.75,TOPA_lift,by="TreeShrew_Peak_rank", all =T)
#
rownames(TS_to_H_0.10) <- TS_to_H_0.10$TreeShrew_Peak_rank
rownames(TS_to_H_0.25) <- TS_to_H_0.25$TreeShrew_Peak_rank
rownames(TS_to_H_0.50) <- TS_to_H_0.50$TreeShrew_Peak_rank
rownames(TS_to_H_0.75) <- TS_to_H_0.75$TreeShrew_Peak_rank
#
TS_to_H_0.10 <- TS_to_H_0.10[TOPA_lift$TreeShrew_Peak_rank,]
TS_to_H_0.25 <- TS_to_H_0.25[TOPA_lift$TreeShrew_Peak_rank,]
TS_to_H_0.50 <- TS_to_H_0.50[TOPA_lift$TreeShrew_Peak_rank,]
TS_to_H_0.75 <- TS_to_H_0.75[TOPA_lift$TreeShrew_Peak_rank,]
#
head(TS_to_H_0.10)
head(TS_to_H_0.25)
head(TS_to_H_0.50)
head(TS_to_H_0.75)
#
head(TOPA_lift)
head(TreeShrew_keep_Final)
#
TreeShrew_keep_FinalX <- TreeShrew_keep_Final
#
TreeShrew_keep_FinalX$TS_to_Human_minMatch0.10 <- grepl("ree",TS_to_H_0.10$V4)
TreeShrew_keep_FinalX$TS_to_Human_minMatch0.25 <- grepl("ree",TS_to_H_0.25$V4)
TreeShrew_keep_FinalX$TS_to_Human_minMatch0.50 <- grepl("ree",TS_to_H_0.50$V4)
TreeShrew_keep_FinalX$TS_to_Human_minMatch0.75 <- grepl("ree",TS_to_H_0.75$V4)
#
TreeShrew_keep_FinalX$TS_to_Mouse_minMatch0.10 <- grepl("ree",TS_to_M_0.10$V4)
TreeShrew_keep_FinalX$TS_to_Mouse_minMatch0.25 <- grepl("ree",TS_to_M_0.25$V4)
TreeShrew_keep_FinalX$TS_to_Mouse_minMatch0.50 <- grepl("ree",TS_to_M_0.50$V4)
TreeShrew_keep_FinalX$TS_to_Mouse_minMatch0.75 <- grepl("ree",TS_to_M_0.75$V4)
#
head(TreeShrew_keep_FinalX)
#
which(colnames(TreeShrew_keep_FinalX)=="ArchR_peakType")
#
x.h <- TreeShrew_keep_FinalX[,c(8,13:16)]
x.m <- TreeShrew_keep_FinalX[,c(8,17:20)]
#
head(x.h)
head(x.m)
#
library(dplyr)
library(tidyr)
library(ggplot2)
##
res <- x.h %>% dplyr::group_by(ArchR_peakType) %>% dplyr::summarise(dplyr::across(dplyr::starts_with("TS_to_Human"),~ mean(.x, na.rm = TRUE),.names = "{.col}_ratio"))
#res <- x.m %>% dplyr::group_by(ArchR_peakType) %>% dplyr::summarise(dplyr::across(dplyr::starts_with("TS_to_Mouse"),~ mean(.x, na.rm = TRUE),.names = "{.col}_ratio"))
res <- as.data.frame(res)
res
res2 <- rbind(res, c("Shuffle", c(218304,204705,161041,113465)/sum(c(218304,204705,161041,113465)))) # 人
#res2 <- rbind(res, c("Shuffle", c(128593,114061,84822,55127)/sum(c(128593,114061,84822,55127)))) # 鼠
colnames(res2) <- colnames(res)
# 转成 data.frame
res2 <- as.data.frame(res2, stringsAsFactors = FALSE)
# 把后面几列转成数值
res2[, -1] <- lapply(res2[, -1], as.numeric)
library(dplyr)
library(tidyr)
library(ggplot2)
# 人
res_long2 <- res2 %>%pivot_longer(cols = -ArchR_peakType,names_to = "minMatch",values_to = "ratio") %>% mutate(minMatch = gsub("_ratio", "", minMatch),minMatch = gsub("TS_to_Human_", "", minMatch),ArchR_peakType = factor(ArchR_peakType,levels = c("Exonic", "Promoter", "Intronic", "Distal", "Shuffle")))
# 鼠
#res_long2 <- res2 %>%pivot_longer(cols = -ArchR_peakType,names_to = "minMatch",values_to = "ratio") %>% mutate(minMatch = gsub("_ratio", "", minMatch),minMatch = gsub("TS_to_Mouse_", "", minMatch),ArchR_peakType = factor(ArchR_peakType,levels = c("Exonic", "Promoter", "Intronic", "Distal", "Shuffle")))
ggplot(res_long2, aes(x = ArchR_peakType, y = ratio, fill = minMatch)) +
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = c(
    "minMatch0.10" = "#C9DBED",
    "minMatch0.25" = "#A9C3DF",
    "minMatch0.50" = "#89ADD2",
    "minMatch0.75" = "#6896C4"
  ))  +theme_bw() +
  labs(x = "ChIPseeker annotation",y = "Proportion TRUE",fill = "minMatch") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
####
#### 生成 树鼩基因组 随机的 位置
#
TS_site.random.bed <- read.csv("/home/sunkeyong/Multispecies_analysis/M_TreeShrew/Analysis_T0808/Peak_analysis_T251220/conservation/TS_site.random.bed",header = F,sep = "\t")
TS_chr <- read.csv("/home/sunkeyong/softeware/ref_tree_shrew/Download_from_NCBI/chr.txt",header = F,sep = "\t")
head(TS_site.random.bed)
head(TS_chr)
TS_site.random.bed2 <- merge(TS_site.random.bed,TS_chr,by="V1")
head(TS_site.random.bed2)
TS_site.random.bed2 <- cbind(TS_site.random.bed2$V2.y,TS_site.random.bed2$V2.x,TS_site.random.bed2$V3)
write.table(TS_site.random.bed2,file = "/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/random_shuffle/TS_site.random.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#
#
##
#
##
table(TreeShrew_keep_Final$ArchR_peakType)
table(TreeShrew_keep_Final$ArchR_peakType)/sum(table(TreeShrew_keep_Final$ArchR_peakType))
#
table(TreeShrew_keep_Final$ChIPseeker_anno)
table(TreeShrew_keep_Final$ChIPseeker_anno)/sum(table(TreeShrew_keep_Final$ChIPseeker_anno))
# 输入频数
anno_counts <- c(Promoter = 33878, Exonic = 126465, Intronic = 580647, Distal = 497446) # ChIPseeker_anno
anno_counts <- c(Promoter = 28725, Exonic = 60688, Intronic = 593920, Distal = 555103) # ArchR_peakType 选用
# 转换为数据框，并按指定顺序设置因子
df <- data.frame(Category = factor(names(anno_counts), levels = c("Promoter", "Exonic", "Intronic", "Distal")),Count = as.numeric(anno_counts))
# 按因子顺序排序（保证扇形顺序）
df <- df[order(df$Category), ]
# 计算百分比和累积和
df$Fraction <- df$Count / sum(df$Count)
df$ymax <- cumsum(df$Fraction)
df$ymin <- c(0, head(df$ymax, n = -1))
df$label_pos <- (df$ymin + df$ymax) / 2
df$Percentage <- paste0(round(df$Fraction * 100, 1), "%")
#
library(ggpubr)
ggdonutchart(df, "Count", label = "Percentage", 
             lab.pos = "in", fill = "Category", 
             color = "white", palette = "Set2")
##

##

##

##

##
head(TOPA_lift)
#
table(TreeShrew_keep_Final$TreeShrew_Peak_rank==TOPA_lift$TreeShrew_Peak_rank)
rownames(TOPA_lift) <- TOPA_lift$TreeShrew_Peak_rank
TOPA_lift2 <- TOPA_lift[TreeShrew_keep_FinalX$TreeShrew_Peak_rank,]
table(TreeShrew_keep_FinalX$TreeShrew_Peak_rank==TOPA_lift2$TreeShrew_Peak_rank)
#
#
TreeShrew_keep_FinalX$TreeShrew_to_Hg38_State <- TOPA_lift2$TreeShrew_to_Hg38_State
TreeShrew_keep_FinalX$TreeShrew_to_Mm10_State <- TOPA_lift2$TreeShrew_to_Mm_State
#
head(TreeShrew_keep_FinalX)
#
write.table(TreeShrew_keep_FinalX,file = "./Human_and_Treeshrew_peak_save/TreeShrew_keep_Final.bed", sep = "\t", col.names = T, row.names = F, quote = F)
#
#
save.image("TreeShrew_TOPA_TE.T0428_correction.RData")
#