library(ChIPseeker)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ArchR)
#
setwd("/home/sunkeyong/Human_em/Cross_species/3_species")
#
HE_peak_raw <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/HE_peaks_by_celltypes.csv",sep = ",",header = T)
head(HE_peak_raw)
tail(HE_peak_raw$Peaks)
dim(HE_peak_raw)
#
HE_peak_raw$X <- paste("Human_Peak_rank",1:dim(HE_peak_raw)[1],sep = "")
head(HE_peak_raw)
#
HE_peak_list <- cbind(HE_peak_raw$Peaks,HE_peak_raw$X)
HE_peak_list <- as.data.frame(HE_peak_list)
HE_peak_list$X <- unlist(strsplit(HE_peak_list$V1,":"))[seq(2,2*dim(HE_peak_list)[1],2)]
HE_peak_list$chr <- unlist(strsplit(HE_peak_list$V1,":"))[seq(1,2*dim(HE_peak_list)[1],2)]
HE_peak_list$start <- unlist(strsplit(HE_peak_list$X,"-"))[seq(1,2*dim(HE_peak_list)[1],2)]
HE_peak_list$end <- unlist(strsplit(HE_peak_list$X,"-"))[seq(2,2*dim(HE_peak_list)[1],2)]
head(HE_peak_list)
#
HE_peak_list2 <- cbind(HE_peak_list$chr,HE_peak_list$start,HE_peak_list$end,HE_peak_list$V2)
HE_peak_list2 <- as.data.frame(HE_peak_list2)
#
HE_peak_list2 <- subset(HE_peak_list2, !grepl("chrM", V1))
HE_peak_list2 <- subset(HE_peak_list2, !grepl("chrY", V1))
#
head(HE_peak_list2)
table(HE_peak_list2$V1)
#
write.table(HE_peak_list2,file = "./human_//HE_peak_list2.raw.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#
#
human_keep <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/human_/HE_peak_list2.rmBlacklist.bed",sep = "\t",header = F)
colnames(human_keep) <- c("chr","snapATAC2_start","end","Human_Peak_rank")
human_keep$ArchR_start <- sprintf("%.0f", as.numeric(human_keep$snapATAC2_start) + 1)
head(human_keep)
#
# 注释peaks
# 我们在这里只需要输入 GRanges的Peak
# 我们在这里只需要输入 geneAnnotation
# 我们在这里只需要输入 genomeAnnotation
#
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ArchR)
#
Hg38_genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Hsapiens.UCSC.hg38)
Hg38_geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, OrgDb = org.Hs.eg.db)
#
#
#
human_keep_GRange <- makeGRangesFromDataFrame(human_keep, seqnames = "chr", start.field = "ArchR_start", end.field = "end")
human_keep_GRange$rank <- human_keep$Human_Peak_rank
human_keep_GRange
#
# 注释
annotated_peaks <- annotatePeakList(
  peaks = human_keep_GRange,
  genomeAnnotation = Hg38_genomeAnnotation,
  geneAnnotation = Hg38_geneAnnotation)
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
table(annotated_peaks.meta$rank==human_keep$Human_Peak_rank)
#
human_keep$UCSC_distToGeneStart <- annotated_peaks.meta$distToGeneStart
human_keep$UCSC_nearestGene <- annotated_peaks.meta$nearestGene
human_keep$UCSC_peakType <- annotated_peaks.meta$peakType
human_keep$UCSC_distToTSS <- annotated_peaks.meta$distToTSS
human_keep$UCSC_GC <- annotated_peaks.meta$GC
#
head(human_keep)
#
addArchRGenome("hg38")
ArchR_hg38_GeneAnnotation <- getGeneAnnotation()
ArchR_hg38_GenomeAnnotation <- getGenomeAnnotation()
#
# 注释
annotated_peaks2 <- annotatePeakList(
  peaks = human_keep_GRange,
  genomeAnnotation = ArchR_hg38_GenomeAnnotation,
  geneAnnotation = ArchR_hg38_GeneAnnotation)
#
annotated_peaks2
#
table(annotated_peaks2$peakType)
max(annotated_peaks2$N)
hist(annotated_peaks2$GC)
which(annotated_peaks2$GC > 0.8)
which(annotated_peaks2$GC < 0.2)
#
annotated_peaks2
#
annotated_peaks2.meta <- as.data.frame(annotated_peaks2)
head(annotated_peaks2.meta)
#
table(annotated_peaks2.meta$rank==human_keep$Human_Peak_rank)
#
human_keep$ArchR_distToGeneStart <- annotated_peaks2.meta$distToGeneStart
human_keep$ArchR_nearestGene <- annotated_peaks2.meta$nearestGene
human_keep$ArchR_peakType <- annotated_peaks2.meta$peakType
human_keep$ArchR_distToTSS <- annotated_peaks2.meta$distToTSS
human_keep$ArchR_GC <- annotated_peaks2.meta$GC
#
head(human_keep)
#
table(human_keep$UCSC_peakType)
table(human_keep$ArchR_peakType)
#
#
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(org.Hs.eg.db)
#
annotated_peaks3 <- annotatePeak(
  human_keep_GRange,
  TxDb = txdb,
  tssRegion = c(-2000, 100),
  annoDb = "org.Hs.eg.db"
)
df <- as.data.frame(annotated_peaks3)
head(df)
table(df$annotation)
plotAnnoPie(annotated_peaks3)
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
human_keep$ChIPseeker_RawAnnotation <- df$annotation
human_keep$ChIPseeker_anno <- df$anno4
#
head(human_keep)
#
table(human_keep$UCSC_peakType)
table(human_keep$ArchR_peakType)
table(human_keep$ChIPseeker_anno)
#
human_keep_Final <- human_keep[which(human_keep$UCSC_GC <= 0.8),]
human_keep_Final <- human_keep_Final[which(human_keep_Final$UCSC_GC >= 0.2),]
dim(human_keep)
dim(human_keep_Final)
#
human_keep_Final <- human_keep_Final[which(human_keep_Final$ArchR_GC <= 0.8),]
human_keep_Final <- human_keep_Final[which(human_keep_Final$ArchR_GC >= 0.2),]
dim(human_keep_Final)
#
hist(human_keep$UCSC_GC)
hist(human_keep_Final$UCSC_GC)
#
head(human_keep_Final)
#
#
human_keep_Final_bed_tmp <- cbind(human_keep_Final$chr,human_keep_Final$snapATAC2_start,human_keep_Final$end,human_keep_Final$Human_Peak_rank)
head(human_keep_Final_bed_tmp)
#
write.table(human_keep_Final_bed_tmp,file = "./human_//human_keep_Final_bed_tmp.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#
#HOPA_TE_anno <- read.csv("./human_/human_keep_Final_bed_tmp.bed.dfam.TEonly_overlap_final.txt",sep = "\t",header = F)
HOPA_TE_anno <- read.csv("./human_/human_keep_Final_bed_tmp.bed.UCSC.TEonly_overlap_final.txt",sep = "\t",header = F)
head(HOPA_TE_anno)
table(HOPA_TE_anno$V8)
table(HOPA_TE_anno$V11)
table(HOPA_TE_anno$V12)
dim(HOPA_TE_anno) #  962353
HOPA_TE_anno$ratio <- HOPA_TE_anno$V13/HOPA_TE_anno$V15
#HOPA_TE_anno2 <- HOPA_TE_anno[which(HOPA_TE_anno$ratio > 0.5),]
HOPA_TE_anno2 <- HOPA_TE_anno[which(HOPA_TE_anno$V13 > 250),]
head(HOPA_TE_anno2)
dim(HOPA_TE_anno2) # 218811
#
# HOPA_TE_anno2 <- cbind(HOPA_TE_anno2$V4,HOPA_TE_anno2$V11,HOPA_TE_anno2$V12,HOPA_TE_anno2$V8)
HOPA_TE_anno2 <- cbind(HOPA_TE_anno2$V4,HOPA_TE_anno2$V11,HOPA_TE_anno2$V12,HOPA_TE_anno2$V8)
head(HOPA_TE_anno2)
colnames(HOPA_TE_anno2) <- c("Human_Peak_rank","repClass","repFamily","repName")
HOPA_TE_anno2 <- as.data.frame(HOPA_TE_anno2)
head(HOPA_TE_anno2)
#
length(unique(HOPA_TE_anno2$rank))
#
head(human_keep_Final)
#
human_keep_Final_TE <- merge(human_keep_Final,HOPA_TE_anno2,by="Human_Peak_rank", all = TRUE)
dim(human_keep_Final_TE)
head(human_keep_Final_TE)
#
dim(human_keep_Final_TE)
dim(HOPA_TE_anno2)
dim(human_keep_Final)
length(unique(HOPA_TE_anno2$Human_Peak_rank))
length(unique(human_keep_Final_TE$Human_Peak_rank))
#
#
#
#
human_keep_Final_TE2 <- human_keep_Final_TE
human_keep_Final_TE2$repClass[is.na(human_keep_Final_TE2$repClass)] <- "other"
human_keep_Final_TE2$tmp_rep <- paste(human_keep_Final_TE2$Human_Peak_rank,human_keep_Final_TE2$repClass,sep = "__")
human_keep_Final_TE3 <- human_keep_Final_TE2[!duplicated(human_keep_Final_TE2$tmp_rep), ]
head(human_keep_Final_TE2)
dim(human_keep_Final_TE3)
dim(human_keep_Final_TE2)
dim(human_keep_Final_TE)
#
table(human_keep_Final_TE2$repClass,human_keep_Final_TE2$ChIPseeker_anno)
table(human_keep_Final_TE2$ChIPseeker_anno)
colSums(table(human_keep_Final_TE2$repClass,human_keep_Final_TE2$ChIPseeker_anno)[-4,])/table(human_keep_Final_TE2$ChIPseeker_anno)
#
#
colSums(table(human_keep_Final_TE3$repClass,human_keep_Final_TE3$ChIPseeker_anno)[-4,])/table(human_keep_Final_TE$ChIPseeker_anno)
#   Distal    Exonic   Intronic  Promoter 
# 0.26871856 0.09270628 0.22134562 0.10984451
#
colSums(table(human_keep_Final_TE2$repClass,human_keep_Final_TE2$ChIPseeker_anno)[-4,])/table(human_keep_Final_TE$ChIPseeker_anno)
#
#  Distal    Exonic  Intronic  Promoter 
# 0.26878129 0.09272919 0.22140720 0.10988865
#
colSums(table(human_keep_Final_TE3$repClass,human_keep_Final_TE3$ArchR_peakType)[-4,])/table(human_keep_Final_TE$ChIPseeker_anno)
#   Distal    Exonic   Intronic  Promoter 
# 0.35675035 0.02739726 0.19243965 0.02881294 
#
colSums(table(human_keep_Final_TE2$repClass,human_keep_Final_TE2$ArchR_peakType)[-4,])/table(human_keep_Final_TE$ChIPseeker_anno)
#
# Distal    Exonic  Intronic  Promoter 
# 0.35681935 0.02739726 0.19250700 0.02882398
#
#
#
#
#
head(HE_peak_raw)
dim(HE_peak_raw)
#
#
colnames(HE_peak_raw)
celltype_list <- colnames(HE_peak_raw)[3:133]
#
i <-130
peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
head(peaks_tmpX)
peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
dim(peaks_tmpX)
peaks_tmpX2 <- subset(human_keep_Final_TE3, Human_Peak_rank %in% peaks_tmpX$Human_Peak_rank )
dim(peaks_tmpX2)
head(peaks_tmpX2)
table(peaks_tmpX2$repClass)/dim(peaks_tmpX)[1]
t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ChIPseeker_anno))/as.numeric(table(peaks_tmpX2$ChIPseeker_anno)))
#
peaks_tmpX3 <- peaks_tmpX2[which(peaks_tmpX2$repClass!="other"),]
length(peaks_tmpX3$Human_Peak_rank)/dim(peaks_tmpX)[1]
length(unique(peaks_tmpX3$Human_Peak_rank))/dim(peaks_tmpX)[1]
#
TE_ratio_ArchR <- c()
TE_ratio_UCSC <- c()
TE_ratio_ChipSeeker <- c()
TE_CT_ratio <- c()
#
for (i in 1:131) {
  peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  peaks_tmpX2 <- subset(human_keep_Final_TE3, Human_Peak_rank %in% peaks_tmpX$Human_Peak_rank )
  TE_ratio_ArchR <- cbind(TE_ratio_ArchR,t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ArchR_peakType))/as.numeric(table(peaks_tmpX2$ArchR_peakType))))
  TE_ratio_UCSC <- cbind(TE_ratio_UCSC,t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$UCSC_peakType))/as.numeric(table(peaks_tmpX2$UCSC_peakType))))
  TE_ratio_ChipSeeker <- cbind(TE_ratio_ChipSeeker,t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ChIPseeker_anno))/as.numeric(table(peaks_tmpX2$ChIPseeker_anno))))
  peaks_tmpX3 <- peaks_tmpX2[which(peaks_tmpX2$repClass!="other"),]
  TE_CT_ratio <- cbind(TE_CT_ratio,length(unique(peaks_tmpX3$Human_Peak_rank))/dim(peaks_tmpX)[1])
}
#
#
#
ID_rank <- read.csv("/home/sunkeyong/Human_em/sub_grep_analysi/Data_analysis/Peak_Validation/liftover/He_evo_rank.csv",sep = "\t",header = F)
#
#
TE_ratio <- TE_ratio_ChipSeeker
#
colnames(TE_ratio) <- paste(rep(celltype_list,each=4),colnames(TE_ratio) ,sep="_")
TE_ratio_Distal <- TE_ratio[,paste(ID_rank$V1,"Distal" ,sep="_")]
TE_ratio_Exonic <- TE_ratio[,paste(ID_rank$V1,"Exonic" ,sep="_")]
TE_ratio_Intronic <- TE_ratio[,paste(ID_rank$V1,"Intronic" ,sep="_")]
TE_ratio_Promoter <- TE_ratio[,paste(ID_rank$V1,"Promoter" ,sep="_")]
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
# save.image("Human_HOPA_TE.RData")
#
median(1-TE_ratio_Distal["other",]) # 
median(1-TE_ratio_Intronic["other",]) # 
median(1-TE_ratio_Exonic["other",]) # 
median(1-TE_ratio_Promoter["other",]) # 
#
head(TE_CT_ratio)
TE_CT_ratio <- as.data.frame(t(TE_CT_ratio))
rownames(TE_CT_ratio) <- celltype_list
TE_CT_ratio$celltype <- celltype_list
head(TE_CT_ratio)
#
#
TE_CT_ratio2 <-TE_CT_ratio
TE_CT_ratio2 <- TE_CT_ratio2[ID_rank$V1,]
head(TE_CT_ratio2)
#
barplot(TE_CT_ratio2$V1)
boxplot(TE_CT_ratio2$V1)
median(TE_CT_ratio2$V1) # 0.0972785
#
#
#
save.image("Human_HOPA_TE.RData")
##
#
#
#
############
head(human_keep_Final_TE)
#
human_keep_Final_TE3 <- human_keep_Final_TE
human_keep_Final_TE3$repClass[is.na(human_keep_Final_TE3$repClass)] <- "None"
human_keep_Final_TE3$repFamily[is.na(human_keep_Final_TE3$repFamily)] <- "None"
human_keep_Final_TE3$repName[is.na(human_keep_Final_TE3$repName)] <- "None"
head(human_keep_Final_TE3)
#
length(human_keep_Final_TE3$Human_Peak_rank)
length(unique(human_keep_Final_TE3$Human_Peak_rank))
#
#
HOPA_lift_mouse <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/Human_to_Mouse/human.classified.bed",sep = "\t",header = F)
HOPA_lift_TS <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/Human_to_TreeShrew/human.classified.bed",sep = "\t",header = F)
#
head(HOPA_lift_mouse)
head(HOPA_lift_TS)
#
dim(HOPA_lift_mouse)
dim(HOPA_lift_TS)
#
table(HOPA_lift_mouse$V4==HOPA_lift_TS$V4)
#
HOPA_lift <- HOPA_lift_TS
colnames(HOPA_lift) <- c("chr","snapATAC2_start","end","Human_Peak_rank","Human_to_TS_State")
HOPA_lift$Human_to_Mm_State <- HOPA_lift_mouse$V5
head(HOPA_lift)
#
table(HOPA_lift$Human_to_TS_State)
table(HOPA_lift$Human_to_Mm_State)
#
table(HOPA_lift$Human_to_TS_State)/sum(table(HOPA_lift$Human_to_TS_State))
# CA_conserved   CA_divergent Human_specific 
# 0.5064041      0.1878526      0.3057433
#
labels <- paste("Human_to_TS", c("CA_conserved","CA_divergent","Human_specific"),sep="-")
labels
pie(table(HOPA_lift$Human_to_TS_State)/sum(table(HOPA_lift$Human_to_TS_State)), 
    labels = labels, col = c("#F72F48","#8E5030","#F6CAAC"))
#
table(HOPA_lift$Human_to_Mm_State)/sum(table(HOPA_lift$Human_to_Mm_State))
# CA_conserved   CA_divergent  Human_specific 
#  0.2160015      0.3417426      0.4422560 
labels <- paste("Human_to_Mm", c("CA_conserved","CA_divergent","Human_specific"),sep="-")
labels
pie(table(HOPA_lift$Human_to_Mm_State)/sum(table(HOPA_lift$Human_to_Mm_State)), 
    labels = labels, col = c("#F72F48","#8E5030","#F6CAAC"))
#
#
#
### 
# Full peak set
table(human_keep_Final_TE3$repClass)/sum(table(human_keep_Final_TE3$repClass))
#      DNA       LINE        LTR       None       SINE 
# 0.01583665 0.08362323 0.08751820 0.78472807 0.02829385 
#
targetX <- HOPA_lift[which(HOPA_lift$Human_to_TS_State=="Human_specific" & HOPA_lift$Human_to_Mm_State=="Human_specific"),]
targetTE <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE$repClass)/sum(table(targetTE$repClass))
#.    DNA       LINE        LTR       None       SINE 
# 0.03935105 0.18185727 0.23008261 0.47237185 0.07633722
targetX <- HOPA_lift[which(HOPA_lift$Human_to_TS_State=="Human_specific"),]
targetTE <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE$repClass)/sum(table(targetTE$repClass))
#.    DNA       LINE        LTR       None       SINE 
# 0.03410886 0.16304107 0.19852702 0.53805514 0.06626791  
targetX <- HOPA_lift[which(HOPA_lift$Human_to_TS_State=="CA_conserved"),]
targetTE <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE$repClass)/sum(table(targetTE$repClass))
#.     DNA       LINE        LTR       None       SINE 
#  0.006739783 0.038396752 0.029828737 0.917459516 0.007575213 
#
#
########
head(HE_peak_raw)
head(HOPA_lift)
#
i <-130
peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
head(peaks_tmpX)
peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
dim(peaks_tmpX)
peaks_tmpX2 <- subset(HOPA_lift, Human_Peak_rank %in% peaks_tmpX$Human_Peak_rank )
table(peaks_tmpX2$Human_to_TS_State)/dim(peaks_tmpX2)[1]
table(peaks_tmpX2$Human_to_Mm_State)/dim(peaks_tmpX2)[1]
#
CA_to_TS_ratio <- c()
CA_to_Mm_ratio <- c()
#
for (i in 1:131) {
  peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
  head(peaks_tmpX)
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  dim(peaks_tmpX)
  peaks_tmpX2 <- subset(HOPA_lift, Human_Peak_rank %in% peaks_tmpX$Human_Peak_rank )
  CA_to_TS_ratio <- cbind(CA_to_TS_ratio,table(peaks_tmpX2$Human_to_TS_State)/dim(peaks_tmpX2)[1])
  CA_to_Mm_ratio <- cbind(CA_to_Mm_ratio,table(peaks_tmpX2$Human_to_Mm_State)/dim(peaks_tmpX2)[1])
  #
}
#
#
head(CA_to_TS_ratio)
colnames(CA_to_TS_ratio) <- celltype_list
CA_to_TS_ratio <- CA_to_TS_ratio[,ID_rank$V1]
pheatmap::pheatmap(CA_to_TS_ratio,cluster_rows = F,cluster_cols = F)
###
head(CA_to_Mm_ratio)
colnames(CA_to_Mm_ratio) <- celltype_list
CA_to_Mm_ratio <- CA_to_Mm_ratio[,ID_rank$V1]
pheatmap::pheatmap(CA_to_Mm_ratio,cluster_rows = F,cluster_cols = F)
#
#
save.image("Human_HOPA_TE.RData")
#
#
#
#
##################### ###################### #
##################### ###################### #
# 分析一下TE repFamily 在不同细胞里面的比例
head(human_keep_Final_TE)
length(human_keep_Final_TE$Human_Peak_rank)
length(unique(human_keep_Final_TE$Human_Peak_rank))
#
human_keep_Final_TE4 <- human_keep_Final_TE
head(human_keep_Final_TE4)
human_keep_Final_TE4$repClass[is.na(human_keep_Final_TE4$repClass)] <- "other"
human_keep_Final_TE4$repFamily[is.na(human_keep_Final_TE4$repFamily)] <- "other"
human_keep_Final_TE4$repName[is.na(human_keep_Final_TE4$repName)] <- "other"
#
table(human_keep_Final_TE4$repClass)
table(human_keep_Final_TE4$repFamily) # 39
table(human_keep_Final_TE4$repName) # 971
#
TE_repFamily_type  <- table(human_keep_Final_TE4$repFamily)
TE_repFamily_type <- as.data.frame(TE_repFamily_type)
# View(TE_repFamily_type)
dim(TE_repFamily_type) # 39 只有39种 就不过滤了
head(TE_repFamily_type)
#
TE_repFamily_Num <- TE_repFamily_type
TE_repFamily_Num$Freq <- 0
head(TE_repFamily_Num)
tail(TE_repFamily_Num)
dim(TE_repFamily_Num)
View(TE_repFamily_Num)
#
TE_repFamily_Num_sub <- TE_repFamily_Num[-27,]
TE_repFamily_Num_all <- TE_repFamily_Num
#
dim(TE_repFamily_Num_sub)
dim(TE_repFamily_Num_all)
#
#
#
i <-130
peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
head(peaks_tmpX)
colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
head(peaks_tmpX)
peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
dim(peaks_tmpX)
head(peaks_tmpX)
peaks_tmpX2 <- subset(human_keep_Final_TE4, Human_Peak_rank %in% peaks_tmpX$Human_Peak_rank )
dim(peaks_tmpX2)
head(peaks_tmpX2)
table(peaks_tmpX2$repFamily)
TRN_all_tmp <- as.data.frame(table(peaks_tmpX2$repFamily)) # 除以 dim(peaks_tmpX)[1]
dim(TRN_all_tmp)
head(TRN_all_tmp)
colnames(TRN_all_tmp) <- c("Var1","Freq_C")
TRN_all_tmpX <- merge(TRN_all_tmp,TE_repFamily_Num,all=T)
head(TRN_all_tmpX)
rownames(TRN_all_tmpX) <- TRN_all_tmpX$Var1
TRN_all_tmpX <- TRN_all_tmpX[TE_repFamily_Num$Var1,]
head(TRN_all_tmpX)
TRN_all_tmpX[, 2][is.na(TRN_all_tmpX[, 2])] <- 0
#
TRN_sub_tmp <- TRN_all_tmpX[rownames(TRN_all_tmpX) != "other", ]
#
head(TRN_all_tmpX)
head(TRN_sub_tmp)
#
TRN_all_tmpX$Freq_C <- TRN_all_tmpX$Freq_C/dim(peaks_tmpX)[1]
TRN_sub_tmp$Freq_C <- TRN_sub_tmp$Freq_C/dim(peaks_tmpX)[1]
#
TE_repFamily_Num_all <- cbind(TE_repFamily_Num_all,TRN_all_tmpX$Freq_C)
TE_repFamily_Num_sub <- cbind(TE_repFamily_Num_sub,TRN_sub_tmp$Freq_C)
#
#
#
TE_repFamily_Num_sub <- TE_repFamily_Num[-27,]
TE_repFamily_Num_all <- TE_repFamily_Num
#
for (i in 1:131) {
  peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  TRN_all_tmp <- as.data.frame(table(peaks_tmpX2$repFamily)) # 除以 dim(peaks_tmpX)[1]
  colnames(TRN_all_tmp) <- c("Var1","Freq_C")
  TRN_all_tmpX <- merge(TRN_all_tmp,TE_repFamily_Num,all=T)
  rownames(TRN_all_tmpX) <- TRN_all_tmpX$Var1
  TRN_all_tmpX[, 2][is.na(TRN_all_tmpX[, 2])] <- 0
  TRN_all_tmpX <- TRN_all_tmpX[as.character(TE_repFamily_Num$Var1),]
  TRN_sub_tmp <- TRN_all_tmpX[rownames(TRN_all_tmpX) != "other", ]
  TRN_all_tmpX$Freq_C <- TRN_all_tmpX$Freq_C/sum(TRN_all_tmpX$Freq_C)
  TRN_sub_tmp$Freq_C <- TRN_sub_tmp$Freq_C/sum(TRN_sub_tmp$Freq_C)
  TE_repFamily_Num_all <- cbind(TE_repFamily_Num_all,TRN_all_tmpX$Freq_C)
  TE_repFamily_Num_sub <- cbind(TE_repFamily_Num_sub,TRN_sub_tmp$Freq_C)
}
#
head(TE_repFamily_Num_all)
head(TE_repFamily_Num_sub)
#
TE_repFamily_Num_all2 <- TE_repFamily_Num_all
TE_repFamily_Num_sub2 <- TE_repFamily_Num_sub
#
rownames(TE_repFamily_Num_all2) <- TE_repFamily_Num_all2$Var1
rownames(TE_repFamily_Num_sub2) <- TE_repFamily_Num_sub2$Var1
#
TE_repFamily_Num_all2 <- TE_repFamily_Num_all2[,-c(1,2)]
TE_repFamily_Num_sub2 <- TE_repFamily_Num_sub2[,-c(1,2)]
#
colnames(TE_repFamily_Num_all2) <- celltype_list
colnames(TE_repFamily_Num_sub2) <- celltype_list
#
TE_repFamily_Num_all2 <- TE_repFamily_Num_all2[,ID_rank$V1]
TE_repFamily_Num_sub2 <- TE_repFamily_Num_sub2[,ID_rank$V1]
#
#
pheatmap::pheatmap(TE_repFamily_Num_all2,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(TE_repFamily_Num_all2[-27,],cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(TE_repFamily_Num_sub2,cluster_rows = F,cluster_cols = F)
#
#############
##############
#
HOPA_lift_TS_sp <- HOPA_lift[which(HOPA_lift$Human_to_TS_State=="Human_specific"),]
#HOPA_lift_TS_sp <- HOPA_lift[which(HOPA_lift$Human_to_TS_State=="Human_specific" & HOPA_lift$Human_to_Mm_State=="Human_specific"),]
#
#
TE_repFamily_Num_SP_sub <- TE_repFamily_Num[-27,]
TE_repFamily_Num_SP_all <- TE_repFamily_Num
#
for (i in 1:131) {
  peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  peaks_tmpX2 <- subset(human_keep_Final_TE4, Human_Peak_rank %in% peaks_tmpX$Human_Peak_rank )
  peaks_tmpX2 <- subset(peaks_tmpX2, Human_Peak_rank %in% HOPA_lift_TS_sp$Human_Peak_rank )
  TRN_all_tmp <- as.data.frame(table(peaks_tmpX2$repFamily)) # 除以 dim(peaks_tmpX)[1]
  colnames(TRN_all_tmp) <- c("Var1","Freq_C")
  TRN_all_tmpX <- merge(TRN_all_tmp,TE_repFamily_Num,all=T)
  rownames(TRN_all_tmpX) <- TRN_all_tmpX$Var1
  TRN_all_tmpX[, 2][is.na(TRN_all_tmpX[, 2])] <- 0
  TRN_all_tmpX <- TRN_all_tmpX[as.character(TE_repFamily_Num$Var1),]
  TRN_sub_tmp <- TRN_all_tmpX[rownames(TRN_all_tmpX) != "other", ]
  TRN_all_tmpX$Freq_C <- TRN_all_tmpX$Freq_C/sum(TRN_all_tmpX$Freq_C)
  TRN_sub_tmp$Freq_C <- TRN_sub_tmp$Freq_C/sum(TRN_sub_tmp$Freq_C)
  TE_repFamily_Num_SP_all <- cbind(TE_repFamily_Num_SP_all,TRN_all_tmpX$Freq_C)
  TE_repFamily_Num_SP_sub <- cbind(TE_repFamily_Num_SP_sub,TRN_sub_tmp$Freq_C)
}
#
head(TE_repFamily_Num_SP_all)
head(TE_repFamily_Num_SP_sub)
#
TE_repFamily_Num_SP_all2 <- TE_repFamily_Num_SP_all
TE_repFamily_Num_SP_sub2 <- TE_repFamily_Num_SP_sub
#
rownames(TE_repFamily_Num_SP_all2) <- TE_repFamily_Num_SP_all2$Var1
rownames(TE_repFamily_Num_SP_sub2) <- TE_repFamily_Num_SP_sub2$Var1
#
TE_repFamily_Num_SP_all2 <- TE_repFamily_Num_SP_all2[,-c(1,2)]
TE_repFamily_Num_SP_sub2 <- TE_repFamily_Num_SP_sub2[,-c(1,2)]
#
colnames(TE_repFamily_Num_SP_all2) <- celltype_list
colnames(TE_repFamily_Num_SP_sub2) <- celltype_list
#
TE_repFamily_Num_SP_all2 <- TE_repFamily_Num_SP_all2[,ID_rank$V1]
TE_repFamily_Num_SP_sub2 <- TE_repFamily_Num_SP_sub2[,ID_rank$V1]
#
#
pheatmap::pheatmap(TE_repFamily_Num_SP_all2,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(TE_repFamily_Num_SP_all2[-27,],cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(TE_repFamily_Num_SP_sub2,cluster_rows = F,cluster_cols = F)
#
pheatmap::pheatmap(TE_repFamily_Num_SP_all2[-27,],cluster_rows = T,cluster_cols = T)
pheatmap::pheatmap(TE_repFamily_Num_all2[-27,],cluster_rows = T,cluster_cols = T)
#
pheatmap::pheatmap(TE_repFamily_Num_SP_all2[-27,],cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(TE_repFamily_Num_all2[-27,],cluster_rows = F,cluster_cols = F)
#
save.image("Human_HOPA_TE.RData")
#
#
#
#
#
#
#
#
#
##################### #
##################### #
##################### #
#  
#  算保守性！！！
#
## 1. 加载包和数据
library(GenomicScores)
availableGScores()
.libPaths("/home/sunkeyong/R/x86_64-pc-linux-gnu-library/4.4")
#
head(human_keep_Final)
dim(human_keep_Final)
#
human_keep_Final_GRange <- makeGRangesFromDataFrame(human_keep_Final, seqnames = "chr", start.field = "ArchR_start", end.field = "end")
human_keep_Final_GRange$rank <- human_keep_Final$Human_Peak_rank
human_keep_Final_GRange
#
#
library(phastCons100way.UCSC.hg38)
library(phastCons30way.UCSC.hg38)
library(phastCons7way.UCSC.hg38)
library(phyloP100way.UCSC.hg38)
#
phast <- phastCons100way.UCSC.hg38
#
peak_cons <- gscores(phast, human_keep_Final_GRange, summaryFun = mean)
peak_cons
#
human_keep_Final$phastCons100way <- peak_cons$default
head(human_keep_Final)
#
Human_random <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/Human_Peak_validation/Human_random_site.bed",header = F,sep = "\t")
head(Human_random)
Human_random_phastCons <- gscores(phast, makeGRangesFromDataFrame(Human_random, seqnames = "V1", start.field = "V2", end.field = "V3"), summaryFun = mean)
#
table(human_keep_Final$ArchR_peakType)
#
phastCons100way_self <- human_keep_Final
dim(phastCons100way_self)
phastCons100way_self <- phastCons100way_self[!is.na(phastCons100way_self[,18]), ]
dim(phastCons100way_self)
#
Human_random_phastCons
Human_random_phastConsX <- as.data.frame(Human_random_phastCons)
head(Human_random_phastConsX)
Human_random_phastConsX <- Human_random_phastConsX[!is.na(Human_random_phastConsX[,6]), ]
dim(Human_random_phastConsX)
#
#
boxplot(phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Exonic"),]$phastCons100way,
        phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Promoter"),]$phastCons100way,
        phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Intronic"),]$phastCons100way,
        phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Distal"),]$phastCons100way,
        Human_random_phastConsX$default,palette = "Set2", 
        outline =F)
#
wilcox.test(phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Distal"),]$phastCons100way,
            Human_random_phastConsX$default)
#
median(phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Exonic"),]$phastCons100way) # 0.2169661
median(phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Promoter"),]$phastCons100way) # 0.1001996
median(phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Intronic"),]$phastCons100way) # 0.05169661
median(phastCons100way_self[which(phastCons100way_self$ArchR_peakType=="Distal"),]$phastCons100way) # 0.0506986
median(Human_random_phastConsX$default) # 0.04311377
##
#
##
############################ ################## #
############################ ################## #
############################ ################## #
# 

############## ####################### ####################### #
###################### ####################### ####################### #
##
#
table(HOPA_lift$Human_to_TS_State)
table(HOPA_lift$Human_to_Mm_State)
#
table(human_keep_Final$ChIPseeker_anno)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "CA_conserved")$Human_Peak_rank)$ChIPseeker_anno)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "CA_divergent")$Human_Peak_rank)$ChIPseeker_anno)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "Human_specific")$Human_Peak_rank)$ChIPseeker_anno)
#
table(human_keep_Final$ArchR_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "CA_conserved")$Human_Peak_rank)$ArchR_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "CA_divergent")$Human_Peak_rank)$ArchR_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "Human_specific")$Human_Peak_rank)$ArchR_peakType)
#
table(human_keep_Final$UCSC_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "CA_conserved")$Human_Peak_rank)$UCSC_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "CA_divergent")$Human_Peak_rank)$UCSC_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_TS_State %in% "Human_specific")$Human_Peak_rank)$UCSC_peakType)
#
#
table(human_keep_Final$ChIPseeker_anno)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "CA_conserved")$Human_Peak_rank)$ChIPseeker_anno)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "CA_divergent")$Human_Peak_rank)$ChIPseeker_anno)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "Human_specific")$Human_Peak_rank)$ChIPseeker_anno)
#
table(human_keep_Final$ArchR_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "CA_conserved")$Human_Peak_rank)$ArchR_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "CA_divergent")$Human_Peak_rank)$ArchR_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "Human_specific")$Human_Peak_rank)$ArchR_peakType)
#
table(human_keep_Final$UCSC_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "CA_conserved")$Human_Peak_rank)$UCSC_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "CA_divergent")$Human_Peak_rank)$UCSC_peakType)
table(subset(human_keep_Final,  Human_Peak_rank %in% subset(HOPA_lift,Human_to_Mm_State %in% "Human_specific")$Human_Peak_rank)$UCSC_peakType)
#
##
############################ ################## #
############################ ################## #
############################ ################## #
# 

#
TE_ratio_CA_conserved <- c()
TE_ratio_CA_divergent <- c()
TE_ratio_Human_specific <- c()
#
i <- 1
peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
head(peaks_tmpX)
peaks_tmpX_A <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_TS_State %in% "CA_conserved")$Human_Peak_rank))
peaks_tmpX_B <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_TS_State %in% "CA_divergent")$Human_Peak_rank))
peaks_tmpX_C <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_TS_State %in% "Human_specific")$Human_Peak_rank))
#
peaks_tmpX_A2 <- table(peaks_tmpX_A$repClass)/dim(peaks_tmpX_A)[1]
peaks_tmpX_B2 <- table(peaks_tmpX_B$repClass)/dim(peaks_tmpX_B)[1]
peaks_tmpX_C2 <- table(peaks_tmpX_C$repClass)/dim(peaks_tmpX_C)[1]
TE_ratio_CA_conserved <- cbind(TE_ratio_CA_conserved,1-peaks_tmpX_A2[4])
TE_ratio_CA_divergent <- cbind(TE_ratio_CA_divergent,1-peaks_tmpX_B2[4])
TE_ratio_Human_specific <- cbind(TE_ratio_Human_specific,1-peaks_tmpX_C2[4])
#
#
TE_ratio_CA_conserved_TS <- c()
TE_ratio_CA_divergent_TS <- c()
TE_ratio_Human_specific_TS <- c()
#
TE_ratio_CA_conserved_mm10 <- c()
TE_ratio_CA_divergent_mm10 <- c()
TE_ratio_Human_specific_mm10 <- c()
#
#
for (i in 1:131) {
  peaks_tmpX <-  HE_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("Human_Peak_rank","celltype")
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  head(peaks_tmpX)
  peaks_tmpX_A <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_TS_State %in% "CA_conserved")$Human_Peak_rank))
  peaks_tmpX_B <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_TS_State %in% "CA_divergent")$Human_Peak_rank))
  peaks_tmpX_C <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_TS_State %in% "Human_specific")$Human_Peak_rank))
  #
  peaks_tmpX_A2 <- table(peaks_tmpX_A$repClass)/dim(peaks_tmpX_A)[1]
  peaks_tmpX_B2 <- table(peaks_tmpX_B$repClass)/dim(peaks_tmpX_B)[1]
  peaks_tmpX_C2 <- table(peaks_tmpX_C$repClass)/dim(peaks_tmpX_C)[1]
  TE_ratio_CA_conserved_TS <- cbind(TE_ratio_CA_conserved_TS,1-peaks_tmpX_A2[4])
  TE_ratio_CA_divergent_TS <- cbind(TE_ratio_CA_divergent_TS,1-peaks_tmpX_B2[4])
  TE_ratio_Human_specific_TS <- cbind(TE_ratio_Human_specific_TS,1-peaks_tmpX_C2[4])
  #
  peaks_tmpX_A <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_Mm_State %in% "CA_conserved")$Human_Peak_rank))
  peaks_tmpX_B <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_Mm_State %in% "CA_divergent")$Human_Peak_rank))
  peaks_tmpX_C <- subset(human_keep_Final_TE3, Human_Peak_rank %in% intersect(peaks_tmpX$Human_Peak_rank,subset(HOPA_lift,Human_to_Mm_State %in% "Human_specific")$Human_Peak_rank))
  #
  peaks_tmpX_A2 <- table(peaks_tmpX_A$repClass)/dim(peaks_tmpX_A)[1]
  peaks_tmpX_B2 <- table(peaks_tmpX_B$repClass)/dim(peaks_tmpX_B)[1]
  peaks_tmpX_C2 <- table(peaks_tmpX_C$repClass)/dim(peaks_tmpX_C)[1]
  TE_ratio_CA_conserved_mm10 <- cbind(TE_ratio_CA_conserved_mm10,1-peaks_tmpX_A2[4])
  TE_ratio_CA_divergent_mm10 <- cbind(TE_ratio_CA_divergent_mm10,1-peaks_tmpX_B2[4])
  TE_ratio_Human_specific_mm10 <- cbind(TE_ratio_Human_specific_mm10,1-peaks_tmpX_C2[4])
}
#
#
boxplot(as.numeric(TE_ratio_CA_conserved_TS),as.numeric(TE_ratio_CA_divergent_TS),as.numeric(TE_ratio_Human_specific_TS),
        as.numeric(TE_ratio_CA_conserved_mm10),as.numeric(TE_ratio_CA_divergent_mm10),as.numeric(TE_ratio_Human_specific_mm10),
        TE_CT_ratio2$V1)
#
HE_TE_ratio_merge <- cbind(as.numeric(TE_ratio_CA_conserved_TS),as.numeric(TE_ratio_CA_divergent_TS),as.numeric(TE_ratio_Human_specific_TS),
                           as.numeric(TE_ratio_CA_conserved_mm10),as.numeric(TE_ratio_CA_divergent_mm10),as.numeric(TE_ratio_Human_specific_mm10),
                           TE_CT_ratio2$V1)
HE_TE_ratio_merge <- as.data.frame(HE_TE_ratio_merge)
colnames(HE_TE_ratio_merge) <- c("CA_conserved_TS","CA_divergent_TS","Human_specific_TS","CA_conserved_mm10","CA_divergent_mm10","Human_specific_mm10","TE_all_ratio")
head(HE_TE_ratio_merge)
#
boxplot(HE_TE_ratio_merge[,1],HE_TE_ratio_merge[,2],HE_TE_ratio_merge[,3],HE_TE_ratio_merge[,4],HE_TE_ratio_merge[,5],HE_TE_ratio_merge[,6],HE_TE_ratio_merge[,7])
#
#
save(HE_TE_ratio_merge,file = "./data_save_For_plot2/HE_TE_ratio_merge.RData")
#
plot(CA_to_TS_ratio[3,],HE_TE_ratio_merge$TE_all_ratio)
plot(CA_to_Mm_ratio[3,],HE_TE_ratio_merge$TE_all_ratio)
#
HE_CA_to_TS_ratio <- CA_to_TS_ratio
HE_CA_to_Mm_ratio <- CA_to_Mm_ratio
#
save(HE_CA_to_TS_ratio,file = "./data_save_For_plot2/HE_CA_to_TS_ratio.RData")
save(HE_CA_to_Mm_ratio,file = "./data_save_For_plot2/HE_CA_to_Mm_ratio.RData")
#
#

#
############################ ################## #
############################ ################## ############################# ################## ############################# ################## #
############################ ################## #
#
head(human_keep_Final)
head(HOPA_lift)
#
table(HOPA_lift$Human_Peak_rank==human_keep_Final$Human_Peak_rank)
rownames(HOPA_lift) <- HOPA_lift$Human_Peak_rank
length(intersect(HOPA_lift$Human_Peak_rank,human_keep_Final$Human_Peak_rank))
length(intersect(human_keep_Final$Human_Peak_rank,HOPA_lift$Human_Peak_rank))
HOPA_lift2 <- HOPA_lift[human_keep_Final$Human_Peak_rank,]
table(HOPA_lift2$Human_Peak_rank==human_keep_Final$Human_Peak_rank)
#
head(human_keep_Final)
human_keep_Final$Human_to_Treeshrew_State <- HOPA_lift2$Human_to_TS_State
human_keep_Final$Human_to_Mouse_State <- HOPA_lift2$Human_to_Mm_State
#
head(human_keep_Final)
#
write.table(human_keep_Final,file = "./Human_and_Treeshrew_peak_save/human_keep_Final.bed", sep = "\t", col.names = T, row.names = F, quote = F)
#
#
#
############################ ################## #
# 
#
table(human_keep_Final$ArchR_peakType)
table(human_keep_Final$ArchR_peakType)/sum(table(human_keep_Final$ArchR_peakType))
# 输入频数
anno_counts <- c(Promoter = 30026, Exonic = 46746, Intronic = 500059, Distal = 439551) # ArchR_peakType 选用
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




