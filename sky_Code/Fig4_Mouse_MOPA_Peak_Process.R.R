library(ChIPseeker)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ArchR)
#
setwd("/home/sunkeyong/Human_em/Cross_species/3_species")
#
Mm_peak_raw <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/MOPA_peaks_by_celltypes.csv",sep = "\t",header = T)
head(Mm_peak_raw)
tail(Mm_peak_raw$Peaks)
dim(Mm_peak_raw)
#
Mm_peak_raw <- cbind(paste("Mouse_Peak_rank",1:dim(Mm_peak_raw)[1],sep = ""),Mm_peak_raw)
head(Mm_peak_raw)
colnames(Mm_peak_raw)[1] <- "Mouse_Peak_rank"
#
#
Mm_peak_list <- cbind(Mm_peak_raw$Peaks,Mm_peak_raw$Mouse_Peak_rank)
Mm_peak_list <- as.data.frame(Mm_peak_list)
Mm_peak_list$X <- unlist(strsplit(Mm_peak_list$V1,":"))[seq(2,2*dim(Mm_peak_list)[1],2)]
Mm_peak_list$chr <- unlist(strsplit(Mm_peak_list$V1,":"))[seq(1,2*dim(Mm_peak_list)[1],2)]
Mm_peak_list$start <- unlist(strsplit(Mm_peak_list$X,"-"))[seq(1,2*dim(Mm_peak_list)[1],2)]
Mm_peak_list$end <- unlist(strsplit(Mm_peak_list$X,"-"))[seq(2,2*dim(Mm_peak_list)[1],2)]
head(Mm_peak_list)
#
Mm_peak_list2 <- cbind(Mm_peak_list$chr,Mm_peak_list$start,Mm_peak_list$end,Mm_peak_list$V2)
Mm_peak_list2 <- as.data.frame(Mm_peak_list2)
head(Mm_peak_list2)
table(Mm_peak_list2$V1)
#
Mm_peak_list2 <- subset(Mm_peak_list2, !grepl("chrM", V1))
Mm_peak_list2 <- subset(Mm_peak_list2, !grepl("chrY", V1))
#
head(Mm_peak_list2)
table(Mm_peak_list2$V1)
#
write.table(Mm_peak_list2,file = "./Mouse_MOPA_/Mm_peak_list2.raw.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#
#
Mouse_keep <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/Mouse_MOPA_/Mm_peak_list2.rmBlacklist.bed",sep = "\t",header = F)
colnames(Mouse_keep) <- c("chr","snapATAC2_start","end","Mouse_Peak_rank")
Mouse_keep$ArchR_start <- sprintf("%.0f", as.numeric(Mouse_keep$snapATAC2_start) + 1)
head(Mouse_keep)
#
# 注释peaks
# 我们在这里只需要输入 GRanges的Peak
# 我们在这里只需要输入 geneAnnotation
# 我们在这里只需要输入 genomeAnnotation
#
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ArchR)
#
Mm10_genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmusculus.UCSC.mm10)
Mm10_geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, OrgDb = org.Mm.eg.db)
#
#
#
Mouse_keep_GRange <- makeGRangesFromDataFrame(Mouse_keep, seqnames = "chr", start.field = "ArchR_start", end.field = "end")
Mouse_keep_GRange$rank <- Mouse_keep$Mouse_Peak_rank
Mouse_keep_GRange
#
# 注释
annotated_peaks <- annotatePeakList(
  peaks = Mouse_keep_GRange,
  genomeAnnotation = Mm10_genomeAnnotation,
  geneAnnotation = Mm10_geneAnnotation)
#
annotated_peaks
#
table(annotated_peaks$peakType)
max(annotated_peaks$N)
hist(annotated_peaks$GC)
which(annotated_peaks$GC > 0.75)
which(annotated_peaks$GC < 0.25)
#
annotated_peaks
#
annotated_peaks.meta <- as.data.frame(annotated_peaks)
head(annotated_peaks.meta)
#
table(annotated_peaks.meta$rank==Mouse_keep$Mouse_Peak_rank)
#
Mouse_keep$UCSC_distToGeneStart <- annotated_peaks.meta$distToGeneStart
Mouse_keep$UCSC_nearestGene <- annotated_peaks.meta$nearestGene
Mouse_keep$UCSC_peakType <- annotated_peaks.meta$peakType
Mouse_keep$UCSC_distToTSS <- annotated_peaks.meta$distToTSS
Mouse_keep$UCSC_GC <- annotated_peaks.meta$GC
#
head(Mouse_keep)
#
addArchRGenome("Mm10")
ArchR_Mm10_GeneAnnotation <- getGeneAnnotation()
ArchR_Mm10_GenomeAnnotation <- getGenomeAnnotation()
#
# 注释
annotated_peaks2 <- annotatePeakList(
  peaks = Mouse_keep_GRange,
  genomeAnnotation = ArchR_Mm10_GenomeAnnotation,
  geneAnnotation = ArchR_Mm10_GeneAnnotation)
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
table(annotated_peaks2.meta$rank==Mouse_keep$Mouse_Peak_rank)
#
Mouse_keep$ArchR_distToGeneStart <- annotated_peaks2.meta$distToGeneStart
Mouse_keep$ArchR_nearestGene <- annotated_peaks2.meta$nearestGene
Mouse_keep$ArchR_peakType <- annotated_peaks2.meta$peakType
Mouse_keep$ArchR_distToTSS <- annotated_peaks2.meta$distToTSS
Mouse_keep$ArchR_GC <- annotated_peaks2.meta$GC
#
head(Mouse_keep)
#
table(Mouse_keep$UCSC_peakType)
table(Mouse_keep$ArchR_peakType)
#
#
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(org.Mm.eg.db)
#
annotated_peaks3 <- annotatePeak(
  Mouse_keep_GRange,
  TxDb = txdb,
  tssRegion = c(-2000, 100),
  annoDb = "org.Mm.eg.db"
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
Mouse_keep$ChIPseeker_RawAnnotation <- df$annotation
Mouse_keep$ChIPseeker_anno <- df$anno4
#
head(Mouse_keep)
#
table(Mouse_keep$UCSC_peakType)
table(Mouse_keep$ArchR_peakType)
table(Mouse_keep$ChIPseeker_anno)
#
Mouse_keep_Final <- Mouse_keep[which(Mouse_keep$UCSC_GC <= 0.75),]
Mouse_keep_Final <- Mouse_keep_Final[which(Mouse_keep_Final$UCSC_GC >= 0.25),]
dim(Mouse_keep)
dim(Mouse_keep_Final)
#
Mouse_keep_Final <- Mouse_keep_Final[which(Mouse_keep_Final$ArchR_GC <= 0.8),]
Mouse_keep_Final <- Mouse_keep_Final[which(Mouse_keep_Final$ArchR_GC >= 0.2),]
dim(Mouse_keep_Final)
#
hist(Mouse_keep$UCSC_GC)
hist(Mouse_keep_Final$UCSC_GC)
#
head(Mouse_keep_Final)
#
#
Mouse_keep_Final_bed_tmp <- cbind(Mouse_keep_Final$chr,Mouse_keep_Final$snapATAC2_start,Mouse_keep_Final$end,Mouse_keep_Final$Mouse_Peak_rank)
head(Mouse_keep_Final_bed_tmp)
#
write.table(Mouse_keep_Final_bed_tmp,file = "./Mouse_MOPA_/Mouse_keep_Final_bed_tmp.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#
MOPA_TE_anno <- read.csv("./Mouse_MOPA_/Mouse_keep_Final_bed_tmp.dfam.TEonly_overlap_final.txt",sep = "\t",header = F)
head(MOPA_TE_anno)
table(MOPA_TE_anno$V8)
table(MOPA_TE_anno$V11)
table(MOPA_TE_anno$V12)
dim(MOPA_TE_anno) # 448612
#
MOPA_TE_anno2 <- MOPA_TE_anno[which(MOPA_TE_anno$V17 > 250),]
head(MOPA_TE_anno2)
dim(MOPA_TE_anno2)
#
MOPA_TE_anno2 <- cbind(MOPA_TE_anno2$V4,MOPA_TE_anno2$V13,MOPA_TE_anno2$V14,MOPA_TE_anno2$V8)
head(MOPA_TE_anno2)
colnames(MOPA_TE_anno2) <- c("Mouse_Peak_rank","repClass","repFamily","repName")
MOPA_TE_anno2 <- as.data.frame(MOPA_TE_anno2)
head(MOPA_TE_anno2)
#
length(unique(MOPA_TE_anno2$rank))
#
head(Mouse_keep_Final)
#
Mouse_keep_Final_TE <- merge(Mouse_keep_Final,MOPA_TE_anno2,by="Mouse_Peak_rank", all = TRUE)
dim(Mouse_keep_Final_TE)
head(Mouse_keep_Final_TE)
#
dim(Mouse_keep_Final_TE)
dim(MOPA_TE_anno2)
dim(Mouse_keep_Final)
length(unique(MOPA_TE_anno2$Mouse_Peak_rank))
length(unique(Mouse_keep_Final_TE$Mouse_Peak_rank))
#
#
#
#
Mouse_keep_Final_TE2 <- Mouse_keep_Final_TE
Mouse_keep_Final_TE2$repClass[is.na(Mouse_keep_Final_TE2$repClass)] <- "other"
Mouse_keep_Final_TE2$tmp_rep <- paste(Mouse_keep_Final_TE2$Mouse_Peak_rank,Mouse_keep_Final_TE2$repClass,sep = "__")
Mouse_keep_Final_TE3 <- Mouse_keep_Final_TE2[!duplicated(Mouse_keep_Final_TE2$tmp_rep), ]
head(Mouse_keep_Final_TE2)
dim(Mouse_keep_Final_TE3)
dim(Mouse_keep_Final_TE2)
dim(Mouse_keep_Final_TE)
#
table(Mouse_keep_Final_TE2$repClass,Mouse_keep_Final_TE2$ChIPseeker_anno)
table(Mouse_keep_Final_TE2$ChIPseeker_anno)
colSums(table(Mouse_keep_Final_TE2$repClass,Mouse_keep_Final_TE2$ChIPseeker_anno)[-4,])/table(Mouse_keep_Final_TE2$ChIPseeker_anno)
#
#
colSums(table(Mouse_keep_Final_TE3$repClass,Mouse_keep_Final_TE3$ChIPseeker_anno)[-4,])/table(Mouse_keep_Final_TE$ChIPseeker_anno)
#   Distal    Exonic   Intronic  Promoter 
# 0.21177890 0.06369888 0.15682374 0.05886291 
#
colSums(table(Mouse_keep_Final_TE2$repClass,Mouse_keep_Final_TE2$ChIPseeker_anno)[-4,])/table(Mouse_keep_Final_TE$ChIPseeker_anno)
#
# Distal    Exonic  Intronic  Promoter 
# 0.21949106 0.06537278 0.16214495 0.0602664
#
colSums(table(Mouse_keep_Final_TE3$repClass,Mouse_keep_Final_TE3$ArchR_peakType)[-4,])/table(Mouse_keep_Final_TE$ChIPseeker_anno)
#   Distal    Exonic   Intronic  Promoter 
#0.24147245 0.01345910 0.14948870 0.02591177
#
colSums(table(Mouse_keep_Final_TE2$repClass,Mouse_keep_Final_TE2$ArchR_peakType)[-4,])/table(Mouse_keep_Final_TE$ChIPseeker_anno)
#
# Distal    Exonic  Intronic  Promoter 
# 0.25005662 0.01413771 0.15444242 0.02655956 
#
#
#
#
#
head(Mm_peak_raw)
dim(Mm_peak_raw)
#
#
colnames(Mm_peak_raw)
#
celltype_list <- colnames(Mm_peak_raw)[3:73]
#
i <-13
peaks_tmpX <-  Mm_peak_raw[c(1,i+2)]
colnames(peaks_tmpX) <- c("Mouse_Peak_rank","celltype")
head(peaks_tmpX)
peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
dim(peaks_tmpX)
peaks_tmpX2 <- subset(Mouse_keep_Final_TE3, Mouse_Peak_rank %in% peaks_tmpX$Mouse_Peak_rank )
dim(peaks_tmpX2)
head(peaks_tmpX2)
table(peaks_tmpX2$repClass)/dim(peaks_tmpX)[1]
t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ChIPseeker_anno))/as.numeric(table(peaks_tmpX2$ChIPseeker_anno)))
#
peaks_tmpX3 <- peaks_tmpX2[which(peaks_tmpX2$repClass!="other"),]
length(peaks_tmpX3$Mouse_Peak_rank)/dim(peaks_tmpX)[1]
length(unique(peaks_tmpX3$Mouse_Peak_rank))/dim(peaks_tmpX)[1]
#
TE_ratio_ArchR <- c()
TE_ratio_UCSC <- c()
TE_ratio_ChipSeeker <- c()
TE_CT_ratio <- c()
#
for (i in 1:71) {
  peaks_tmpX <-  Mm_peak_raw[c(1,i+2)]
  colnames(peaks_tmpX) <- c("Mouse_Peak_rank","celltype")
  peaks_tmpX <- peaks_tmpX[which(peaks_tmpX$celltype=="True"),]
  peaks_tmpX2 <- subset(Mouse_keep_Final_TE3, Mouse_Peak_rank %in% peaks_tmpX$Mouse_Peak_rank )
  TE_ratio_ArchR <- cbind(TE_ratio_ArchR,t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ArchR_peakType))/as.numeric(table(peaks_tmpX2$ArchR_peakType))))
  TE_ratio_UCSC <- cbind(TE_ratio_UCSC,t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$UCSC_peakType))/as.numeric(table(peaks_tmpX2$UCSC_peakType))))
  TE_ratio_ChipSeeker <- cbind(TE_ratio_ChipSeeker,t(t(table(peaks_tmpX2$repClass,peaks_tmpX2$ChIPseeker_anno))/as.numeric(table(peaks_tmpX2$ChIPseeker_anno))))
  peaks_tmpX3 <- peaks_tmpX2[which(peaks_tmpX2$repClass!="other"),]
  TE_CT_ratio <- cbind(TE_CT_ratio,length(unique(peaks_tmpX3$Mouse_Peak_rank))/dim(peaks_tmpX)[1])
}
#
#
#
#
#
TE_ratio <- TE_ratio_ArchR
#
colnames(TE_ratio) <- paste(rep(celltype_list,each=4),colnames(TE_ratio) ,sep="_")
TE_ratio_Distal <- TE_ratio[,paste(celltype_list,"Distal" ,sep="_")]
TE_ratio_Exonic <- TE_ratio[,paste(celltype_list,"Exonic" ,sep="_")]
TE_ratio_Intronic <- TE_ratio[,paste(celltype_list,"Intronic" ,sep="_")]
TE_ratio_Promoter <- TE_ratio[,paste(celltype_list,"Promoter" ,sep="_")]
#
barplot(1-TE_ratio_Distal["other",])
barplot(1-TE_ratio_Exonic["other",])
barplot(1-TE_ratio_Intronic["other",])
barplot(1-TE_ratio_Promoter["other",])
#
#View(t(TE_ratio_Distal))
#
boxplot(1-TE_ratio_Distal["other",],
        1-TE_ratio_Intronic["other",],
        1-TE_ratio_Exonic["other",],
        1-TE_ratio_Promoter["other",])
#
wilcox.test(1-TE_ratio_Exonic["other",],
            1-TE_ratio_Promoter["other",])
#
# save.image("Mouse_MOPA_TE.RData")
#
median(1-TE_ratio_Distal["other",]) # 0.1254227 
median(1-TE_ratio_Intronic["other",]) # 0.07368003
median(1-TE_ratio_Exonic["other",]) # 0.01096823
median(1-TE_ratio_Promoter["other",]) # 0.01094645
#
head(TE_CT_ratio)
TE_CT_ratio <- as.data.frame(t(TE_CT_ratio))
rownames(TE_CT_ratio) <- celltype_list
TE_CT_ratio$celltype <- celltype_list
head(TE_CT_ratio)
#
#
TE_CT_ratio2 <-TE_CT_ratio
TE_CT_ratio2 <- TE_CT_ratio2[celltype_list,]
head(TE_CT_ratio2)
#
View(TE_CT_ratio2)
#
barplot(TE_CT_ratio2$V1)
boxplot(TE_CT_ratio2$V1)
median(TE_CT_ratio2$V1) # 0.06332897
#
#
#
save.image("Mouse_MOPA_TE.RData")
##
#

#

#

#
#

#

#

############################## #
#
MOPA_TE <- read.csv("/home/sunkeyong/Human_em/Cross_species/Mouse_MOPA_NumC71_Downsample/MOPA_TE.txt",sep = "\t",header = F)
MOPA_TE
#
barplot(MOPA_TE$V1,col=MOPA_TE$V4)
#
# 定义分组顺序（按你给出的顺序）
group_order <- c("immune", "epithelial", "endothelium", "glia", "PNS",
                 "Cardiac_muscle_lineages", "Myocytes", "Mesoderm", "CNS_Neurons")

# 定义颜色（与顺序对应）
group_colors <- c("immune" = "#FDCDA1",
                  "epithelial" = "#00BFF1",
                  "endothelium" = "#F4CAE1",
                  "glia" = "#C9DBEC",
                  "PNS" = "#66BD61",
                  "Cardiac_muscle_lineages" = "#F6E8C1",
                  "Myocytes" = "#FF33C1",
                  "Mesoderm" = "#B3EE31",
                  "CNS_Neurons" = "#FF7F01")

# 将数据框的 V3 转换为因子，按指定顺序排列
MOPA_TE$V3 <- factor(MOPA_TE$V3, levels = group_order)
#
library(ggplot2)

ggplot(MOPA_TE, aes(x = V3, y = V1, fill = V3)) +
  geom_boxplot(outlier.shape = NA) +                 # 箱线图，隐藏离群点
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, shape = 16) +  # 实心点
  scale_fill_manual(values = group_colors) +         # 应用你指定的颜色
  scale_x_discrete(drop = FALSE) +                   # 保留因子所有水平（即使无数据）
  labs(x = "Group (V3)", y = "Value (V1)", 
       title = "Boxplot with all solid points") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 倾斜标签避免重叠
