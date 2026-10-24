#
setwd("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek")
#
#
hpeak2_alls <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/raw_back_up/HE_DS_snapATAC2_pseudoCell_macs3_q0.001_peaks_by_celltypes.csv",sep = ",",header = T)
head(hpeak2_alls)
hpeak2_alls[1:4,1:4]
rownames(hpeak2_alls) <- hpeak2_alls$Peaks
hpeak2_alls <- hpeak2_alls[,-1]
hpeak2_alls <- hpeak2_alls[,-1]
hpeak2_alls[1:4,1:4]
#
cell_pick <- c("C10.1","C10.3","C3.1","C10.8","C3.3",
               "C10.5","C10.14","C10.11","C10.2","C10.6",
               "C10.7","C10.9","C10.10")
#
hpeak2_alls_tmp <- hpeak2_alls[,cell_pick]
hpeak2_alls[1:4,1:4]
#
HE_peak_raw_HAR_num <- as.data.frame(lapply(hpeak2_alls_tmp, function(x) as.numeric(x == "True")))
HE_peak_raw_HAR_num[1:4,1:4]
rownames(HE_peak_raw_HAR_num) <- rownames(hpeak2_alls)
HE_peak_raw_HAR_num[1:4,1:4]
#
HE_peak_raw_HAR_num$x <- rowSums(HE_peak_raw_HAR_num)
HE_peak_raw_HAR_num[1:4,1:4]
#
HE_peak_raw_HAR_num2 <- HE_peak_raw_HAR_num[which(HE_peak_raw_HAR_num$x > 0),]
HE_peak_raw_HAR_num2[1:4,]
dim(HE_peak_raw_HAR_num2)
#
HE_immune_peak_raw <- rownames(HE_peak_raw_HAR_num2)
write.table(HE_immune_peak_raw,file = "/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/human_immune_peak_raw.bed",
            sep = "\t", col.names = F, row.names = F, quote = F)
#
#
#
tpeak2_alls <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/raw_back_up/TS_DS_snapATAC2_pseudoCell_macs3_q0.001_peaks_by_celltypes.csv",sep = ",",header = T)
head(tpeak2_alls)
tpeak2_alls[1:4,1:4]
rownames(tpeak2_alls) <- tpeak2_alls$Peaks
tpeak2_alls <- tpeak2_alls[,-1]
tpeak2_alls <- tpeak2_alls[,-1]
tpeak2_alls[1:4,1:4]
#
cell_pick <- c("C24.21","C24.4","C24.8","C24.12","C24.14","C24.10","C24.1", "C24.22","C24.23","C24.18","C24.3","C24.17","C24.11","C24.15","C24.16")
#
tpeak2_alls_tmp <- tpeak2_alls[,cell_pick]
tpeak2_alls[1:4,1:4]
tpeak2_alls_tmp[1:4,1:4]
#
TS_peak_raw_HAR_num <- as.data.frame(lapply(tpeak2_alls_tmp, function(x) as.numeric(x == "True")))
TS_peak_raw_HAR_num[1:4,1:4]
rownames(TS_peak_raw_HAR_num) <- rownames(tpeak2_alls)
TS_peak_raw_HAR_num[1:4,1:4]
#
TS_peak_raw_HAR_num$x <- rowSums(TS_peak_raw_HAR_num)
#
TS_peak_raw_HAR_num2 <- TS_peak_raw_HAR_num[which(TS_peak_raw_HAR_num$x > 0),]
TS_peak_raw_HAR_num2[1:4,]
head(TS_peak_raw_HAR_num2)
dim(TS_peak_raw_HAR_num2)
#
TS_immune_peak_raw <- rownames(TS_peak_raw_HAR_num2)
write.table(TS_immune_peak_raw,file = "/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/Treeshrew_immune_peak_raw.bed",
            sep = "\t", col.names = F, row.names = F, quote = F)
#
#
#
######### ######### ######### #########  # 
######### ######### ######### #########  # 
#
TS_bed <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/Human_and_Treeshrew_peak_save/TreeShrew_keep_Final.bed",header = T,sep = "\t")
head(TS_bed)
dim(TS_bed)
HE_bed <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/Human_and_Treeshrew_peak_save/human_keep_Final.bed",header = T,sep = "\t")
head(HE_bed)
dim(HE_bed)
#
#
TS_bed$PEAKS <- paste(paste(TS_bed$chr,TS_bed$snapATAC2_start,sep=":"),TS_bed$end,sep = "-")
TS_bedX <- subset(TS_bed, PEAKS %in% rownames(TS_peak_raw_HAR_num2))
length(rownames(TS_peak_raw_HAR_num2))
dim(TS_bedX) # 270659
table(TS_bedX$TreeShrew_to_Hg38_State)
# CA_conserved CA_divergent  TS_specific 
# 141484        64162        65013 
# 141484 + 64162
# 205646
# 
#
head(HE_bed)
dim(HE_bed)
HE_bed$PEAKS <- paste(paste(HE_bed$chr,HE_bed$snapATAC2_start,sep=":"),HE_bed$end,sep = "-")
head(HE_bed$PEAKS)
HE_bedX <- subset(HE_bed, PEAKS %in% rownames(HE_peak_raw_HAR_num2))
length(rownames(HE_peak_raw_HAR_num2))
dim(244562) # 244562
table(HE_bedX$Human_to_Treeshrew_State)
# CA_conserved   CA_divergent Human_specific 
# 131446          39429          73687 
#
head(HE_bedX)
#
#
# 
save(HE_bedX,file = "Human_immune_peak_full_anno.RData")
save(TS_bedX,file = "Treeshrew_immune_peak_full_anno.RData")
#
######### ######### ######### #########  # 
######### ######### ######### #########  # 
######### ######### ######### #########  # 
TS_HS_con_list <- read.csv("/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/TreeShrew_self/immune_analysis_minimatch0.1_0.002_TE_version/TS_Human_conserved_min010/TS_to_Human.conserved_activity.bestPairs.Min0.10.recip10.tsv",header = F,sep = "\t")
#
head(TS_HS_con_list)
dim(TS_HS_con_list)
TS_HS_con_list2 <- TS_HS_con_list[which(TS_HS_con_list$V9 > 250),]
dim(TS_HS_con_list2)
#
head(TS_HS_con_list2)
#
TS_HS_con_list2 <- TS_HS_con_list2[!duplicated(TS_HS_con_list2$V4), ]
TS_HS_con_list2 <- TS_HS_con_list2[!duplicated(TS_HS_con_list2$V8), ]
#
dim(TS_HS_con_list2)
#
#
#
HE_bed_keep <- subset(HE_bed, Human_Peak_rank %in% TS_HS_con_list2$V8)
TS_bed_keep <- subset(TS_bed, TreeShrew_Peak_rank %in% TS_HS_con_list2$V4)
#
head(HE_bed_keep)
head(TS_bed_keep)
#
dim(HE_bed_keep)
dim(TS_bed_keep)
#
#
HE_bed_keep$PEAKS <- paste(paste(HE_bed_keep$chr,HE_bed_keep$snapATAC2_start,sep=":"),HE_bed_keep$end,sep = "-")
TS_bed_keep$PEAKS <- paste(paste(TS_bed_keep$chr,TS_bed_keep$snapATAC2_start,sep=":"),TS_bed_keep$end,sep = "-")
#
# HE_bed_keep <- subset(HE_bed_keep, PEAKS %in% TS_HS_con_list2$V8)
# TS_bed_keep <- subset(TS_bed_keep, PEAKS %in% TS_HS_con_list2$V4)
#
#
head(HE_bed_keep)
head(TS_bed_keep)
#
dim(HE_bed_keep)
dim(TS_bed_keep)
#
HE_immune_peak_raw
TS_immune_peak_raw
#
length(HE_immune_peak_raw)
length(TS_immune_peak_raw)
#
#
TS_bed_keep2 <- subset(TS_bed_keep, PEAKS %in% TS_immune_peak_raw)
head(TS_bed_keep2)
#
head(TS_HS_con_list2)
#
HE_keepX <- subset(TS_HS_con_list2, V4 %in% TS_bed_keep2$TreeShrew_Peak_rank)
#
#
head(HE_bed_keep)
HE_bed_keep2 <- subset(HE_bed_keep, Human_Peak_rank %in% HE_keepX$V8)
#
head(HE_bed_keep2)
head(TS_bed_keep2)
#
dim(HE_bed_keep2)
dim(TS_bed_keep2)
#
#
Keep_con_list <- HE_keepX
#
dim(Keep_con_list)
head(Keep_con_list)
#
head(HE_bed_keep2)
head(TS_bed_keep2)
#
dim(HE_bed_keep2)
dim(TS_bed_keep2)
#
#
rownames(HE_bed_keep2) <- HE_bed_keep2$Human_Peak_rank
rownames(TS_bed_keep2) <- TS_bed_keep2$TreeShrew_Peak_rank
#
HE_bed_keep2X <- HE_bed_keep2[Keep_con_list$V8,]
TS_bed_keep2X <- TS_bed_keep2[Keep_con_list$V4,]
#
Keep_con_list2 <- Keep_con_list
#
Keep_con_list2$V1 <- TS_bed_keep2X$chr
Keep_con_list2$V2 <- TS_bed_keep2X$snapATAC2_start
Keep_con_list2$V3 <- TS_bed_keep2X$end
#
Keep_con_list2$V5 <- HE_bed_keep2X$chr
Keep_con_list2$V6 <- HE_bed_keep2X$snapATAC2_start
Keep_con_list2$V7 <- HE_bed_keep2X$end
#
head(Keep_con_list2)
#
colnames(Keep_con_list2) <- c("TS_chr","TS_start","TS_end","TreeShrew_Peak_rank","HE_chr","HE_start","HE_end","Human_Peak_rank","overlap_length")
#
head(Keep_con_list2)
#
dim(Keep_con_list2)
Keep_con_list2$orth_index <- paste("orth_index",1:108744,sep="")
head(Keep_con_list2)
dim(Keep_con_list2)
#
#
orth_hg38 <- Keep_con_list2[,c(5,6,7,10)]
orth_TS <- Keep_con_list2[,c(1,2,3,10)]
#
head(orth_hg38)
head(orth_TS)
#
write.table(orth_hg38,file = "/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/orth.human.raw.bed",
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(orth_TS,file = "/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/orth.ts.raw.bed",
            sep = "\t", col.names = F, row.names = F, quote = F)
#
#
write.table(Keep_con_list2,file = "/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step1_conserved_peak_seek/immune_human_TS_Keep_con_list2.bed",
            sep = "\t", col.names = T, row.names = F, quote = F)
#

#

#



#


#





#
barplot(c(270659,205646,108744)) #immune_all_seqorth_peakorth.num.pdf
#