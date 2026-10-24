library(data.table)
library(ggplot2)
#
ratio_dt <- as.data.table(CA_to_Hg38_ratio, keep.rownames = "category")
ratio_dt <- as.data.table(CA_to_Mm_ratio, keep.rownames = "category")
#
ratio_long <- melt(
  ratio_dt,
  id.vars = "category",
  variable.name = "subclass",
  value.name = "ratio"
)
#
head(ratio_long)
#
ratio_long[, category := factor(
  category,
  levels = c("CA_conserved", "CA_divergent", "TS_specific")
)]
#
ratio_long[, subclass := factor(subclass, levels = colnames(CA_to_Mm_ratio))]
#
my_cols <- c(
  "CA_conserved" = "#EAA0A0",
  "CA_divergent" = "#46ADD1",
  "TS_specific" = "#2F6EAF"
)
#
p1 <- ggplot(ratio_long, aes(x = subclass, y = ratio, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_flip() +
  scale_fill_manual(values = my_cols) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = NULL, y = "Fraction", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.position = "top"
  )
#
p1
#
#Treeshrew_to_hg38_ratio <- ratio_long
#Treeshrew_to_Mm10_ratio <- ratio_long
#
Treeshrew_to_hg38_ratio <- as.data.frame(Treeshrew_to_hg38_ratio)
Treeshrew_to_hg38_ratio$class <- unlist(strsplit(as.character(Treeshrew_to_hg38_ratio$subclass),".",fixed=T))[seq(1,2*dim(Treeshrew_to_hg38_ratio)[1],2)]
Treeshrew_to_Mm10_ratio <- as.data.frame(Treeshrew_to_Mm10_ratio)
Treeshrew_to_Mm10_ratio$class <- unlist(strsplit(as.character(Treeshrew_to_Mm10_ratio$subclass),".",fixed=T))[seq(1,2*dim(Treeshrew_to_Mm10_ratio)[1],2)]
##
wilcox.test(Treeshrew_to_hg38_ratio[which(Treeshrew_to_hg38_ratio$class=="C5" & Treeshrew_to_hg38_ratio$category=="TS_specific"),]$ratio,
            Treeshrew_to_hg38_ratio[which(Treeshrew_to_hg38_ratio$class=="C24" & Treeshrew_to_hg38_ratio$category=="TS_specific"),]$ratio)
#
wilcox.test(Treeshrew_to_Mm10_ratio[which(Treeshrew_to_Mm10_ratio$class=="C5" & Treeshrew_to_Mm10_ratio$category=="TS_specific"),]$ratio,
            Treeshrew_to_Mm10_ratio[which(Treeshrew_to_Mm10_ratio$class=="C24" & Treeshrew_to_Mm10_ratio$category=="TS_specific"),]$ratio)
#
plot(Treeshrew_to_hg38_ratio[which(Treeshrew_to_hg38_ratio$class=="C5" & Treeshrew_to_hg38_ratio$category=="TS_specific"),]$ratio,
     Treeshrew_to_Mm10_ratio[which(Treeshrew_to_Mm10_ratio$class=="C5" & Treeshrew_to_Mm10_ratio$category=="TS_specific"),]$ratio)
plot(Treeshrew_to_hg38_ratio[which(Treeshrew_to_hg38_ratio$class=="C24" & Treeshrew_to_hg38_ratio$category=="TS_specific"),]$ratio,
     Treeshrew_to_Mm10_ratio[which(Treeshrew_to_Mm10_ratio$class=="C24" & Treeshrew_to_Mm10_ratio$category=="TS_specific"),]$ratio)
#
#
#save(Treeshrew_to_hg38_ratio,file = "./data_save_For_plot/Treeshrew_to_hg38_ratio.RData")
#save(Treeshrew_to_Mm10_ratio,file = "./data_save_For_plot/Treeshrew_to_Mm10_ratio.RData")
#
plot(Treeshrew_to_hg38_ratio$ratio,Treeshrew_to_Mm10_ratio$ratio)
## 都是 W = 0, p-value = 1.851e-14
#
#
#
library(ggplot2)
#
TS.co.for.Main <- c("#DDAD3B","#00BFFF","#9970AB","#00FF7F",## 1-4
                    "#FDCDAC","#3288BD","#32CD32","#48D1CC",## 5-8
                    "#FF33CC","#B15928","#6495ED","#D6604D",## 9-12
                    "#0000FF","#C2A5CF","#EE82EE","#FF7F00",## 13-16
                    "#9E0142","#B3CDE3","#87CEFA","#00FFFF",## 17-20
                    "#CDE8B9","#F4CAE4","#0099FF","#DE77AE")## 21-24
#
names(TS.co.for.Main) <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24")
#
df_hs <- subset(Treeshrew_to_hg38_ratio, category == "TS_specific")
head(df_hs)
#
df_hs$class <- factor(
  df_hs$class,
  levels = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24")
)
#
ggplot(df_hs, aes(x = class, y = ratio, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = TS.co.for.Main) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "Human-specific ratio") +
  theme(legend.position = "none")
#
#
#
#
###
####### ####### #######. #
# TE 含量 弄 TreeShrew_to_Hg38
# Full peak set
table(TOPA_lift$TreeShrew_to_Hg38_State)
#
table(TreeShrew_keep_Final_TE3$repClass)/sum(table(TreeShrew_keep_Final_TE3$repClass))
#
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Hg38_State=="TS_specific"),]
targetTE1 <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE1$repClass)/sum(table(targetTE1$repClass))
#.    DNA       LINE        LTR       None       SINE 
# 0.03410886 0.16304107 0.19852702 0.53805514 0.06626791  
#
table(TreeShrew_keep_Final_TE3$repClass)
table(targetTE1$repClass)
#
#
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Hg38_State=="CA_conserved"),]
targetTE2 <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE2$repClass)/sum(table(targetTE2$repClass))
#.     DNA       LINE        LTR       None       SINE 
#  0.006739783 0.038396752 0.029828737 0.917459516 0.007575213 
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Hg38_State=="CA_divergent"),]
targetTE3 <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE3$repClass)/sum(table(targetTE3$repClass))
#.     DNA       LINE        LTR       None       SINE 
# 0.01061612 0.07626639 0.06233502 0.82845037 0.02233209 
#
HE_to_TS_TE_ratio <- cbind(table(TreeShrew_keep_Final_TE3$repClass)/sum(table(TreeShrew_keep_Final_TE3$repClass)),
                           table(targetTE2$repClass)/sum(table(targetTE2$repClass)),
                           table(targetTE3$repClass)/sum(table(targetTE3$repClass)),
                           table(targetTE1$repClass)/sum(table(targetTE1$repClass)))
HE_to_TS_TE_ratio <- as.data.frame(HE_to_TS_TE_ratio)
colnames(HE_to_TS_TE_ratio) <- c("Full set","CA_conserved","CA_divergent","TS_specific_to_Hg38")
HE_to_TS_TE_ratio
HE_to_TS_TE_ratio2 <- HE_to_TS_TE_ratio[-4,]
HE_to_TS_TE_ratio2
#
#     Full set     CA_conserved  CA_divergent       TS_specific_to_Hg38
# DNA  0.006840732  0.005916022  0.008096147         0.006777595
# LINE 0.054372190  0.023462404  0.049837857         0.117033517
# LTR  0.034995769  0.020597765  0.034547307         0.061866409
# SINE 0.010250601  0.002773135  0.008323416         0.026565723
#
library(data.table)
library(ggplot2)
#
#
# 转成长表
df <- as.data.table(HE_to_TS_TE_ratio2, keep.rownames = "TE_class")

df_long <- melt(
  df,
  id.vars = "TE_class",
  variable.name = "Category",
  value.name = "Ratio"
)

# 顺序
df_long$TE_class <- factor(
  df_long$TE_class,
  levels = c("DNA", "LINE", "LTR", "SINE")
)

df_long$Category <- factor(
  df_long$Category,
  levels = c("Full set", "CA_conserved", "CA_divergent", "TS_specific_to_Hg38")
)

# 颜色
TE_cols <- c(
  "LTR"  = "#EB5A2E",
  "LINE" = "#F5BC5B",
  "SINE" = "#AFACDA",
  "DNA"  = "#5D3A92"
)

# 画气泡图
p <- ggplot(df_long, aes(x = Category, y = TE_class)) +
  geom_point(
    aes(size = Ratio, fill = TE_class),
    shape = 21, color = "black", stroke = 0.4
  ) +
  scale_fill_manual(values = TE_cols) +
  scale_size(range = c(4, 16)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = NULL, size = "Ratio", fill = "TE class") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
#
p
#
Treeshrew_to_hg38_TE <- df_long
#
#
#########
# 弄小鼠的
####### ####### #######. #
# Full peak set
table(TOPA_lift$TreeShrew_to_Mm_State)
#
table(TreeShrew_keep_Final_TE3$repClass)/sum(table(TreeShrew_keep_Final_TE3$repClass))
#
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Mm_State=="TS_specific"),]
targetTE1 <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE1$repClass)/sum(table(targetTE1$repClass))
#.    DNA       LINE        LTR       None       SINE 
# 0.02794986 0.14283013 0.16375082 0.61035913 0.05511005 
table(TreeShrew_keep_Final_TE3$repClass)
table(targetTE1$repClass)
#
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Mm_State=="CA_conserved"),]
targetTE2 <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE2$repClass)/sum(table(targetTE2$repClass))
#.     DNA       LINE        LTR       None       SINE 
#  0.004254350 0.023503690 0.015445932 0.952865082 0.003930947
targetX <- TOPA_lift[which(TOPA_lift$TreeShrew_to_Mm_State=="CA_divergent"),]
targetTE3 <- subset(TreeShrew_keep_Final_TE3, TreeShrew_Peak_rank %in% targetX$TreeShrew_Peak_rank)
table(targetTE3$repClass)/sum(table(targetTE3$repClass))
#.     DNA       LINE        LTR       None       SINE 
# 0.007479595 0.044992731 0.034406714 0.904135658 0.008985303  
#
HE_to_Mm_TE_ratio <- cbind(table(TreeShrew_keep_Final_TE3$repClass)/sum(table(TreeShrew_keep_Final_TE3$repClass)),
                           table(targetTE2$repClass)/sum(table(targetTE2$repClass)),
                           table(targetTE3$repClass)/sum(table(targetTE3$repClass)),
                           table(targetTE1$repClass)/sum(table(targetTE1$repClass)))
HE_to_Mm_TE_ratio <- as.data.frame(HE_to_Mm_TE_ratio)
colnames(HE_to_Mm_TE_ratio) <- c("Full set","CA_conserved","CA_divergent","TS_specific_to_Mm10")
HE_to_Mm_TE_ratio
HE_to_Mm_TE_ratio2 <- HE_to_Mm_TE_ratio[-4,]
HE_to_Mm_TE_ratio2
#
#      Full set    CA_conserved CA_divergent       TS_specific_to_Mm10
# DNA  0.006840732  0.004071729  0.005944547          0.00829862
# LINE 0.054372190  0.013256361  0.024566242          0.08579338
# LTR  0.034995769  0.008812784  0.018230295          0.05369283
# SINE 0.010250601  0.001464149  0.002660993          0.01768806
#
library(data.table)
library(ggplot2)
#
#
# 转成长表
df <- as.data.table(HE_to_Mm_TE_ratio2, keep.rownames = "TE_class")

df_long <- melt(
  df,
  id.vars = "TE_class",
  variable.name = "Category",
  value.name = "Ratio"
)

# 顺序
df_long$TE_class <- factor(
  df_long$TE_class,
  levels = c("DNA", "LINE", "LTR", "SINE")
)

df_long$Category <- factor(
  df_long$Category,
  levels = c("Full set", "CA_conserved", "CA_divergent", "TS_specific_to_Mm10")
)

# 颜色
TE_cols <- c(
  "LTR"  = "#EB5A2E",
  "LINE" = "#F5BC5B",
  "SINE" = "#AFACDA",
  "DNA"  = "#5D3A92"
)

# 画气泡图
p <- ggplot(df_long, aes(x = Category, y = TE_class)) +
  geom_point(
    aes(size = Ratio, fill = TE_class),
    shape = 21, color = "black", stroke = 0.4
  ) +
  scale_fill_manual(values = TE_cols) +
  scale_size(range = c(4, 16)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = NULL, size = "Ratio", fill = "TE class") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
#
p
#
#
TreeShrew_to_Mm_TE <- df_long
#
#
Treeshrew_to_hg38_TE
TreeShrew_to_Mm_TE
#
save(Treeshrew_to_hg38_TE,file = "./data_save_For_plot/Treeshrew_to_hg38_TE.RData")
save(TreeShrew_to_Mm_TE,file = "./data_save_For_plot/TreeShrew_to_Mm_TE.RData")
#
#
############ ########
#
head(TE_CT_ratio2)
#
#
ID_rank <- read.csv("/home/sunkeyong/Multispecies_analysis/M_TreeShrew/Analysis_T0808/pseudo_bulk_analysis/bed_dedu/ArchR_showing/GWAS_celltype_rank.csv",sep = "\t")
head(ID_rank)
#
#
table(ID_rank$Subtype_ID_Num==TE_CT_ratio2$celltype)
#
TE_CT_ratio2$color <- ID_rank$rank_color
head(TE_CT_ratio2)
#
barplot(TE_CT_ratio2$V1,col=TE_CT_ratio2$color)
#
TE_CT_ratio2$class <- unlist(strsplit(TE_CT_ratio2$celltype,".",fixed=T))[seq(1,2*197,2)]
head(TE_CT_ratio2)
#
library(ggplot2)
#
TE_CT_ratio2$class <- factor(
  TE_CT_ratio2$class,
  levels = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24")
)
#
ggplot(TE_CT_ratio2, aes(x = class, y = V1, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8, color = "black") +
  scale_fill_manual(values = TS.co.for.Main) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "TE ratio") +
  theme(legend.position = "none")
#
head(TE_CT_ratio2)
C5 <- TE_CT_ratio2[which(TE_CT_ratio2$class=="C5"),]$V1
C24 <- TE_CT_ratio2[which(TE_CT_ratio2$class=="C24"),]$V1
#
wilcox.test(C5,C24)
wilcox.test(C24,C5) # p-value = 6.49e-11
#
#
TE_CT_ratio2$V1 #
#
ratio_long2X <- Treeshrew_to_hg38_ratio[which(Treeshrew_to_hg38_ratio$category=="TS_specific"),]
#
plot(TE_CT_ratio2$V1,ratio_long2X$ratio)
#
#
Treeshrew_TE_ratio <- TE_CT_ratio2
save(Treeshrew_TE_ratio,file = "Treeshrew_TE_ratio.RData")
#
#






























#
all_counts <- c(DNA=8472, LINE=67338, LTR=43341, `None-TE`=1106618, SINE=12695)
# sp_counts  <- c(DNA=5333, LINE=55134, LTR=34505, `None-TE`=536298, SINE=11367) # TS -to- mouse
sp_counts  <- c(DNA=1991, LINE=34380, LTR=18174, `None-TE`=231413, SINE=7804) # TS -to- hg38
#
non_sp_counts <- all_counts - sp_counts

calc_or <- function(te_name) {
  a <- sp_counts[te_name]
  b <- sum(sp_counts) - a
  c <- non_sp_counts[te_name]
  d <- sum(non_sp_counts) - c
  
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  rownames(mat) <- c("HumanSpec", "NonHumanSpec")
  colnames(mat) <- c(te_name, paste0("not_", te_name))
  
  ft <- fisher.test(mat)
  
  data.frame(
    TE = te_name,
    OR = unname(ft$estimate),
    p.value = ft$p.value,
    CI.low = ft$conf.int[1],
    CI.high = ft$conf.int[2]
  )
}

res <- do.call(rbind, lapply(names(sp_counts), calc_or))
res
res
# # TS -to- mouse
# TE        OR      p.value    CI.low   CI.high
# 1     DNA 1.5800023 3.903657e-94 1.5113790 1.6520502
# 2    LINE 4.4874034 0.000000e+00 4.3993305 4.5780898
# 3     LTR 3.7689180 0.000000e+00 3.6812475 3.8599416
# 4 None-TE 0.2255563 0.000000e+00 0.2223792 0.2288064
# 5    SINE 8.0608459 0.000000e+00 7.6150005 8.5414866
#
# TS -to- hg38
# TE        OR   p.value    CI.low   CI.high
# 1     DNA 0.9878509 0.6445661 0.9388313 1.0390832
# 2    LINE 3.6670676 0.0000000 3.6091297 3.7247422
# 3     LTR 2.4095032 0.0000000 2.3624269 2.4572584
# 4 None-TE 0.2947095 0.0000000 0.2912906 0.2982351
# 5    SINE 5.2440071 0.0000000 5.0580082 5.4375005
























#
#
load("./data_save_For_plot/human_TE_ratio.RData")
load("./data_save_For_plot/Human_to_Mm_ratio.RData")
load("./data_save_For_plot/Human_to_Mm_TE.RData")
load("./data_save_For_plot/Human_to_TS_ratio.RData")
load("./data_save_For_plot/Human_to_TS_TE.RData")
#
load("./data_save_For_plot/Treeshrew_TE_ratio.RData")
load("./data_save_For_plot/Treeshrew_to_hg38_ratio.RData")
load("./data_save_For_plot/Treeshrew_to_hg38_TE.RData")
load("./data_save_For_plot/TreeShrew_to_Mm_TE.RData")
load("./data_save_For_plot/Treeshrew_to_Mm10_ratio.RData")
#
##### ##### #####   ##### ##### ##### ##### ##### ##### 
#
Human_to_TS_ratio
Treeshrew_to_hg38_ratio
Human_to_TS_ratio$species <- "Human"
Treeshrew_to_hg38_ratio$species <- "Treeshrew"
Human_to_TS_ratio$subclass <- paste("Human",Human_to_TS_ratio$subclass)
Treeshrew_to_hg38_ratio$subclass <- paste("Treeshrew",Treeshrew_to_hg38_ratio$subclass)
Human_to_TS_ratio$class <- paste("Human",Human_to_TS_ratio$class)
Treeshrew_to_hg38_ratio$class <- paste("Treeshrew",Treeshrew_to_hg38_ratio$class)
#
median(Treeshrew_to_hg38_ratio$ratio)
median(Human_to_TS_ratio$ratio)
mean(Treeshrew_to_hg38_ratio$ratio)
mean(Human_to_TS_ratio$ratio)
#
#
HT_ratio <- rbind(Human_to_TS_ratio,Treeshrew_to_hg38_ratio)
head(HT_ratio)
table(HT_ratio$category)
HT_ratio <- HT_ratio[which(HT_ratio$category=="Human_specific" |HT_ratio$category=="TS_specific" ),]
HT_ratio
head(HT_ratio)
#
#
HT_ratio$class <- factor(
  HT_ratio$class,
  levels = c(c("Human C1","Human C2","Human C3","Human C4","Human C5","Human C6","Human C7","Human C8","Human C9","Human C10","Human C11"),
             c("Treeshrew C1","Treeshrew C2","Treeshrew C3","Treeshrew C4","Treeshrew C5","Treeshrew C6","Treeshrew C7","Treeshrew C8",
               "Treeshrew C9","Treeshrew C10","Treeshrew C11","Treeshrew C12","Treeshrew C13","Treeshrew C14","Treeshrew C15","Treeshrew C16",
               "Treeshrew C17","Treeshrew C18","Treeshrew C19","Treeshrew C20","Treeshrew C21","Treeshrew C22","Treeshrew C23","Treeshrew C24"))
)
#
# HE_11color,TS.co.for.Main
#
p1 <- ggplot(subset(HT_ratio, species %in% "Human"), aes(x = class, y = ratio, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = HE_11color) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "Human-specific ratio") +
  theme(legend.position = "none")
#
p2 <- ggplot(subset(HT_ratio, species %in% "Treeshrew"), aes(x = class, y = ratio, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = TS.co.for.Main) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "Treeshrew-specific ratio") +
  theme(legend.position = "none")
#
p1 + p2
#
#
#
##### ##### ###
Human_to_Mm_ratio
Treeshrew_to_Mm10_ratio
#
Human_to_Mm_ratio
Treeshrew_to_Mm10_ratio
Human_to_Mm_ratio$species <- "Human"
Treeshrew_to_Mm10_ratio$species <- "Treeshrew"
Human_to_Mm_ratio$subclass <- paste("Human",Human_to_Mm_ratio$subclass)
Treeshrew_to_Mm10_ratio$subclass <- paste("Treeshrew",Treeshrew_to_Mm10_ratio$subclass)
Human_to_Mm_ratio$class <- paste("Human",Human_to_Mm_ratio$class)
Treeshrew_to_Mm10_ratio$class <- paste("Treeshrew",Treeshrew_to_Mm10_ratio$class)
#
median(Treeshrew_to_Mm10_ratio$ratio)
median(Human_to_Mm_ratio$ratio)
mean(Treeshrew_to_Mm10_ratio$ratio)
mean(Human_to_Mm_ratio$ratio)
#
#
HT_ratio <- rbind(Human_to_Mm_ratio,Treeshrew_to_Mm10_ratio)
head(HT_ratio)
table(HT_ratio$category)
HT_ratio <- HT_ratio[which(HT_ratio$category=="Human_specific" |HT_ratio$category=="TS_specific" ),]
HT_ratio
head(HT_ratio)
#
#
HT_ratio$class <- factor(
  HT_ratio$class,
  levels = c(c("Human C1","Human C2","Human C3","Human C4","Human C5","Human C6","Human C7","Human C8","Human C9","Human C10","Human C11"),
             c("Treeshrew C1","Treeshrew C2","Treeshrew C3","Treeshrew C4","Treeshrew C5","Treeshrew C6","Treeshrew C7","Treeshrew C8",
               "Treeshrew C9","Treeshrew C10","Treeshrew C11","Treeshrew C12","Treeshrew C13","Treeshrew C14","Treeshrew C15","Treeshrew C16",
               "Treeshrew C17","Treeshrew C18","Treeshrew C19","Treeshrew C20","Treeshrew C21","Treeshrew C22","Treeshrew C23","Treeshrew C24"))
)
#
#
p1 <- ggplot(subset(HT_ratio, species %in% "Human"), aes(x = class, y = ratio, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = HE_11color) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "Human-specific ratio") +
  theme(legend.position = "none")
#
p2 <- ggplot(subset(HT_ratio, species %in% "Treeshrew"), aes(x = class, y = ratio, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = TS.co.for.Main) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "Treeshrew-specific ratio") +
  theme(legend.position = "none")
#
p1 + p2
#
#
##### ##### #####   ##### ##### ##### ##### ##### ##### 
##### ##### ###
Human_to_Mm_TE
Human_to_TS_TE
Treeshrew_to_hg38_TE
TreeShrew_to_Mm_TE
# 整理数据
df <- data.frame(
  row.names = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"),
  X1 = c(0.015836646, 0.083623234, 0.087518201, 0.028293849,     0.006840732, 0.054372190, 0.034995769, 0.010250601),
  X2 = c(0.006739783, 0.038396752, 0.029828737, 0.007575213,     0.005916022, 0.023462404, 0.020597765, 0.002773135),
  X3 = c(0.010616123, 0.076266393, 0.062335023, 0.02233209,      0.008096147, 0.049837857, 0.034547307, 0.008323416),
  X4 = c(0.034108861, 0.163041065, 0.198527022, 0.066267909,     0.006777595, 0.117033517, 0.061866409, 0.026565723),
  X5 = c(0.00425435, 0.02350369, 0.015445932, 0.003930947,       0.004071729, 0.013256361, 0.008812784, 0.001464149),
  X6 = c(0.007479595, 0.044992731, 0.034406714, 0.008985303,     0.005944547, 0.024566242, 0.018230295, 0.002660993),
  X7 = c(0.027949862, 0.142830132, 0.16375082, 0.055110053,      0.008298620, 0.085793379, 0.053692831, 0.017688057)
)
library(tidyr)
library(ggplot2)
library(tibble)
#
df_long <- df %>%
  rownames_to_column(var = "L") %>%
  pivot_longer(cols = starts_with("X"), names_to = "X", values_to = "value")

# 将 L 和 X 转为因子，保持原始顺序（可选）
df_long$L <- factor(df_long$L, levels = rev(rownames(df)))  # 反转 y 轴使 L1 在上方
df_long$X <- factor(df_long$X, levels = colnames(df))
#
p <- ggplot(df_long, aes(x = X, y = L)) +
  geom_point(aes(size = value, color = value), alpha = 0.8) +
  scale_size_continuous(range = c(2, 12), name = "Value") +   # 控制气泡最小/最大尺寸
  scale_color_viridis_c(option = "plasma", name = "Value") +  # 颜色渐变
  theme_minimal(base_size = 12) +
  labs(x = "X variable", y = "L variable", title = "Bubble plot of matrix values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank())

print(p)
#
#
#
#
##### ##### #####   ##### ##### ##### ##### ##### ##### 
#
human_TE_ratio
Treeshrew_TE_ratio
#
head(human_TE_ratio)
head(Treeshrew_TE_ratio)
#
human_TE_ratio$species <- "Human"
Treeshrew_TE_ratio$species <- "Treeshrew"
human_TE_ratio$celltype <- paste("Human",human_TE_ratio$celltype)
Treeshrew_TE_ratio$celltype <- paste("Treeshrew",Treeshrew_TE_ratio$celltype)
human_TE_ratio$class <- paste("Human",human_TE_ratio$class)
Treeshrew_TE_ratio$class <- paste("Treeshrew",Treeshrew_TE_ratio$class)
#
head(human_TE_ratio)
head(Treeshrew_TE_ratio)
#
HT_ratio <- rbind(human_TE_ratio,Treeshrew_TE_ratio)
head(HT_ratio)
#
#
HT_ratio$class <- factor(
  HT_ratio$class,
  levels = c(c("Human C1","Human C2","Human C3","Human C4","Human C5","Human C6","Human C7","Human C8","Human C9","Human C10","Human C11"),
             c("Treeshrew C1","Treeshrew C2","Treeshrew C3","Treeshrew C4","Treeshrew C5","Treeshrew C6","Treeshrew C7","Treeshrew C8",
               "Treeshrew C9","Treeshrew C10","Treeshrew C11","Treeshrew C12","Treeshrew C13","Treeshrew C14","Treeshrew C15","Treeshrew C16",
               "Treeshrew C17","Treeshrew C18","Treeshrew C19","Treeshrew C20","Treeshrew C21","Treeshrew C22","Treeshrew C23","Treeshrew C24"))
)
#
# HE_11color TS.co.for.Main 
#
p1 <- ggplot(subset(HT_ratio, species %in% "Human"), aes(x = class, y = V1, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = HE_11color) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "Human-specific ratio") +
  theme(legend.position = "none")
#
p2 <- ggplot(subset(HT_ratio, species %in% "Treeshrew"), aes(x = class, y = V1, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = TS.co.for.Main) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "Treeshrew-specific ratio") +
  theme(legend.position = "none")
#
p1 + p2
#
#
#

#

#

#

#

HE_11color <- c("#FDCDA1","#00BFF1","#9970A1","#9999F9",## 1-4
                "#F4CAE1","#C9DBEC","#66BD61","#F6E8C1",## 5-8
                "#FF33C1","#B3EE31","#FF7F01") ## 9-11




############################## #
############################## #
############################## #
#
library(data.table)
#			
# 输入数据 人 to 树鼩
full_counts <- c(Distal   = 439551,Exonic   = 46746,Intronic = 500059,Promoter = 30026)
cat_counts <- list(CA_conserved = c(Distal   = 199847,Exonic   = 29961,Intronic = 265466,Promoter = 19426),
                   CA_divergent = c(Distal   = 79824,Exonic   = 11427,Intronic = 94576,Promoter = 5103),
                   specific = c(Distal   = 159880,Exonic   = 5358,Intronic = 140017,Promoter = 5497))
#
# 输入数据 人 to 小鼠
full_counts <- c(Distal   = 439551,Exonic   = 46746,Intronic = 500059,Promoter = 30026)
cat_counts <- list(CA_conserved = c(Distal   = 76278,Exonic   = 15849,Intronic = 111742,Promoter = 15671),
                   CA_divergent = c(Distal   = 137719,Exonic   = 23144,Intronic = 180148,Promoter = 6330),
                   specific = c(Distal   = 225554,Exonic   = 7753,Intronic = 208169,Promoter = 8025))
#
# 输入数据 树鼩 to 人
table(TreeShrew_keep_FinalX$ArchR_peakType,TreeShrew_keep_FinalX$TreeShrew_to_Hg38_State)
#
full_counts <- c(Distal   = 555103,Exonic   = 60688,Intronic = 593920,Promoter = 28725)
cat_counts <- list(CA_conserved = c(Distal   = 229234,Exonic   = 31313,Intronic = 259810,Promoter = 15138),
                   CA_divergent = c(Distal   = 185717,Exonic   = 20811,Intronic = 195777,Promoter = 6895),
                   specific = c(Distal   = 140152,Exonic   = 8564,Intronic = 138333,Promoter = 6692))
#
# 输入数据 树鼩 to 小鼠
table(TreeShrew_keep_FinalX$ArchR_peakType,TreeShrew_keep_FinalX$TreeShrew_to_Mm10_State)
#
full_counts <- c(Distal   = 555103,Exonic   = 60688,Intronic = 593920,Promoter = 28725)
cat_counts <- list(CA_conserved = c(Distal   = 83561,Exonic   = 17468,Intronic = 104489,Promoter = 9624),
                   CA_divergent = c(Distal   = 167887,Exonic   = 25463,Intronic = 181134,Promoter = 6198),
                   specific = c(Distal   = 303655,Exonic   = 17757,Intronic = 308297,Promoter = 12903))

#
# 计算 OR 的函数：某个类别 vs 其余类别；某个context vs 其他context
calc_or_one <- function(cat_name, context_name, full_counts, cat_counts) {
  cat_vec <- cat_counts[[cat_name]]
  rest_vec <- full_counts - cat_vec
  
  a <- cat_vec[context_name]
  b <- sum(cat_vec) - a
  c <- rest_vec[context_name]
  d <- sum(rest_vec) - c
  
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  rownames(mat) <- c(cat_name, "Others")
  colnames(mat) <- c(context_name, paste0("not_", context_name))
  
  ft <- fisher.test(mat)
  
  data.table(
    Category = cat_name,
    Context = context_name,
    a = as.numeric(a),
    b = as.numeric(b),
    c = as.numeric(c),
    d = as.numeric(d),
    OR = unname(ft$estimate),
    CI_low = ft$conf.int[1],
    CI_high = ft$conf.int[2],
    p_value = ft$p.value
  )
}

# 循环所有类别和context
or_df <- rbindlist(
  lapply(names(cat_counts), function(cat_name) {
    rbindlist(
      lapply(names(full_counts), function(context_name) {
        calc_or_one(cat_name, context_name, full_counts, cat_counts)
      })
    )
  })
)

# 调整顺序，方便后面画图
or_df[, Category := factor(Category, levels = c("CA_conserved", "CA_divergent", "specific"))]
or_df[, Context := factor(Context, levels = c("Promoter", "Exonic", "Intronic", "Distal"))]
#
# 额外加一列，标记富集还是缺失
or_df[, Direction := ifelse(OR > 1, "Enriched", "Depleted")]
#
# 查看结果
print(or_df)
#
library(ggplot2)
#
ggplot(or_df, aes(x = OR, y = Context, color = Direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, size = 0.6) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("Enriched" = "#E64B35", "Depleted" = "#4DBBD5")) +
  facet_wrap(~ Category, ncol = 1) +
  theme_classic(base_size = 12) +
  labs(x = "Odds ratio", y = NULL, color = NULL)
#
#
#
#
######################################################################################################################################################################################################
###################################################################################################
#
load("./data_save_For_plot2/HE_CA_to_Mm_ratio.RData")
load("./data_save_For_plot2/HE_CA_to_TS_ratio.RData")
load("./data_save_For_plot2/HE_TE_ratio_merge.RData")
load("./data_save_For_plot2/TS_CA_to_Hg38_ratio.RData")
load("./data_save_For_plot2/TS_CA_to_Mm_ratio.RData")
load("./data_save_For_plot2/TS_TE_ratio_merge.RData")
#
HE_ID_rank <- read.csv("/home/sunkeyong/Human_em/sub_grep_analysi/Data_analysis/GWAS_by_snap_Peaks_q0.01/subtype.color.txt",sep = "\t")
TS_ID_rank <- read.csv("/home/sunkeyong/Multispecies_analysis/M_TreeShrew/Analysis_T0808/pseudo_bulk_analysis/bed_dedu/ArchR_showing/GWAS_celltype_rank.csv",sep = "\t")
#
plot(HE_TE_ratio_merge$TE_all_ratio,HE_CA_to_TS_ratio[3,],col=HE_ID_rank$Major_ID_color,type="p",pch = 16)
plot(HE_TE_ratio_merge$TE_all_ratio,HE_CA_to_Mm_ratio[3,],col=HE_ID_rank$Major_ID_color,type="p",pch = 16)
#
plot(TS_TE_ratio_merge$TE_all_ratio,TS_CA_to_Hg38_ratio[3,],col=TS_ID_rank$rank_color,type="p",pch = 16)
plot(TS_TE_ratio_merge$TE_all_ratio,TS_CA_to_Mm_ratio[3,],col=TS_ID_rank$rank_color,type="p",pch = 16)
#
cor.test(HE_TE_ratio_merge$TE_all_ratio,HE_CA_to_TS_ratio[3,],col=HE_ID_rank$Major_ID_color,method = c("spearman")) # 0.9741205 
cor.test(HE_TE_ratio_merge$TE_all_ratio,HE_CA_to_Mm_ratio[3,],col=HE_ID_rank$Major_ID_color,method = c("spearman")) # 0.9787754
#
cor.test(TS_TE_ratio_merge$TE_all_ratio,TS_CA_to_Hg38_ratio[3,],col=TS_ID_rank$rank_color,method = c("spearman")) # 0.7674314
cor.test(TS_TE_ratio_merge$TE_all_ratio,TS_CA_to_Mm_ratio[3,],col=TS_ID_rank$rank_color,method = c("spearman")) # 0.7770845
#
#
#
#
head(HE_TE_ratio_merge)
#
boxplot(HE_TE_ratio_merge[,1],HE_TE_ratio_merge[,2],HE_TE_ratio_merge[,3],HE_TE_ratio_merge[,4],HE_TE_ratio_merge[,5],HE_TE_ratio_merge[,6],HE_TE_ratio_merge[,7],
        TS_TE_ratio_merge[,1],TS_TE_ratio_merge[,2],TS_TE_ratio_merge[,3],TS_TE_ratio_merge[,4],TS_TE_ratio_merge[,5],TS_TE_ratio_merge[,6],TS_TE_ratio_merge[,7],
        outline=F)
#
boxplot(HE_TE_ratio_merge[,7],
        HE_TE_ratio_merge[,1],HE_TE_ratio_merge[,2],HE_TE_ratio_merge[,3],HE_TE_ratio_merge[,4],HE_TE_ratio_merge[,5],HE_TE_ratio_merge[,6],
        outline=F,
        col=c("gray","#76EEB0","#76EEB0","#76EEB0","#5AB000","#5AB000","#5AB000"))
#
boxplot(TS_TE_ratio_merge[,7],
        TS_TE_ratio_merge[,1],TS_TE_ratio_merge[,2],TS_TE_ratio_merge[,3],TS_TE_ratio_merge[,4],TS_TE_ratio_merge[,5],TS_TE_ratio_merge[,6],
        outline=F,
        col=c("gray","#F0DFF9","#F0DFF9","#F0DFF9","#FF99C9","#FF99C9","#FF99C9"))
#
median(HE_TE_ratio_merge[,7])
boxplot(HE_TE_ratio_merge[,7])
#        
#


#



#


#



#


#


#
