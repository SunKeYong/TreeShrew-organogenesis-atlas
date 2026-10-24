library(data.table)
library(ggplot2)
#
#ratio_dt <- as.data.table(CA_to_Mm_ratio, keep.rownames = "category")
ratio_dt <- as.data.table(CA_to_TS_ratio, keep.rownames = "category")
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
  levels = c("CA_conserved", "CA_divergent", "Human_specific")
)]
#
ratio_long[, subclass := factor(subclass, levels = colnames(CA_to_TS_ratio))]
#
my_cols <- c(
  "CA_conserved" = "#EAA0A0",
  "CA_divergent" = "#8E5030",
  "Human_specific" = "#F6CAAC"
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
#Human_to_TS_ratio <- ratio_long
#Human_to_Mm_ratio <- ratio_long
#
Human_to_TS_ratio <- as.data.frame(Human_to_TS_ratio)
Human_to_TS_ratio$class <- unlist(strsplit(as.character(Human_to_TS_ratio$subclass),".",fixed=T))[seq(1,2*dim(Human_to_TS_ratio)[1],2)]
Human_to_Mm_ratio <- as.data.frame(Human_to_Mm_ratio)
Human_to_Mm_ratio$class <- unlist(strsplit(as.character(Human_to_Mm_ratio$subclass),".",fixed=T))[seq(1,2*dim(Human_to_Mm_ratio)[1],2)]
##
wilcox.test(Human_to_TS_ratio[which(Human_to_TS_ratio$class=="C1" & Human_to_TS_ratio$category=="Human_specific"),]$ratio,
            Human_to_TS_ratio[which(Human_to_TS_ratio$class=="C10" & Human_to_TS_ratio$category=="Human_specific"),]$ratio)
#
wilcox.test(Human_to_Mm_ratio[which(Human_to_Mm_ratio$class=="C1" & Human_to_Mm_ratio$category=="Human_specific"),]$ratio,
            Human_to_Mm_ratio[which(Human_to_Mm_ratio$class=="C10" & Human_to_Mm_ratio$category=="Human_specific"),]$ratio)
#
plot(Human_to_TS_ratio[which(Human_to_TS_ratio$class=="C1" & Human_to_TS_ratio$category=="Human_specific"),]$ratio,
     Human_to_Mm_ratio[which(Human_to_Mm_ratio$class=="C1" & Human_to_Mm_ratio$category=="Human_specific"),]$ratio)
plot(Human_to_TS_ratio[which(Human_to_TS_ratio$class=="C10" & Human_to_TS_ratio$category=="Human_specific"),]$ratio,
     Human_to_Mm_ratio[which(Human_to_Mm_ratio$class=="C10" & Human_to_Mm_ratio$category=="Human_specific"),]$ratio)
#
#
#save(Human_to_TS_ratio,file = "./data_save_For_plot/Human_to_TS_ratio.RData")
#save(Human_to_Mm_ratio,file = "./data_save_For_plot/Human_to_Mm_ratio.RData")
#
plot(Human_to_TS_ratio$ratio,Human_to_Mm_ratio$ratio)
## 小鼠的 p-value = 2.068e-10 # # 树鼩的 p-value = 2.068e-10
#
#
#
library(ggplot2)
#
HE_11color <- c("#FDCDA1","#00BFF1","#9970A1","#9999F9",## 1-4
                "#F4CAE1","#C9DBEC","#66BD61","#F6E8C1",## 5-8
                "#FF33C1","#B3EE31","#FF7F01") ## 9-11
#
names(HE_11color) <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11")
#
df_hs <- subset(Human_to_Mm_ratio, category == "Human_specific")
head(df_hs)
#
df_hs$class <- factor(
  df_hs$class,
  levels = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11")
)
#
ggplot(df_hs, aes(x = class, y = ratio, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = HE_11color) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "Human-specific ratio") +
  theme(legend.position = "none")
#
#
#
#
###
####### ####### #######. #
# TE 含量 弄 human to TS
# Full peak set
table(HOPA_lift$Human_to_TS_State)
#
table(human_keep_Final_TE3$repClass)/sum(table(human_keep_Final_TE3$repClass))
#
targetX <- HOPA_lift[which(HOPA_lift$Human_to_TS_State=="Human_specific"),]
targetTE1 <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE1$repClass)/sum(table(targetTE1$repClass))
#.    DNA       LINE        LTR       None       SINE 
# 0.03410886 0.16304107 0.19852702 0.53805514 0.06626791  
targetX <- HOPA_lift[which(HOPA_lift$Human_to_TS_State=="CA_conserved"),]
targetTE2 <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE2$repClass)/sum(table(targetTE2$repClass))
#.     DNA       LINE        LTR       None       SINE 
#  0.006739783 0.038396752 0.029828737 0.917459516 0.007575213 
targetX <- HOPA_lift[which(HOPA_lift$Human_to_TS_State=="CA_divergent"),]
targetTE3 <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE3$repClass)/sum(table(targetTE3$repClass))
#.     DNA       LINE        LTR       None       SINE 
# 0.01061612 0.07626639 0.06233502 0.82845037 0.02233209 
#
HE_to_TS_TE_ratio <- cbind(table(human_keep_Final_TE3$repClass)/sum(table(human_keep_Final_TE3$repClass)),
                           table(targetTE2$repClass)/sum(table(targetTE2$repClass)),
                           table(targetTE3$repClass)/sum(table(targetTE3$repClass)),
                           table(targetTE1$repClass)/sum(table(targetTE1$repClass)))
HE_to_TS_TE_ratio <- as.data.frame(HE_to_TS_TE_ratio)
colnames(HE_to_TS_TE_ratio) <- c("Full set","CA_conserved","CA_divergent","Human_specific_to_TS")
HE_to_TS_TE_ratio
HE_to_TS_TE_ratio2 <- HE_to_TS_TE_ratio[-4,]
HE_to_TS_TE_ratio2
#
#       Full set  CA_conserved  CA_divergent     Human_specific_to_TS
# DNA  0.01583665  0.006739783   0.01061612           0.03410886
# LINE 0.08362323  0.038396752   0.07626639           0.16304107
# LTR  0.08751820  0.029828737   0.06233502           0.19852702
# SINE 0.02829385  0.007575213   0.02233209           0.06626791
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
  levels = c("Full set", "CA_conserved", "CA_divergent", "Human_specific_to_TS")
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
Human_to_TS_TE <- df_long
#
#
#########
# 弄小鼠的
####### ####### #######. #
# Full peak set
table(HOPA_lift$Human_to_Mm_State)
#
table(human_keep_Final_TE3$repClass)/sum(table(human_keep_Final_TE3$repClass))
#
targetX <- HOPA_lift[which(HOPA_lift$Human_to_Mm_State=="Human_specific"),]
targetTE1 <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE1$repClass)/sum(table(targetTE1$repClass))
#.    DNA       LINE        LTR       None       SINE 
# 0.02794986 0.14283013 0.16375082 0.61035913 0.05511005 
targetX <- HOPA_lift[which(HOPA_lift$Human_to_Mm_State=="CA_conserved"),]
targetTE2 <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE2$repClass)/sum(table(targetTE2$repClass))
#.     DNA       LINE        LTR       None       SINE 
#  0.004254350 0.023503690 0.015445932 0.952865082 0.003930947
targetX <- HOPA_lift[which(HOPA_lift$Human_to_Mm_State=="CA_divergent"),]
targetTE3 <- subset(human_keep_Final_TE3, Human_Peak_rank %in% targetX$Human_Peak_rank)
table(targetTE3$repClass)/sum(table(targetTE3$repClass))
#.     DNA       LINE        LTR       None       SINE 
# 0.007479595 0.044992731 0.034406714 0.904135658 0.008985303  
#
HE_to_Mm_TE_ratio <- cbind(table(human_keep_Final_TE3$repClass)/sum(table(human_keep_Final_TE3$repClass)),
                           table(targetTE2$repClass)/sum(table(targetTE2$repClass)),
                           table(targetTE3$repClass)/sum(table(targetTE3$repClass)),
                           table(targetTE1$repClass)/sum(table(targetTE1$repClass)))
HE_to_Mm_TE_ratio <- as.data.frame(HE_to_Mm_TE_ratio)
colnames(HE_to_Mm_TE_ratio) <- c("Full set","CA_conserved","CA_divergent","Human_specific_to_TS")
HE_to_Mm_TE_ratio
HE_to_Mm_TE_ratio2 <- HE_to_Mm_TE_ratio[-4,]
HE_to_Mm_TE_ratio2
#
#       Full set CA_conserved CA_divergent Human_specific_to_TS
# DNA  0.01583665  0.004254350  0.007479595           0.02794986
# LINE 0.08362323  0.023503690  0.044992731           0.14283013
# LTR  0.08751820  0.015445932  0.034406714           0.16375082
# SINE 0.02829385  0.003930947  0.008985303           0.05511005
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
  levels = c("Full set", "CA_conserved", "CA_divergent", "Human_specific_to_TS")
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
Human_to_Mm_TE <- df_long
#
#
Human_to_Mm_TE
Human_to_TS_TE
#
save(Human_to_TS_TE,file = "./data_save_For_plot/Human_to_TS_TE.RData")
save(Human_to_Mm_TE,file = "./data_save_For_plot/Human_to_Mm_TE.RData")
#
#
############ ########
#
head(TE_CT_ratio2)
#
ID_rank <- read.csv("/home/sunkeyong/Human_em/sub_grep_analysi/Data_analysis/GWAS_by_snap_Peaks_q0.01/subtype.color.txt",sep = "\t")
head(ID_rank)
#
table(ID_rank$Subtype_ID_Num==TE_CT_ratio2$celltype)
#
TE_CT_ratio2$color <- ID_rank$Major_ID_color
head(TE_CT_ratio2)
#
barplot(TE_CT_ratio2$V1,col=TE_CT_ratio2$color)
#
TE_CT_ratio2$class <- unlist(strsplit(TE_CT_ratio2$celltype,".",fixed=T))[seq(1,2*131,2)]
head(TE_CT_ratio2)
#
library(ggplot2)
#
TE_CT_ratio2$class <- factor(
  TE_CT_ratio2$class,
  levels = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11")
)
#
ggplot(TE_CT_ratio2, aes(x = class, y = V1, fill = class)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8, color = "black") +
  scale_fill_manual(values = HE_11color) +
  theme_classic(base_size = 12) +
  labs(x = "Class", y = "TE ratio") +
  theme(legend.position = "none")
#
head(TE_CT_ratio2)
C1 <- TE_CT_ratio2[which(TE_CT_ratio2$class=="C1"),]$V1
C10 <- TE_CT_ratio2[which(TE_CT_ratio2$class=="C10"),]$V1
#
wilcox.test(C1,C10)
wilcox.test(C10,C1)
#
#
TE_CT_ratio2$V1
#
ratio_long2X <- ratio_long2[which(ratio_long2$category=="Human_specific"),]
#
plot(TE_CT_ratio2$V1,ratio_long2X$ratio)
#
#
human_TE_ratio <- TE_CT_ratio2
save(human_TE_ratio,file = "human_TE_ratio.RData")
#
