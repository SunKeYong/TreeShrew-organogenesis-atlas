#
setwd("/home/sunkeyong/Human_em/Cross_species/immune_hs_ts_reanalysis2/step3_peak_plot")
#
library(ArchR)
# 
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
# 画 motif
# 
load("/home/sunkeyong/Human_em/Cross_species/Immune_Human/ArchR/Human_groupMotif.RData")
load("/home/sunkeyong/Human_em/Cross_species/Immune_Human/ArchR/Human_plotVarDev.RData")
#
load("/home/sunkeyong/Human_em/Cross_species/Immune_TreeShrew/ArchR/TS_groupMotif.RData")
load("/home/sunkeyong/Human_em/Cross_species/Immune_TreeShrew/ArchR/TS_plotVarDev.RData")
#
Human_plotVarDev
TS_plotVarDev
#
colnames(Human_groupMotif)
colnames(Human_groupMotif) <- c("B_cells",
                                "Macrophages",
                                "CLP",
                                "cDCs",
                                "Erythrocyte",
                                "HSPCs",
                                "HSPCs (Mpo+)",
                                "Kupffer cells",
                                "Mast cells",
                                "Megakaryocyte",
                                "Microglia",
                                "Monocytes",
                                "NK cells",
                                "Osteoclasts",
                                "pDCs",
                                "T cells",
                                "Terminal megakaryocyte")
#
HT_keep_celltypes <- c("B_cells",
                       "Macrophages",
                       "CLP",
                       "Erythrocyte",
                       "HSPCs",
                       "HSPCs (Mpo+)",
                       "Kupffer cells",
                       "Mast cells",
                       "Megakaryocyte",
                       "Microglia",
                       "Monocytes",
                       "NK cells",
                       "T cells")
#
colnames(TS_groupMotif) <- c("Activated T cells",
                             "B cells",
                             "Macrophages",
                             "CLP",
                             "Common myeloid progenitors",
                             "cDC",
                             "Eosinophils",
                             "Erythrocyte",
                             "Granulocytes",
                             "HSPCs",
                             "Kupffer cells",
                             "Mast cells",
                             "Megakaryocyte-erythroid progenitor",
                             "Megakaryocyte",
                             "Microglia",
                             "Monocytes",
                             "HSPCs (Mpo+)",
                             "NK cells",
                             "Osteoclasts",
                             "pDC",
                             "Pre-B cells",
                             "Pre-T cells",
                             "T cells")
TS_keep_celltypes <- c("Pre-B cells",
                       "Macrophages",
                       "CLP",
                       "Erythrocyte",
                       "HSPCs",
                       "HSPCs (Mpo+)",
                       "Kupffer cells",
                       "Mast cells",
                       "Megakaryocyte",
                       "Microglia",
                       "Monocytes",
                       "NK cells",
                       "T cells")
#
TS_groupMotif
#
co_features <- rownames(TS_groupMotif)[1:869]
#
#
TS_groupMotif_keep <- TS_groupMotif[co_features,TS_keep_celltypes]
head(TS_groupMotif_keep)
#
Human_groupMotif_keep <- Human_groupMotif[co_features,HT_keep_celltypes]
head(Human_groupMotif_keep)
#
plot(Human_groupMotif_keep[1,],TS_groupMotif_keep[1,])
cor(Human_groupMotif_keep[1,],TS_groupMotif_keep[1,])
#
cor_tmp <-c()
cor_Motif <-c()
cor_pValue <-c()
#
for (i in 1:869) {
  cor_tmp <- cor.test(Human_groupMotif_keep[i,],TS_groupMotif_keep[i,],method = c("pearson"))
  cor_Motif <- c(cor_Motif,as.numeric(cor_tmp$estimate))
  cor_pValue <- c(cor_pValue,as.numeric(cor_tmp$p.value))
}
#
hist(cor_Motif)
cor_Motif <- as.data.frame(cor_Motif)
rownames(cor_Motif) <- as.character(rownames(Human_groupMotif_keep))
#
head(cor_Motif)
#
#
cor_data <- cor_Motif$cor_Motif
cor_data <- as.data.frame(cor_data)
head(cor_data)
#
ggplot(cor_data, aes(x=cor_data)) +
  geom_density(color="red", size=1.5, adjust=1.5)
#
#
#
#
################################
################################
################################
cor_x <- cbind(TS_groupMotif_keep,Human_groupMotif_keep)
#
# cor_x <- cbind(TS_groupMotif_keep[var_motif,],Human_groupMotif_keep[var_motif,])
#
dim(cor_x)
dim(TS_groupMotif_keep)
dim(Human_groupMotif_keep)
head(cor_x)
#
#
cor_x <- cor(cor_x)[1:13,14:26]
min(cor_x)
max(cor_x)
#
pheatmap(cor_x,color =  colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
#
rank_id_TS <- c("HSPCs (Mpo+)",
                "Monocytes",
                "Microglia",
                "Kupffer cells",
                "Macrophages",
                "Mast cells",
                "Erythrocyte",
                "Megakaryocyte",
                "HSPCs",
                "CLP",
                "Pre-B cells",
                "NK cells",
                "T cells")
rank_id_HS <-  c("HSPCs (Mpo+)",
                 "Monocytes",
                 "Microglia",
                 "Kupffer cells",
                 "Macrophages",
                 "Mast cells",
                 "Erythrocyte",
                 "Megakaryocyte",
                 "HSPCs",
                 "CLP",
                 "B_cells",
                 "NK cells",
                 "T cells")
#
pheatmap(cor_x[rank_id_TS,rank_id_HS],
         cluster_rows = F,
         cluster_cols = F,
         color =  colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
#
#
#######
#######
#######
#######
#######
cor_data <- cor(t(rbind(TS_groupMotif_keep,Human_groupMotif_keep)),method = c("pearson"))
cor_data <- as.data.frame(cor_data)
dim(cor_data)
cor_data <- cor_data[1:869,870:1738]
#
#
head(cor_data)
##########
boxplot(diag(as.matrix(cor_data)))
boxplot(as.vector(cor_data))
###########
cor_random <- as.data.frame(as.numeric(cor_data[lower.tri(cor_data) | upper.tri(cor_data)]))
colnames(cor_random) <- "cor_random"
dim(cor_random)
head(cor_random)
#
ggplot(cor_random, aes(x=cor_random)) +
  geom_density(color="red", size=1.5, adjust=1.5)
#
###########
cor_diag <- as.data.frame(as.numeric(diag(as.matrix(cor_data))))
head(cor_diag)
colnames(cor_diag) <- "cor_diag"
#
ggplot(cor_diag, aes(x=cor_diag)) +
  geom_density(color="red", size=1.5, adjust=1.5)
##
#
###########
#
#
cor_estimate <- cbind(rownames(Human_groupMotif_keep),cor_Motif,cor_pValue)
head(cor_estimate)
plot(cor_estimate$cor_Motif,-log10(cor_estimate$cor_pValue))
#
cor_estimate$FDR <- p.adjust(cor_estimate$cor_pValue, method = "BH")
#
plot(cor_estimate$cor_Motif,-log10(cor_estimate$FDR))
#
cor_estimate_high <- cor_estimate[which(cor_estimate$FDR < 0.01),]
dim(cor_estimate_high) # 444
max(cor_estimate_high$FDR)
min(cor_estimate_high$cor_Motif) # 0.7250703
#
#
#
#
library(ggplot2)
library(dplyr)

# 1. 读取两组数据
orth <- cor_diag$cor_diag
shuf <- cor_random$cor_random
#
# 去掉 NA
orth <- orth[!is.na(orth)]
shuf <- shuf[!is.na(shuf)]

# 2. 合成长表，允许两组长度不一样
df_long <- bind_rows(
  data.frame(Group = "Orthologous", Correlation = orth),
  data.frame(Group = "Shuffled", Correlation = shuf)
)

df_long$Group <- factor(df_long$Group, levels = c("Orthologous", "Shuffled"))

# 3. Wilcoxon rank-sum test / Mann-Whitney U test
wilcox_res <- wilcox.test(orth, shuf)
#
p_value <- wilcox_res$p.value
#
p_label <- ifelse(
  p_value < 1e-15,
  "Mann-Whitney U test: P < 1e-15",
  paste0("Mann-Whitney U test: P = ", signif(p_value, 3))
)
#
# 4. 计算 median
median_df <- df_long %>%
  group_by(Group) %>%
  summarise(
    median_cor = median(Correlation, na.rm = TRUE),
    .groups = "drop"
  )

# 5. 画 density distribution
p <- ggplot(df_long, aes(x = Correlation, fill = Group, color = Group)) +
  geom_density(alpha = 0.35, linewidth = 1.1, adjust = 1.2) +
  geom_vline(
    xintercept = 0.7250703,
    linetype = "dashed",
    linewidth = 0.8,
    color = "black")+
  scale_fill_manual(values = c(
    "Orthologous" = "#4DBBD5",
    "Shuffled" = "#E64B35"
  )) +
  scale_color_manual(values = c(
    "Orthologous" = "#0072B5",
    "Shuffled" = "#BC3C29"
  )) +
  labs(
    x = "Pearson correlation",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = p_label,
    hjust = 1.05,
    vjust = 1.5,
    size = 4.2
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.78, 0.82),
    legend.background = element_blank()
  )
#
print(p)
#
#
############## v ####### #
############## v ####### #
############## v ####### #
############## v ####### #
############## v ####### #
############## v ####### #
############## v ####### #
############## v ####### #
############## v ####### #
plot(TS_groupMotif_keep["deviations__ENSG00000102145_LINE2081_GATA1_D_N7",],
     Human_groupMotif_keep["deviations__ENSG00000102145_LINE2081_GATA1_D_N7",],cex = 1,type="p", pch = 19, col=tmp_color)
plot(TS_groupMotif_keep["deviations__ENSG00000164330_LINE251_EBF1_D_N2",],
     Human_groupMotif_keep["deviations__ENSG00000164330_LINE251_EBF1_D_N2",],cex = 1,type="p", pch = 19, col=tmp_color)
#
plot(TS_groupMotif_keep["deviations__ENSG00000196092_LINE3170_PAX5_D_N8",],
     Human_groupMotif_keep["deviations__ENSG00000196092_LINE3170_PAX5_D_N8",],cex = 1,type="p", pch = 19, col=tmp_color)
plot(TS_groupMotif_keep["deviations__ENSG00000066336_LINE1813_SPI1_D_N5",],
     Human_groupMotif_keep["deviations__ENSG00000066336_LINE1813_SPI1_D_N5",],cex = 1,type="p", pch = 19, col=tmp_color)
plot(TS_groupMotif_keep["deviations__ENSG00000092067_LINE336_CEBPE_D_N1",],
     Human_groupMotif_keep["deviations__ENSG00000092067_LINE336_CEBPE_D_N1",],cex = 1,type="p", pch = 19, col=tmp_color)
plot(TS_groupMotif_keep["deviations__ENSG00000119138_LINE843_KLF9_I",],
     Human_groupMotif_keep["deviations__ENSG00000119138_LINE843_KLF9_I",],cex = 1,type="p", pch = 19, col=tmp_color)
plot(TS_groupMotif_keep["deviations__ENSG00000139083_LINE1875_ETV6_D_N1",],
     Human_groupMotif_keep["deviations__ENSG00000139083_LINE1875_ETV6_D_N1",],cex = 1,type="p", pch = 19, col=tmp_color)
plot(TS_groupMotif_keep["deviations__ENSG00000124813_LINE3275_RUNX2_D_N2",],
     Human_groupMotif_keep["deviations__ENSG00000124813_LINE3275_RUNX2_D_N2",],cex = 1,type="p", pch = 19, col=tmp_color)
#
#
plot(TS_groupMotif_keep["deviations__ENSG00000128272_LINE397_ATF4_D",],
     Human_groupMotif_keep["deviations__ENSG00000128272_LINE397_ATF4_D",],cex = 1,type="p", pch = 19, col=tmp_color)
#
plot(TS_groupMotif_keep["deviations__ENSG00000134954_LINE1864_ETS1_D_N5",],
     Human_groupMotif_keep["deviations__ENSG00000134954_LINE1864_ETS1_D_N5",],cex = 1,type="p", pch = 19, col=tmp_color)
plot(TS_groupMotif_keep["deviations__ENSG00000151702_LINE1884_FLI1_D_N2",],
     Human_groupMotif_keep["deviations__ENSG00000151702_LINE1884_FLI1_D_N2",],cex = 1,type="p", pch = 19, col=tmp_color)
plot(TS_groupMotif_keep["deviations__ENSG00000125347_LINE2729_IRF1_D_N2",],
     Human_groupMotif_keep["deviations__ENSG00000125347_LINE2729_IRF1_D_N2",],cex = 1,type="p", pch = 19, col=tmp_color)
#    
HEi_color <- as.data.frame(cbind(colnames(Human_groupMotif_keep),
                   c("#FDCDAC","#00BFFF","#9970AB","#32CD32","#FF33CC","#B15928","#DE77AE","#B3CDE3",
                     "#3375EF","#C2A5CF","#EE82EE","#9E0142","#CDE8B9")))
colnames(HEi_color) <- c("celltypes","color")
rownames(HEi_color) <- HEi_color$celltypes
#
tmp_color <- HEi_color[colnames(Human_groupMotif_keep),]$color
#
rank_id_HS
barplot(rep(1,13),col=HEi_color[rank_id_HS,]$color, space = 0)
#
#
#
#
conserved_motif <- cor_estimate_high$`rownames(Human_groupMotif_keep)`
#
pheatmap::pheatmap(Human_groupMotif_keep[conserved_motif,],scale = "row",cluster_cols = F)
pheatmap::pheatmap(TS_groupMotif_keep[conserved_motif,],scale = "row",cluster_cols = F)
#
#
p <- pheatmap::pheatmap(TS_groupMotif_keep[conserved_motif,],scale = "row",cluster_cols = T)
sorted_row_names <- rownames(TS_groupMotif_keep[conserved_motif,])[p$tree_row$order]
colnames(TS_groupMotif_keep)[p$tree_col$order]
colnames(Human_groupMotif_keep)[p$tree_col$order]
#
#
motifheatmap_rank_id_TS <- c("Macrophages","Kupffer cells","Microglia",  "Monocytes",  
                             "Pre-B cells","Erythrocyte","Megakaryocyte",
                             "HSPCs (Mpo+)","Mast cells",
                             "CLP","HSPCs","NK cells","T cells"  )
motifheatmap_rank_id_HS <-  c("Macrophages","Kupffer cells","Microglia",  "Monocytes",  
                              "B_cells", "Erythrocyte","Megakaryocyte","HSPCs (Mpo+)","Mast cells",
                              "CLP","HSPCs","NK cells","T cells")
#
pheatmap::pheatmap(Human_groupMotif_keep[sorted_row_names,motifheatmap_rank_id_HS],
                   scale = "row",
                   show_rownames = F,
                   color = paletteContinuous("blueYellow"),
                   cluster_cols = F,cluster_rows = F)
#
pheatmap::pheatmap(TS_groupMotif_keep[sorted_row_names,motifheatmap_rank_id_TS],
                   scale = "row",
                   show_rownames = F,
                   color = paletteContinuous("blueYellow"),
                   cluster_cols = F,cluster_rows = F)
#
#