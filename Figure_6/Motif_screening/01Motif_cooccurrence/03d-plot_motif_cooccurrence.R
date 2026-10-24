
#C3.1
TFTF3 <- read.csv("~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict/C3.1/finemo/heart_tissue_VISTA_HsToTS_05.filter/counts/TS_C3.1_VISTA_HsToTS_05.filter_motif_presence_matrix.csv",row.names = 1)
#C3.2
TFTF3 <- read.csv("~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict/C3.2/finemo/heart_tissue_VISTA_HsToTS_05.filter/counts/TS_C3.2_VISTA_HsToTS_05.filter_motif_presence_matrix.csv",row.names = 1)
#C6.2
TFTF3 <- read.csv("~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/chrombpnet_model_predicting/finemo/heart_tissue_VISTA/counts/Hs_C6.2_heart_tissue_VISTA_motif_presence_matrix.csv",row.names = 1)

TFTF3
pheatmap::pheatmap(TFTF3)
pheatmap::pheatmap(cor(TFTF3),display_numbers = TRUE)
pheatmap::pheatmap(cor(TFTF3[,TFTF_name]),display_numbers = TRUE)


C3.1_peak <- rownames(TFTF3)
C3.2_peak <- rownames(TFTF3)
C6.2_peak <- rownames(TFTF3)
length(intersect(C3.1_peak,C3.2_peak))
co_peak <- intersect(intersect(C3.1_peak,C3.2_peak),C6.2_peak)
length(co_peak)

#1计算每个motif出现的频率
motif_totals <- colSums(TFTF3)
top6_sums <- sort(motif_totals, decreasing = TRUE)[1:6]
top6_sums
C3.1_sums <- sort(motif_totals, decreasing = TRUE)[TFTF_name]
C3.1_sums_df <- as.data.frame(C3.1_sums)
rownames(C3.1_sums_df) <- c("GATA4","MEF2D","TEAD","NFI","KLF12","TBX20")
C3.1_sums_df
#
#1计算每个motif出现的频率
motif_totals <- colSums(TFTF3)
top6_sums <- sort(motif_totals, decreasing = TRUE)[1:6]
top6_sums
C3.2_sums <- sort(motif_totals, decreasing = TRUE)[TFTF_name]
C3.2_sums_df <- as.data.frame(C3.2_sums)
rownames(C3.2_sums_df) <- c("GATA4","MEF2D","TEAD","NFI","KLF12","TBX20")
C3.2_sums_df
#
#1计算每个motif出现的频率
motif_totals <- colSums(TFTF3)
top6_sums <- sort(motif_totals, decreasing = TRUE)[1:6]
top6_sums
C6.2_sums <- sort(motif_totals, decreasing = TRUE)[TFTF_name]
C6.2_sums_df <- as.data.frame(C6.2_sums)
rownames(C6.2_sums_df) <- c("GATA4","KLF12","MEF2D","TEAD","TBX20","NFI")
C6.2_sums_df
#
C3.1_sums_df <- C3.1_sums_df[c("GATA4","KLF12","MEF2D","TEAD","TBX20","NFI"), ,drop=FALSE]
C3.2_sums_df <- C3.2_sums_df[c("GATA4","KLF12","MEF2D","TEAD","TBX20","NFI"), ,drop=FALSE]
merge_df <- cbind(C6.2_sums_df,C3.1_sums_df, C3.2_sums_df)
#
pdf(file.path(file_dir,"HsC6.2_TSC3.1_TSC3.2_top6motif_frequency_plot_hit.pdf"),width = 8, height = 8)
par(mar = c(10, 4, 4, 2)) #c(bottom, left, top, right)
mat <- as.matrix(t(merge_df))  # 转置motif为列，细胞类型为行
# 绘制并排柱状图
barplot(mat,
        beside = TRUE,                 # 并排显示
        col = c("#0977a3","#027342c3", "#e99b0a"),  # 设置三种颜色
        legend.text = rownames(mat),   # 添加图例（细胞类型名）
        args.legend = list(x = "topright", bty = "n"), # 图例位置
        main = "Motif Frequency Comparison",
        xlab = "Motif", ylab = "Frequency",
        ylim = c(0, max(merge_df) * 1.1) # 留出顶部空间
)
dev.off()
#2
# 计算每个peak的motif总数
peak_totals <- rowSums(TFTF3[co_peak, TFTF_name])
# 按总数排序，取前20个
#top20_names <- names(sort(peak_totals, decreasing = TRUE))[1:min(20, length(peak_totals))]
sorted_names <- names(sort(peak_totals, decreasing = TRUE))
# 提取对应的矩阵子集
#top20_peaks_mat <- as.matrix(TFTF3[top20_names, TFTF_name, drop = FALSE])
sorted_mat <- as.matrix(TFTF3[sorted_names,TFTF_name,drop = FALSE])
pheatmap::pheatmap(sorted_mat,filename=file.path(file_dir,"Hs_C6.2_VISTA_HsToTS_05.filter_motif_presence.sorted.pdf"))
#
sorted_mat[sorted_mat > 1] <- 1
pheatmap::pheatmap(sorted_mat,filename=file.path(file_dir,"Hs_C6.2_VISTA_HsToTS_05.filter_motif_presence.binarized.sorted.pdf"))
#
colors <- c("#2E86C1", "#A93226", "#28B463", "#D4AC0D",
                    "#6C3483", "#138D75", "#CA6F1E", "#7D3C98")
barplot(t(as.matrix(top20_peaks_mat)), #转置
        names.arg = top20_names,
        cex.names = 0.7,
        las = 2,
        col = colors,
        xlab = "Peak ID",
        ylab = "Motif数量",
        main = "Motif数量最多的前20个Peaks",
        legend.text = colnames(top20_peaks_mat),
        args.legend = list(x = "topright", bty = "n", cex = 0.8),
        ylim = c(0, max(rowSums(top20_peaks_mat)) * 1.1))  # 扩展y轴范围

top20_names_C6.2 <- top20_names
top20_names_C3.1 <- top20_names
top20_names_C3.2 <- top20_names
intersect(top20_names_C3.1,top20_names_C3.2)
intersect(intersect(top20_names_C3.1,top20_names_C3.2),top20_names_C6.2)


#C3.1
TF_cooccurrence <- read.csv("~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict/C3.1/finemo/heart_tissue_VISTA_HsToTS_05.filter/counts/TS_C3.1_VISTA_HsToTS_05.filter_motif_cooccurrence_matrix.csv",row.names = 1)
file_dir <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict/C3.1/vista"
#C3.2
TF_cooccurrence <- read.csv("~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict/C3.2/finemo/heart_tissue_VISTA_HsToTS_05.filter/counts/TS_C3.2_VISTA_HsToTS_05.filter_motif_cooccurrence_matrix.csv",row.names = 1)
file_dir <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict/C3.2/vista"
#C6.2
TF_cooccurrence <- read.csv("~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/chrombpnet_model_predicting/finemo/heart_tissue_VISTA/counts/Hs_C6.2_heart_tissue_VISTA_motif_cooccurrence_matrix.csv",row.names = 1)
file_dir <- "~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/chrombpnet_model_predicting/vista"
#

#统一TFname排序
#C3.1
TFTF_name <- c('pos_patterns.pattern_3',#GATA4
               'pos_patterns.pattern_6',#KLF12
               'pos_patterns.pattern_2',#MEF2D
               'pos_patterns.pattern_5',#TEAD1
               'pos_patterns.pattern_4',#TBX20
               'pos_patterns.pattern_1')#NFIA

#C3.2
TFTF_name <- c('pos_patterns.pattern_5',#GATA4
               'pos_patterns.pattern_3',#KLF12
               'pos_patterns.pattern_2',#MEF2D
               'pos_patterns.pattern_4',#TEAD1
               'pos_patterns.pattern_6',#TBX20
               'pos_patterns.pattern_1')#NFIC

#C6.2
TFTF_name <- c('pos_patterns.pattern_2',#GATA4
               'pos_patterns.pattern_5',#KLF12
               'pos_patterns.pattern_1',#MEF2D
               'pos_patterns.pattern_6',#TEAD3
               'pos_patterns.pattern_3',#TBX20
               'pos_patterns.pattern_4')#NFIC

# 将对角线设为0
diag(TF_cooccurrence) <- 0
TF_cooccurrence
max(TF_cooccurrence)
pheatmap::pheatmap(TF_cooccurrence)
pheatmap::pheatmap(TF_cooccurrence[TFTF_name,TFTF_name], cluster_rows = TRUE, cluster_cols = TRUE,breaks=seq(0, 93, length.out = 101),colorRampPalette(c("#F7FBFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594"))(100))
pheatmap::pheatmap(TF_cooccurrence[TFTF_name,TFTF_name], cluster_rows = TRUE, cluster_cols = TRUE,breaks=seq(0, 93, length.out = 101),colorRampPalette(c("#F7FBFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594"))(100),filename=file.path(file_dir,"Hs_C6.2_heart_tissue_VISTA_motif_cooccurrence_matrix.pdf"))
#
pheatmap::pheatmap(TF_cooccurrence[TFTF_name,TFTF_name], cluster_rows = FALSE, cluster_cols = FALSE,breaks=seq(0, 93, length.out = 101),colorRampPalette(c("#F7FBFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594"))(100))
pheatmap::pheatmap(TF_cooccurrence[TFTF_name,TFTF_name], cluster_rows = FALSE, cluster_cols = FALSE,breaks=seq(0, 93, length.out = 101),colorRampPalette(c("#F7FBFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594"))(100),filename=file.path(file_dir,"Hs_C6.2_heart_tissue_VISTA_motif_cooccurrence_matrix_nocluster.pdf"))

