#Refer to https://github.com/LlorensLab/zamboni_et_al/blob/main/code/explore/ml_interpretations.R
setwd("~/chrombpnet/chrombpnet_use_snakefile/qc_metrics/treeshrew")
library(jsonlite)
library(dplyr)
library(ggplot2)
library(rvest)

# extracts qc metrics from chrombpnet json output file
extract_metrics_to_df <- function(json_file_path) {
  json_data <- fromJSON(json_file_path)
  
  spearmanr <- json_data$counts_metrics$peaks$spearmanr
  pearsonr <- json_data$counts_metrics$peaks$pearsonr
  mse <- json_data$counts_metrics$peaks$mse
  median_jsd <- json_data$profile_metrics$peaks$median_jsd
  median_norm_jsd <- json_data$profile_metrics$peaks$median_norm_jsd
  
  # Create a data frame
  df <- data.frame(
    spearmanr = spearmanr,
    pearsonr = pearsonr,
    mse = mse,
    median_jsd = median_jsd,
    median_norm_jsd = median_norm_jsd
  )
  
  return(df)
}


# qc metrics for each validation fold

#-----C1.1-C24.3--24_representative_celltypes-------------------------------------
json_C1.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C1Adrenal_cells/output/models"
json_C2.5_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C2PNS/output/models"
json_C3.2_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/models"
json_C4.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C4NAG/output/models"
json_C5.16_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C5CNS_neurons/output/models"
json_C6.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C6Definitive_erythroid_cells/output/models"
json_C7.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C7Endothelial_cells/output/models"
json_C8.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C8Epithelial_cells/output/models"
json_C9.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C9Gastrointestinal_tract/output/models"
json_C10.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C10Hepatocytes/output/models"
json_C11.2_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C11Intermediate_neuronal_progenitors/output/models"
json_C12.4_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C12Kidney_cells/output/models"
json_C13.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C13Lung_cells/output/models"
json_C14.7_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C14Mesoderm/output/models"
json_C15.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C15Muscle_cells/output/models"
json_C16.5_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C16Neural_crest_PNS_glia/output/models"
json_C17.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C17Neural_crest_PNS_neurons/output/models"
json_C18.20_path <-"~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C18Neuroectoderm_and_glia/output/models"
json_C19.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C19Olfactory_sensory_neurons/output/models"
json_C20.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C20Oligodendrocytes/output/models"
json_C21.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C21Otic_cells/output/models"
json_C22.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C22Primitive_erythroid_cells/output/models"
json_C23.6_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C23Retinal_cells/output/models"
json_C24.3_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C24white_blood_cells/output/models"
#
jsons <- c("json_C1.1_path", "json_C2.5_path", "json_C3.2_path", "json_C4.1_path",
           "json_C5.16_path", "json_C6.1_path", "json_C7.1_path", "json_C8.1_path",
           "json_C9.1_path", "json_C10.1_path", "json_C11.2_path","json_C12.4_path",
           "json_C13.1_path", "json_C14.7_path", "json_C15.1_path", "json_C16.5_path",
           "json_C17.1_path", "json_C18.20_path", "json_C19.1_path", "json_C20.1_path",
           "json_C21.1_path", "json_C22.1_path", "json_C23.6_path", "json_C24.3_path")
celltypes <- c("C1.1","C2.5","C3.2","C4.1","C5.16","C6.1","C7.1","C8.1","C9.1","C10.1","C11.2",
  "C12.4","C13.1","C14.7","C15.1","C16.5","C17.1","C18.20","C19.1","C20.1","C21.1","C22.1","C23.6","C24.3")

names(jsons)<- celltypes
jsons
metrics <- list()
for(i in names(jsons)){
  for(j in c("fold_0", "fold_1", "fold_2", "fold_3", "fold_4")){
    cross_validation <- paste0(i, "_", j)
    json_file_path <- paste0(get(jsons[i]), "/",j, "/", i, "/evaluation/chrombpnet_metrics.json")
    metrics[[cross_validation]] <- extract_metrics_to_df(json_file_path)
    metrics[[cross_validation]]$celltype <- i
    metrics[[cross_validation]]$fold <- j
  }
}

metrics
#
performance <- do.call("rbind", metrics)
#
summary_performance <- performance %>%
  group_by(celltype) %>%
  summarise(
    mean_pearsonr = mean(pearsonr),
    std_pearsonr = sd(pearsonr),
    mean_spearmanr = mean(spearmanr),
    std_spearmanr = sd(spearmanr)
  )
write.csv(performance, file = "TS_24celltype_performance.csv")
write.csv(summary_performance, file = "TS_24celltype_summary_performance.csv")
#
performance <- read.csv("TS_24celltype_performance.csv")
summary_performance <- read.csv("TS_24celltype_summary_performance.csv")
head(performance)
head(summary_performance)
# Barplot
performance$celltype <- factor(performance$celltype,
                                        levels = c("C1.1","C2.5","C3.2","C4.1","C5.16","C6.1","C7.1","C8.1","C9.1","C10.1","C11.2","C12.4","C13.1","C14.7","C15.1","C16.5","C17.1","C18.20","C19.1","C20.1","C21.1","C22.1","C23.6","C24.3")) 
summary_performance$celltype <- factor(summary_performance$celltype,
                                        levels = c("C1.1","C2.5","C3.2","C4.1","C5.16","C6.1","C7.1","C8.1","C9.1","C10.1","C11.2","C12.4","C13.1","C14.7","C15.1","C16.5","C17.1","C18.20","C19.1","C20.1","C21.1","C22.1","C23.6","C24.3")) 
TS_24color <- c("#FFFF00","#00BFFF","#9970AB","#00FF7F","#FDCDAC",
"#3288BD","#32CD32","#48D1CC","#FF33CC","#B15928",
"#6495ED","#D6604D","#0000FF","#C2A5CF","#EE82EE",
"#FF7F00","#9E0142","#B3CDE3","#87CEFA","#00FFFF",
"#CDE8B9","#F4CAE4","#0099FF","#DE77AE")
#boxplot-pearsonr------------------------------------
p <- ggplot(performance, aes(x = celltype, y = pearsonr, fill = celltype)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    aes(color = celltype),
    width = 0.15,
    size = 1.8,
    alpha = 0.7
  ) +
  scale_fill_manual(values = TS_24color) +
  scale_color_manual(values = TS_24color) +
  coord_cartesian(ylim = c(0.55, 0.85)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
ggsave("TS_24celltype_pearson_boxplot_0.55_0.85.pdf", 
       plot = p,
       width = 10,    # 宽度（英寸）
       height = 4,    # 高度（英寸）
       units = "in")  # 单位
#boxplot-spearmanr------------------------------------
p <- ggplot(performance, aes(x = celltype, y = spearmanr, fill = celltype)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    aes(color = celltype),
    width = 0.15,
    size = 1.8,
    alpha = 0.7
  ) +
  scale_fill_manual(values = TS_24color) +
  scale_color_manual(values = TS_24color) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
ggsave("TS_24celltype_spearman_boxplot.pdf", 
       plot = p,
       width = 10,    # 宽度（英寸）
       height = 6,    # 高度（英寸）
       units = "in")  # 单位
#boxplot-median_JSD------------------------------------
p <- ggplot(performance, aes(x = celltype, y = median_jsd, fill = celltype)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    aes(color = celltype),
    width = 0.15,
    size = 1.8,
    alpha = 0.7
  ) +
  scale_fill_manual(values = TS_24color) +
  scale_color_manual(values = TS_24color) +
  coord_cartesian(ylim = c(0.35, 0.65)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
ggsave("TS_24celltype_median_jsd_boxplot_0.35_0.65.pdf", 
       plot = p,
       width = 10,    # 宽度（英寸）
       height = 4,    # 高度（英寸）
       units = "in")  # 单位

#-----TS_All_cell_types-------------------------------------------
major_celltypes <- c("TS_C1Adrenal_cells","TS_C2PNS","TS_C3Cardiomyocytes","TS_C4NAG","TS_C5CNS_neurons","TS_C6Definitive_erythroid_cells",
  "TS_C7Endothelial_cells","TS_C8Epithelial_cells","TS_C9Gastrointestinal_tract","TS_C10Hepatocytes","TS_C11Intermediate_neuronal_progenitors",
  "TS_C12Kidney_cells","TS_C13Lung_cells","TS_C14Mesoderm","TS_C15Muscle_cells","TS_C16Neural_crest_PNS_glia","TS_C17Neural_crest_PNS_neurons","TS_C18Neuroectoderm_and_glia",
  "TS_C19Olfactory_sensory_neurons","TS_C20Oligodendrocytes","TS_C21Otic_cells","TS_C22Primitive_erythroid_cells","TS_C23Retinal_cells","TS_C24white_blood_cells")
#major_celltypes <- c("Hs_C1CNS_neurons")
#main_path <- paste(path_to_cross_valid, "Hs_C1CNS_neurons", "output/models","fold_0",sep="/")
#list.dirs(main_path, full.names = FALSE, recursive = FALSE)
metrics <- list()
path_to_cross_valid <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew"
for(i in major_celltypes){
  for(j in c("fold_0","fold_1","fold_2","fold_3","fold_4")){
    main_path <- paste(path_to_cross_valid, i, "output/models",j,sep="/")
    # 获取所有子文件夹（具体的细胞文件夹）
    sub_cells <- list.dirs(main_path, full.names = FALSE, recursive = FALSE)
    for(cell in sub_cells) {
      cell_path <- paste(main_path, cell, sep="/")
      cross_validation <- paste0(cell, "_", j)
      json_file_path <- paste0(cell_path, "/evaluation/chrombpnet_metrics.json")
      metrics[[cross_validation]] <- extract_metrics_to_df(json_file_path)
      metrics[[cross_validation]]$celltype <- cell
      metrics[[cross_validation]]$fold <- j
    }
  }
}
#add metrics from other paths
major_celltypes <- c("TS_C2PNS","TS_C5CNS_neurons","TS_C6Definitive_erythroid_cells",
  "TS_C7Endothelial_cells","TS_C8Epithelial_cells","TS_C9Gastrointestinal_tract","TS_C11Intermediate_neuronal_progenitors",
  "TS_C12Kidney_cells","TS_C14Mesodem","TS_C15Muscle_cells","TS_C16Neural_crest_PNS_glia","TS_C17Neural_crest_PNS_neurons","TS_C18Neuroectoderm_and_glia",
  "TS_C23Retinal_cells","TS_C24white_blood_cells")
path_to_cross_valid <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/Download_metrics"
for(i in major_celltypes){
  for(j in c("fold_0","fold_1","fold_2","fold_3","fold_4")){
    main_path <- paste(path_to_cross_valid, i, "models",j,sep="/")
    # 获取所有子文件夹（具体的细胞文件夹）
    sub_cells <- list.dirs(main_path, full.names = FALSE, recursive = FALSE)
    for(cell in sub_cells) {
      cell_path <- paste(main_path, cell, sep="/")
      cross_validation <- paste0(cell, "_", j)
      json_file_path <- paste0(cell_path, "/evaluation/chrombpnet_metrics.json")
      metrics[[cross_validation]] <- extract_metrics_to_df(json_file_path)
      metrics[[cross_validation]]$celltype <- cell
      metrics[[cross_validation]]$fold <- j
    }
  }
}
#
metrics
length(metrics)

performance <- do.call("rbind", metrics)
performance
performance[which(performance$pearsonr < 0.6), ]#None
head(performance[order(performance$pearsonr), ])
head(performance[order(performance$spearmanr), ])
head(performance[order(performance$median_jsd), ])

summary_performance <- performance %>%
  group_by(celltype) %>%
  summarise(
    mean_pearsonr = mean(pearsonr),
    std_pearsonr = sd(pearsonr),
    mean_spearmanr = mean(spearmanr),
    std_spearmanr = sd(spearmanr),
    mean_jsd=mean(median_jsd),
    std_jsd=sd(median_jsd)
  )
head(summary_performance)
write.csv(performance, file = "~/chrombpnet/chrombpnet_use_snakefile/qc_metrics/treeshrew/TS_all_celltype_performance.csv")
write.csv(summary_performance, file = "~/chrombpnet/chrombpnet_use_snakefile/qc_metrics/treeshrew/TS_all_celltype_summary_performance.csv")
#
performance <- read.csv("TS_all_celltype_performance.csv")
summary_performance <- read.csv("TS_all_celltype_summary_performance.csv")
head(performance)
head(summary_performance)
#
#最大值和对应细胞类型
summary_performance[which.max(summary_performance$mean_pearsonr), c("celltype", "mean_pearsonr")]
#最小值和对应细胞类型
summary_performance[which.min(summary_performance$mean_pearsonr), c("celltype", "mean_pearsonr")]

library(ggplot2)
library(tidyr)
library(dplyr)
# 转成长表格式，方便 ggplot 作图
summary_performance_long <- summary_performance %>%
  select(mean_spearmanr, mean_pearsonr,mean_jsd) %>%
  pivot_longer(
    cols = everything(),
    names_to = "metric",
    values_to = "correlation"
  )
head(summary_performance_long)
p <- ggplot(summary_performance_long, aes(x = metric, y = correlation)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Correlation",
    title = "ChromBPNet prediction correlation"
  )
p
ggsave(
  filename = "~/chrombpnet/chrombpnet_use_snakefile/qc_metrics/treeshrew/TS_all_celltype_correlation_plot.pdf",
  plot = p,
  width = 8,
  height = 6
)

#---check modisco reports------------------------------
#-----C1.1-C24.3---------------------------------------
# check modisco reports
modisco_C1.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C1Adrenal_cells/output/predict"
modisco_C2.5_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C2PNS/output/predict"
modisco_C3.2_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict"
modisco_C4.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C4NAG/output/predict"
modisco_C5.16_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C5CNS_neurons/output/predict"
modisco_C6.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C6Definitive_erythroid_cells/output/predict"
modisco_C7.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C7Endothelial_cells/output/predict"
modisco_C8.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C8Epithelial_cells/output/predict"
modisco_C9.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C9Gastrointestinal_tract/output/predict"
modisco_C10.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C10Hepatocytes/output/predict"
modisco_C11.2_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C11Intermediate_neuronal_progenitors/output/predict"
modisco_C12.4_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C12Kidney_cells/output/predict"
modisco_C13.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C13Lung_cells/output/predict"
modisco_C14.7_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C14Mesoderm/output/predict"
modisco_C15.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C15Muscle_cells/output/predict"
modisco_C16.5_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C16Neural_crest_PNS_glia/output/predict"
modisco_C17.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C17Neural_crest_PNS_neurons/output/predict"
modisco_C18.20_path <-"~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C18Neuroectoderm_and_glia/output/predict"
modisco_C19.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C19Olfactory_sensory_neurons/output/predict"
modisco_C20.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C20Oligodendrocytes/output/predict"
modisco_C21.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C21Otic_cells/output/predict"
modisco_C22.1_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C22Primitive_erythroid_cells/output/predict"
modisco_C23.6_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C23Retinal_cells/output/predict"
modisco_C24.3_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C24white_blood_cells/output/predict"
modiscos <- c("modisco_C1.1_path", "modisco_C2.5_path", "modisco_C3.2_path", "modisco_C4.1_path",
           "modisco_C5.16_path", "modisco_C6.1_path", "modisco_C7.1_path", "modisco_C8.1_path",
           "modisco_C9.1_path", "modisco_C10.1_path", "modisco_C11.2_path","modisco_C12.4_path",
           "modisco_C13.1_path", "modisco_C14.7_path", "modisco_C15.1_path", "modisco_C16.5_path",
           "modisco_C17.1_path", "modisco_C18.20_path", "modisco_C19.1_path", "modisco_C20.1_path",
           "modisco_C21.1_path", "modisco_C22.1_path", "modisco_C23.6_path", "modisco_C24.3_path")
celltypes <- c("C1.1","C2.5","C3.2","C4.1","C5.16","C6.1","C7.1","C8.1","C9.1","C10.1","C11.2",
  "C12.4","C13.1","C14.7","C15.1","C16.5","C17.1","C18.20","C19.1","C20.1","C21.1","C22.1","C23.6","C24.3")
names(modiscos)<- celltypes
modiscos
#
modisco_report <- list()
for(i in names(modiscos)){
  html_file <- paste0(get(modiscos[i]), "/",i, "/tfmodisco/report/motifs.html")
  html_content <- read_html(html_file)
  modisco_report[[i]] <- html_content %>% html_nodes("table") %>% .[[1]] %>% html_table()
  modisco_report[[i]]$celltype <- i
  modisco_report[[i]]$match0 <- sub("[.]", "_", modisco_report[[i]]$match0)
  modisco_report[[i]]$match1 <- sub("[.]", "_", modisco_report[[i]]$match1)
  modisco_report[[i]]$match2 <- sub("[.]", "_", modisco_report[[i]]$match2)
  modisco_report[[i]]$tf0 <- sapply(strsplit(modisco_report[[i]]$match0, "_"), "[[", 1) %>% stringr::str_to_title()
  modisco_report[[i]]$tf1 <- sapply(strsplit(modisco_report[[i]]$match1, "_"), "[[", 1) %>% stringr::str_to_title()
  modisco_report[[i]]$tf2 <- sapply(strsplit(modisco_report[[i]]$match2, "_"), "[[", 1) %>% stringr::str_to_title()
  #modisco_report[[i]] <- replace_human_genes_with_mouse(modisco_report[[i]], human_to_mouse_dict)
}
#
modisco_report
modisco_report_df_all <- do.call("rbind", modisco_report)
unique(modisco_report_df_all$celltype)
#
# save modisco_report_df_all and manually curate to order the TF matches consistently across celltypes
write.csv(modisco_report_df_all,"modisco_report_df_all.csv")
modisco_report_df_all <- read.csv("modisco_report_df_all.csv")
#
celltype_order <- c("C1.1","C2.5","C3.2","C4.1","C5.16","C6.1","C7.1","C8.1","C9.1","C10.1","C11.2",
  "C12.4","C13.1","C14.7","C15.1","C16.5","C17.1","C18.20","C19.1","C20.1","C21.1","C22.1","C23.6","C24.3")
motifs_freq_singles <- modisco_report_df_all %>%
  group_by(celltype) %>%
  mutate(
    total_seqlets_in_celltype = sum(num_seqlets)
  ) %>%
  group_by(celltype, tf0) %>%
  summarise(
    motif_seqlets = sum(num_seqlets),
    total_seqlets_in_celltype = first(total_seqlets_in_celltype),
    freq_withpattern = motif_seqlets / total_seqlets_in_celltype
  ) %>%
  top_n(10, freq_withpattern) %>%
  ungroup()%>%
  # 按定义的顺序排序
  mutate(
    celltype = factor(celltype, levels = celltype_order)
  ) %>%
  arrange(celltype, desc(freq_withpattern))


motifs_freq_singles$round_freq <- round(motifs_freq_singles$freq_withpattern, 3)
write.csv(motifs_freq_singles,"motifs_freq_singles.csv")

motifs_freq_wide <- motifs_freq_singles[, c("celltype", "tf0", "round_freq")] %>% #转换为宽格式:celltype x tf
  tidyr::spread(key = tf0, value = round_freq, fill = 0)
celltype <- motifs_freq_wide$celltype
rownames(motifs_freq_wide) <- celltype
motifs_freq_wide <- motifs_freq_wide[,c(2:ncol(motifs_freq_wide))]
###
fun_color_range <- colorRampPalette(c("#2D558E", "grey99","#8A267C"))  # Create color generating function
hm_colors <- fun_color_range(100) 

# 转置后，tf x celltype.Figure 2d按行（每个tf在各个细胞类型之间）标准化
pheatmap::pheatmap(t(motifs_freq_wide), cluster_rows=T, cluster_cols=F,
                   scale="row", clustering_method="ward.D2",
                   fontsize=6, color = hm_colors, 
                   border_color = NA)#filename="Hs_11celltype_motif_frequency.pdf"

###新添加，想要先转换为矩阵画热图260205###
motifs_freq_wide_mat <- as.matrix(motifs_freq_wide)
rownames(motifs_freq_wide_mat)<-celltype
#scale(按列标准化，默认)tf在各个细胞类型之间标准化
motifs_freq_wide_mat_scaled <- scale(motifs_freq_wide_mat)
#转置矩阵,tf x celltype
motifs_freq_wide_mat_scaled_T <- t(motifs_freq_wide_mat_scaled)
all_colnames <- colnames(motifs_freq_wide_mat_scaled_T)
#sorted_cols <- str_sort(all_colnames, numeric = TRUE)
#mat <- mat[, sorted_cols]
motifs_freq_wide_mat_scaled_T <- motifs_freq_wide_mat_scaled_T[, all_colnames]
#
# 找出每一行（每一个 tf）最大值所在的列索引。apply(mat, 1, FUN)，1表示 对每一行执行 FUN
max_col_index <- apply(motifs_freq_wide_mat_scaled_T, 1, which.max)
# 根据最大值出现的列位置进行排序；如果位置相同，则按该位置的数值大小排序
motif_order <- order(max_col_index, -apply(motifs_freq_wide_mat_scaled_T, 1, max))
# 重新排列矩阵的行
motifs_freq_wide_mat_scaled_T <- motifs_freq_wide_mat_scaled_T[motif_order, ]
#
pheatmap::pheatmap(motifs_freq_wide_mat_scaled_T, cluster_rows=F, cluster_cols=F,
                   fontsize=6, color = hm_colors, 
                   border_color = NA,filename="TS_24celltype_top10motif_frequency.pdf")#filename="Hs_11celltype_motif_frequency.pdf"
