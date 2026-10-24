# plot observed and predicted accessibility across cell type models
# Figure 2b
library(Seurat)
library(Signac)
#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
library(GenomicRanges)
#library(EnsDb.Mmusculus.v79)
#library(BSgenome.Mmusculus.UCSC.mm10)
#library(DropletUtils)
#library(JASPAR2020)
#library(TFBSTools)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(gridExtra)

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
#library(DoubletFinder)
#library(scDblFinder)

#library(TCseq)
#library(gprofiler2)
#library(UpSetR)

library(Matrix)
library(patchwork)
library(parallel)

library(rvest)
library(jsonlite)



regions <- c("chr1-86640024-86641024",
"chr1-92688664-92689664",
"chr1-133532669-133533669",
"chr4-75350950-75351950",
"chr4-111028549-111029549",
"chr5-30668489-30669489",
"chr7-36437497-36438497",
"chr7-93907809-93908809",
"chr8-38991250-38992250",
"chr9-63918377-63919377",
"chr11-25125633-25126633",
"chr12-79875081-79876081",
"chr13-43530745-43531745",
"chr14-21944419-21945419",
"chr14-47854661-47855661",
"chr14-78196359-78197359",
"chr15-71663934-71664934",
"chr16-57604635-57605635",
"chr16-68751956-68752956",
"chr17-1717534-1718534",
"chr17-6087538-6088538",
"chr17-37609327-37610327",
"chr18-57398475-57399475",
"chr19-1079732-1080732",
"chr21-7551466-7552466",
"chr22-5571631-5572631",
"chr22-22053868-22054868",
"chr22-50197736-50198736",
"chr4-64736287-64737287",
"chr5-30668503-30669503",
"chr5-30734960-30735960",
"chr7-36437902-36438902",
"chr8-38991253-38992253",
"chr8-51697262-51698262",
"chr8-53216905-53217905",
"chr8-62449563-62450563",
"chr9-63918357-63919357",
"chr10-56364696-56365696",
"chr10-58033717-58034717",
"chr12-71069241-71070241",
"chr14-77626365-77627365",
"chr14-85593270-85594270",
"chr15-34612072-34613072",
"chr15-47582755-47583755",
"chr16-57247556-57248556",
"chr16-57254798-57255798",
"chr16-57604610-57605610",
"chr16-68751962-68752962",
"chr17-1717547-1718547",
"chr17-37609343-37610343",
"chr18-36128171-36129171",
"chr18-57398482-57399482",
"chr19-1079746-1080746",
"chr19-46641645-46642645",
"chr22-50197757-50198757")
names(regions) <- seq_along(regions)
length(regions)


#all celltype在C3.1_C3.2_marker_peak的预测结果
clusters <- c("C1.1","C2.1","C2.2","C2.3","C2.4","C2.5","C2.6","C2.7","C3.1","C3.2","C4.1","C5.1","C5.2","C5.3","C5.4","C5.5","C5.6","C5.7","C5.8","C5.9","C5.10","C5.11","C5.12","C5.13","C5.14","C5.15","C5.16","C5.17",
  "C5.18","C5.19","C5.20","C5.21","C5.22","C5.23","C5.24","C5.25","C5.26","C5.27","C6.1","C6.2","C7.1","C7.2","C7.3","C7.4","C7.5","C7.6","C7.7","C7.8","C7.9","C7.10",
  "C8.1","C8.2","C8.3","C8.4","C8.5","C8.6","C8.7","C8.8","C8.9","C8.10","C8.11","C8.12","C8.13","C8.14","C8.15","C8.16","C8.17","C8.18","C8.19","C8.20","C8.21",
  "C8.22","C8.23","C8.24","C8.25","C8.26","C8.27","C8.28","C9.1","C9.2","C9.3","C9.4","C10.1","C11.1","C11.2","C11.3","C11.4","C12.1","C12.2","C12.3","C12.4","C12.5","C12.6","C12.7","C12.8","C12.9","C12.10","C12.11","C12.12","C13.1",
  "C14.1","C14.2","C14.3","C14.4","C14.5","C14.6","C14.7","C14.8","C14.9","C14.10","C14.11","C14.12","C14.13","C14.14","C14.15","C14.16","C14.17","C14.18","C14.19","C14.20","C14.21","C14.22","C14.23","C14.24",
  "C15.1","C15.2","C16.1","C16.2","C16.3","C16.4","C16.5","C16.6","C16.7","C17.1","C17.2","C17.3","C17.4","C18.1","C18.2","C18.3","C18.4","C18.5","C18.6","C18.7","C18.8","C18.9","C18.10","C18.11","C18.12","C18.13","C18.14","C18.15","C18.16","C18.17","C18.18","C18.19","C18.20",
  "C19.1","C20.1","C21.1","C22.1","C23.1","C23.2","C23.3","C23.4","C23.5","C23.6","C23.7","C23.8","C23.9","C23.10","C23.11","C23.12","C23.13",
  "C24.1","C24.2","C24.3","C24.4","C24.5","C24.6","C24.7","C24.8","C24.9","C24.10","C24.11","C24.12","C24.13","C24.14","C24.15","C24.16","C24.17","C24.18","C24.19","C24.20","C24.21","C24.22","C24.23")
length(clusters)
# 1. 定义大类与子类的对应关系列表
mapping_list <- list(
  "TS_C1Adrenal_cells" = c("C1.1"),
  "TS_C2PNS"           = c("C2.1", "C2.2", "C2.3", "C2.4", "C2.5","C2.6","C2.7"),
  "TS_C3Cardiomyocytes" = c("C3.1", "C3.2"),
  "TS_C4NAG"            = c("C4.1"),
  "TS_C5CNS_neurons"    = c("C5.1", "C5.2", "C5.3", "C5.4", "C5.5", "C5.6", "C5.7", "C5.8", "C5.9", "C5.10", "C5.11","C5.12","C5.13","C5.14","C5.15","C5.16","C5.17","C5.18","C5.19","C5.20","C5.21","C5.22","C5.23","C5.24","C5.25","C5.26","C5.27" ),
  "TS_C6Definitive_erythroid_cells" = c("C6.1", "C6.2"),
  "TS_C7Endothelial_cells" = c("C7.1", "C7.2", "C7.3", "C7.4", "C7.5", "C7.6", "C7.7", "C7.8", "C7.9", "C7.10"),
  "TS_C8Epithelial_cells" = c("C8.1", "C8.2", "C8.3", "C8.4", "C8.5", "C8.6", "C8.7", "C8.8", "C8.9", "C8.10", "C8.11", "C8.12", "C8.13", "C8.14", "C8.15", "C8.16", "C8.17", "C8.18", "C8.19", "C8.20", "C8.21", "C8.22", "C8.23", "C8.24", "C8.25", "C8.26", "C8.27", "C8.28" ),
  "TS_C9Gastrointestinal_tract" = c("C9.1", "C9.2", "C9.3", "C9.4"),
  "TS_C10Hepatocytes" = c("C10.1"),
  "TS_C11Intermediate_neuronal_progenitors" = c("C11.1", "C11.2", "C11.3", "C11.4"),
  "TS_C12Kidney_cells" = c("C12.1", "C12.2", "C12.3", "C12.4", "C12.5", "C12.6", "C12.7", "C12.8", "C12.9", "C12.10", "C12.11", "C12.12" ),
  "TS_C13Lung_cells" = c("C13.1"),
  "TS_C14Mesoderm" = c("C14.1", "C14.2", "C14.3", "C14.4", "C14.5", "C14.6", "C14.7", "C14.8", "C14.9", "C14.10", "C14.11", "C14.12", "C14.13", "C14.14", "C14.15", "C14.16", "C14.17", "C14.18", "C14.19", "C14.20", "C14.21", "C14.22", "C14.23", "C14.24" ),
  "TS_C15Muscle_cells" = c("C15.1", "C15.2"),
  "TS_C16Neural_crest_PNS_glia" = c("C16.1", "C16.2", "C16.3", "C16.4", "C16.5", "C16.6", "C16.7"),
  "TS_C17Neural_crest_PNS_neurons" = c("C17.1", "C17.2", "C17.3", "C17.4" ),
  "TS_C18Neuroectoderm_and_glia" = c("C18.1", "C18.2", "C18.3", "C18.4", "C18.5", "C18.6", "C18.7", "C18.8", "C18.9", "C18.10", "C18.11", "C18.12", "C18.13", "C18.14", "C18.15", "C18.16", "C18.17","C18.18", "C18.19", "C18.20" ),
  "TS_C19Olfactory_sensory_neurons" = c("C19.1" ),
  "TS_C20Oligodendrocytes" = c("C20.1" ),
  "TS_C21Otic_cells" = c("C21.1" ),
  "TS_C22Primitive_erythroid_cells" = c("C22.1" ),
  "TS_C23Retinal_cells" = c("C23.1", "C23.2", "C23.3", "C23.4", "C23.5", "C23.6", "C23.7", "C23.8", "C23.9", "C23.10", "C23.11", "C23.12", "C23.13" ),
  "TS_C24white_blood_cells" = c("C24.1", "C24.2", "C24.3", "C24.4", "C24.5", "C24.6", "C24.7", "C24.8", "C24.9", "C24.10", "C24.11", "C24.12", "C24.13", "C24.14", "C24.15", "C24.16", "C24.17", "C24.18", "C24.19", "C24.20", "C24.21", "C24.22", "C24.23" )
)
# 2. 将列表转换为长表格 (Data Frame)
mapping_df <- do.call(rbind, lapply(names(mapping_list), function(x) {
  data.frame(group = x, id = mapping_list[[x]])
}))
head(mapping_df)
dim(mapping_df)
# 3. 定义路径模板
template_predicted <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/%s/output/predict/%s/pre_bw/fold_0/C3.1_C3.2_marker_peak_filter_with_motif_pred_bw_chrombpnet_nobias.bw"
# 4. 一次性生成所有路径
predicted <- sprintf(template_predicted, mapping_df$group, mapping_df$id)
predicted <- as.list(predicted)
names(predicted) <- clusters
print(head(predicted, 10))
#
template_observed <- "~/Treeshrew/FileSave/Bigwig_files_C1.5W_per_celltype/%s.bw"
observed <- sprintf(template_observed, mapping_df$id)
observed <- as.list(observed)
names(observed) <- clusters
print(head(observed, 10))

#1获得每个细胞在每个区域的信号数据
#不向上下游延伸，直接算1000bppeak的信号
p_predicted_getdata <- list()
for(i in names(regions)){
  p_predicted_getdata[[i]] <- BigwigTrack_getdata(
    regions[i],
    observed, #predicted/observed
    extend.upstream = 0,
    extend.downstream = 0,
    type = "line",
    y_label = i,
    bigwig.scale = "common")
}
length(p_predicted_getdata)
head(p_predicted_getdata[[1]]$all.data)
#         position  score   bw    bin
# 1     86640024   0.12   C1.1    43
# 2     86640025   0.15   C1.1    43
# ...
# (每个 bw 有 1000 行，n 个 bw → n * 1000 行)

#2计算每个区域的平均信号值
# 初始化空数据框
avg_signal_df <- data.frame()
# 遍历每个区域的结果
for (region_name in names(p_predicted_getdata)) {
  region_data <- p_predicted_getdata[[region_name]]$all.data  # 取出该区域的所有碱基数据
  # 按细胞类型（bw）分组，计算平均信号
  library(dplyr)
  region_avg <- region_data %>%
    group_by(bw) %>%
    summarise(avg_signal = mean(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(region = region_name)  # 添加区域列
  # 合并到总数据框
  avg_signal_df <- rbind(avg_signal_df, region_avg)
}
# 查看结果结构
head(avg_signal_df)
#write.csv(avg_signal_df,"~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/data/region_for_predict/C3.1_C3.2_marker_peak_filter_with_motif_avg_signal_predicted.csv",row.names = F)
#write.csv(avg_signal_df,"~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/data/region_for_predict/C3.1_C3.2_marker_peak_filter_with_motif_avg_signal_observed.csv",row.names = F)
#avg_signal_df <- read.csv("~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/data/region_for_predict/C3.1_C3.2_marker_peak_filter_with_motif_avg_signal_predicted.csv")
avg_signal_df <- read.csv("~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/data/region_for_predict/C3.1_C3.2_marker_peak_filter_with_motif_avg_signal_observed.csv")
#3绘制气泡图
library(ggplot2)
p <- ggplot(avg_signal_df, aes(x = region, y = bw, size = avg_signal)) +
  geom_point(alpha = 0.7, color = "steelblue") +  # 可改用 color 映射其他变量
  scale_size_continuous(range = c(1, 10)) +       # 调整气泡大小范围
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 旋转X轴标签
    axis.title = element_text(size = 12)
  ) +
  labs(
    title = "Average signal per cell type across regions",
    x = "Region",
    y = "Cell type",
    size = "Avg signal"
  )

print(p)

#4针对每一个Peak，看在C3.2里面的 预测/观测 值是不是超过95%的模型
# 假设 avg_signal_df 已包含所有区域各细胞类型的平均信号
# 列名：bw, avg_signal, region
library(dplyr)

# 1. 计算每个区域跨细胞类型的 95% 分位数 和 最大值
#region_stats <- avg_signal_df %>%
 # group_by(region) %>%
  #summarise(
   # region_threshold = quantile(avg_signal, probs = 0.95, na.rm = TRUE),
    #max_signal = max(avg_signal, na.rm = TRUE),
    #.groups = "drop"
  #)

# 2.1 筛选出 C3.1 或 C3.2 满足条件的区域
#selected_regions <- avg_signal_df %>%
  #filter(bw %in% c("C3.1", "C3.2")) %>%
  #left_join(region_stats, by = "region") %>%
  #大于 95% 阈值且等于该区域的最大值
  #filter(avg_signal > region_threshold & avg_signal == max_signal) %>%
  #filter(avg_signal > region_threshold) %>%
  #pull(region) %>%
  #unique()

#2.2 将C3.1和C3.2的marker分开筛选，要求：C3.1的marker peak在C3.1的score中超过95%阈值，C3.2的marker peak在C3.2的score中超过95%阈值
#感觉不同细胞的predict score不能相互比较，有测序深度的影响，后面只关注observe score选择top的
#selected_regions_C3.1 <- avg_signal_df[which(avg_signal_df$region<=28),] %>%
 # filter(bw %in% c("C3.1")) %>%
  #left_join(region_stats, by = "region") %>%
  #大于 95% 阈值且等于该区域的最大值
  #filter(avg_signal > region_threshold & avg_signal == max_signal) %>%
  #filter(avg_signal > region_threshold) %>%
  #pull(region) %>%
  #unique()
# 查看结果
#print(selected_regions_C3.1)
  
#selected_regions_C3.2 <- avg_signal_df[which(avg_signal_df$region>28),] %>%
 # filter(bw %in% c("C3.2")) %>%
  #left_join(region_stats, by = "region") %>%
  #大于 95% 阈值且等于该区域的最大值
  #filter(avg_signal > region_threshold & avg_signal == max_signal) %>%
  #filter(avg_signal > region_threshold) %>%
  #unique()
# 查看结果
#print(selected_regions_C3.2)

#将selected_region按信号强度排序，分别保留前5个
#sorted_all <- avg_signal_df[which(avg_signal_df$region %in% selected_regions_C3.1), ] %>%
 # filter(bw == "C3.1") %>%
  #arrange(desc(avg_signal))
#selected_regions_C3.1_Top5 <- head(sorted_all$region, 5)
#selected_regions_C3.1_Top5

#Final: avg_signal_df using observed
sorted_all <- avg_signal_df[which(avg_signal_df$region>28),] %>% #C3.2 的 marker 区域
  filter(bw == "C3.2") %>%
  arrange(desc(avg_signal))
selected_regions_C3.2_Top5 <- head(sorted_all$region, 5)
selected_regions_C3.2_Top5 #49 30 55 50 46
selected_regions_C3.2_Top10 <- head(sorted_all$region, 10)
selected_regions_C3.2_Top10 #49 30 55 50 46 33 35 37 54 44
selected_regions_C3.2_Top15 <- head(sorted_all$region, 15)
selected_regions_C3.2_Top15 #49 30 55 50 46 33 35 37 54 44 53 29 31 51 47
selected_regions_C3.2_Top20 <- head(sorted_all$region, 20)
selected_regions_C3.2_Top20 #49 30 55 50 46 33 35 37 54 44 53 29 31 51 47 38 40 41 43 48

#4-2绘制气泡图
library(ggplot2)
p2 <- ggplot(avg_signal_df[avg_signal_df$region %in% c(selected_regions_C3.2_Top10), ], aes(x = factor(region), y = bw, size = avg_signal)) +
  geom_point(alpha = 0.7, color = "steelblue") +  # 可改用 color 映射其他变量
  scale_size_continuous(range = c(1, 10)) +       # 调整气泡大小范围
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 旋转X轴标签
    axis.title = element_text(size = 12)
  ) +
  labs(
    title = "Average signal per cell type across regions",
    x = "Region",
    y = "Cell type",
    size = "Avg signal"
  )

print(p2)
