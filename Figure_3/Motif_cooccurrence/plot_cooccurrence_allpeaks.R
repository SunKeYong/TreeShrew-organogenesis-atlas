setwd("~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict/C3.2/finemo/counts")
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(tidyr)
#-----------------
library(rvest)
modisco_C3.2_path <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict"
modiscos <- c("modisco_C3.2_path")
celltypes <- c("C3.2")
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
modisco_report$"C3.2"
class(modisco_report$"C3.2")
#
modisco_report_df_all <- do.call("rbind", modisco_report)
#
unique(modisco_report_df_all$celltype)
#---建立 pattern → tf0 对照表-------------
motif_annotation <- modisco_report[["C3.2"]] %>%
  transmute(
    name = pattern,
    motif_label = tf0
  )
head(motif_annotation)
#检查共现矩阵的名称能否匹配 pattern
head(rownames(co_mat))
head(motif_annotation$name)
setdiff(
  rownames(co_mat),
  motif_annotation$name
)
setdiff(
  motif_annotation$name,
  rownames(co_mat)
)
#--------------------------------------------------
# 1. 假设 co_mat 是 motif 共现矩阵
#--------------------------------------------------
co_mat <- read.csv("~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C3Cardiomyocytes/output/predict/C3.2/finemo/counts/TS_C3.2_Allpeaks_motif_cooccurrence_matrix.csv",
  row.names=1,
  check.names=FALSE)
co_mat <- as.matrix(co_mat)
head(co_mat)
str(co_mat)
mode(co_mat) #numeric
#
nrow(co_mat) == ncol(co_mat) #方阵
identical(rownames(co_mat), colnames(co_mat)) #行名和列名是否完全一致
isSymmetric(co_mat) #检查矩阵是否对称:motif A 和 motif B 的共现次数，应该等于 motif B 和 motif A
#对角线是每个 motif 自身出现的 peak 数量,保存后清零
motif_frequency <- diag(co_mat)
diag(co_mat) <- 0

#--------------------------------------------------
# 2. 从共现矩阵提取上三角，转成长表
#--------------------------------------------------
index <- upper.tri(co_mat, diag = FALSE) #得到一个逻辑矩阵：True-属于矩阵上三角，需要保留,FALSE-不保留，diag = FALSE：不包含对角线

edges <- data.frame(
  motif1 = rownames(co_mat)[row(co_mat)[index]],
  motif2 = colnames(co_mat)[col(co_mat)[index]],
  weight = co_mat[index],
  stringsAsFactors = FALSE
)

head(edges)
dim(edges)
#--删除没有共现的 motif pair-----------------------------
edges <- edges %>%
  filter(
    !is.na(weight),
    weight > 0
  )

dim(edges)
summary(edges$weight)
#--------------------------------------------------
# 3. 过滤太弱的边：数值或者排序
#--------------------------------------------------
# Method 1
#edges_plot <- edges %>%
 # filter(weight >= 50)
#head(edges_plot)
#dim(edges_plot)
# Method 2
edges_plot <- edges %>%
  slice_max(
     order_by = weight,
     n = 20,
     with_ties = FALSE
   )
# Method 3
#weight_cutoff <- quantile(
 # edges$weight,
  #0.70,
  #na.rm = TRUE
#)
#weight_cutoff
#edges_plot <- edges %>%
 # filter(weight >= weight_cutoff)

dim(edges_plot)
head(edges_plot)
#--------------------------------------------------
# 4. 构建节点信息
#    节点大小这里用“与其他motif总共现权重”
#     添加 tf0 标签（motif_label）
#--------------------------------------------------
nodes <- data.frame(
  name = rownames(co_mat),
  # motif出现在多少个peak中
  motif_frequency = motif_frequency,
  # 该motif与其他motif的总共现次数
  total_cooccurrence = rowSums(co_mat),
  stringsAsFactors = FALSE
)%>%
  left_join(
    motif_annotation,
    by = "name"
  ) %>%
  mutate(
    # tf0 缺失时，仍然显示原来的 pattern
    motif_label = ifelse(
      is.na(motif_label) | motif_label == "",
      name,
      motif_label
    )
  )
head(nodes)
# 保留出现在 edges_plot 中的 motif--------------------------------------------------
motifs_in_edges <- unique(
  c(
    edges_plot$motif1,
    edges_plot$motif2
  )
)

nodes_plot <- nodes %>%
  filter(name %in% motifs_in_edges)
#--------------------------------------------------
# 5. 构建网络对象 igraph 对象
# 把边表和节点表转换为 igraph 网络对象，会做三件事：
#（1）创建节点。vertices = nodes 中第一列必须是节点名称。nodes 中其他列会变成节点属性，比如motif_frequency
# (2)创建连线。d = edges_plot 的前两列会被解释为 from,to。directed = FALSE是无向网络
# (3)保存连线属性。边表第三列 weight成为边的属性
#--------------------------------------------------
g <- graph_from_data_frame(
  d = edges_plot[, c("motif1", "motif2", "weight")],
  vertices = nodes_plot,
  directed = FALSE
)
g
#--------------------------------------------------
# 6. 作图
#整体逻辑是：
  #g：网络对象；
  #圆圈：motif；
  #连线：两个 motif 在同一个 peak 中共同出现；
  #线宽/线颜色深浅：motif pair 的共现频率；
  #圆圈大小：单个 motif 的出现频率。
#--------------------------------------------------

set.seed(123)

p <- ggraph(
  g,
  layout = "fr"
) +

  # 线的颜色深浅表示共现频率
  geom_edge_link(
    aes(colour = weight),
    width = 1,
    alpha = 0.85
  ) +

  # 浅灰色表示低频，黑色表示高频
  scale_edge_colour_gradient(
    low = "grey90",
    high = "grey10",
    name = "Co-occurrence frequency",
    guide = guide_edge_colourbar(
      title.position = "top",
      barwidth = 1,
      barheight = 5
    )
  ) +

  # 圆圈大小表示motif出现频率
  geom_node_point(
    aes(size = motif_frequency),
    shape = 21,
    fill = "#F4A582",
    colour = "black",
    stroke = 0.8
  ) +

  scale_size_continuous(
    range = c(6, 16),
    name = "Motif frequency"
  ) +

  geom_node_text(
    aes(label = motif_label),
    repel = TRUE,
    size = 4
  ) +

  theme_void() +

  theme(
    legend.position = "right"
  )

p

ggsave(
  "motif_cooccurrence_network.pdf",
  plot = p,
  width = 8,
  height = 7
)
