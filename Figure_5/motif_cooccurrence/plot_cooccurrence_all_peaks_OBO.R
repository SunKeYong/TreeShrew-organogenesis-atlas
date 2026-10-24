library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rvest)
library(stringr)


# ============================================================
# 0. 基本参数
# ============================================================
#species <- "human"
#base_predict_dir <- "~/chrombpnet/chrombpnet_use_snakefile/Hs_C10white_blood_cells/output/predict"
#output_dir <- "~/chrombpnet/Hs_TS_immune_chrombpnet/motif_cooccurrence"
# 需要处理的细胞类型
#celltypes <- c(
#  "C10.1","C10.2","C10.3","C10.4","C10.5","C10.6","C10.7","C10.8","C10.9","C10.10","C10.11","C10.12","C10.13","C10.14"
#)

species <- "human"
base_predict_dir <- "~/chrombpnet/chrombpnet_use_snakefile/Hs_C3Erythroid/output/predict"
output_dir <- "~/chrombpnet/Hs_TS_immune_chrombpnet/motif_cooccurrence"
# 需要处理的细胞类型
celltypes <- c(
  "C3.1","C3.3"
)

#species <- "treeshrew"
#base_predict_dir <- "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C24white_blood_cells/output/predict"
#output_dir <- "~/chrombpnet/Hs_TS_immune_chrombpnet/motif_cooccurrence"
# 需要处理的细胞类型
#celltypes <- c(
 # "C24.1","C24.2","C24.4","C24.5","C24.6","C24.7","C24.8","C24.9","C24.10","C24.11","C24.12","C24.13","C24.14","C24.15","C24.16","C24.17","C24.18","C24.19","C24.20","C24.21","C24.22","C24.23"
#)
# 每个细胞类型保留共现频率最高的多少条边
top_n_edges <- 20


# ============================================================
# 1. 读取每个细胞类型的 TF-MoDISco motif annotation
# ============================================================

modisco_report <- list()

for (celltype in celltypes) {
  html_file <- file.path(base_predict_dir,"tfmodisco",celltype,"report","motifs.html")
  message("Reading TF-MoDISco report: ", celltype)
  # 检查文件是否存在
  if (!file.exists(html_file)) {
    warning("motifs.html does not exist: ", html_file)
    next
  }
  html_content <- read_html(html_file)
  tmp <- html_content %>% html_nodes("table") %>% .[[1]] %>% html_table()
  tmp$celltype <- celltype
  # 将 match0/match1/match2 中的 "." 替换为 "_"
  tmp$match0 <- sub("[.]", "_", tmp$match0)
  tmp$match1 <- sub("[.]", "_", tmp$match1)
  tmp$match2 <- sub("[.]", "_", tmp$match2)
  # 提取 TF 名称
  tmp$tf0 <- sapply(strsplit(tmp$match0, "_"),"[[",1) %>% str_to_title()
  tmp$tf1 <- sapply(strsplit(tmp$match1, "_"),"[[",1) %>% str_to_title()
  tmp$tf2 <- sapply(strsplit(tmp$match2, "_"),"[[",1) %>% str_to_title()
  modisco_report[[celltype]] <- tmp
}


# 合并所有细胞类型的 modisco annotation
modisco_report_df_all <- bind_rows(modisco_report)
unique(modisco_report_df_all$celltype)


# ============================================================
# 2. 定义单个细胞类型的 motif 共现网络分析函数
# ============================================================

plot_motif_cooccurrence_network <- function(
    species,
    celltype,
    base_predict_dir,
    output_dir,
    modisco_report,
    top_n_edges = 20,
    output_filename = NULL
) {

  message("\n========================================")
  message("Processing: ", celltype)
  message("========================================")

  # ----------------------------------------------------------
  # 2.1 检查该细胞类型是否存在 TF-MoDISco annotation
  # ----------------------------------------------------------

  if (!celltype %in% names(modisco_report)) {
    warning("No TF-MoDISco annotation found for ",celltype,". Skip.")
    return(NULL)
  }

  # ----------------------------------------------------------
  # 2.2 建立 pattern → tf0 对照表
  # ----------------------------------------------------------

  motif_annotation <- modisco_report[[celltype]] %>%
    transmute(
      name = pattern,
      motif_label = tf0
    )

  # ----------------------------------------------------------
  # 2.3 读取 motif 共现矩阵
  # ----------------------------------------------------------

  co_mat_dir <- file.path(output_dir,species)

  co_mat_file <- file.path(
    co_mat_dir,
    paste0(species,"_",celltype,"_Allpeaks_motif_cooccurrence_matrix.csv")
  )

  message("Reading co-occurrence matrix:")
  message(co_mat_file)

  if (!file.exists(co_mat_file)) {
    warning(
      "Co-occurrence matrix does not exist: ",
      co_mat_file
    )
    return(NULL)
  }

  co_mat <- read.csv(
    co_mat_file,
    row.names = 1,
    check.names = FALSE
  )

  co_mat <- as.matrix(co_mat)


  # ----------------------------------------------------------
  # 2.4 基本检查
  # ----------------------------------------------------------

  message("Matrix dimension: ",nrow(co_mat), " x ",ncol(co_mat))

  if (nrow(co_mat) != ncol(co_mat)) {
    stop(
      celltype,": co-occurrence matrix is not square."
    )
  }

  if (!identical(rownames(co_mat), colnames(co_mat))) {
    stop(
      celltype,": rownames and colnames are different."
    )
  }

  if (!isSymmetric(co_mat)) {
    stop(
      celltype,": co-occurrence matrix is not symmetric."
    )
  }

  # ----------------------------------------------------------
  # 2.5 检查 motif annotation 和 co-occurrence matrix 是否匹配
  # ----------------------------------------------------------

  missing_annotation <- setdiff(
    rownames(co_mat),
    motif_annotation$name
  )

  missing_matrix <- setdiff(
    motif_annotation$name,
    rownames(co_mat)
  )

  if (length(missing_annotation) > 0) {

    message(
      "Patterns in co_mat but missing from motif annotation:"
    )

    print(missing_annotation)
  }

  if (length(missing_matrix) > 0) {

    message(
      "Patterns in motif annotation but missing from co_mat:"
    )

    print(missing_matrix)
  }

  # ----------------------------------------------------------
  # 2.6 保存 motif frequency
  # 对角线 = motif 出现的 peak 数
  # ----------------------------------------------------------

  motif_frequency <- diag(co_mat)

  # 清零对角线，避免影响 motif pair 共现分析
  diag(co_mat) <- 0


  # ----------------------------------------------------------
  # 2.7 提取矩阵上三角，转换成长表
  # ----------------------------------------------------------

  index <- upper.tri(co_mat,diag = FALSE)

  edges <- data.frame(
    motif1 = rownames(co_mat)[row(co_mat)[index]],
    motif2 = colnames(co_mat)[col(co_mat)[index]],
    weight = co_mat[index],
    stringsAsFactors = FALSE
  )

  # 删除没有共现的 motif pair
  edges <- edges %>%
    filter(
      !is.na(weight),
      weight > 0
    )

  message(
    "Number of motif pairs with co-occurrence > 0: ",
    nrow(edges)
  )

  # ----------------------------------------------------------
  # 2.8 选取共现频率最高的 top N motif pair
  # ----------------------------------------------------------

  if (nrow(edges) == 0) {

    warning(
      celltype,
      ": no motif pairs with weight > 0."
    )
    return(NULL)
  }

  edges_plot <- edges %>%
    slice_max(
      order_by = weight,
      n = top_n_edges,
      with_ties = FALSE
    )

  message(
    "Number of edges used for plotting: ",
    nrow(edges_plot)
  )

  # ----------------------------------------------------------
  # 2.9 构建节点信息
  # ----------------------------------------------------------

  nodes <- data.frame(
    name = rownames(co_mat),

    # motif 出现在多少个 peak 中
    motif_frequency = motif_frequency,
    # motif 与其他 motif 的总共现次数
    total_cooccurrence = rowSums(co_mat),
    stringsAsFactors = FALSE

  ) %>%

    left_join(
      motif_annotation,
      by = "name"
    ) %>%
    mutate(
      # 如果没有 tf0 annotation，
      # 则用原始 pattern 名作为标签
      motif_label = ifelse(
        is.na(motif_label) | motif_label == "",
        name,
        motif_label
      )
    )


  # ----------------------------------------------------------
  # 2.10 只保留真正出现在 edges_plot 中的 motif
  # ----------------------------------------------------------

  motifs_in_edges <- unique(
    c(
      edges_plot$motif1,
      edges_plot$motif2
    )
  )


  nodes_plot <- nodes %>%
    filter(name %in% motifs_in_edges)


  # ----------------------------------------------------------
  # 2.11 构建 igraph 网络
  # ----------------------------------------------------------

  g <- graph_from_data_frame(
    d = edges_plot[,c("motif1","motif2","weight")],
    vertices = nodes_plot,
    directed = FALSE
  )


  # ----------------------------------------------------------
  # 2.12 绘制网络
  # ----------------------------------------------------------

  # 为了让不同细胞类型使用相同的随机布局逻辑
  set.seed(123)

  p <- ggraph(
    g,
    layout = "fr"
  ) +

    # 边：颜色深浅表示 motif pair 共现频率

    geom_edge_link(
      aes(colour = weight),
      width = 1,
      alpha = 0.85
    ) +

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

    # 节点：大小表示 motif frequency

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

    # motif annotation

    geom_node_text(
      aes(label = motif_label),
      repel = TRUE,
      size = 4
    ) +

    # 标题

    ggtitle(
      paste0(
        celltype,
        " motif co-occurrence"
      )
    ) +

    theme_void() +

    theme(
      legend.position = "right",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )


  # ----------------------------------------------------------
  # 2.13 输出文件名
  # ----------------------------------------------------------
  plot_dir <- file.path(output_dir,species,"network_results")

  if (is.null(output_filename)) {

    output_filename <- file.path(
      plot_dir,
      paste0(species,"_",celltype,"_motif_cooccurrence_network.pdf")
    )
  }

  # ----------------------------------------------------------
  # 2.14 保存图片
  # ----------------------------------------------------------

  ggsave(
    filename = output_filename,
    plot = p,
    width = 8,
    height = 7
  )

  message("Saved: ",output_filename)


  # ----------------------------------------------------------
  # 2.15 返回结果
  # ----------------------------------------------------------

  return(
    list(
      celltype = celltype,
      co_mat = co_mat,
      motif_annotation = motif_annotation,
      motif_frequency = motif_frequency,
      edges = edges,
      edges_plot = edges_plot,
      nodes = nodes,
      nodes_plot = nodes_plot,
      graph = g,
      plot = p
    )
  )
}

# ============================================================
# 3. 批量处理所有细胞类型
# ============================================================

network_results <- list()


for (celltype in celltypes) {
  network_results[[celltype]] <-
    plot_motif_cooccurrence_network(
      species=species,
      celltype = celltype,
      base_predict_dir = base_predict_dir,
      output_dir = output_dir,
      modisco_report = modisco_report,
      top_n_edges = top_n_edges
    )
}


# ============================================================
# 4. 查看结果
# ============================================================

names(network_results)

# 例如查看 C3.1
network_results[["C3.1"]]$plot

# 查看 C3.1 最强的 motif pair
network_results[["C3.1"]]$edges_plot

# 查看 C3.1 节点信息
network_results[["C3.1"]]$nodes_plot