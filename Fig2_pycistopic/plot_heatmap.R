# 1. 安装并加载包
# if (!require("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(stringr)  # 用于自然排序 
library(RColorBrewer)
# 2. 读取数据
df <- read.csv("~/Treeshrew/pycistopic_260418/TS_Average_Topic_Subtype_ID_Num.csv",row.names=1)
mat <- as.matrix(df)
all_colnames <- colnames(mat)
sorted_cols <- str_sort(all_colnames, numeric = TRUE)
mat <- mat[, sorted_cols]
#
# 找出每一行（每一个 Topic）最大值所在的列索引
max_col_index <- apply(mat, 1, which.max)
# 根据最大值出现的列位置进行排序；如果位置相同，则按该位置的数值大小排序
topic_order <- order(max_col_index, -apply(mat, 1, max))
# 重新排列矩阵的行
mat <- mat[topic_order, ]

# 3. 颜色映射修改
# 方式 A：自动映射（由 ComplexHeatmap 根据数据最大/最小值自动决定）
# 方式 B：如果你想要更鲜艳的渐变，可以手动定义最大值
max_val <- max(mat, na.rm = TRUE)
min_val <- min(mat, na.rm = TRUE)
# 使用 colorRamp2 定义一个覆盖全量程的颜色梯度
# 比如：最小值是白色/蓝色，中间值是橙色，最大值（14）是深红色
#col_fun = colorRamp2(c(min_val, (min_val + max_val)/2, max_val), 
 #                    c("#362E93","#4EAB54", "#EEDF4A"))
col_fun = colorRamp2(c(min_val, (min_val + max_val)/2, max_val), 
                     c("#FFFFFF","#FDAE61", "#D73027"))
col_fun = colorRamp2(c(min_val, (min_val + max_val)/2, max_val), 
                     c("#FFFFE5", "#FDAE61", "#A50F15"))
col_fun = colorRamp2(c(min_val, (min_val + max_val)/2, max_val), 
                     c("#2166AC", "white", "#B2182B"))

# 4.提取细胞分类信息
# 具体的细胞类型名称
cell_types <- colnames(mat)
# 自动提取大类：匹配 "C" 后面直到第一个点 "." 之前的数字
# 例如 "C1.1" -> "C1", "C10.11" -> "C10"
major_types <- str_extract(cell_types, "C\\d+")
# A. 为大类分配颜色 (如 C1, C2, C3...)
TS.co.for.Main <- c("#DDAD3B","#00BFFF","#9970AB","#00FF7F",## 1-4
                    "#FDCDAC","#3288BD","#32CD32","#48D1CC",## 5-8
                    "#FF33CC","#B15928","#6495ED","#D5004D",## 9-12
                    "#0000FF","#C2A5CF","#EE82EE","#FF7F00",## 13-16
                    "#9E0142","#B3CDE3","#87CEFA","#00FFFF",## 17-20
                    "#CDE8B9","#F4CAE4","#0099FF","#DE77AE")## 21-24
#
# 构造成带有名字的颜色向量（这就是 ComplexHeatmap 要求的格式）
major_levels <- str_sort(unique(major_types), numeric = TRUE)
major_colors <- setNames(TS.co.for.Main, major_levels) #给颜色向量加名字，生成一个“类别 → 颜色”的对应关系


# B. 为具体类型分配颜色 (C1.1, C1.2...)
# 如果类型太多，建议使用渐变色板，或者直接随机生成
#brewer.pal(12, "Paired")会从 RColorBrewer 的 "Paired" 调色板里取 12 个基础颜色
#colorRampPalette(brewer.pal(12, "Paired"))(length(all_types))会基于这 12 个颜色做插值，生成和 all_types 数量一样多的颜色
#smaple()会把这 50 个颜色随机打乱顺序.为了把颜色顺序随机打乱，避免相邻 cell type 总是拿到调色板中很相近的颜色。
set.seed(42)
all_types <- unique(cell_types)
type_colors <- setNames(
  sample(colorRampPalette(brewer.pal(12, "Paired"))(length(all_types))), 
  all_types
)

# 5. 创建顶部的多层注释 (Annotation) ---
top_ha <- HeatmapAnnotation(
  MajorType = major_types,         # 第一层：大类
  CellType = cell_types,           # 第二层：具体类型
  col = list(
    MajorType = major_colors, 
    CellType = type_colors
  ),
  annotation_legend_param = list(
    MajorType = list(
      title = "Major Groups",
      at = major_levels,        # 强制顺序：C1, C2, C3...
      labels = major_levels,
      ncol = 2                  # 分成两列，防止太长
    )
  ),
  show_legend = c(TRUE, FALSE),    # 通常只显示大类的图例，具体类型的图例太长可以隐藏
  annotation_name_side = "left",
  simple_anno_size = unit(0.4, "cm") # 调整 color bar 的高度
)

# 4. 绘制热图
#pdf("Topic_Heatmap_R.pdf", width = 8, height = 10)
ht <- Heatmap(mat, 
        name = "Mean topic probability", 
        col = col_fun,                # 使用我们定义的颜色映射
        top_annotation = top_ha,       # 添加刚才定义的注释条
        cluster_rows = FALSE,          # 对 Topic 进行聚类
        cluster_columns = FALSE,       # 对细胞类型进行聚类
        show_row_names = FALSE, 
        show_column_names = FALSE,
        column_title = NULL,
        row_title = "Topics",
        use_raster = TRUE,
        raster_quality = 2,
        raster_device = "png",
        #rect_gp = gpar(col = "white", lwd = 0.5), # 给方块加个细白边，更好看
        row_names_gp = gpar(fontsize = 8),        # 调整字体大小
        column_names_gp = gpar(fontsize = 10))
pdf(
    "~/Treeshrew/pycistopic_260418/TS_Average_Topic_Subtype.pdf",
    width = 20,
    height = 10,
    useDingbats = FALSE
)

draw(
    ht,
    annotation_legend_side = "right",
    heatmap_legend_side = "right"
)
dev.off()
write.csv(
  mat,
  file = "~/Treeshrew/pycistopic_260418/TS_Average_Topic_Subtype_ID_Num_sortedForplot.csv",
  quote = FALSE
)


#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################



df <- read.csv("~/Treeshrew/pycistopic_260418/TS_Topic_pseudocell.csv",row.names=1)
mat <- as.matrix(df)
all_colnames <- colnames(mat)
sorted_cols <- str_sort(all_colnames, numeric = TRUE)#numeric = TRUE 是对整个字符串中所有数字片段做自然排序
mat <- mat[, sorted_cols]
#
# 找出每一行（每一个 Topic）最大值所在的列索引
max_col_index <- apply(mat, 1, which.max)
# 根据最大值出现的列位置进行排序；如果位置相同，则按该位置的数值大小排序
topic_order <- order(max_col_index, -apply(mat, 1, max))
# 重新排列矩阵的行
mat <- mat[topic_order, ]

# 3. 颜色映射修改
# 方式 A：自动映射（由 ComplexHeatmap 根据数据最大/最小值自动决定）
# 方式 B：如果你想要更鲜艳的渐变，可以手动定义最大值
max_val <- max(mat, na.rm = TRUE)
min_val <- min(mat, na.rm = TRUE)
# 使用 colorRamp2 定义一个覆盖全量程的颜色梯度
# 比如：最小值是白色/蓝色，中间值是橙色，最大值（14）是深红色
col_fun = colorRamp2(c(min_val, (min_val + max_val)/2, max_val), 
                     c("#362E93","#4EAB54", "#EEDF4A"))
col_fun = colorRamp2(c(min_val, (min_val + max_val)/2, max_val), 
                     c("#FFFFFF","#FDAE61", "#D73027"))
col_fun = colorRamp2(c(min_val, (min_val + max_val)/2, max_val), 
                     c("#FFFFE5", "#FDAE61", "#A50F15"))
library(viridisLite)
col_fun <- colorRamp2(
  seq(min_val, max_val, length.out = 256),
  viridis(256)
)#生成 256 个颜色，用于构建连续的热图颜色梯度

# 4.提取细胞分类信息
# 具体的细胞类型名称
cell_id <- colnames(mat)
cell_types = str_split_i(cell_id, "__", 1)
# 自动提取大类：匹配 "C" 后面直到第一个点 "." 之前的数字
# 例如 "C1.1" -> "C1", "C10.11" -> "C10"
major_types <- str_extract(cell_types, "C\\d+")
# A. 为大类分配颜色 (如 C1, C2, C3...)
TS.co.for.Main <- c("#DDAD3B","#00BFFF","#9970AB","#00FF7F",## 1-4
                    "#FDCDAC","#3288BD","#32CD32","#48D1CC",## 5-8
                    "#FF33CC","#B15928","#6495ED","#D5004D",## 9-12
                    "#0000FF","#C2A5CF","#EE82EE","#FF7F00",## 13-16
                    "#9E0142","#B3CDE3","#87CEFA","#00FFFF",## 17-20
                    "#CDE8B9","#F4CAE4","#0099FF","#DE77AE")## 21-24
#
# 构造成带有名字的颜色向量（这就是 ComplexHeatmap 要求的格式）
major_levels <- str_sort(unique(major_types), numeric = TRUE)
major_colors <- setNames(TS.co.for.Main, major_levels)

# B. 为具体类型分配颜色 (C1.1, C1.2...)
# 如果类型太多，建议使用渐变色板，或者直接随机生成
set.seed(42)
all_types <- unique(cell_types)
type_colors <- setNames(
  sample(colorRampPalette(brewer.pal(12, "Paired"))(length(all_types))), 
  all_types
)
#library(ArchR)
#pal <- paletteContinuous("solarExtra") #把 solarExtra 调色板中的一段颜色，按顺序分配给 all_types 里的每个 cell type
#idx <- round(
 # seq(1, length(pal), length.out = length(all_types))
#)
#type_colors <- setNames(
 # pal[idx],
  #all_types
#)

# 5. 创建顶部的多层注释 (Annotation) ---
top_ha <- HeatmapAnnotation(
  MajorType = major_types,         # 第一层：大类
  CellType = cell_types,           # 第二层：具体类型
  col = list(
    MajorType = major_colors, 
    CellType = type_colors
  ),
  annotation_legend_param = list(
    MajorType = list(
      title = "Major Groups",
      at = major_levels,        # 强制顺序：C1, C2, C3...
      labels = major_levels,
      ncol = 2                  # 分成两列，防止太长
    )
  ),
  show_legend = c(TRUE, FALSE),    # 通常只显示大类的图例，具体类型的图例太长可以隐藏
  annotation_name_side = "left",
  simple_anno_size = unit(0.4, "cm") # 调整 color bar 的高度
)

# 4. 绘制热图
#pdf("Topic_Heatmap_R.pdf", width = 8, height = 10)
ht <- Heatmap(mat, 
        name = "Topic probability", 
        col = col_fun,                # 使用我们定义的颜色映射
        top_annotation = top_ha,       # 添加刚才定义的注释条
        cluster_rows = FALSE,          # 对 Topic 进行聚类
        cluster_columns = FALSE,       # 对细胞类型进行聚类
        show_row_names = FALSE, 
        show_column_names = FALSE,
        column_title = NULL,
        row_title = "Topics",
        use_raster = TRUE,
        raster_quality = 4,
        raster_device = "png",
        #rect_gp = gpar(col = "white", lwd = 0.5), # 给方块加个细白边，更好看
        row_names_gp = gpar(fontsize = 8),        # 调整字体大小
        column_names_gp = gpar(fontsize = 10))
pdf(
    "~/Treeshrew/pycistopic_260418/TS_Topic_pseudocell.pdf",
    width = 20,
    height = 10,
    useDingbats = FALSE
)
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()
write.csv(
  mat,
  file = "~/Treeshrew/pycistopic_260418/TS_Topic_pseudocell_sortedForplot.csv",
  quote = FALSE
)
