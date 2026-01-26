

cat("========================================\n")
cat("E vs MVL 综合模块分析\n")
cat("结合单细胞和空间数据的模块通讯分析\n")
cat("========================================\n\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# 0. 加载必要的库
# =============================================================================
cat("0. 加载必要的库...\n")

required_packages <- c("dplyr", "NeuronChat", "ggplot2", "reshape2", "stats", 
                       "ComplexHeatmap", "circlize", "RColorBrewer", "scales", 
                       "ggrepel", "patchwork", "grid", "gridExtra", "viridis")

missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  cat("⚠ 以下包未安装:", paste(missing_packages, collapse=", "), "\n")
  cat("  正在尝试安装...\n")
  for (pkg in missing_packages) {
    tryCatch({
      install.packages(pkg, dependencies = TRUE)
      cat("✓", pkg, "安装成功\n")
    }, error = function(e) {
      cat("✗", pkg, "安装失败:", e$message, "\n")
    })
  }
}

for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "已加载\n")
  } else {
    stop("✗ 无法加载包: ", pkg)
  }
}

# =============================================================================
# 1. 加载 NeuronChat 分析结果
# =============================================================================
cat("\n1. 加载 NeuronChat 分析结果...\n")

neuronchat_path <- NULL
possible_paths <- c(
  "neuronchat_analysis_results_merged/neuronchat_analysis_result.rds",
  "../neuronchat_analysis_results_merged/neuronchat_analysis_result.rds",
  "../../neuronchat_analysis_results_merged/neuronchat_analysis_result.rds",
  "../NeuronChat-main/neuronchat_analysis_results_merged/neuronchat_analysis_result.rds",
  "D:/r/NeuronChat-main/NeuronChat-main/neuronchat_analysis_results_merged/neuronchat_analysis_result.rds",
  "D:/r/data/neuronchat/NeuronChat-main/neuronchat_analysis_results_merged/neuronchat_analysis_result.rds"
)

for (path in possible_paths) {
  if (file.exists(path)) {
    neuronchat_path <- normalizePath(path, winslash = "/")
    cat("✓ 找到文件:", neuronchat_path, "\n")
    break
  }
}

if (is.null(neuronchat_path)) {
  stop("✗ 未找到 NeuronChat 结果文件\n",
       "  请检查文件路径或修改脚本中的路径")
}

tryCatch({
  x <- readRDS(neuronchat_path)
  cat("✓ NeuronChat 对象已加载\n")
  cat("  对象类型:", class(x), "\n")
  cat("  细胞类型数量:", length(unique(x@idents)), "\n")
  cat("  总细胞数:", length(x@idents), "\n")
}, error = function(e) {
  stop("✗ 加载 NeuronChat 对象失败: ", e$message)
})

# 提取单细胞数据
single_cell_data <- x@data
single_cell_meta <- data.frame(
  cell_id = colnames(single_cell_data),
  cell_type = as.character(x@idents),
  stringsAsFactors = FALSE
)

# 检查是否有 anno_level_3_merged（如果存在，使用它）
if ("anno_level_3_merged" %in% names(x@meta)) {
  single_cell_meta$anno_level_3_merged <- as.character(x@meta$anno_level_3_merged)
  cat("✓ 找到 anno_level_3_merged，将使用它作为细胞类型\n")
  # 使用 anno_level_3_merged 作为主要细胞类型
  single_cell_meta$cell_type <- single_cell_meta$anno_level_3_merged
} else {
  cat("⚠ 未找到 anno_level_3_merged，使用默认的 idents\n")
}

cat("✓ 单细胞数据已提取\n")
cat("  - 表达矩阵维度:", dim(single_cell_data), "\n")
cat("  - 总细胞数:", ncol(single_cell_data), "\n")
cat("  - 总基因数:", nrow(single_cell_data), "\n")
cat("  - 细胞类型数:", length(unique(single_cell_meta$cell_type)), "\n")

# =============================================================================
# 2. 加载空间数据
# =============================================================================
cat("\n2. 加载空间数据...\n")

spatial_data_path <- "D:/R/data/neuronchat/spatial neuronchat/rawdata/all_cells_coordinates.csv"

if (!file.exists(spatial_data_path)) {
  alternative_paths <- c(
    "all_cells_coordinates.csv",
    "../all_cells_coordinates.csv",
    "../../all_cells_coordinates.csv",
    "D:/r/data/neurochat/spatial neuronchat/raw data/all_cells_coordinates.csv"
  )
  
  spatial_data_path <- NULL
  for (alt_path in alternative_paths) {
    if (file.exists(alt_path)) {
      spatial_data_path <- normalizePath(alt_path, winslash = "/")
      cat("  使用备选路径:", spatial_data_path, "\n")
      break
    }
  }
  
  if (is.null(spatial_data_path)) {
    stop("✗ 未找到空间数据文件\n")
  }
}

tryCatch({
  spatial_data <- read.csv(spatial_data_path, stringsAsFactors = FALSE, check.names = FALSE)
  cat("✓ 空间数据已加载\n")
  cat("  - 数据维度:", nrow(spatial_data), "行,", ncol(spatial_data), "列\n")
  cat("  - 列名:", paste(colnames(spatial_data), collapse=", "), "\n")
}, error = function(e) {
  stop("✗ 加载空间数据失败: ", e$message)
})

# 检查必需的列
required_cols <- c("Assigned_Cell_Type", "ROI_Region")
if (!all(required_cols %in% colnames(spatial_data))) {
  stop("✗ 空间数据缺少必需的列: ", paste(setdiff(required_cols, colnames(spatial_data)), collapse=", "),
       "\n  可用列:", paste(colnames(spatial_data), collapse=", "))
}

# 检查坐标列
if (!all(c("X_Coordinate", "Y_Coordinate") %in% colnames(spatial_data))) {
  # 尝试其他可能的坐标列名
  coord_cols <- grep("coord|position|location", colnames(spatial_data), ignore.case = TRUE, value = TRUE)
  if (length(coord_cols) >= 2) {
    spatial_data$X_Coordinate <- spatial_data[[coord_cols[1]]]
    spatial_data$Y_Coordinate <- spatial_data[[coord_cols[2]]]
    cat("  ✓ 使用坐标列:", paste(coord_cols, collapse=", "), "\n")
  } else {
    stop("✗ 未找到坐标列\n")
  }
}

# 检查ROI_Region的值
unique_regions <- unique(spatial_data$ROI_Region)
cat("\nROI_Region的唯一值:", paste(unique_regions, collapse=", "), "\n")

# 分离E和MVL区域的数据
e_region_data <- spatial_data[spatial_data$ROI_Region == "e", ]
mvl_region_data <- spatial_data[spatial_data$ROI_Region == "mvl", ]

cat("\n区域统计:\n")
cat("  - E区域细胞数:", nrow(e_region_data), "\n")
cat("  - MVL区域细胞数:", nrow(mvl_region_data), "\n")
cat("  - 其他区域细胞数:", nrow(spatial_data) - nrow(e_region_data) - nrow(mvl_region_data), "\n")

# =============================================================================
# 3. 基于空间数据定义E和MVL模块（结合空间比例和细胞类型比例）
# =============================================================================
cat("\n3. 基于空间数据定义E和MVL模块...\n")

# 3.1 计算E和MVL区域的细胞类型比例
cat("3.1 计算区域细胞类型比例...\n")

e_celltype_counts <- table(e_region_data$Assigned_Cell_Type)
mvl_celltype_counts <- table(mvl_region_data$Assigned_Cell_Type)

e_celltype_props <- prop.table(e_celltype_counts)
mvl_celltype_props <- prop.table(mvl_celltype_counts)

cat("  - E区域细胞类型数:", length(e_celltype_props), "\n")
cat("  - MVL区域细胞类型数:", length(mvl_celltype_props), "\n")

# 3.2 计算细胞类型的富集度（E vs MVL）
cat("\n3.2 计算细胞类型富集度...\n")

# 获取所有细胞类型
all_celltypes <- unique(c(names(e_celltype_props), names(mvl_celltype_props)))

# 计算富集度（E中的比例 / MVL中的比例）
enrichment_scores <- data.frame(
  CellType = all_celltypes,
  E_Proportion = ifelse(all_celltypes %in% names(e_celltype_props), 
                        e_celltype_props[all_celltypes], 0),
  MVL_Proportion = ifelse(all_celltypes %in% names(mvl_celltype_props), 
                          mvl_celltype_props[all_celltypes], 0),
  stringsAsFactors = FALSE
)

# 避免除零
enrichment_scores$MVL_Proportion[enrichment_scores$MVL_Proportion == 0] <- 0.0001
enrichment_scores$E_Proportion[enrichment_scores$E_Proportion == 0] <- 0.0001

# 计算富集度（对数比值比）
enrichment_scores$E_Enrichment <- log2(enrichment_scores$E_Proportion / enrichment_scores$MVL_Proportion)
enrichment_scores$MVL_Enrichment <- log2(enrichment_scores$MVL_Proportion / enrichment_scores$E_Proportion)

# 按富集度排序
enrichment_scores <- enrichment_scores[order(-enrichment_scores$E_Enrichment), ]

cat("\nE区域富集的前10个细胞类型:\n")
print(head(enrichment_scores[order(-enrichment_scores$E_Enrichment), ], 10))

cat("\nMVL区域富集的前10个细胞类型:\n")
print(head(enrichment_scores[order(-enrichment_scores$MVL_Enrichment), ], 10))

# 3.3 基于富集度和统计显著性定义模块
cat("\n3.3 基于富集度定义模块...\n")

# 统计显著性检验（Fisher精确检验）
enrichment_results <- data.frame(
  CellType = character(),
  E_Count = integer(),
  MVL_Count = integer(),
  E_Total = integer(),
  MVL_Total = integer(),
  P_Value = numeric(),
  E_Enrichment = numeric(),
  stringsAsFactors = FALSE
)

for (celltype in all_celltypes) {
  e_count <- sum(e_region_data$Assigned_Cell_Type == celltype)
  mvl_count <- sum(mvl_region_data$Assigned_Cell_Type == celltype)
  e_total <- nrow(e_region_data)
  mvl_total <- nrow(mvl_region_data)
  
  if (e_count + mvl_count > 0) {
    # Fisher精确检验
    contingency_table <- matrix(
      c(e_count, mvl_count, e_total - e_count, mvl_total - mvl_count),
      nrow = 2, ncol = 2
    )
    
    fisher_test <- fisher.test(contingency_table)
    
    e_prop <- ifelse(e_total > 0, e_count / e_total, 0)
    mvl_prop <- ifelse(mvl_total > 0, mvl_count / mvl_total, 0)
    
    e_enrichment <- ifelse(e_prop > 0 & mvl_prop > 0, 
                           log2(e_prop / mvl_prop), 
                           ifelse(e_prop > 0, 10, -10))
    
    enrichment_results <- rbind(enrichment_results, data.frame(
      CellType = celltype,
      E_Count = e_count,
      MVL_Count = mvl_count,
      E_Total = e_total,
      MVL_Total = mvl_total,
      P_Value = fisher_test$p.value,
      E_Enrichment = e_enrichment,
      stringsAsFactors = FALSE
    ))
  }
}

# FDR校正
enrichment_results$P_Value_Adj <- p.adjust(enrichment_results$P_Value, method = "BH")

# 定义E模块：显著富集在E区域且富集度>1（即E中比例是MVL中的2倍以上）
e_module_cell_types <- enrichment_results[
  enrichment_results$P_Value_Adj < 0.05 & 
    enrichment_results$E_Enrichment > 1 & 
    enrichment_results$E_Count >= 10,  # 至少10个细胞
  "CellType"
]

# 定义MVL模块：显著富集在MVL区域且富集度>1
mvl_module_cell_types <- enrichment_results[
  enrichment_results$P_Value_Adj < 0.05 & 
    enrichment_results$E_Enrichment < -1 &  # MVL富集意味着E_Enrichment < -1
    enrichment_results$MVL_Count >= 10,
  "CellType"
]

cat("\n定义的模块:\n")
cat("  - E模块细胞类型数:", length(e_module_cell_types), "\n")
if (length(e_module_cell_types) > 0) {
  cat("  - E模块类型:", paste(e_module_cell_types, collapse=", "), "\n")
}

cat("  - MVL模块细胞类型数:", length(mvl_module_cell_types), "\n")
if (length(mvl_module_cell_types) > 0) {
  cat("  - MVL模块类型:", paste(mvl_module_cell_types, collapse=", "), "\n")
}

# 如果模块定义失败，使用备用方法（基于前N个富集的类型）
if (length(e_module_cell_types) < 2) {
  cat("\n⚠ E模块定义过少，使用备用方法（前5个E富集的类型）\n")
  e_module_cell_types <- head(enrichment_results[order(-enrichment_results$E_Enrichment), "CellType"], 5)
  cat("  - E模块类型（备用）:", paste(e_module_cell_types, collapse=", "), "\n")
}

if (length(mvl_module_cell_types) < 2) {
  cat("\n⚠ MVL模块定义过少，使用备用方法（前5个MVL富集的类型）\n")
  mvl_module_cell_types <- head(enrichment_results[order(enrichment_results$E_Enrichment), "CellType"], 5)
  cat("  - MVL模块类型（备用）:", paste(mvl_module_cell_types, collapse=", "), "\n")
}

# =============================================================================
# 4. 创建模块级别的细胞分配（匹配单细胞数据）
# =============================================================================
cat("\n4. 创建模块级别的细胞分配...\n")

# 匹配细胞类型名称（NeuronChat中使用anno_level_3_merged，空间数据中使用Assigned_Cell_Type）
neuronchat_celltypes <- unique(single_cell_meta$cell_type)
spatial_celltypes <- unique(spatial_data$Assigned_Cell_Type)

# 完全匹配
matched_types <- intersect(neuronchat_celltypes, spatial_celltypes)
cat("  - 完全匹配的细胞类型数:", length(matched_types), "\n")

# 匹配E和MVL模块的细胞类型
e_module_matched <- intersect(e_module_cell_types, matched_types)
mvl_module_matched <- intersect(mvl_module_cell_types, matched_types)

cat("  - E模块匹配的细胞类型数:", length(e_module_matched), "\n")
if (length(e_module_matched) > 0) {
  cat("  - E模块匹配类型:", paste(e_module_matched, collapse=", "), "\n")
}

cat("  - MVL模块匹配的细胞类型数:", length(mvl_module_matched), "\n")
if (length(mvl_module_matched) > 0) {
  cat("  - MVL模块匹配类型:", paste(mvl_module_matched, collapse=", "), "\n")
}

# 分配细胞到模块
single_cell_meta$module <- "Other"
single_cell_meta$module[single_cell_meta$cell_type %in% e_module_matched] <- "E"
single_cell_meta$module[single_cell_meta$cell_type %in% mvl_module_matched] <- "MVL"

# 统计模块细胞数
e_cells <- sum(single_cell_meta$module == "E")
mvl_cells <- sum(single_cell_meta$module == "MVL")
other_cells <- sum(single_cell_meta$module == "Other")

cat("\n模块细胞分布:\n")
cat("  - E模块细胞数:", e_cells, "(", round(e_cells/nrow(single_cell_meta)*100, 2), "%)\n")
cat("  - MVL模块细胞数:", mvl_cells, "(", round(mvl_cells/nrow(single_cell_meta)*100, 2), "%)\n")
cat("  - 其他细胞数:", other_cells, "(", round(other_cells/nrow(single_cell_meta)*100, 2), "%)\n")
cat("  - 总细胞数:", nrow(single_cell_meta), "\n")

# =============================================================================
# 5. 创建模块级别的NeuronChat对象（仅E和MVL模块）
# =============================================================================
cat("\n5. 创建模块级别的NeuronChat对象...\n")

# 过滤数据，仅包含E和MVL模块细胞
module_cells <- single_cell_meta$module %in% c("E", "MVL")
module_data <- single_cell_data[, module_cells]
module_meta <- single_cell_meta[module_cells, ]

cat("  - 模块数据维度:", dim(module_data), "\n")
cat("  - E细胞数:", sum(module_meta$module == "E"), "\n")
cat("  - MVL细胞数:", sum(module_meta$module == "MVL"), "\n")

# 创建模块级别的NeuronChat对象
tryCatch({
  module_neuronchat <- NeuronChat::createNeuronChat(
    object = module_data,
    meta = module_meta,
    group.by = module_meta$module,
    DB = "mouse"
  )
  cat("✓ 模块级别NeuronChat对象创建成功\n")
}, error = function(e) {
  cat("✗ 创建模块级别NeuronChat对象失败:", e$message, "\n")
  module_neuronchat <- NULL
})

# =============================================================================
# 6. 运行模块级别的通讯分析
# =============================================================================
cat("\n6. 运行模块级别的通讯分析...\n")

if (!is.null(module_neuronchat)) {
  tryCatch({
    # 运行NeuronChat分析（增加置换次数以提高统计功效）
    module_neuronchat <- NeuronChat::run_NeuronChat(
      object = module_neuronchat,
      N = 0,
      M = 100,  # 置换次数
      fdr = 0.05,
      K = 0.5
    )
    cat("✓ 模块级别通讯分析完成\n")
    
    # 提取模块级别网络
    module_net <- NeuronChat::net_aggregation(module_neuronchat@net, method = "weight")
    cat("✓ 模块级别网络已提取\n")
    cat("  - 网络维度:", dim(module_net), "\n")
    cat("  - 网络细胞类型:", rownames(module_net), "\n")
    
  }, error = function(e) {
    cat("✗ 运行模块级别分析失败:", e$message, "\n")
    module_net <- NULL
  })
} else {
  cat("✗ 跳过模块级别分析（对象创建失败）\n")
  module_net <- NULL
}

# =============================================================================
# 7. 分析模块内部连接和模块间连接
# =============================================================================
cat("\n7. 分析模块连接情况...\n")

if (!is.null(module_net) && "E" %in% rownames(module_net) && "MVL" %in% rownames(module_net)) {
  
  # 7.1 提取连接强度（归一化以避免细胞数量影响）
  cat("7.1 提取连接强度（归一化）...\n")
  
  # 原始连接强度
  e_internal_raw <- ifelse("E" %in% rownames(module_net) && "E" %in% colnames(module_net),
                           module_net["E", "E"], 0)
  mvl_internal_raw <- ifelse("MVL" %in% rownames(module_net) && "MVL" %in% colnames(module_net),
                             module_net["MVL", "MVL"], 0)
  e_to_mvl_raw <- ifelse("E" %in% rownames(module_net) && "MVL" %in% colnames(module_net),
                         module_net["E", "MVL"], 0)
  mvl_to_e_raw <- ifelse("MVL" %in% rownames(module_net) && "E" %in% colnames(module_net),
                         module_net["MVL", "E"], 0)
  
  # 归一化：除以细胞数量（避免细胞数量影响）
  e_cell_count <- sum(module_meta$module == "E")
  mvl_cell_count <- sum(module_meta$module == "MVL")
  
  # 归一化强度（每1000个细胞的强度）
  e_internal_normalized <- e_internal_raw / e_cell_count * 1000
  mvl_internal_normalized <- mvl_internal_raw / mvl_cell_count * 1000
  e_to_mvl_normalized <- e_to_mvl_raw / sqrt(e_cell_count * mvl_cell_count) * 1000
  mvl_to_e_normalized <- mvl_to_e_raw / sqrt(e_cell_count * mvl_cell_count) * 1000
  
  cat("  - E内部连接强度（原始）:", round(e_internal_raw, 6), "\n")
  cat("  - E内部连接强度（归一化，每1000细胞）:", round(e_internal_normalized, 6), "\n")
  cat("  - MVL内部连接强度（原始）:", round(mvl_internal_raw, 6), "\n")
  cat("  - MVL内部连接强度（归一化，每1000细胞）:", round(mvl_internal_normalized, 6), "\n")
  cat("  - E→MVL连接强度（原始）:", round(e_to_mvl_raw, 6), "\n")
  cat("  - E→MVL连接强度（归一化）:", round(e_to_mvl_normalized, 6), "\n")
  cat("  - MVL→E连接强度（原始）:", round(mvl_to_e_raw, 6), "\n")
  cat("  - MVL→E连接强度（归一化）:", round(mvl_to_e_normalized, 6), "\n")
  
  # 7.2 分析LR对
  cat("\n7.2 分析LR对...\n")
  
  if (!is.null(module_neuronchat) && !is.null(module_neuronchat@net)) {
    lr_pairs_module <- names(module_neuronchat@net)
    cat("  - 模块级别LR对数:", length(lr_pairs_module), "\n")
    
    # 提取LR对强度（按连接类型分类）
    e_internal_lr_pairs <- list()
    mvl_internal_lr_pairs <- list()
    e_to_mvl_lr_pairs <- list()
    mvl_to_e_lr_pairs <- list()
    
    for (lr_pair in lr_pairs_module) {
      if (lr_pair %in% names(module_neuronchat@net)) {
        lr_net <- module_neuronchat@net[[lr_pair]]
        
        # E内部
        if ("E" %in% rownames(lr_net) && "E" %in% colnames(lr_net)) {
          strength <- lr_net["E", "E"]
          if (strength > 0) {
            e_internal_lr_pairs[[lr_pair]] <- strength / e_cell_count * 1000  # 归一化
          }
        }
        
        # MVL内部
        if ("MVL" %in% rownames(lr_net) && "MVL" %in% colnames(lr_net)) {
          strength <- lr_net["MVL", "MVL"]
          if (strength > 0) {
            mvl_internal_lr_pairs[[lr_pair]] <- strength / mvl_cell_count * 1000  # 归一化
          }
        }
        
        # E→MVL
        if ("E" %in% rownames(lr_net) && "MVL" %in% colnames(lr_net)) {
          strength <- lr_net["E", "MVL"]
          if (strength > 0) {
            e_to_mvl_lr_pairs[[lr_pair]] <- strength / sqrt(e_cell_count * mvl_cell_count) * 1000  # 归一化
          }
        }
        
        # MVL→E
        if ("MVL" %in% rownames(lr_net) && "E" %in% colnames(lr_net)) {
          strength <- lr_net["MVL", "E"]
          if (strength > 0) {
            mvl_to_e_lr_pairs[[lr_pair]] <- strength / sqrt(e_cell_count * mvl_cell_count) * 1000  # 归一化
          }
        }
      }
    }
    
    cat("  - E内部LR对数:", length(e_internal_lr_pairs), "\n")
    cat("  - MVL内部LR对数:", length(mvl_internal_lr_pairs), "\n")
    cat("  - E→MVL LR对数:", length(e_to_mvl_lr_pairs), "\n")
    cat("  - MVL→E LR对数:", length(mvl_to_e_lr_pairs), "\n")
    
    # 计算LR密度（活跃LR对 / 总LR对）
    total_lr_pairs <- length(lr_pairs_module)
    e_internal_lr_density <- length(e_internal_lr_pairs) / total_lr_pairs
    mvl_internal_lr_density <- length(mvl_internal_lr_pairs) / total_lr_pairs
    cross_lr_density <- (length(e_to_mvl_lr_pairs) + length(mvl_to_e_lr_pairs)) / total_lr_pairs
    
    cat("\n  LR密度（活跃LR对占比）:\n")
    cat("    - E内部LR密度:", round(e_internal_lr_density, 4), "\n")
    cat("    - MVL内部LR密度:", round(mvl_internal_lr_density, 4), "\n")
    cat("    - 模块间LR密度:", round(cross_lr_density, 4), "\n")
    
  } else {
    cat("⚠ 无法访问LR对数据\n")
    e_internal_lr_pairs <- list()
    mvl_internal_lr_pairs <- list()
    e_to_mvl_lr_pairs <- list()
    mvl_to_e_lr_pairs <- list()
    e_internal_lr_density <- 0
    mvl_internal_lr_density <- 0
    cross_lr_density <- 0
  }
  
  # 7.3 分析信息流（in/out）
  cat("\n7.3 分析信息流（in/out）...\n")
  
  # E模块的信息流
  e_outgoing <- e_internal_raw + e_to_mvl_raw  # E发出的信息
  e_incoming <- e_internal_raw + mvl_to_e_raw  # E接收的信息
  
  # MVL模块的信息流
  mvl_outgoing <- mvl_internal_raw + mvl_to_e_raw  # MVL发出的信息
  mvl_incoming <- mvl_internal_raw + e_to_mvl_raw  # MVL接收的信息
  
  # 归一化信息流
  e_outgoing_norm <- e_outgoing / e_cell_count * 1000
  e_incoming_norm <- e_incoming / e_cell_count * 1000
  mvl_outgoing_norm <- mvl_outgoing / mvl_cell_count * 1000
  mvl_incoming_norm <- mvl_incoming / mvl_cell_count * 1000
  
  cat("  E模块信息流:\n")
  cat("    - 发出（outgoing）:", round(e_outgoing, 6), "（归一化:", round(e_outgoing_norm, 6), "）\n")
  cat("    - 接收（incoming）:", round(e_incoming, 6), "（归一化:", round(e_incoming_norm, 6), "）\n")
  cat("    - 净信息流（outgoing - incoming）:", round(e_outgoing - e_incoming, 6), "\n")
  
  cat("  MVL模块信息流:\n")
  cat("    - 发出（outgoing）:", round(mvl_outgoing, 6), "（归一化:", round(mvl_outgoing_norm, 6), "）\n")
  cat("    - 接收（incoming）:", round(mvl_incoming, 6), "（归一化:", round(mvl_incoming_norm, 6), "）\n")
  cat("    - 净信息流（outgoing - incoming）:", round(mvl_outgoing - mvl_incoming, 6), "\n")
  
  # 7.4 分析模块间连接的平衡性
  cat("\n7.4 分析模块间连接平衡性...\n")
  
  total_cross_strength <- e_to_mvl_raw + mvl_to_e_raw
  if (total_cross_strength > 0) {
    e_to_mvl_prop <- e_to_mvl_raw / total_cross_strength
    mvl_to_e_prop <- mvl_to_e_raw / total_cross_strength
    
    # 平衡性指标（对称性）
    balance_index <- 1 - abs(e_to_mvl_prop - mvl_to_e_prop)  # 0-1，1表示完全平衡
    asymmetry_index <- abs(e_to_mvl_prop - mvl_to_e_prop)  # 0-1，0表示完全平衡
    
    # 相对差异
    if (mvl_to_e_raw > 0) {
      relative_diff <- (e_to_mvl_raw - mvl_to_e_raw) / mvl_to_e_raw * 100
    } else {
      relative_diff <- Inf
    }
    
    cat("  - E→MVL占比:", round(e_to_mvl_prop, 4), "(", round(e_to_mvl_prop*100, 2), "%)\n")
    cat("  - MVL→E占比:", round(mvl_to_e_prop, 4), "(", round(mvl_to_e_prop*100, 2), "%)\n")
    cat("  - 平衡性指数（1=完全平衡）:", round(balance_index, 4), "\n")
    cat("  - 不对称性指数（0=完全平衡）:", round(asymmetry_index, 4), "\n")
    cat("  - 相对差异:", round(relative_diff, 2), "%\n")
    
    # 统计显著性检验（二项检验）
    if (total_cross_strength > 0) {
      # 将强度转换为计数（用于统计检验）
      e_to_mvl_count <- round(e_to_mvl_raw * 1000)
      mvl_to_e_count <- round(mvl_to_e_raw * 1000)
      
      if (e_to_mvl_count + mvl_to_e_count > 0) {
        binom_test <- binom.test(e_to_mvl_count, e_to_mvl_count + mvl_to_e_count, p = 0.5)
        cat("  - 二项检验p值:", format(binom_test$p.value, scientific = TRUE), "\n")
        cat("  - 显著不对称:", binom_test$p.value < 0.05, "\n")
      }
    }
  }
  
  # 存储结果
  module_analysis_results <- list(
    # 原始强度
    e_internal_raw = e_internal_raw,
    mvl_internal_raw = mvl_internal_raw,
    e_to_mvl_raw = e_to_mvl_raw,
    mvl_to_e_raw = mvl_to_e_raw,
    
    # 归一化强度
    e_internal_normalized = e_internal_normalized,
    mvl_internal_normalized = mvl_internal_normalized,
    e_to_mvl_normalized = e_to_mvl_normalized,
    mvl_to_e_normalized = mvl_to_e_normalized,
    
    # LR对
    e_internal_lr_pairs = e_internal_lr_pairs,
    mvl_internal_lr_pairs = mvl_internal_lr_pairs,
    e_to_mvl_lr_pairs = e_to_mvl_lr_pairs,
    mvl_to_e_lr_pairs = mvl_to_e_lr_pairs,
    
    # LR密度
    e_internal_lr_density = e_internal_lr_density,
    mvl_internal_lr_density = mvl_internal_lr_density,
    cross_lr_density = cross_lr_density,
    
    # 信息流
    e_outgoing = e_outgoing,
    e_incoming = e_incoming,
    mvl_outgoing = mvl_outgoing,
    mvl_incoming = mvl_incoming,
    e_outgoing_norm = e_outgoing_norm,
    e_incoming_norm = e_incoming_norm,
    mvl_outgoing_norm = mvl_outgoing_norm,
    mvl_incoming_norm = mvl_incoming_norm,
    
    # 平衡性
    balance_index = balance_index,
    asymmetry_index = asymmetry_index,
    relative_diff = relative_diff,
    
    # 细胞数量
    e_cell_count = e_cell_count,
    mvl_cell_count = mvl_cell_count
  )
  
} else {
  cat("✗ 模块网络不可用\n")
  module_analysis_results <- NULL
}

# =============================================================================
# 8. 创建输出目录并保存结果
# =============================================================================
cat("\n8. 创建输出目录并保存结果...\n")

output_dir <- "E_MVL_integrated_module_analysis_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

plots_dir <- file.path(output_dir, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

cat("✓ 输出目录已创建:", output_dir, "\n")

# 保存模块定义结果
if (exists("enrichment_results")) {
  write.csv(enrichment_results, 
            file.path(output_dir, "cell_type_enrichment_analysis.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("✓ 细胞类型富集分析结果已保存\n")
}

# 保存模块连接分析结果
if (!is.null(module_analysis_results)) {
  # 汇总统计
  summary_stats <- data.frame(
    Metric = c(
      "E_Cell_Count", "MVL_Cell_Count",
      "E_Internal_Strength_Raw", "MVL_Internal_Strength_Raw",
      "E_Internal_Strength_Normalized", "MVL_Internal_Strength_Normalized",
      "E_to_MVL_Strength_Raw", "MVL_to_E_Strength_Raw",
      "E_to_MVL_Strength_Normalized", "MVL_to_E_Strength_Normalized",
      "E_Internal_LR_Count", "MVL_Internal_LR_Count",
      "E_to_MVL_LR_Count", "MVL_to_E_LR_Count",
      "E_Internal_LR_Density", "MVL_Internal_LR_Density", "Cross_LR_Density",
      "E_Outgoing_Info", "E_Incoming_Info", "E_Net_Info_Flow",
      "MVL_Outgoing_Info", "MVL_Incoming_Info", "MVL_Net_Info_Flow",
      "Balance_Index", "Asymmetry_Index", "Relative_Difference"
    ),
    Value = c(
      module_analysis_results$e_cell_count,
      module_analysis_results$mvl_cell_count,
      module_analysis_results$e_internal_raw,
      module_analysis_results$mvl_internal_raw,
      module_analysis_results$e_internal_normalized,
      module_analysis_results$mvl_internal_normalized,
      module_analysis_results$e_to_mvl_raw,
      module_analysis_results$mvl_to_e_raw,
      module_analysis_results$e_to_mvl_normalized,
      module_analysis_results$mvl_to_e_normalized,
      length(module_analysis_results$e_internal_lr_pairs),
      length(module_analysis_results$mvl_internal_lr_pairs),
      length(module_analysis_results$e_to_mvl_lr_pairs),
      length(module_analysis_results$mvl_to_e_lr_pairs),
      module_analysis_results$e_internal_lr_density,
      module_analysis_results$mvl_internal_lr_density,
      module_analysis_results$cross_lr_density,
      module_analysis_results$e_outgoing,
      module_analysis_results$e_incoming,
      module_analysis_results$e_outgoing - module_analysis_results$e_incoming,
      module_analysis_results$mvl_outgoing,
      module_analysis_results$mvl_incoming,
      module_analysis_results$mvl_outgoing - module_analysis_results$mvl_incoming,
      module_analysis_results$balance_index,
      module_analysis_results$asymmetry_index,
      module_analysis_results$relative_diff
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_stats,
            file.path(output_dir, "module_connection_summary.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  
  # 保存LR对详细信息
  if (length(module_analysis_results$e_internal_lr_pairs) > 0) {
    e_internal_df <- data.frame(
      LR_Pair = names(module_analysis_results$e_internal_lr_pairs),
      Strength_Normalized = unlist(module_analysis_results$e_internal_lr_pairs),
      Direction = "E_Internal",
      stringsAsFactors = FALSE
    )
    write.csv(e_internal_df,
              file.path(output_dir, "e_internal_lr_pairs.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
  
  if (length(module_analysis_results$mvl_internal_lr_pairs) > 0) {
    mvl_internal_df <- data.frame(
      LR_Pair = names(module_analysis_results$mvl_internal_lr_pairs),
      Strength_Normalized = unlist(module_analysis_results$mvl_internal_lr_pairs),
      Direction = "MVL_Internal",
      stringsAsFactors = FALSE
    )
    write.csv(mvl_internal_df,
              file.path(output_dir, "mvl_internal_lr_pairs.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
  
  if (length(module_analysis_results$e_to_mvl_lr_pairs) > 0) {
    e_to_mvl_df <- data.frame(
      LR_Pair = names(module_analysis_results$e_to_mvl_lr_pairs),
      Strength_Normalized = unlist(module_analysis_results$e_to_mvl_lr_pairs),
      Direction = "E_to_MVL",
      stringsAsFactors = FALSE
    )
    write.csv(e_to_mvl_df,
              file.path(output_dir, "e_to_mvl_lr_pairs.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
  
  if (length(module_analysis_results$mvl_to_e_lr_pairs) > 0) {
    mvl_to_e_df <- data.frame(
      LR_Pair = names(module_analysis_results$mvl_to_e_lr_pairs),
      Strength_Normalized = unlist(module_analysis_results$mvl_to_e_lr_pairs),
      Direction = "MVL_to_E",
      stringsAsFactors = FALSE
    )
    write.csv(mvl_to_e_df,
              file.path(output_dir, "mvl_to_e_lr_pairs.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
  
  # 保存模块网络
  if (!is.null(module_net)) {
    write.csv(module_net,
              file.path(output_dir, "module_level_network.csv"),
              fileEncoding = "UTF-8")
  }
  
  cat("✓ 结果已保存\n")
}

# =============================================================================
