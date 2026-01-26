
# 创建鸡同源鼠基因Seurat对象

library(Seurat)
library(dplyr)

cat("=== Creating Chicken Ortholog Mouse Seurat Object ===\n")
cat("Starting analysis at:", Sys.time(), "\n")

# 1. 加载原始鸡单细胞数据
cat("\n1. Loading original chicken single-cell data...\n")
chicken_seurat <- readRDS("chicken_mn_Ex_Inh_seurat.rds")
cat("✓ Original chicken data loaded:", ncol(chicken_seurat), "cells,", nrow(chicken_seurat), "genes\n")

# 2. 加载转换后的表达矩阵
cat("\n2. Loading converted expression matrix...\n")
if (file.exists("chicken_to_mouse_converted_matrix.rds")) {
  converted_expr_matrix <- readRDS("chicken_to_mouse_converted_matrix.rds")
  cat("✓ Converted expression matrix loaded:", nrow(converted_expr_matrix), "genes,", ncol(converted_expr_matrix), "cells\n")
} else {
  stop("✗ Converted expression matrix not found: chicken_to_mouse_converted_matrix.rds")
}

# 3. 检查数据一致性
cat("\n3. Checking data consistency...\n")

# 检查细胞数量是否一致
if (ncol(converted_expr_matrix) != ncol(chicken_seurat)) {
  cat("⚠ Warning: Cell count mismatch\n")
  cat("  Original cells:", ncol(chicken_seurat), "\n")
  cat("  Converted cells:", ncol(converted_expr_matrix), "\n")
  
  # 如果细胞数量不匹配，尝试匹配细胞ID
  original_cells <- colnames(chicken_seurat)
  converted_cells <- colnames(converted_expr_matrix)
  
  common_cells <- intersect(original_cells, converted_cells)
  cat("  Common cells:", length(common_cells), "\n")
  
  if (length(common_cells) > 0) {
    cat("  Using common cells only\n")
    chicken_seurat <- chicken_seurat[, common_cells]
    converted_expr_matrix <- converted_expr_matrix[, common_cells]
  } else {
    stop("✗ No common cells found between original and converted data")
  }
} else {
  cat("✓ Cell count matches\n")
}

# 检查细胞ID是否一致
if (!all(colnames(chicken_seurat) == colnames(converted_expr_matrix))) {
  cat("⚠ Warning: Cell IDs don't match exactly\n")
  cat("  Reordering converted matrix to match original...\n")
  converted_expr_matrix <- converted_expr_matrix[, colnames(chicken_seurat)]
}

cat("✓ Data consistency check completed\n")
cat("  Final cells:", ncol(chicken_seurat), "\n")
cat("  Final genes:", nrow(converted_expr_matrix), "\n")

# 4. 创建新的Seurat对象
cat("\n4. Creating new Seurat object...\n")

# 创建新的Seurat对象，使用转换后的表达矩阵
chicken_orth_mice_seurat <- CreateSeuratObject(
  counts = converted_expr_matrix,
  project = "ChickenOrthologMouse",
  assay = "RNA"
)

cat("✓ New Seurat object created\n")
cat("  Cells:", ncol(chicken_orth_mice_seurat), "\n")
cat("  Genes:", nrow(chicken_orth_mice_seurat), "\n")

# 4.1 标准化数据以填充data层
cat("\n4.1 Normalizing data to fill data layer...\n")

# 首先进行基本标准化
chicken_orth_mice_seurat <- NormalizeData(chicken_orth_mice_seurat, verbose = FALSE)
cat("✓ Basic normalization completed\n")

# 4.2 使用SCTransform进行高级标准化
cat("\n4.2 Performing SCTransform normalization...\n")

# 检查可用的回归变量
available_vars <- colnames(chicken_orth_mice_seurat@meta.data)
cat("Available metadata variables:", paste(available_vars, collapse = ", "), "\n")

# 选择可用的回归变量
regress_vars <- c()
if ("nCount_RNA" %in% available_vars) {
  regress_vars <- c(regress_vars, "nCount_RNA")
}
if ("percent_mt" %in% available_vars) {
  regress_vars <- c(regress_vars, "percent_mt")
}
if ("nFeature_RNA" %in% available_vars) {
  regress_vars <- c(regress_vars, "nFeature_RNA")
}

cat("Variables to regress:", paste(regress_vars, collapse = ", "), "\n")

tryCatch({
  if (length(regress_vars) > 0) {
    chicken_orth_mice_seurat <- SCTransform(chicken_orth_mice_seurat, 
                                           vars.to.regress = regress_vars,
                                           verbose = FALSE)
  } else {
    chicken_orth_mice_seurat <- SCTransform(chicken_orth_mice_seurat, 
                                           verbose = FALSE)
  }
  cat("✓ SCTransform normalization completed\n")
}, error = function(e) {
  cat("⚠ SCTransform failed:", e$message, "\n")
  cat("  Continuing with basic normalization only\n")
})

# 检查层是否已填充
available_layers <- Layers(chicken_orth_mice_seurat)
cat("Available layers after normalization:", paste(available_layers, collapse = ", "), "\n")

# 检查SCTransform结果（Seurat v5兼容）
if ("scale.data" %in% available_layers) {
  cat("✓ SCTransform completed successfully (scale.data layer created)\n")
  cat("✓ This is the correct behavior for Seurat v5\n")
} else if ("SCT" %in% available_layers) {
  cat("✓ SCT layer created successfully (Seurat v4 behavior)\n")
} else {
  cat("⚠ SCTransform may have failed - no scale.data or SCT layer found\n")
}

# 5. 复制原始元数据
cat("\n5. Copying metadata from original Seurat object...\n")

# 复制所有元数据
chicken_orth_mice_seurat@meta.data <- chicken_seurat@meta.data

# 添加基因转换信息
chicken_orth_mice_seurat@misc$gene_conversion_info <- list(
  original_genes = nrow(chicken_seurat),
  converted_genes = nrow(converted_expr_matrix),
  conversion_method = "weighted_average",
  conversion_date = Sys.time()
)

cat("✓ Metadata copied and gene conversion info added\n")

# 6. 检查细胞类型信息
cat("\n6. Checking cell type information...\n")

# 检查细胞类型列
if ("anno_level_3_final" %in% colnames(chicken_orth_mice_seurat@meta.data)) {
  cell_types <- chicken_orth_mice_seurat@meta.data$anno_level_3_final
  cat("✓ Found anno_level_3_final column\n")
} else if ("anno_level_3" %in% colnames(chicken_orth_mice_seurat@meta.data)) {
  # 如果anno_level_3_final不存在，从anno_level_3创建它
  cat("⚠ anno_level_3_final not found, creating from anno_level_3\n")
  chicken_orth_mice_seurat@meta.data$anno_level_3_final <- chicken_orth_mice_seurat@meta.data$anno_level_3
  cell_types <- chicken_orth_mice_seurat@meta.data$anno_level_3_final
  cat("✓ Created anno_level_3_final column from anno_level_3\n")
} else {
  cat("⚠ No cell type columns found\n")
}

if (exists("cell_types")) {
  cat("  Cell types:", length(unique(cell_types)), "\n")
  cat("  Sample cell types:", paste(head(unique(cell_types), 5), collapse = ", "), "\n")
  
  # 检查细胞类型分布
  cell_type_counts <- table(cell_types)
  cat("  Cell type distribution:\n")
  for (i in 1:min(5, length(cell_type_counts))) {
    cat("    ", names(cell_type_counts)[i], ":", cell_type_counts[i], "\n")
  }
}

# 7. 添加基本QC指标
cat("\n7. Adding basic QC metrics...\n")

# 计算基本QC指标
chicken_orth_mice_seurat$nFeature_RNA <- colSums(GetAssayData(chicken_orth_mice_seurat, layer = "counts") > 0)
chicken_orth_mice_seurat$nCount_RNA <- colSums(GetAssayData(chicken_orth_mice_seurat, layer = "counts"))
chicken_orth_mice_seurat$percent_mt <- PercentageFeatureSet(chicken_orth_mice_seurat, pattern = "^MT-")

cat("✓ Basic QC metrics added\n")
cat("  Features per cell range:", range(chicken_orth_mice_seurat$nFeature_RNA), "\n")
cat("  Counts per cell range:", range(chicken_orth_mice_seurat$nCount_RNA), "\n")
cat("  MT percentage range:", range(chicken_orth_mice_seurat$percent_mt), "\n")

# 8. 检查基因表达情况
cat("\n8. Checking gene expression...\n")

# 检查可用的层
available_layers <- Layers(chicken_orth_mice_seurat)
cat("Available layers:", paste(available_layers, collapse = ", "), "\n")

# 使用counts层进行基因表达统计（因为data层可能为空）
if ("counts" %in% available_layers) {
  cat("Using counts layer for gene expression analysis\n")
  expr_data <- GetAssayData(chicken_orth_mice_seurat, layer = "counts")
} else if ("data" %in% available_layers) {
  cat("Using data layer for gene expression analysis\n")
  expr_data <- GetAssayData(chicken_orth_mice_seurat, layer = "data")
} else {
  cat("⚠ No suitable layer found for expression analysis\n")
  expr_data <- NULL
}

if (!is.null(expr_data)) {
  # 计算基因表达统计
  gene_means <- rowMeans(expr_data)
  gene_expressed <- rowSums(expr_data > 0)
  gene_expression_rate <- gene_expressed / ncol(expr_data)
  
  cat("✓ Gene expression statistics:\n")
  cat("  Mean expression range:", round(range(gene_means), 4), "\n")
  cat("  Mean expression average:", round(mean(gene_means), 4), "\n")
  cat("  Expression rate range:", round(range(gene_expression_rate), 4), "\n")
  cat("  Expression rate average:", round(mean(gene_expression_rate), 4), "\n")
  
  # 检查NeuronChat相关基因
  neuronchat_genes <- c("Nrxn1", "Nrxn3", "Nlgn1", "Nlgn3", "Grm1", "Grm2", "Grm3", "Grm4", "Grm5")
  existing_neuronchat_genes <- intersect(rownames(chicken_orth_mice_seurat), neuronchat_genes)
  cat("  NeuronChat genes present:", length(existing_neuronchat_genes), "out of", length(neuronchat_genes), "\n")
  
  if (length(existing_neuronchat_genes) > 0) {
    cat("  Sample NeuronChat gene expression:\n")
    for (gene in existing_neuronchat_genes[1:min(3, length(existing_neuronchat_genes))]) {
      mean_expr <- gene_means[gene]
      expr_rate <- gene_expression_rate[gene]
      cat("    ", gene, ": mean =", round(mean_expr, 4), ", rate =", round(expr_rate, 4), "\n")
    }
  }
} else {
  cat("⚠ Cannot perform gene expression analysis - no suitable data layer\n")
}

# 9. 保存结果
cat("\n9. Saving results...\n")

# 保存Seurat对象
saveRDS(chicken_orth_mice_seurat, "chicken_orth_mice_new.rds")
cat("✓ Seurat object saved: chicken_orth_mice_seurat.rds\n")

# 保存转换信息
conversion_info <- data.frame(
  metric = c("original_genes", "converted_genes", "cells", "conversion_method", "conversion_date"),
  value = c(nrow(chicken_seurat), nrow(converted_expr_matrix), ncol(chicken_orth_mice_seurat), 
            "weighted_average", as.character(Sys.time())),
  stringsAsFactors = FALSE
)

write.csv(conversion_info, "chicken_orth_mice_conversion_info.csv", row.names = FALSE)
cat("✓ Conversion info saved: chicken_orth_mice_conversion_info.csv\n")

# 10. 输出最终统计信息
cat("\n=== Final Statistics ===\n")
cat("Analysis completed at:", Sys.time(), "\n")
cat("Original genes:", nrow(chicken_seurat), "\n")
cat("Converted genes:", nrow(converted_expr_matrix), "\n")
cat("Cells:", ncol(chicken_orth_mice_seurat), "\n")
cat("Cell types:", length(unique(cell_types)), "\n")
cat("Conversion rate:", round(nrow(converted_expr_matrix)/nrow(chicken_seurat)*100, 2), "%\n")

cat("\nFiles created:\n")
cat("- chicken_orth_mice_seurat.rds: Seurat object with converted genes\n")
cat("- chicken_orth_mice_conversion_info.csv: Conversion information\n")

cat("\n=== Seurat Object Creation Complete ===\n")
cat("The Seurat object is ready for NeuronChat analysis.\n")
cat("Gene conversion used weighted average method for duplicate genes.\n")
