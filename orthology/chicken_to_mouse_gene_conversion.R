# Chicken to Mouse Gene ID Conversion Script
# 鸡基因ID转换为鼠基因ID的脚本

library(Seurat)
library(dplyr)

cat("=== Chicken to Mouse Gene ID Conversion ===\n")
cat("Starting analysis at:", Sys.time(), "\n")

# 1. 加载鸡单细胞数据
cat("\n1. Loading chicken single-cell data...\n")
chicken_seurat <- readRDS("chicken_mn_Ex_Inh_seurat.rds")
cat("✓ Chicken data loaded:", ncol(chicken_seurat), "cells,", nrow(chicken_seurat), "genes\n")

# 2. 检查基因ID格式
cat("\n2. Checking gene ID formats...\n")
chicken_genes <- rownames(chicken_seurat)
cat("Total genes in chicken data:", length(chicken_genes), "\n")
cat("Sample chicken genes (first 10):", paste(head(chicken_genes, 10), collapse = ", "), "\n")

# 检查基因ID格式
if (all(grepl("^ENSGALG", chicken_genes))) {
  cat("✓ Gene IDs are chicken Ensembl IDs (ENSGALG format)\n")
} else if (all(grepl("^ENSG", chicken_genes))) {
  cat("✓ Gene IDs are Ensembl IDs (ENSG format)\n")
} else {
  cat("? Gene IDs appear to be gene symbols or mixed format\n")
  cat("Sample formats:", paste(unique(substr(chicken_genes[1:min(20, length(chicken_genes))], 1, 8)), collapse = ", "), "\n")
}

# 3. 加载同源关系数据
cat("\n3. Loading homolog relationship data...\n")

# 检查OrthoFinder同源文件
orthofinder_file <- "G.gallus__v__M.musculus.tsv"
if (!file.exists(orthofinder_file)) {
  stop("✗ OrthoFinder homolog file not found: ", orthofinder_file)
}

cat("✓ Found OrthoFinder homolog file:", orthofinder_file, "\n")

# 加载OrthoFinder同源关系数据
homolog_data <- read.delim(orthofinder_file, stringsAsFactors = FALSE)
cat("✓ OrthoFinder data loaded:", nrow(homolog_data), "orthogroups\n")

# 检查文件结构
cat("OrthoFinder file columns:", paste(colnames(homolog_data), collapse = ", "), "\n")

# 4. 处理OrthoFinder同源关系
cat("\n4. Processing OrthoFinder homolog relationships...\n")

# 创建鸡Ensembl ID到鼠Ensembl ID的映射
chicken_to_mouse <- data.frame(
  chicken_ensembl = character(),
  mouse_ensembl = character(),
  stringsAsFactors = FALSE
)

# 处理每个同源组
for (i in 1:nrow(homolog_data)) {
  if (i %% 1000 == 0) {
    cat("Processed", i, "orthogroups...\n")
  }
  
  chicken_genes_list <- strsplit(homolog_data$G.gallus[i], ", ")[[1]]
  mouse_genes_list <- strsplit(homolog_data$M.musculus[i], ", ")[[1]]
  
  # 清理基因ID（移除版本号）
  chicken_genes_list <- gsub("\\.[0-9]+$", "", chicken_genes_list)
  mouse_genes_list <- gsub("\\.[0-9]+$", "", mouse_genes_list)
  
  # 创建所有可能的鸡-鼠基因对
  for (chicken_gene in chicken_genes_list) {
    for (mouse_gene in mouse_genes_list) {
      chicken_to_mouse <- rbind(chicken_to_mouse, data.frame(
        chicken_ensembl = chicken_gene,
        mouse_ensembl = mouse_gene,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("✓ OrthoFinder mapping created:", nrow(chicken_to_mouse), "gene pairs\n")

# 检查是否有重复的鸡基因
duplicated_chicken <- sum(duplicated(chicken_to_mouse$chicken_ensembl))
if (duplicated_chicken > 0) {
  cat("⚠ Warning:", duplicated_chicken, "chicken genes have multiple mouse homologs\n")
  cat("Sample duplicated chicken genes:", paste(head(chicken_to_mouse$chicken_ensembl[duplicated(chicken_to_mouse$chicken_ensembl)], 3), collapse = ", "), "\n")
  
  # 处理重复：保留第一个匹配（或可以基于其他标准选择最佳匹配）
  cat("Resolving duplicates by keeping first match...\n")
  chicken_to_mouse <- chicken_to_mouse[!duplicated(chicken_to_mouse$chicken_ensembl), ]
  cat("✓ After removing duplicates:", nrow(chicken_to_mouse), "pairs\n")
}

# 5. 检查鸡基因在单细胞数据中的存在情况
cat("\n5. Checking gene presence in single-cell data...\n")

# 清理单细胞数据中的基因ID（移除版本号，如果有的话）
chicken_genes_clean <- gsub("\\.[0-9]+$", "", chicken_genes)

# 找到在单细胞数据中存在的鸡基因
available_chicken_genes <- intersect(chicken_genes_clean, chicken_to_mouse$chicken_ensembl)
cat("✓ Chicken genes with homologs in single-cell data:", length(available_chicken_genes), "\n")

if (length(available_chicken_genes) == 0) {
  cat("✗ No chicken genes found in single-cell data that have homologs\n")
  cat("This might indicate a gene ID format mismatch\n")
  
  # 尝试部分匹配
  cat("\nTrying partial matching...\n")
  chicken_ensembl_prefix <- unique(substr(chicken_to_mouse$chicken_ensembl, 1, 15))
  chicken_data_prefix <- unique(substr(chicken_genes, 1, 15))
  
  cat("Sample chicken Ensembl IDs from OrthoFinder:", paste(head(chicken_ensembl_prefix, 3), collapse = ", "), "\n")
  cat("Sample chicken Ensembl IDs from data:", paste(head(chicken_data_prefix, 3), collapse = ", "), "\n")
  
  # 检查是否有部分匹配
  partial_matches <- intersect(chicken_ensembl_prefix, chicken_data_prefix)
  if (length(partial_matches) > 0) {
    cat("✓ Found partial matches:", length(partial_matches), "\n")
    cat("This suggests gene ID format differences\n")
  }
  
  stop("Cannot proceed without matching gene IDs")
}

# 6. 创建基因ID转换映射
cat("\n6. Creating gene ID conversion mapping...\n")

# 筛选出在单细胞数据中存在的同源关系
valid_homologs <- chicken_to_mouse[chicken_to_mouse$chicken_ensembl %in% available_chicken_genes, ]
cat("✓ Valid homolog pairs:", nrow(valid_homologs), "\n")

# 创建原始基因ID到清理后基因ID的映射
original_to_clean <- setNames(chicken_genes_clean, chicken_genes)

# 为valid_homologs添加原始基因ID列
valid_homologs$original_chicken_ensembl <- names(original_to_clean)[match(valid_homologs$chicken_ensembl, original_to_clean)]

# 检查是否有重复的鸡基因（在valid_homologs中应该已经处理过了）
duplicated_chicken_valid <- sum(duplicated(valid_homologs$chicken_ensembl))
if (duplicated_chicken_valid > 0) {
  cat("⚠ Warning:", duplicated_chicken_valid, "chicken genes still have multiple homologs\n")
  # 保留第一个匹配
  valid_homologs <- valid_homologs[!duplicated(valid_homologs$chicken_ensembl), ]
  cat("✓ After removing duplicates:", nrow(valid_homologs), "pairs\n")
}

# 7. 加载鼠基因信息并转换为基因符号
cat("\n7. Loading mouse gene information and converting to gene symbols...\n")

# 检查是否有鼠基因注释文件
mouse_gene_files <- c("mouse_gene_info.csv", "mouse_gene_annotation.csv", "Mus_musculus.GRCm38.genes.gtf")
mouse_gene_info <- NULL

# 首先尝试使用biomaRt获取完整的基因信息
if (requireNamespace("biomaRt", quietly = TRUE)) {
  cat("Using biomaRt to get mouse gene symbols...\n")
  tryCatch({
    mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    mouse_genes <- unique(valid_homologs$mouse_ensembl)
    
    # 批量查询基因符号
    gene_info <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = mouse_genes,
      mart = mart
    )
    
    mouse_gene_info <- gene_info
    colnames(mouse_gene_info) <- c("ensembl_id", "gene_symbol")
    cat("✓ Retrieved", nrow(mouse_gene_info), "gene symbols from Ensembl\n")
  }, error = function(e) {
    cat("✗ Failed to retrieve gene symbols from Ensembl:", e$message, "\n")
    mouse_gene_info <- NULL
  })
}

# 如果biomaRt失败，尝试本地文件
if (is.null(mouse_gene_info)) {
  for (file in mouse_gene_files) {
    if (file.exists(file)) {
      if (grepl("\\.csv$", file)) {
        mouse_gene_info <- read.csv(file, stringsAsFactors = FALSE)
        # 检查列名并调整
        if ("ensembl_gene_id" %in% colnames(mouse_gene_info)) {
          mouse_gene_info <- mouse_gene_info[, c("ensembl_gene_id", "gene_symbol")]
          colnames(mouse_gene_info) <- c("ensembl_id", "gene_symbol")
        }
      } else if (grepl("\\.gtf$", file)) {
        # 对于GTF文件，我们需要解析基因ID和基因符号
        cat("GTF file found, parsing gene information...\n")
        gtf_lines <- readLines(file)
        gene_lines <- gtf_lines[grepl("gene_id", gtf_lines)]
        
        # 提取基因ID和基因符号
        gene_ids <- gsub(".*gene_id \"([^\"]+)\".*", "\\1", gene_lines)
        gene_symbols <- gsub(".*gene_name \"([^\"]+)\".*", "\\1", gene_lines)
        
        mouse_gene_info <- data.frame(
          ensembl_id = gene_ids,
          gene_symbol = gene_symbols,
          stringsAsFactors = FALSE
        )
      }
      cat("✓ Found mouse gene info file:", file, "\n")
      break
    }
  }
}

if (is.null(mouse_gene_info)) {
  cat("⚠ Cannot get gene symbols, will use Ensembl IDs as gene names\n")
  # 创建简单的映射，使用Ensembl ID作为基因符号
  mouse_gene_info <- data.frame(
    ensembl_id = unique(valid_homologs$mouse_ensembl),
    gene_symbol = unique(valid_homologs$mouse_ensembl),
    stringsAsFactors = FALSE
  )
} else {
  cat("✓ Mouse gene info loaded:", nrow(mouse_gene_info), "entries\n")
  cat("Mouse gene info columns:", paste(colnames(mouse_gene_info), collapse = ", "), "\n")
}

# 8. 转换基因ID
cat("\n8. Converting gene IDs...\n")

# 将鼠Ensembl ID转换为基因符号
mouse_ensembl_to_symbol <- setNames(mouse_gene_info$gene_symbol, mouse_gene_info$ensembl_id)

# 创建最终的鸡基因到鼠基因符号的映射
converted_genes <- mouse_ensembl_to_symbol[valid_homologs$mouse_ensembl]
names(converted_genes) <- valid_homologs$chicken_ensembl

# 移除NA值
converted_genes <- converted_genes[!is.na(converted_genes)]

cat("✓ Gene conversion mapping created for", length(converted_genes), "genes\n")
cat("Sample conversions:\n")
sample_conversions <- head(converted_genes, 5)
for (i in 1:length(sample_conversions)) {
  cat("  ", names(sample_conversions)[i], "->", sample_conversions[i], "\n")
}

# 9. 创建转换后的表达矩阵
cat("\n9. Creating converted expression matrix...\n")

# 获取原始表达矩阵
original_expr_matrix <- GetAssayData(chicken_seurat, slot = "data")

# 筛选出有同源关系的基因（使用原始基因ID）
genes_to_keep <- valid_homologs$original_chicken_ensembl

# 检查基因数量并确保保持矩阵结构
if (length(genes_to_keep) == 1) {
  # 如果只有一个基因，需要特殊处理以保持矩阵结构
  cat("⚠ Warning: Only 1 gene found, creating single-row matrix\n")
  filtered_expr_matrix <- original_expr_matrix[genes_to_keep, , drop = FALSE]
} else {
  filtered_expr_matrix <- original_expr_matrix[genes_to_keep, ]
}

# 重命名基因为鼠基因符号
# 确保converted_genes的索引与valid_homologs$chicken_ensembl匹配
if (length(genes_to_keep) == 1) {
  # 单基因情况
  rownames(filtered_expr_matrix) <- converted_genes[1]
} else {
  # 多基因情况
  rownames(filtered_expr_matrix) <- converted_genes[valid_homologs$chicken_ensembl]
}

cat("✓ Converted expression matrix created\n")
cat("  Original genes:", nrow(original_expr_matrix), "\n")
cat("  Converted genes:", nrow(filtered_expr_matrix), "\n")
cat("  Cells:", ncol(filtered_expr_matrix), "\n")

# 10. 检查转换结果
cat("\n10. Checking conversion results...\n")

# 检查是否有重复的鼠基因符号
duplicated_mouse <- sum(duplicated(rownames(filtered_expr_matrix)))
if (duplicated_mouse > 0) {
  cat("⚠ Warning:", duplicated_mouse, "mouse gene symbols are duplicated\n")
  cat("Sample duplicated genes:", paste(head(rownames(filtered_expr_matrix)[duplicated(rownames(filtered_expr_matrix))], 3), collapse = ", "), "\n")
  
  # 处理重复：保留表达量最高的
  cat("Resolving duplicates by keeping highest expressed genes...\n")
  
  # 计算每个基因的平均表达量
  gene_means <- rowMeans(filtered_expr_matrix)
  
  # 为重复基因选择表达量最高的
  unique_genes <- unique(rownames(filtered_expr_matrix))
  final_genes <- c()
  
  for (gene in unique_genes) {
    gene_indices <- which(rownames(filtered_expr_matrix) == gene)
    if (length(gene_indices) > 1) {
      # 选择表达量最高的
      best_index <- gene_indices[which.max(gene_means[gene_indices])]
      final_genes <- c(final_genes, best_index)
    } else {
      final_genes <- c(final_genes, gene_indices)
    }
  }
  
  # 确保保持矩阵结构
  if (length(final_genes) == 1) {
    filtered_expr_matrix <- filtered_expr_matrix[final_genes, , drop = FALSE]
  } else {
    filtered_expr_matrix <- filtered_expr_matrix[final_genes, ]
  }
  cat("✓ After resolving duplicates:", nrow(filtered_expr_matrix), "genes\n")
}

# 11. 保存结果
cat("\n11. Saving results...\n")

# 保存转换后的表达矩阵
saveRDS(filtered_expr_matrix, "chicken_to_mouse_converted_matrix.rds")
cat("✓ Converted expression matrix saved: chicken_to_mouse_converted_matrix.rds\n")

# 保存基因转换映射
gene_conversion_df <- data.frame(
  original_chicken_ensembl = valid_homologs$original_chicken_ensembl,
  chicken_ensembl = valid_homologs$chicken_ensembl,
  mouse_ensembl = valid_homologs$mouse_ensembl,
  mouse_symbol = converted_genes[valid_homologs$chicken_ensembl],
  stringsAsFactors = FALSE
)
write.csv(gene_conversion_df, "gene_conversion_mapping.csv", row.names = FALSE)
cat("✓ Gene conversion mapping saved: gene_conversion_mapping.csv\n")

# 保存转换统计信息
conversion_stats <- data.frame(
  metric = c("Original_genes", "Genes_with_homologs", "Converted_genes", "Final_genes", "Cells"),
  value = c(length(chicken_genes), length(available_chicken_genes), length(converted_genes), nrow(filtered_expr_matrix), ncol(filtered_expr_matrix))
)
write.csv(conversion_stats, "conversion_statistics.csv", row.names = FALSE)
cat("✓ Conversion statistics saved: conversion_statistics.csv\n")

# 12. 输出最终统计信息
cat("\n=== Final Statistics ===\n")
cat("Original chicken genes:", length(chicken_genes), "\n")
cat("Chicken genes with homologs:", length(available_chicken_genes), "\n")
cat("Successfully converted genes:", nrow(filtered_expr_matrix), "\n")
cat("Cells:", ncol(filtered_expr_matrix), "\n")
cat("Conversion rate:", round(nrow(filtered_expr_matrix)/length(chicken_genes)*100, 2), "%\n")

cat("\nSample converted genes:\n")
sample_final_genes <- head(rownames(filtered_expr_matrix), 10)
for (i in 1:length(sample_final_genes)) {
  cat("  ", sample_final_genes[i], "\n")
}

cat("\nFiles created:\n")
cat("- chicken_to_mouse_converted_matrix.rds: Converted expression matrix\n")
cat("- gene_conversion_mapping.csv: Gene ID conversion mapping\n")
cat("- conversion_statistics.csv: Conversion statistics\n")

cat("\n=== Conversion Complete ===\n")
cat("Analysis completed at:", Sys.time(), "\n")
