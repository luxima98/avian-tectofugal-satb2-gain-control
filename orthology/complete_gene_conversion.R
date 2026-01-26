
# 完整的鸡基因ID转换为鼠基因ID的脚本（使用加权平均处理重复基因）

library(Seurat)
library(dplyr)

# 检查并安装必要的包
required_packages <- c("Seurat", "dplyr")
optional_packages <- c("biomaRt")

# 检查必需包
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package ", pkg, " is not installed. Please install it first.")
  }
}

# 检查可选包
biomaRt_available <- requireNamespace("biomaRt", quietly = TRUE)
if (biomaRt_available) {
  cat("✓ biomaRt package available\n")
} else {
  cat("⚠ biomaRt package not available, will use Ensembl IDs as gene symbols\n")
}

cat("=== Complete Chicken to Mouse Gene ID Conversion (Weighted Average) ===\n")
cat("Starting analysis at:", Sys.time(), "\n")

# 1. 加载鸡单细胞数据
cat("\n1. Loading chicken single-cell data...\n")
chicken_seurat <- readRDS("chicken_mn_Ex_Inh_seurat.rds")
cat("✓ Chicken data loaded:", ncol(chicken_seurat), "cells,", nrow(chicken_seurat), "genes\n")

# 2. 检查基因ID格式
cat("\n2. Checking gene ID formats...\n")
chicken_genes <- rownames(chicken_seurat)
cat("Total genes in chicken data:", length(chicken_genes), "\n")
cat("Sample chicken genes:", paste(head(chicken_genes, 5), collapse = ", "), "\n")

# 3. 检查并使用现成的同源关系数据
cat("\n3. Loading existing homolog data...\n")

# 优先使用现成的基因转换映射文件
gene_conversion_file <- "gene_conversion_mapping.csv"
homologs_file <- "chicken_homologs.csv"

# 检查是否有现成的基因转换映射
if (file.exists(gene_conversion_file)) {
  cat("✓ Found existing gene conversion mapping:", gene_conversion_file, "\n")
  chicken_to_mouse <- read.csv(gene_conversion_file, stringsAsFactors = FALSE)
  cat("✓ Loaded", nrow(chicken_to_mouse), "gene pairs from existing mapping\n")
  
  # 检查列名并标准化
  if ("original_chicken_ensembl" %in% colnames(chicken_to_mouse)) {
    colnames(chicken_to_mouse)[colnames(chicken_to_mouse) == "original_chicken_ensembl"] <- "chicken_ensembl"
  }
  
  # 确保mouse_symbol列存在
  if (!"mouse_symbol" %in% colnames(chicken_to_mouse)) {
    cat("⚠ Warning: mouse_symbol column not found in gene conversion mapping\n")
  } else {
    cat("✓ Found mouse_symbol column in gene conversion mapping\n")
  }
  
} else if (file.exists(homologs_file)) {
  cat("✓ Found existing homolog data:", homologs_file, "\n")
  homologs_data <- read.csv(homologs_file, stringsAsFactors = FALSE)
  
  # 筛选有效的同源关系
  valid_homologs <- homologs_data[homologs_data$found_homolog == TRUE, ]
  
  chicken_to_mouse <- data.frame(
    chicken_ensembl = valid_homologs$chicken_ensembl,
    mouse_ensembl = valid_homologs$mouse_ensembl,
    stringsAsFactors = FALSE
  )
  cat("✓ Loaded", nrow(chicken_to_mouse), "gene pairs from homolog data\n")
  
} else {
  # 如果都没有，则使用OrthoFinder数据
  cat("⚠ No existing mapping found, using OrthoFinder data...\n")
  
  orthofinder_file <- "G.gallus__v__M.musculus.tsv"
  if (!file.exists(orthofinder_file)) {
    stop("No gene mapping files found. Please ensure at least one of the following exists:\n",
         "- gene_conversion_mapping.csv\n",
         "- chicken_homologs.csv\n", 
         "- G.gallus__v__M.musculus.tsv")
  }
  
  # 读取OrthoFinder结果
  homolog_data <- read.delim(orthofinder_file, stringsAsFactors = FALSE)
  cat("✓ OrthoFinder data loaded:", nrow(homolog_data), "orthogroups\n")
  
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
    gene_pairs <- expand.grid(
      chicken_ensembl = chicken_genes_list,
      mouse_ensembl = mouse_genes_list,
      stringsAsFactors = FALSE
    )
    
    # 批量添加到数据框
    chicken_to_mouse <- rbind(chicken_to_mouse, gene_pairs)
  }
  
  cat("✓ OrthoFinder mapping created:", nrow(chicken_to_mouse), "gene pairs\n")
}

# 5. 检查鸡基因在单细胞数据中的存在情况
cat("\n5. Checking gene presence in single-cell data...\n")

# 清理单细胞数据中的基因ID（移除版本号，如果有的话）
chicken_genes_clean <- gsub("\\.[0-9]+$", "", chicken_genes)

# 找到在单细胞数据中存在的鸡基因
available_chicken_genes <- intersect(chicken_genes_clean, chicken_to_mouse$chicken_ensembl)
cat("✓ Chicken genes with homologs in single-cell data:", length(available_chicken_genes), "\n")

if (length(available_chicken_genes) == 0) {
  stop("✗ No chicken genes found in single-cell data that have homologs")
}

# 6. 创建基因ID转换映射
cat("\n6. Creating gene ID conversion mapping...\n")

# 筛选出在单细胞数据中存在的同源关系
valid_homologs <- chicken_to_mouse[chicken_to_mouse$chicken_ensembl %in% available_chicken_genes, ]
cat("✓ Valid homolog pairs:", nrow(valid_homologs), "\n")

# 调试信息：检查valid_homologs的列名
cat("✓ Valid homologs columns:", paste(colnames(valid_homologs), collapse = ", "), "\n")
cat("✓ Has mouse_symbol column:", "mouse_symbol" %in% colnames(valid_homologs), "\n")

# 数据验证
if (nrow(valid_homologs) == 0) {
  stop("✗ No valid homolog pairs found after filtering. Please check gene ID formats and mapping data.")
}

# 创建原始基因ID到清理后基因ID的映射
original_to_clean <- setNames(chicken_genes_clean, chicken_genes)

# 为valid_homologs添加原始基因ID列
valid_homologs$original_chicken_ensembl <- names(original_to_clean)[match(valid_homologs$chicken_ensembl, original_to_clean)]

# 检查是否有重复的鸡基因
duplicated_chicken <- sum(duplicated(valid_homologs$chicken_ensembl))
if (duplicated_chicken > 0) {
  cat("⚠ Warning:", duplicated_chicken, "chicken genes have multiple mouse homologs\n")
  
  # 使用最佳匹配策略：基于同源组大小和基因ID相似性选择最佳匹配
  cat("Using best-match strategy for duplicate genes...\n")
  
  # 计算每个鸡基因的匹配得分
  valid_homologs$match_score <- 0
  
  # 得分规则：
  # 1. 同源组中基因数量越少，得分越高（更特异）
  # 2. 基因ID相似性越高，得分越高
  for (i in 1:nrow(valid_homologs)) {
    chicken_gene <- valid_homologs$chicken_ensembl[i]
    mouse_gene <- valid_homologs$mouse_ensembl[i]
    
    # 找到这个鸡基因所在的所有同源组
    chicken_groups <- which(sapply(1:nrow(homolog_data), function(j) {
      chicken_genes_in_group <- strsplit(homolog_data$G.gallus[j], ", ")[[1]]
      chicken_genes_in_group <- gsub("\\.[0-9]+$", "", chicken_genes_in_group)
      chicken_gene %in% chicken_genes_in_group
    }))
    
    # 计算同源组特异性得分（基因数量越少得分越高）
    group_specificity <- 1 / length(chicken_groups)
    
    # 计算基因ID相似性得分（简单的字符串相似性）
    gene_similarity <- 1 - adist(chicken_gene, mouse_gene) / max(nchar(chicken_gene), nchar(mouse_gene))
    
    # 综合得分
    valid_homologs$match_score[i] <- group_specificity * 0.6 + gene_similarity * 0.4
  }
  
  # 选择每个鸡基因的最佳匹配
  valid_homologs <- valid_homologs %>%
    group_by(chicken_ensembl) %>%
    slice_max(match_score, n = 1) %>%
    ungroup()
  
  cat("✓ After best-match selection:", nrow(valid_homologs), "pairs\n")
  cat("✓ Match score range:", round(range(valid_homologs$match_score), 3), "\n")
}

# 7. 转换鼠Ensembl ID为基因符号
cat("\n7. Converting mouse Ensembl IDs to gene symbols...\n")

# 优先使用现成的鼠基因信息文件
mouse_gene_info_file <- "mouse_gene_info.csv"
mouse_gene_info <- NULL

if (file.exists(mouse_gene_info_file)) {
  cat("✓ Found existing mouse gene info:", mouse_gene_info_file, "\n")
  mouse_gene_info <- read.csv(mouse_gene_info_file, stringsAsFactors = FALSE)
  
  # 检查列名并标准化
  if ("ensembl_gene_id" %in% colnames(mouse_gene_info)) {
    colnames(mouse_gene_info)[colnames(mouse_gene_info) == "ensembl_gene_id"] <- "ensembl_id"
  }
  if ("external_gene_name" %in% colnames(mouse_gene_info)) {
    colnames(mouse_gene_info)[colnames(mouse_gene_info) == "external_gene_name"] <- "gene_symbol"
  }
  
  cat("✓ Loaded", nrow(mouse_gene_info), "mouse gene symbols from existing file\n")
  
} else {
  # 如果没有现成文件，尝试使用biomaRt
  cat("⚠ No existing mouse gene info found, trying biomaRt...\n")
  
  if (biomaRt_available) {
    cat("Using biomaRt to get mouse gene symbols...\n")
    tryCatch({
      # 使用更稳定的Ensembl镜像
      mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      mouse_genes <- unique(valid_homologs$mouse_ensembl)
      
      # 分批查询基因符号（避免一次性查询太多）
      batch_size <- 1000
      gene_info_list <- list()
      
      for (i in seq(1, length(mouse_genes), batch_size)) {
        end_idx <- min(i + batch_size - 1, length(mouse_genes))
        batch_genes <- mouse_genes[i:end_idx]
        
        cat("Processing batch", ceiling(i/batch_size), "of", ceiling(length(mouse_genes)/batch_size), "\n")
        
        batch_info <- biomaRt::getBM(
          attributes = c("ensembl_gene_id", "external_gene_name"),
          filters = "ensembl_gene_id",
          values = batch_genes,
          mart = mart
        )
        
        gene_info_list[[length(gene_info_list) + 1]] <- batch_info
      }
      
      # 合并所有批次的结果
      gene_info <- do.call(rbind, gene_info_list)
      
      mouse_gene_info <- gene_info
      colnames(mouse_gene_info) <- c("ensembl_id", "gene_symbol")
      cat("✓ Retrieved", nrow(mouse_gene_info), "gene symbols from Ensembl\n")
    }, error = function(e) {
      cat("✗ Failed to retrieve gene symbols from Ensembl:", e$message, "\n")
      cat("This might be due to network issues or Ensembl server problems\n")
      mouse_gene_info <- NULL
    })
  }
  
  if (is.null(mouse_gene_info)) {
    cat("⚠ Cannot get gene symbols, will use Ensembl IDs as gene names\n")
    # 创建简单的映射，使用Ensembl ID作为基因符号
    mouse_gene_info <- data.frame(
      ensembl_id = unique(valid_homologs$mouse_ensembl),
      gene_symbol = unique(valid_homologs$mouse_ensembl),
      stringsAsFactors = FALSE
    )
  }
}

# 8. 创建最终的基因转换映射
cat("\n8. Creating final gene conversion mapping...\n")

# 检查是否已经有现成的mouse_symbol列
if ("mouse_symbol" %in% colnames(valid_homologs)) {
  cat("✓ Found existing mouse_symbol column in mapping data\n")
  # 直接使用现成的mouse_symbol列
  converted_genes <- setNames(valid_homologs$mouse_symbol, valid_homologs$chicken_ensembl)
  # 移除NA值
  converted_genes <- converted_genes[!is.na(converted_genes)]
  cat("✓ Using existing gene symbols for", length(converted_genes), "genes\n")
} else {
  # 如果没有现成的mouse_symbol列，则从mouse_gene_info转换
  cat("Converting mouse Ensembl IDs to gene symbols...\n")
  mouse_ensembl_to_symbol <- setNames(mouse_gene_info$gene_symbol, mouse_gene_info$ensembl_id)
  # 创建最终的鸡基因到鼠基因符号的映射
  converted_genes <- mouse_ensembl_to_symbol[valid_homologs$mouse_ensembl]
  names(converted_genes) <- valid_homologs$chicken_ensembl
  # 移除NA值
  converted_genes <- converted_genes[!is.na(converted_genes)]
  cat("✓ Gene conversion mapping created for", length(converted_genes), "genes\n")
}
cat("Sample conversions:\n")
sample_conversions <- head(converted_genes, 5)
for (i in 1:length(sample_conversions)) {
  cat("  ", names(sample_conversions)[i], "->", sample_conversions[i], "\n")
}

# 9. 创建转换后的表达矩阵
cat("\n9. Creating converted expression matrix...\n")

# 修改获取矩阵的部分，使用counts slot
if (packageVersion("SeuratObject") >= "5.0.0") {
  original_expr_matrix <- GetAssayData(chicken_seurat, layer = "counts")
} else {
  original_expr_matrix <- GetAssayData(chicken_seurat, slot = "counts")
}
cat("✓ Using raw counts matrix for conversion\n")

# 筛选出有同源关系的基因（使用原始基因ID）
genes_to_keep <- valid_homologs$original_chicken_ensembl
filtered_expr_matrix <- original_expr_matrix[genes_to_keep, ]

# 重命名基因为鼠基因符号
rownames(filtered_expr_matrix) <- converted_genes[valid_homologs$chicken_ensembl]

cat("✓ Converted expression matrix created\n")
cat("  Original genes:", nrow(original_expr_matrix), "\n")
cat("  Converted genes:", nrow(filtered_expr_matrix), "\n")
cat("  Cells:", ncol(filtered_expr_matrix), "\n")

# 10. 检查转换结果并处理重复基因
cat("\n10. Checking conversion results and handling duplicates...\n")

# 检查是否有重复的鼠基因符号
duplicated_mouse <- sum(duplicated(rownames(filtered_expr_matrix)))
if (duplicated_mouse > 0) {
  cat("⚠ Warning:", duplicated_mouse, "mouse gene symbols are duplicated\n")
  cat("Sample duplicated genes:", paste(head(rownames(filtered_expr_matrix)[duplicated(rownames(filtered_expr_matrix))], 3), collapse = ", "), "\n")
  
  # 处理重复：使用加权平均
  cat("Resolving duplicates using weighted average...\n")
  
  # 获取所有唯一的鼠基因符号，过滤掉NA
  unique_genes <- unique(rownames(filtered_expr_matrix))
  # 过滤掉NA基因名
  valid_genes <- unique_genes[!is.na(unique_genes) & unique_genes != "NA"]
  cat("Processing", length(valid_genes), "unique mouse genes (filtered out", length(unique_genes) - length(valid_genes), "NA genes)...\n")
  
  # 创建新的表达矩阵
  final_expr_matrix <- matrix(0, nrow = length(valid_genes), ncol = ncol(filtered_expr_matrix))
  rownames(final_expr_matrix) <- valid_genes
  colnames(final_expr_matrix) <- colnames(filtered_expr_matrix)
  
  # 处理每个唯一的鼠基因
  processed_count <- 0
  for (gene in valid_genes) {
    processed_count <- processed_count + 1
    if (processed_count %% 1000 == 0) {
      cat("Processed", processed_count, "genes...\n")
    }
    
    gene_indices <- which(rownames(filtered_expr_matrix) == gene)
    
    if (length(gene_indices) == 0) {
      cat("  Gene", gene, "not found in matrix\n")
      next
    } else if (length(gene_indices) == 1) {
      # 只有一个鸡基因映射到这个鼠基因
      final_expr_matrix[gene, ] <- filtered_expr_matrix[gene_indices, ]
    } else {
      # 多个鸡基因映射到这个鼠基因，使用加权平均
      cat("  Gene", gene, "has", length(gene_indices), "chicken homologs\n")
      
      # 计算每个鸡基因的表达量权重（基于平均表达量）
      gene_expr_means <- rowMeans(filtered_expr_matrix[gene_indices, , drop = FALSE])
      
      # 避免除零错误
      if (sum(gene_expr_means) > 0) {
        weights <- gene_expr_means / sum(gene_expr_means)
      } else {
        weights <- rep(1/length(gene_indices), length(gene_indices))
      }
      
      # 计算加权平均表达量
      # 确保weights和gene_indices长度匹配
      if (length(weights) != length(gene_indices)) {
        cat("    Warning: weights length (", length(weights), ") != gene_indices length (", length(gene_indices), ")\n")
        weights <- rep(1/length(gene_indices), length(gene_indices))
      }
      
      # 获取基因表达矩阵子集
      gene_expr_subset <- filtered_expr_matrix[gene_indices, , drop = FALSE]
      
      # 确保weights是列向量，与基因表达矩阵的行数匹配
      if (length(weights) != nrow(gene_expr_subset)) {
        cat("    Error: weights length (", length(weights), ") != gene subset rows (", nrow(gene_expr_subset), ")\n")
        next
      }
      
      # 计算加权平均表达量
      weighted_expr <- colSums(gene_expr_subset * weights)
      final_expr_matrix[gene, ] <- weighted_expr
      
      # 输出权重信息（仅对前几个基因）
      if (processed_count <= 5) {
        cat("    Weights:", paste(round(weights, 3), collapse = ", "), "\n")
        cat("    Original means:", paste(round(gene_expr_means, 3), collapse = ", "), "\n")
        cat("    Final mean:", round(mean(weighted_expr), 3), "\n")
      }
    }
  }
  
  # 更新表达矩阵
  filtered_expr_matrix <- final_expr_matrix
  cat("✓ After resolving duplicates with weighted average:", nrow(filtered_expr_matrix), "genes\n")
  
  # 统计重复处理情况
  duplicate_stats <- table(table(rownames(filtered_expr_matrix)))
  cat("Duplicate resolution statistics:\n")
  for (i in 1:length(duplicate_stats)) {
    cat("  ", names(duplicate_stats)[i], "chicken genes ->", duplicate_stats[i], "mouse genes\n")
  }
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

# 如果有匹配得分，也保存
if ("match_score" %in% colnames(valid_homologs)) {
  gene_conversion_df$match_score <- valid_homologs$match_score
}
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
cat("- chicken_to_mouse_converted_matrix.rds: Converted expression matrix (weighted average)\n")
cat("- gene_conversion_mapping.csv: Gene ID conversion mapping\n")
cat("- conversion_statistics.csv: Conversion statistics\n")

cat("\nMethod used for duplicate genes: Weighted Average\n")
cat("- Multiple chicken genes mapping to same mouse gene are averaged\n")
cat("- Weights are based on expression levels of chicken genes\n")
cat("- This preserves more expression information than max selection\n")

# 12. 检查neurongene_lig_mice基因覆盖情况
cat("\n12. Checking neurongene_lig_mice gene coverage...\n")

# 检查neurongene_lig_mice文件是否存在
neurongene_file <- "neurongene_lig_mice.csv"
if (file.exists(neurongene_file)) {
  cat("✓ Found neurongene_lig_mice file:", neurongene_file, "\n")
  
  # 读取neurongene_lig_mice基因列表
  neurongene_lig_mice <- read.csv(neurongene_file, stringsAsFactors = FALSE)
  neurongene_genes <- neurongene_lig_mice$gene_symbol
  
  cat("✓ Loaded", length(neurongene_genes), "neurongene_lig_mice genes\n")
  cat("Sample genes:", paste(head(neurongene_genes, 5), collapse = ", "), "\n")
  
  # 检查这些基因在转换后的矩阵中是否存在（使用大小写模糊匹配）
  converted_genes <- rownames(filtered_expr_matrix)
  
  # 创建大小写不敏感的匹配
  neurongene_upper <- toupper(neurongene_genes)
  converted_upper <- toupper(converted_genes)
  
  # 使用大小写不敏感匹配
  neurongene_found_indices <- match(neurongene_upper, converted_upper)
  neurongene_found <- neurongene_genes[!is.na(neurongene_found_indices)]
  neurongene_missing <- neurongene_genes[is.na(neurongene_found_indices)]
  
  # 获取实际匹配的基因名称（保持原始大小写）
  if (length(neurongene_found) > 0) {
    actual_matched_genes <- converted_genes[neurongene_found_indices[!is.na(neurongene_found_indices)]]
    cat("✓ Case-insensitive matching enabled\n")
    cat("✓ Sample actual matches:", paste(head(actual_matched_genes, 5), collapse = ", "), "\n")
  }
  
  cat("\n=== Neurongene Coverage Analysis ===\n")
  cat("Total neurongene_lig_mice genes:", length(neurongene_genes), "\n")
  cat("Found in converted matrix:", length(neurongene_found), "\n")
  cat("Missing from converted matrix:", length(neurongene_missing), "\n")
  cat("Coverage rate:", round(length(neurongene_found)/length(neurongene_genes)*100, 2), "%\n")
  
  if (length(neurongene_found) > 0) {
    cat("\n✓ Found genes in converted matrix:\n")
    cat(paste(head(neurongene_found, 10), collapse = ", "), "\n")
    if (length(neurongene_found) > 10) {
      cat("... and", length(neurongene_found) - 10, "more\n")
    }
  }
  
  if (length(neurongene_missing) > 0) {
    cat("\n⚠ Missing genes from converted matrix:\n")
    cat(paste(head(neurongene_missing, 10), collapse = ", "), "\n")
    if (length(neurongene_missing) > 10) {
      cat("... and", length(neurongene_missing) - 10, "more\n")
    }
    
    # 分析缺失基因的原因（使用大小写不敏感匹配）
    cat("\nAnalyzing missing genes...\n")
    
    # 检查缺失基因是否在原始鸡数据中（大小写不敏感）
    chicken_genes_upper <- toupper(rownames(chicken_seurat))
    missing_upper <- toupper(neurongene_missing)
    
    missing_in_chicken_indices <- match(missing_upper, chicken_genes_upper)
    missing_in_chicken <- neurongene_missing[is.na(missing_in_chicken_indices)]
    missing_in_mapping <- neurongene_missing[!is.na(missing_in_chicken_indices)]
    
    cat("Missing genes not in original chicken data:", length(missing_in_chicken), "\n")
    cat("Missing genes in chicken data but not in mapping:", length(missing_in_mapping), "\n")
    
    if (length(missing_in_chicken) > 0) {
      cat("Genes not in original chicken data:\n")
      cat(paste(head(missing_in_chicken, 5), collapse = ", "), "\n")
    }
    
    if (length(missing_in_mapping) > 0) {
      cat("Genes in chicken data but not in gene mapping:\n")
      cat(paste(head(missing_in_mapping, 5), collapse = ", "), "\n")
    }
  }
  
  # 保存覆盖情况报告（使用大小写不敏感匹配）
  # 重新计算匹配状态（大小写不敏感）
  neurongene_found_status <- !is.na(match(toupper(neurongene_genes), toupper(converted_genes)))
  neurongene_in_chicken_status <- !is.na(match(toupper(neurongene_genes), toupper(rownames(chicken_seurat))))
  
  # 获取匹配的mouse_ensembl_id（大小写不敏感）
  matched_indices <- match(toupper(neurongene_genes), toupper(gene_conversion_df$mouse_symbol))
  matched_mouse_ensembl <- ifelse(!is.na(matched_indices), 
                                  gene_conversion_df$mouse_ensembl[matched_indices], 
                                  NA)
  
  coverage_report <- data.frame(
    gene_symbol = neurongene_genes,
    found_in_converted_matrix = neurongene_found_status,
    found_in_original_chicken = neurongene_in_chicken_status,
    mouse_ensembl_id = matched_mouse_ensembl,
    stringsAsFactors = FALSE
  )
  
  write.csv(coverage_report, "neurongene_lig_mice_coverage_report.csv", row.names = FALSE)
  cat("\n✓ Coverage report saved: neurongene_lig_mice_coverage_report.csv\n")
  
} else {
  cat("⚠ neurongene_lig_mice.csv file not found, skipping coverage analysis\n")
}

# 13. 创建新的Seurat对象并进行SCT标准化
cat("\n13. Creating new Seurat object and performing SCT normalization...\n")

# 创建新的Seurat对象，使用转换后的表达矩阵
chicken_orth_mice_seurat <- CreateSeuratObject(
  counts = filtered_expr_matrix,
  project = "chicken_orth_mice",
  min.cells = 3,
  min.features = 200
)

cat("✓ New Seurat object created with converted mouse genes\n")
cat("  Genes:", nrow(chicken_orth_mice_seurat), "\n")
cat("  Cells:", ncol(chicken_orth_mice_seurat), "\n")

# 添加细胞元数据（从原始对象复制）
if ("anno_level_3" %in% colnames(chicken_seurat@meta.data)) {
  chicken_orth_mice_seurat@meta.data$anno_level_3 <- chicken_seurat@meta.data$anno_level_3
  cat("✓ Added anno_level_3 metadata\n")
}

# 执行SCT标准化
cat("Performing SCTransform normalization...\n")
chicken_orth_mice_seurat <- SCTransform(
  chicken_orth_mice_seurat,
  method = "glmGamPoi",
  vars.to.regress = NULL,
  verbose = TRUE
)

cat("✓ SCTransform normalization completed\n")

# 保存SCT标准化后的Seurat对象
saveRDS(chicken_orth_mice_seurat, "chicken_orth_mice_seurat_2.rds")
cat("✓ SCT-normalized Seurat object saved: chicken_orth_mice_seurat.rds\n")

# 保存SCT转换统计信息
sct_stats <- data.frame(
  metric = c("Original_genes", "Converted_genes", "Cells", "SCT_completed", "Normalization_method"),
  value = c(length(chicken_genes), nrow(filtered_expr_matrix), ncol(filtered_expr_matrix), "Yes", "SCTransform")
)
write.csv(sct_stats, "sct_normalization_statistics.csv", row.names = FALSE)
cat("✓ SCT normalization statistics saved: sct_normalization_statistics.csv\n")

cat("\n=== SCT Normalization Complete ===\n")
cat("✓ Gene conversion completed using raw counts\n")
cat("✓ SCTransform normalization applied\n")
cat("✓ Final Seurat object ready for downstream analysis\n")

cat("\n=== Conversion Complete ===\n")
cat("Analysis completed at:", Sys.time(), "\n")

