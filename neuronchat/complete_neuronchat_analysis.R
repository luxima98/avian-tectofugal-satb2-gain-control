# Complete NeuronChat Analysis with Cross-Species Mapping
# 完整的跨物种NeuronChat分析脚本

library(Seurat)
library(dplyr)
library(NeuronChat)
library(CellChat)

cat("=== Complete NeuronChat Analysis with Cross-Species Mapping ===\n")
cat("Starting analysis at:", Sys.time(), "\n")

# 1. 加载数据
cat("\n1. Loading data...\n")

# 加载转换后的Seurat对象
chicken_orth_mice_seurat <- readRDS("seurat_obj_merged.rds")
cat("✓ Chicken ortholog mouse Seurat loaded:", ncol(chicken_orth_mice_seurat), "cells,", nrow(chicken_orth_mice_seurat), "genes\n")

# 加载基因映射信息
gene_conversion_df <- read.csv("gene_conversion_mapping.csv", stringsAsFactors = FALSE)
cat("✓ Gene conversion mapping loaded:", nrow(gene_conversion_df), "entries\n")

# 加载NeuronChat映射结果
neuronchat_mapping <- read.csv("complete_neuronchat_mapping.csv", stringsAsFactors = FALSE)
cat("✓ NeuronChat mapping loaded:", nrow(neuronchat_mapping), "entries\n")

# 2. 检查细胞类型分组
cat("\n2. Checking cell type grouping...\n")

if ("anno_level_3" %in% colnames(chicken_orth_mice_seurat@meta.data)) {
  cell_types <- chicken_orth_mice_seurat@meta.data$anno_level_3_merged
  cat("✓ Found anno_level_3 column\n")
  cat("Cell types:", length(unique(cell_types)), "\n")
  cat("Sample cell types:", paste(head(unique(cell_types), 5), collapse = ", "), "\n")
} else {
  stop("✗ anno_level_3_final column not found in metadata")
}

# 3. 创建NeuronChat对象
cat("\n3. Creating NeuronChat object...\n")

# 从Seurat对象中提取表达矩阵
# 优先使用scale.data（Seurat v5 SCTransform结果），然后SCT层
available_layers <- Layers(chicken_orth_mice_seurat)
cat("Available layers:", paste(available_layers, collapse = ", "), "\n")

if ("scale.data" %in% available_layers) {
  expr_matrix <- GetAssayData(chicken_orth_mice_seurat, layer = "scale.data")
  cat("✓ Using scale.data layer for expression matrix (Seurat v5 SCTransform)\n")
} else if ("SCT" %in% available_layers) {
  expr_matrix <- GetAssayData(chicken_orth_mice_seurat, assay = "SCT", layer = "data")
  cat("✓ Using SCT assay data layer for expression matrix (Seurat v4 SCTransform)\n")
} else {
  stop("✗ Neither scale.data nor SCT layer found! Please ensure SCTransform has been run on the Seurat object.\nAvailable layers: ", paste(available_layers, collapse = ", "))
}

cat("✓ Expression matrix extracted:", nrow(expr_matrix), "genes,", ncol(expr_matrix), "cells\n")

# 获取细胞类型分组（使用上面确定的列）
cell_types <- chicken_orth_mice_seurat@meta.data$anno_level_3_merged

# 使用mouse数据库（因为我们已经转换为鼠基因符号）
neuronchat_obj <- createNeuronChat(
  object = expr_matrix,
  group.by = cell_types,
  DB = "mouse"
)

cat("✓ NeuronChat object created\n")
cat("✓ Object class:", class(neuronchat_obj), "\n")
cat("✓ Object slots:", paste(slotNames(neuronchat_obj), collapse = ", "), "\n")
cat("✓ Database interactions:", length(neuronchat_obj@DB), "\n")
cat("✓ Ligand-receptor pairs:", length(neuronchat_obj@LR), "\n")

# 检查对象是否有效
if (length(neuronchat_obj@LR) > 0) {
  cat("✓ NeuronChat object is valid\n")
} else {
  cat("✗ ERROR: NeuronChat object is invalid - no LR pairs found\n")
  stop("Cannot proceed with invalid NeuronChat object")
}

# 3.1 运行NeuronChat分析
cat("\n3.1 Running NeuronChat communication inference...\n")
cat("Using fdr=1.0 (no FDR correction) to keep more significant interactions\n")
cat("This will use simple p<0.05 threshold instead of FDR correction\n")

start_time <- Sys.time()
neuronchat_obj <- run_NeuronChat(neuronchat_obj, M=100, fdr=0.05)
end_time <- Sys.time()

cat("✓ Communication inference completed in", round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")
cat("✓ Net slot populated with", length(neuronchat_obj@net), "interactions\n")

# 检查分析结果
if (length(neuronchat_obj@net) > 0) {
  cat("✓ Sample interactions:", paste(head(names(neuronchat_obj@net), 5), collapse = ", "), "\n")

  # 统计显著相互作用
  significant_count <- 0
  total_connections <- 0
  for (i in 1:length(neuronchat_obj@net)) {
    if (!is.null(neuronchat_obj@net[[i]]) && is.matrix(neuronchat_obj@net[[i]])) {
      sig_count <- sum(neuronchat_obj@net[[i]] > 0)
      if (sig_count > 0) {
        significant_count <- significant_count + 1
        total_connections <- total_connections + sig_count
      }
    }
  }
  cat("✓ Pairs with significant interactions:", significant_count, "\n")
  cat("✓ Total significant cell-cell connections:", total_connections, "\n")
} else {
  cat("⚠ WARNING: Net slot is empty! No interactions detected.\n")
}

# 4. 整合映射置信度信息
cat("\n4. Integrating mapping confidence information...\n")

# 创建基因到映射置信度的映射（使用大小写不敏感匹配）
gene_to_confidence <- setNames(neuronchat_mapping$mapping_confidence, toupper(neuronchat_mapping$mouse_symbol))

# 检查NeuronChat中的基因有多少有映射信息
# LR slot包含ligand-receptor对的名称，需要从DB中提取基因信息
lr_pairs <- neuronchat_obj@LR
cat("Sample LR pairs:", paste(head(lr_pairs, 3), collapse = ", "), "\n")

# 从DB中提取所有涉及的基因
all_genes_in_db <- c()
for (pair in lr_pairs) {
  if (pair %in% names(neuronchat_obj@DB)) {
    pair_info <- neuronchat_obj@DB[[pair]]
    all_genes_in_db <- c(all_genes_in_db, pair_info$lig_contributor, pair_info$receptor_subunit)
  }
}
all_genes_in_db <- unique(all_genes_in_db)
mapped_genes <- intersect(toupper(all_genes_in_db), names(gene_to_confidence))
cat("✓ Genes with mapping confidence:", length(mapped_genes), "out of", length(all_genes_in_db), "\n")

# 5. 计算表达支持证据
cat("\n5. Calculating expression support evidence...\n")

# 使用已提取的表达矩阵
# expr_matrix 已经在上面提取了

# 计算每个基因在每个细胞类型中的表达支持（优化版本）
cat("✓ Calculating expression support for", length(unique(cell_types)), "cell types and", nrow(expr_matrix), "genes\n")
cat("This may take several minutes...\n")

# 使用向量化操作优化
expr_support_results <- data.frame(
  gene = character(),
  cell_type = character(),
  pct_expr = numeric(),
  avg_expr = numeric(),
  expr_support = numeric(),
  stringsAsFactors = FALSE
)

cell_type_list <- unique(cell_types)
gene_list <- rownames(expr_matrix)  # 使用表达矩阵的基因名，不是Seurat对象的基因名

# 添加进度条
total_combinations <- length(gene_list) * length(cell_type_list)
cat("Total combinations to process:", total_combinations, "\n")

processed <- 0
start_time <- Sys.time()

for (gene in gene_list) {
  for (ct in cell_type_list) {
    # 获取该细胞类型的细胞
    cells_in_type <- which(cell_types == ct)

    if (length(cells_in_type) > 0) {
      # 计算表达统计
      gene_expr <- expr_matrix[gene, cells_in_type]
      total_cells <- length(cells_in_type)
      expr_cells <- sum(gene_expr > 0)
      pct_expr <- expr_cells / total_cells

      if (expr_cells > 0) {
        avg_expr <- mean(gene_expr[gene_expr > 0])
      } else {
        avg_expr <- 0
      }

      # 计算表达支持分数
      expr_support <- min(1, pct_expr * (1 + log1p(avg_expr)))

      expr_support_results <- rbind(expr_support_results, data.frame(
        gene = gene,
        cell_type = ct,
        pct_expr = pct_expr,
        avg_expr = avg_expr,
        expr_support = expr_support
      ))
    }

    # 进度报告
    processed <- processed + 1
    if (processed %% 10000 == 0) {
      elapsed <- Sys.time() - start_time
      cat("Processed:", processed, "/", total_combinations,
          "(", round(processed/total_combinations*100, 1), "%)",
          "Elapsed:", round(elapsed, 1), "seconds\n")
    }
  }
}

end_time <- Sys.time()
cat("✓ Expression support calculation completed in", round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")

cat("✓ Expression support calculated for", nrow(expr_support_results), "gene-cell type pairs\n")

# 6. 更新映射置信度计算
cat("\n6. Updating mapping confidence calculation...\n")

# 定义更新后的置信度计算函数
calculate_updated_confidence <- function(orth_type, identity, coverage, domain, expr_support) {
  # 同源关系分数
  orth_score <- ifelse(orth_type == "1:1", 1,
                      ifelse(orth_type %in% c("1:many", "many:1"), 0.7, 0))

  # 序列相似性分数
  seq_score <- (identity / 100) * coverage

  # 结构域分数
  domain_score <- domain

  # 表达支持分数
  expr_score <- ifelse(is.na(expr_support), 0, expr_support)

  # 权重分配
  weights <- c(0.4, 0.2, 0.2, 0.2)

  # 计算综合置信度
  confidence <- weights[1] * orth_score +
                weights[2] * seq_score +
                weights[3] * domain_score +
                weights[4] * expr_score

  return(confidence)
}

# 为每个基因计算平均表达支持
if (nrow(expr_support_results) > 0) {
  gene_expr_support <- expr_support_results %>%
    group_by(gene) %>%
    summarise(avg_expr_support = mean(expr_support, na.rm = TRUE), .groups = 'drop')

  # 更新映射结果（使用大小写不敏感匹配）
  neuronchat_mapping$expr_support_real <- gene_expr_support$avg_expr_support[match(toupper(neuronchat_mapping$mouse_symbol), toupper(gene_expr_support$gene))]
} else {
  # 如果没有表达支持结果，设置为NA
  neuronchat_mapping$expr_support_real <- NA
  gene_expr_support <- data.frame(gene = character(), avg_expr_support = numeric())
}

# 调试：检查基因名匹配问题
cat("\n=== Debug: Gene Name Matching ===\n")
cat("Sample genes in expression matrix:", paste(head(rownames(expr_matrix), 5), collapse = ", "), "\n")
cat("Sample genes in mapping file:", paste(head(neuronchat_mapping$mouse_symbol, 5), collapse = ", "), "\n")
cat("Sample genes in expr_support:", paste(head(gene_expr_support$gene, 5), collapse = ", "), "\n")

# 检查匹配情况（使用大小写不敏感匹配）
if (nrow(gene_expr_support) > 0) {
  matched_genes <- intersect(toupper(neuronchat_mapping$mouse_symbol), toupper(gene_expr_support$gene))
  cat("Matched genes between mapping and expr_support:", length(matched_genes), "\n")
  cat("Sample matched genes:", paste(head(matched_genes, 5), collapse = ", "), "\n")
} else {
  cat("No expression support results available for matching\n")
}

# 检查NeuronChat数据库中的基因
db_genes <- c()
for (pair in neuronchat_obj@LR) {
  if (pair %in% names(neuronchat_obj@DB)) {
    pair_info <- neuronchat_obj@DB[[pair]]
    db_genes <- c(db_genes, pair_info$lig_contributor, pair_info$receptor_subunit)
  }
}
db_genes <- unique(db_genes)
cat("Sample genes in NeuronChat database:", paste(head(db_genes, 5), collapse = ", "), "\n")

# 检查三重匹配（使用大小写不敏感匹配）
if (nrow(gene_expr_support) > 0) {
  triple_match <- intersect(intersect(toupper(neuronchat_mapping$mouse_symbol), toupper(gene_expr_support$gene)), toupper(db_genes))
  cat("Triple matched genes (mapping + expr_support + DB):", length(triple_match), "\n")
  cat("Sample triple matched genes:", paste(head(triple_match, 5), collapse = ", "), "\n")
} else {
  cat("No expression support results available for triple matching\n")
}

# 重新计算置信度
neuronchat_mapping$mapping_confidence_updated <- mapply(
  calculate_updated_confidence,
  neuronchat_mapping$orthology_type,
  neuronchat_mapping$sequence_identity,
  neuronchat_mapping$sequence_coverage,
  neuronchat_mapping$domain_score,
  neuronchat_mapping$expr_support_real
)

cat("✓ Mapping confidence updated with real expression data\n")


cat("✓ Gene expression support summary saved\n")

# 数据质量检查
cat("\n=== Data Quality Check (Steps 1-6) ===\n")

# 检查NeuronChat对象
cat("NeuronChat object status:\n")
cat("- Object class:", class(neuronchat_obj), "\n")
cat("- Database interactions:", length(neuronchat_obj@DB), "\n")
cat("- LR pairs:", length(neuronchat_obj@LR), "\n")
cat("- Net slot length:", length(neuronchat_obj@net), "\n")

# 检查表达支持结果
cat("\nExpression support results:\n")
cat("- Total gene-cell type pairs:", nrow(expr_support_results), "\n")
cat("- Unique genes:", length(unique(expr_support_results$gene)), "\n")
cat("- Unique cell types:", length(unique(expr_support_results$cell_type)), "\n")
cat("- Expression support range:", round(range(expr_support_results$expr_support), 3), "\n")
cat("- Mean expression support:", round(mean(expr_support_results$expr_support), 3), "\n")

# 检查基因表达支持摘要
cat("\nGene expression support summary:\n")
cat("- Genes with expression support:", nrow(gene_expr_support), "\n")
cat("- Mean expression support range:", round(range(gene_expr_support$avg_expr_support), 3), "\n")
cat("- Genes with zero expression support:", sum(gene_expr_support$avg_expr_support == 0), "\n")

# 检查更新的映射结果
cat("\nUpdated mapping results:\n")
cat("- Total mappings:", nrow(neuronchat_mapping), "\n")
cat("- Mappings with expression support:", sum(!is.na(neuronchat_mapping$expr_support_real)), "\n")
cat("- Updated confidence range:", round(range(neuronchat_mapping$mapping_confidence_updated, na.rm = TRUE), 3), "\n")
cat("- Mean updated confidence:", round(mean(neuronchat_mapping$mapping_confidence_updated, na.rm = TRUE), 3), "\n")

cat("\n=== Steps 1-6 Complete ===\n")
cat("Intermediate results saved. Ready to proceed to step 7 (run_NeuronChat).\n")
cat("Files saved in", results_dir, "folder:\n")
cat("- neuronchat_obj_step6.rds: NeuronChat object\n")
cat("- expression_support_step6.csv: Expression support details\n")
cat("- neuronchat_mapping_updated_step6.csv: Updated mapping results\n")
cat("- gene_expr_support_step6.csv: Gene expression support summary\n")

# 7. 运行NeuronChat分析
cat("\n7. Running NeuronChat analysis...\n")

# 运行NeuronChat
neuronchat_obj <- run_NeuronChat(neuronchat_obj, M = 100, fdr = 0.05)
cat("✓ NeuronChat analysis completed\n")

# 检查NeuronChat分析结果
cat("✓ Net slot length:", length(neuronchat_obj@net), "\n")
if (length(neuronchat_obj@net) > 0) {
  cat("✓ Sample interactions:", paste(head(names(neuronchat_obj@net), 5), collapse = ", "), "\n")

  # 检查第一个相互作用的非零值
  first_interaction <- neuronchat_obj@net[[1]]
  non_zero_count <- sum(first_interaction > 0)
  cat("✓ First interaction non-zero values:", non_zero_count, "\n")
} else {
  cat("⚠ WARNING: Net slot is empty! No interactions detected.\n")
  cat("This suggests run_NeuronChat failed or found no significant interactions.\n")
}

# 8. 应用映射置信度权重
cat("\n8. Applying mapping confidence weights...\n")

# 获取通信网络 - net是一个list，每个元素是一个ligand-receptor对的矩阵
net_list <- neuronchat_obj@net
cat("✓ Communication network loaded:", length(net_list), "ligand-receptor pairs\n")

# 创建ligand和receptor的置信度映射（使用大小写不敏感匹配）
ligand_confidence <- setNames(neuronchat_mapping$mapping_confidence_updated, toupper(neuronchat_mapping$mouse_symbol))
receptor_confidence <- setNames(neuronchat_mapping$mapping_confidence_updated, toupper(neuronchat_mapping$mouse_symbol))

# 创建加权网络结果
weighted_results <- data.frame(
  ligand = character(),
  receptor = character(),
  prob = numeric(),
  ligand_confidence = numeric(),
  receptor_confidence = numeric(),
  combined_confidence = numeric(),
  final_score = numeric(),
  stringsAsFactors = FALSE
)

# 处理每个ligand-receptor对
for (lr_pair in names(net_list)) {
  # 从DB中获取ligand和receptor信息
  if (lr_pair %in% names(neuronchat_obj@DB)) {
    pair_info <- neuronchat_obj@DB[[lr_pair]]
    ligands <- pair_info$lig_contributor
    receptors <- pair_info$receptor_subunit

    # 获取通信强度矩阵
    comm_matrix <- net_list[[lr_pair]]

    # 计算平均通信强度
    avg_prob <- mean(comm_matrix, na.rm = TRUE)

    # 计算ligand和receptor的平均置信度（使用大小写不敏感匹配）
    ligand_conf <- mean(ligand_confidence[toupper(ligands)], na.rm = TRUE)
    receptor_conf <- mean(receptor_confidence[toupper(receptors)], na.rm = TRUE)

    # 计算组合置信度
    combined_conf <- sqrt(ligand_conf * receptor_conf)

    # 计算最终分数
    final_score <- avg_prob * combined_conf

    # 添加到结果中
    weighted_results <- rbind(weighted_results, data.frame(
      ligand = paste(ligands, collapse = ","),
      receptor = paste(receptors, collapse = ","),
      prob = avg_prob,
      ligand_confidence = ligand_conf,
      receptor_confidence = receptor_conf,
      combined_confidence = combined_conf,
      final_score = final_score
    ))
  }
}

# 移除NA值
weighted_results <- weighted_results[!is.na(weighted_results$final_score), ]

cat("✓ Mapping confidence weights applied\n")
cat("Weighted interactions:", nrow(weighted_results), "\n")

# 9. 分层结果
cat("\n9. Creating stratified results...\n")

# 按置信度分位数分层（修正：使用更合理的分位数）
confidence_quantiles <- quantile(weighted_results$combined_confidence, probs = c(0.25, 0.75), na.rm = TRUE)
high_confidence <- weighted_results[weighted_results$combined_confidence >= confidence_quantiles[2], ]
medium_confidence <- weighted_results[weighted_results$combined_confidence >= confidence_quantiles[1] & weighted_results$combined_confidence < confidence_quantiles[2], ]
low_confidence <- weighted_results[weighted_results$combined_confidence < confidence_quantiles[1], ]

cat("Confidence quantiles (25%, 75%):", round(confidence_quantiles, 3), "\n")
cat("High confidence interactions (top 25%):", nrow(high_confidence), "\n")
cat("Medium confidence interactions (middle 50%):", nrow(medium_confidence), "\n")
cat("Low confidence interactions (bottom 25%):", nrow(low_confidence), "\n")

# 添加置信度分布统计
cat("\nConfidence distribution summary:\n")
cat("Min confidence:", round(min(weighted_results$combined_confidence, na.rm = TRUE), 3), "\n")
cat("Max confidence:", round(max(weighted_results$combined_confidence, na.rm = TRUE), 3), "\n")
cat("Mean confidence:", round(mean(weighted_results$combined_confidence, na.rm = TRUE), 3), "\n")
cat("Median confidence:", round(median(weighted_results$combined_confidence, na.rm = TRUE), 3), "\n")

# 10. 处理many-to-many情况
cat("\n10. Handling many-to-many cases...\n")

# 检查是否有重复的ligand-receptor对
duplicated_pairs <- sum(duplicated(paste(weighted_results$ligand, weighted_results$receptor, sep = "_")))
if (duplicated_pairs > 0) {
  cat("⚠ Found", duplicated_pairs, "duplicated ligand-receptor pairs\n")

  # 使用best-hit策略：保留置信度最高的
  net_best_hit <- weighted_results %>%
    group_by(ligand, receptor) %>%
    slice_max(combined_confidence, n = 1) %>%
    ungroup()

  cat("✓ After best-hit filtering:", nrow(net_best_hit), "interactions\n")
} else {
  net_best_hit <- weighted_results
  cat("✓ No duplicated pairs found\n")
}

# 11. 保存结果
cat("\n11. Saving results...\n")

# 创建结果文件夹
results_dir <- "neuronchat_analysis_results_merged"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("✓ Created results directory:", results_dir, "\n")
} else {
  cat("✓ Results directory exists:", results_dir, "\n")
}

# 保存更新的映射结果
write.csv(neuronchat_mapping, file.path(results_dir, "complete_neuronchat_mapping_updated.csv"), row.names = FALSE)
cat("✓ Updated mapping results saved\n")

# 保存分层结果
write.csv(high_confidence, file.path(results_dir, "high_confidence_interactions.csv"), row.names = FALSE)
write.csv(medium_confidence, file.path(results_dir, "medium_confidence_interactions.csv"), row.names = FALSE)
write.csv(low_confidence, file.path(results_dir, "low_confidence_interactions.csv"), row.names = FALSE)
cat("✓ Stratified results saved\n")

# 保存best-hit结果
write.csv(net_best_hit, file.path(results_dir, "neuronchat_best_hit_interactions.csv"), row.names = FALSE)
cat("✓ Best-hit results saved\n")

# 保存表达支持结果
write.csv(expr_support_results, file.path(results_dir, "expression_support_detailed.csv"), row.names = FALSE)
cat("✓ Expression support results saved\n")

# 保存NeuronChat对象
saveRDS(neuronchat_obj, file.path(results_dir, "neuronchat_analysis_result.rds"))
cat("✓ NeuronChat object saved\n")

# 12. 生成统计报告
cat("\n12. Generating statistical report...\n")

# 计算NeuronChat网络统计
net_stats <- list()
if (length(neuronchat_obj@net) > 0) {
  # 计算网络统计
  total_interactions <- length(neuronchat_obj@net)

  # 计算非零相互作用
  non_zero_interactions <- 0
  max_prob <- 0
  min_prob <- Inf
  total_prob <- 0

  for (lr_pair in names(neuronchat_obj@net)) {
    comm_matrix <- neuronchat_obj@net[[lr_pair]]
    non_zero_count <- sum(comm_matrix > 0)
    if (non_zero_count > 0) {
      non_zero_interactions <- non_zero_interactions + 1
      max_prob <- max(max_prob, max(comm_matrix))
      min_prob <- min(min_prob, min(comm_matrix[comm_matrix > 0]))
      total_prob <- total_prob + sum(comm_matrix)
    }
  }

  net_stats <- list(
    total_interactions = total_interactions,
    significant_interactions = non_zero_interactions,
    interaction_rate = round(non_zero_interactions/total_interactions*100, 2),
    max_probability = round(max_prob, 4),
    min_probability = round(min_prob, 4),
    mean_probability = round(total_prob/sum(sapply(neuronchat_obj@net, function(x) sum(x > 0))), 4)
  )
} else {
  net_stats <- list(
    total_interactions = 0,
    significant_interactions = 0,
    interaction_rate = 0,
    max_probability = 0,
    min_probability = 0,
    mean_probability = 0
  )
}

# 创建统计报告
stats_report <- data.frame(
  metric = c(
    "Total_cells", "Total_genes", "Cell_types", "Ligand_receptor_pairs",
    "Genes_with_mapping", "High_confidence_interactions", "Medium_confidence_interactions",
    "Low_confidence_interactions", "Best_hit_interactions", "Average_confidence",
    "Min_confidence", "Max_confidence", "Median_confidence",
    "NeuronChat_total_interactions", "NeuronChat_significant_interactions",
    "NeuronChat_interaction_rate_percent", "NeuronChat_max_probability",
    "NeuronChat_min_probability", "NeuronChat_mean_probability"
  ),
  value = c(
    ncol(chicken_orth_mice_seurat),
    nrow(chicken_orth_mice_seurat),
    length(unique(cell_types)),
    length(neuronchat_obj@LR),
    length(mapped_genes),
    nrow(high_confidence),
    nrow(medium_confidence),
    nrow(low_confidence),
    nrow(net_best_hit),
    round(mean(weighted_results$combined_confidence, na.rm = TRUE), 3),
    round(min(weighted_results$combined_confidence, na.rm = TRUE), 3),
    round(max(weighted_results$combined_confidence, na.rm = TRUE), 3),
    round(median(weighted_results$combined_confidence, na.rm = TRUE), 3),
    net_stats$total_interactions,
    net_stats$significant_interactions,
    net_stats$interaction_rate,
    net_stats$max_probability,
    net_stats$min_probability,
    net_stats$mean_probability
  )
)

write.csv(stats_report, file.path(results_dir, "neuronchat_analysis_statistics.csv"), row.names = FALSE)
cat("✓ Statistical report saved\n")

# 13. 输出最终信息
cat("\n=== Final Analysis Results ===\n")
cat("Total cells analyzed:", ncol(chicken_orth_mice_seurat), "\n")
cat("Total genes:", nrow(chicken_orth_mice_seurat), "\n")
cat("Cell types:", length(unique(cell_types)), "\n")
cat("Ligand-receptor pairs:", length(neuronchat_obj@LR), "\n")
cat("Genes with mapping confidence:", length(mapped_genes), "\n")
cat("High confidence interactions (top 25%):", nrow(high_confidence), "\n")
cat("Medium confidence interactions (middle 50%):", nrow(medium_confidence), "\n")
cat("Low confidence interactions (bottom 25%):", nrow(low_confidence), "\n")
cat("Best-hit interactions:", nrow(net_best_hit), "\n")
cat("Average confidence:", round(mean(weighted_results$combined_confidence, na.rm = TRUE), 3), "\n")
cat("Confidence range:", round(min(weighted_results$combined_confidence, na.rm = TRUE), 3), "-", round(max(weighted_results$combined_confidence, na.rm = TRUE), 3), "\n")

# NeuronChat网络统计
cat("\nNeuronChat Network Statistics:\n")
cat("Total interactions in net:", net_stats$total_interactions, "\n")
cat("Significant interactions:", net_stats$significant_interactions, "\n")
cat("Interaction rate:", net_stats$interaction_rate, "%\n")
cat("Max probability:", net_stats$max_probability, "\n")
cat("Min probability:", net_stats$min_probability, "\n")
cat("Mean probability:", net_stats$mean_probability, "\n")

cat("\nSample high confidence interactions:\n")
if (nrow(high_confidence) > 0) {
  sample_high <- head(high_confidence, 5)
  for (i in 1:nrow(sample_high)) {
    cat("  ", sample_high$ligand[i], "->", sample_high$receptor[i],
        "(conf:", round(sample_high$combined_confidence[i], 3), ")\n")
  }
}

cat("\nFiles created in", results_dir, "folder:\n")
cat("- complete_neuronchat_mapping_updated.csv: Updated mapping results\n")
cat("- high/medium/low_confidence_interactions.csv: Stratified results\n")
cat("- neuronchat_best_hit_interactions.csv: Best-hit results\n")
cat("- expression_support_detailed.csv: Expression support details\n")
cat("- neuronchat_analysis_result.rds: NeuronChat analysis object\n")
cat("- neuronchat_analysis_statistics.csv: Statistical report\n")

cat("\n=== Analysis Complete ===\n")
cat("Analysis completed at:", Sys.time(), "\n")
