#!/usr/bin/env Rscript
# Run Cell Type Analysis
# 运行细胞类型分析

cat("=== Running Cell Type Analysis ===\n")

# Load required libraries
library(Seurat)
library(NeuronChat)
library(dplyr)

# Set random seed
set.seed(123)

# =============================================================================
# 1. Load Data
# =============================================================================
cat("\n1. Loading data...\n")

# Load seurat_obj from RDS file
seurat_obj_path <- "neuron_subset_working.rds"
if (!file.exists(seurat_obj_path)) {
  stop("neuron_subset_working.rds file not found.")
}

seurat_obj <- readRDS(seurat_obj_path)
cat("✓ Seurat object loaded\n")

# Get groups and cell types
groups <- unique(seurat_obj$group)

# Find cell type column
possible_cell_type_cols <- c("neuron_subtype", "cell_type", "celltype", "CellType", "cell.types", "annotation", "cluster")
cell_type_col <- NULL

for (col in possible_cell_type_cols) {
  if (col %in% colnames(seurat_obj@meta.data)) {
    cell_type_col <- col
    break
  }
}

if (is.null(cell_type_col)) {
  stop("No cell type column found. Available columns: ", paste(colnames(seurat_obj@meta.data), collapse = ", "))
}

cell_types <- unique(seurat_obj@meta.data[[cell_type_col]])

cat("Groups:", paste(groups, collapse = ", "), "\n")
cat("Cell type column:", cell_type_col, "\n")
cat("Cell types:", paste(cell_types, collapse = ", "), "\n")

# =============================================================================
# 2. Create Group-Specific NeuronChat Objects
# =============================================================================
cat("\n2. Creating group-specific NeuronChat objects...\n")

group_neuronchat_objects <- list()

for (group in groups) {
  cat("Processing group:", group, "\n")
  
  # Get cells for this group
  group_cells <- which(seurat_obj$group == group)
  group_seurat <- subset(seurat_obj, cells = group_cells)
  
  cat("- Cells in group:", ncol(group_seurat), "\n")
  
  # Extract expression and metadata
  expr_matrix <- GetAssayData(group_seurat, slot = "data")
  cell_meta <- data.frame(
    cell_id = colnames(expr_matrix),
    neuron_subtype = as.character(group_seurat@meta.data[[cell_type_col]]),
    group = group,
    stringsAsFactors = FALSE
  )
  
  # Create NeuronChat object
  tryCatch({
    group_neuronchat_obj <- createNeuronChat(
      object = expr_matrix,
      meta = cell_meta,
      group.by = cell_meta$neuron_subtype,
      DB = 'mouse'
    )
    
    # Run NeuronChat analysis
    group_neuronchat_obj <- run_NeuronChat(
      object = group_neuronchat_obj,
      M = 100,
      fdr = 0.05
    )
    
    group_neuronchat_objects[[group]] <- group_neuronchat_obj
    
    # Check results
    net <- group_neuronchat_obj@net
    total_strength <- sum(sapply(net, function(x) sum(x, na.rm = TRUE)))
    sig_interactions <- sum(sapply(net, function(x) sum(x > 0, na.rm = TRUE)))
    
    cat("✓ Group", group, "analyzed\n")
    cat("- Total strength:", total_strength, "\n")
    cat("- Significant interactions:", sig_interactions, "\n")
    cat("- LR pairs with activity:", sum(sapply(net, function(x) sum(x > 0, na.rm = TRUE)) > 0), "\n")
    
  }, error = function(e) {
    cat("Error analyzing group", group, ":", e$message, "\n")
    group_neuronchat_objects[[group]] <- NULL
  })
}

# =============================================================================
# 3. Analyze Cell Type Involvement
# =============================================================================
cat("\n3. Analyzing cell type involvement...\n")

if (length(group_neuronchat_objects) >= 2) {
  group1 <- groups[1]
  group2 <- groups[2]
  
  net1 <- group_neuronchat_objects[[group1]]@net
  net2 <- group_neuronchat_objects[[group2]]@net
  
  cat("Group 1 (", group1, ") network summary:\n")
  cat("- LR pairs:", length(net1), "\n")
  cat("- Total strength:", sum(sapply(net1, function(x) sum(x, na.rm = TRUE))), "\n")
  cat("- Active LR pairs:", sum(sapply(net1, function(x) sum(x > 0, na.rm = TRUE)) > 0), "\n")
  
  cat("Group 2 (", group2, ") network summary:\n")
  cat("- LR pairs:", length(net2), "\n")
  cat("- Total strength:", sum(sapply(net2, function(x) sum(x, na.rm = TRUE))), "\n")
  cat("- Active LR pairs:", sum(sapply(net2, function(x) sum(x > 0, na.rm = TRUE)) > 0), "\n")
  
  # Find common LR pairs
  common_lr_pairs <- intersect(names(net1), names(net2))
  cat("Common LR pairs:", length(common_lr_pairs), "\n")
  
  # Analyze cell type involvement
  cell_type_involvement <- data.frame()
  
  if (length(common_lr_pairs) > 0) {
    cat("Analyzing cell type involvement for common LR pairs...\n")
    
    for (lr_pair in common_lr_pairs[1:min(10, length(common_lr_pairs))]) {
      cat("Processing LR pair:", lr_pair, "\n")
      
      matrix1 <- net1[[lr_pair]]
      matrix2 <- net2[[lr_pair]]
      
      cat("- Matrix1 dimensions:", dim(matrix1), "\n")
      cat("- Matrix2 dimensions:", dim(matrix2), "\n")
      cat("- Matrix1 non-zero entries:", sum(matrix1 > 0, na.rm = TRUE), "\n")
      cat("- Matrix2 non-zero entries:", sum(matrix2 > 0, na.rm = TRUE), "\n")
      
      # Find cell type pairs with activity
      for (sender in rownames(matrix1)) {
        for (receiver in colnames(matrix1)) {
          if (sender %in% rownames(matrix2) && receiver %in% colnames(matrix2)) {
            val1 <- matrix1[sender, receiver]
            val2 <- matrix2[sender, receiver]
            
            if (!is.na(val1) && !is.na(val2) && (val1 > 0 || val2 > 0)) {
              change_ratio <- val2 / (val1 + 1e-10)
              
              cell_type_involvement <- rbind(cell_type_involvement, data.frame(
                lr_pair = lr_pair,
                sender_cell_type = sender,
                receiver_cell_type = receiver,
                group1_strength = val1,
                group2_strength = val2,
                change_ratio = change_ratio,
                log2_change = log2(change_ratio),
                stringsAsFactors = FALSE
              ))
              
              cat("  - Added:", sender, "->", receiver, "val1:", val1, "val2:", val2, "\n")
            }
          }
        }
      }
    }
  }
  
  # If no cell type involvement data, create basic summary
  if (nrow(cell_type_involvement) == 0) {
    cat("No cell type involvement data found. Creating basic summary...\n")
    
    # Create basic cell type summary
    cell_type_summary <- data.frame(
      cell_type = cell_types,
      group1_cells = sapply(cell_types, function(ct) sum(seurat_obj@meta.data[[cell_type_col]] == ct & seurat_obj$group == group1)),
      group2_cells = sapply(cell_types, function(ct) sum(seurat_obj@meta.data[[cell_type_col]] == ct & seurat_obj$group == group2)),
      stringsAsFactors = FALSE
    )
    
    # Calculate ratios
    cell_type_summary$cell_ratio <- cell_type_summary$group2_cells / (cell_type_summary$group1_cells + 1e-10)
    cell_type_summary$log2_ratio <- log2(cell_type_summary$cell_ratio)
    
    cell_type_involvement <- cell_type_summary
    cat("Created basic cell type summary with", nrow(cell_type_involvement), "entries.\n")
  }
  
} else {
  cat("Not enough groups analyzed. Creating basic cell type summary...\n")
  
  # Create basic cell type summary
  cell_type_summary <- data.frame(
    cell_type = cell_types,
    group1_cells = sapply(cell_types, function(ct) sum(seurat_obj@meta.data[[cell_type_col]] == ct & seurat_obj$group == groups[1])),
    group2_cells = sapply(cell_types, function(ct) sum(seurat_obj@meta.data[[cell_type_col]] == ct & seurat_obj$group == groups[2])),
    stringsAsFactors = FALSE
  )
  
  # Calculate ratios
  cell_type_summary$cell_ratio <- cell_type_summary$group2_cells / (cell_type_summary$group1_cells + 1e-10)
  cell_type_summary$log2_ratio <- log2(cell_type_summary$cell_ratio)
  
  cell_type_involvement <- cell_type_summary
}

# =============================================================================
# 4. Save Results
# =============================================================================
cat("\n4. Saving results...\n")

# Save cell type involvement
output_file <- file.path("mice hipp exercise", "neuronchat_complete_analysis", "data", "cell_type_involvement.csv")
write.csv(cell_type_involvement, output_file, 
          row.names = FALSE, fileEncoding = "UTF-8")

cat("✓ Cell type involvement data saved\n")
cat("Entries found:", nrow(cell_type_involvement), "\n")

# Print summary
if (nrow(cell_type_involvement) > 0) {
  cat("\nCell type involvement summary:\n")
  print(head(cell_type_involvement, 10))
} else {
  cat("No cell type involvement data found.\n")
}

cat("==================\n")
