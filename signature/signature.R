# These are steps to run the Single Cell Signature Score Explorer

# Single-Cell Signature Scorer
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer")
# Merging raw data
# Extract each layer explicitly
data_4k  <- GetAssayData(pbmc12k.v4, assay = "RNA", slot = "data.pbmc_4k")
data_8k  <- GetAssayData(pbmc12k.v4, assay = "RNA", slot = "data.pbmc_8k")
data_gdt <- GetAssayData(part1.v4, assay = "RNA", slot = "data.gdt")

# Combine into a single matrix (matching by gene names)
pbmc_12k.data <- RowMergeSparseMatrices(data_4k, data_8k)
pbmc_13k.data <- RowMergeSparseMatrices(pbmc_12k.data, data_gdt)

# Convert to dense matrix and transpose (cells as rows, genes as columns)
expr_12k <- t(as.matrix(pbmc_12k.data))
expr_13k <- t(as.matrix(pbmc_13k.data))

# Clean cell names: replace characters that may crash the downstream tool
rownames(expr_12k) <- gsub("[/\\\\:\\-]", "_", rownames(expr_12k))
rownames(expr_13k) <- gsub("[/\\\\:\\-]", "_", rownames(expr_13k))

expr_12k_df <- data.frame(id = rownames(expr_12k), expr_12k)
expr_13k_df <- data.frame(id = rownames(expr_13k), expr_13k)

# Write to .tsv with 'id' as the first column
write.table(
  expr_12k_df,
  file = "all_data/pbmc12k_expression_rna_data.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  expr_13k_df,
  file = "all_data/pbmc13k_expression_rna_data.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Double click on Single-Cell Signature Scorer exe file to run
file.copy(from = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/results/PNAS_REP_pbmc13k_expression.tsv",
          to = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger/all_scores/PNAS_REP_pbmc13k_expression.tsv")
file.copy(from = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/results/PNAS_REP_pbmc12k_expression.tsv",
          to = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger/all_scores/PNAS_REP_pbmc12k_expression.tsv")

# Single-Cell Signature Merger
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger")

#==================================================
# AFTER merging to create pbmc_12k and pbmc_13k
pbmc_12k <- merge(pbmc_4k, y = pbmc_8k)
pbmc_13k <- merge(pbmc_12k, y = gdt_donor1)

# Run preprocessing and tSNE for visualization only
pbmc_12k <- NormalizeData(pbmc_12k)
pbmc_12k <- FindVariableFeatures(pbmc_12k)
pbmc_12k <- ScaleData(pbmc_12k)
pbmc_12k <- RunPCA(pbmc_12k)
pbmc_12k <- RunTSNE(pbmc_12k, dims = 1:20)

pbmc_13k <- NormalizeData(pbmc_13k)
pbmc_13k <- FindVariableFeatures(pbmc_13k)
pbmc_13k <- ScaleData(pbmc_13k)
pbmc_13k <- RunPCA(pbmc_13k)
pbmc_13k <- RunTSNE(pbmc_13k, dims = 1:20)

DimPlot(pbmc_12k, reduction = "tsne", group.by = "orig.ident") +
  ggtitle("tSNE plot of PBMC 12k no integration")
DimPlot(pbmc_13k, reduction = "tsne", group.by = "orig.ident") +
  ggtitle("tSNE plot of PBMC 13k no integration")

# Export coordinates
tsne_12k <- Embeddings(pbmc_12k, "tsne")
tsne_13k <- Embeddings(pbmc_13k, "tsne")

# Make sure cell barcodes match expression file (clean if needed)
rownames(tsne_12k) <- gsub("[/\\\\:\\-]", "_", rownames(tsne_12k))
rownames(tsne_13k) <- gsub("[/\\\\:\\-]", "_", rownames(tsne_13k))

# Save as required for SCSE
write.table(data.frame(CellID = rownames(tsne_12k), tsne_12k),
            "all_tsne/My_tsne_coordinate_12k.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(CellID = rownames(tsne_13k), tsne_13k),
            "all_tsne/My_tsne_coordinate_13k.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# =======================================
# Single-Cell Signature Combiner
setwd("C:/Users/vnmt1/Downloads/pnas-rep/rep-paper/signature_score_results")
signature_score_12k <- fread("pbmc_12k_signature_score.csv")
signature_score_13k <- fread("pbmc_13k_signature_score.csv")
signature_score_12k_raw <- fread("pbmc12k_signature_score_raw.csv")
signature_score_13k_raw <- fread("pbmc13k_signature_score_raw.csv")
# Count for 12k dataset
sum(signature_score_12k$CombinedScore > 0.35) # 242
sum(signature_score_12k_raw$CombinedScore > 0.35) # 296
# Count for 13k dataset
sum(signature_score_13k$CombinedScore > 0.35) # 265
sum(signature_score_13k_raw$CombinedScore > 0.35) # 889

FeaturePlot(pbmc_12k, features = c("CD3E", "CD8A", "CD4", "TRGC1", "TRGC2", "TRDC"))
