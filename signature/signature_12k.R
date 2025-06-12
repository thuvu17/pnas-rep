# Preprocessing
pbmc12k.list <- list(pbmc_4k, pbmc_8k)
pbmc12k.list <- lapply(pbmc12k.list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  gene_counts <- LayerData(obj[["RNA"]], layer = "counts")
  gene_filter <- Matrix::rowSums(gene_counts > 0) >= 3
  obj <- obj[gene_filter, ]
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt >= 0.05)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize")
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  return(obj)
})

# Integration using Seurat v4
pbmc12k.features <- SelectIntegrationFeatures(pbmc12k.list)
pbmc12k.anchors <- FindIntegrationAnchors(pbmc12k.list, anchor.features = pbmc12k.features)
pbmc_12k.v4 <- IntegrateData(pbmc12k.anchors)

# Dimensionality reduction and clustering
DefaultAssay(pbmc_12k.v4) <- "integrated"
pbmc_12k.v4 <- ScaleData(pbmc_12k.v4)
pbmc_12k.v4 <- RunPCA(pbmc_12k.v4, npcs = 30)
pbmc_12k.v4 <- RunTSNE(pbmc_12k.v4, dims = 1:30)

# ===============================
# On integrated$data
# Single-Cell Signature Scorer
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/data")
write.table(t(as.matrix(pbmc_12k.v4[["integrated"]]$data)), "pbmc12k_expression_integrated_data.tsv", sep="\t", row.names = TRUE, col.names=NA)

file.copy(from = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/results/PNAS_REP_gdt_donor1_expression.tsv",
          to = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger/scores/PNAS_REP_gdt_donor1_expression.tsv")


# Single-Cell Signature Merger
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger")
tsne_coords <- Embeddings(pbmc_12k.v4, reduction = "tsne")
tsne_df <- data.frame(CellName = rownames(tsne_coords), tsne_coords)
write.table(tsne_df,"all_tsne/pbmc12k_tsne_coordinate.tsv", sep="\t", row.names = FALSE)


# Single-Cell Signature Combiner
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureCombiner/results")
pbmc12k_integrated_data_signature_score <- fread("pbmc12k_signature_score_integrated_data.csv")
sum(pbmc12k_integrated_data_signature_score$CombinedScore > 0.35) # 525
length(Cells(pbmc_12k.v4)) # 12721

# ============================================================================
# On RNA$data
# Single-Cell Signature Scorer
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/data")
write.table(t(as.matrix(pbmc_12k.data)), "pbmc12k_expression_rna_data.tsv", sep="\t", row.names = TRUE, col.names=NA)

# Single-Cell Signature Merger
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger")
tsne_coords <- Embeddings(pbmc_12k.v4, reduction = "tsne")
tsne_df <- data.frame(CellName = rownames(tsne_coords), tsne_coords)
write.table(tsne_df,"tsne/pbmc12k_tsne_coordinate.tsv", sep="\t", row.names = FALSE)


# Single-Cell Signature Combiner
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureCombiner/results")
signature_score <- fread("pbmc12k_signature_score_rna_data.csv")
sum(signature_score$CombinedScore > 0.35) # 992
length(Cells(part1.v4)) # 14247


print("CD3E" %in% rownames(pbmc_12k.v4[["integrated"]])) # FALSE
print("CD3D" %in% rownames(pbmc_12k.v4[["integrated"]])) # TRUE
print("TRDC" %in% rownames(pbmc_12k.v4[["integrated"]])) # TRUE
print("TRGC1" %in% rownames(pbmc_12k.v4[["integrated"]])) # TRUE
print("TRGC2" %in% rownames(pbmc_12k.v4[["integrated"]])) # TRUE
print("CD8A" %in% rownames(pbmc_12k.v4[["integrated"]])) # TRUE
print("CD8B" %in% rownames(pbmc_12k.v4[["integrated"]])) # TRUE
