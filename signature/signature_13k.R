# Run Single-Cell Signature on PBMC 13k

# On integrated$data
# Single-Cell Signature Scorer
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/data")
write.table(t(as.matrix(part1.v4[["integrated"]]$data)), "pbmc13k_expression.tsv", sep="\t", row.names = TRUE, col.names=NA)

file.copy(from = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/results/PNAS_REP_gdt_donor1_expression.tsv",
          to = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger/scores/PNAS_REP_gdt_donor1_expression.tsv")


# Single-Cell Signature Merger
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger")
tsne_coords <- Embeddings(part1.v4, reduction = "tsne")
tsne_df <- data.frame(CellName = rownames(tsne_coords), tsne_coords)
write.table(tsne_df,"tsne/pbmc13k_tsne_coordinate.tsv", sep="\t", row.names = FALSE)


# Single-Cell Signature Combiner
signature_score <- fread("pbmc13k_signature_score_integrated_data.csv")
sum(signature_score$CombinedScore > 0.35) # 1242
length(Cells(part1.v4)) # 1526

DimPlot(part1.v4, reduction = "tsne", group.by = "orig.ident") + ggtitle("part 1 Seurat v4")

# ============================================================================
# On RNA$data
# Single-Cell Signature Scorer
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/data")
write.table(t(as.matrix(pbmc_13k.data)), "pbmc13k_expression.tsv", sep="\t", row.names = TRUE, col.names=NA)

file.copy(from = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/results/PNAS_REP_gdt_donor1_expression.tsv",
          to = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger/scores/PNAS_REP_gdt_donor1_expression.tsv")


# Single-Cell Signature Merger
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger")
tsne_coords <- Embeddings(part1.v4, reduction = "tsne")
tsne_df <- data.frame(CellName = rownames(tsne_coords), tsne_coords)
write.table(tsne_df,"tsne/pbmc13k_tsne_coordinate.tsv", sep="\t", row.names = FALSE)


# Single-Cell Signature Combiner
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureCombiner/results")
signature_score <- fread("pbmc13k_signature_score_rna_data.csv")
sum(signature_score$CombinedScore > 0.35) # 2139
length(Cells(part1.v4)) # 14247


print("CD3E" %in% rownames(part1.v4[["integrated"]])) # FALSE
print("CD3D" %in% rownames(part1.v4[["integrated"]])) # TRUE
print("TRDC" %in% rownames(part1.v4[["integrated"]])) # TRUE
print("TRGC1" %in% rownames(part1.v4[["integrated"]])) # TRUE
print("TRGC2" %in% rownames(part1.v4[["integrated"]])) # TRUE
print("CD8A" %in% rownames(part1.v4[["integrated"]])) # TRUE
print("CD8B" %in% rownames(part1.v4[["integrated"]])) # TRUE