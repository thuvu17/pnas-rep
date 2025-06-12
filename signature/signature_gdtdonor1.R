# Run Single-Cell Signature on purified donor 1 GDT

# Single-Cell Signature Scorer
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/data")
write.table(t(as.matrix(gdt_donor1[["RNA"]]$counts)), "gdt_donor1_expression.tsv", sep="\t", row.names = TRUE, col.names=NA)

file.copy(from = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureScorer/results/PNAS_REP_gdt_donor1_expression.tsv",
          to = "C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger/scores/PNAS_REP_gdt_donor1_expression.tsv")


# Single-Cell Signature Merger
setwd("C:/Users/vnmt1/Downloads/SingleCellSignatureExplorer_20240613/SingleCellSignatureMerger")
tsne_coords <- Embeddings(gdt_donor1, reduction = "tsne")
tsne_df <- data.frame(CellName = rownames(tsne_coords), tsne_coords)
write.table(tsne_df,"tsne/gdt_donor1_tsne_coordinate.tsv", sep="\t", row.names = FALSE)


# Single-Cell Signature Combiner
signature_score <- fread("gdt_donor1_signature_score.csv")
sum(signature_score$CombinedScore > 0.35) # 1242
length(Cells(gdt_donor1)) # 1526

