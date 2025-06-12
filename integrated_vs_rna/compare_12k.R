# Matrix subtraction: integrated - rna
pbmc_4k.delta_12k <- pbmc_4k.int_12k - pbmc_4k.rna_12k
pbmc_8k.delta_12k <- pbmc_8k.int_12k - pbmc_8k.rna_12k
pbmc_4k.mean_delta_12k <- Matrix::rowMeans(pbmc_4k.delta_12k)
pbmc_8k.mean_delta_12k <- Matrix::rowMeans(pbmc_8k.delta_12k)
# Top 30 genes with highest mean shift
top_shift_genes_4k_12k <- head(sort(pbmc_4k.mean_delta_12k, decreasing = TRUE), 30)
top_shift_genes_8k_12k <- head(sort(pbmc_8k.mean_delta_12k, decreasing = TRUE), 30)

# Barplot
barplot(top_shift_genes_4k_12k,
        las = 2, cex.names = 0.7,
        main = "Top 30 Genes by Mean Δ (Integrated - RNA) in PBMC 4k(12k)",
        ylab = "Mean Expression Change",
        col = "steelblue")

barplot(top_shift_genes_8k_12k,
        las = 2, cex.names = 0.7,
        main = "Top 30 Genes by Mean Δ (Integrated - RNA) in PBMC 8k(12k)",
        ylab = "Mean Expression Change",
        col = "red")

# Imputed zeros in integrated
zero_imputed_4k_12k <- sum(pbmc_4k.rna_12k == 0 & pbmc_4k.int_12k != 0)
zero_imputed_8k_12k <- sum(pbmc_8k.rna_12k == 0 & pbmc_8k.int_12k != 0)
percent_imputed_4k_12k <- zero_imputed_4k_12k / length(pbmc_4k.rna_12k) * 100
percent_imputed_8k_12k <- zero_imputed_8k_12k / length(pbmc_8k.rna_12k) * 100
zero_imputed_4k_12k # 4600014
zero_imputed_8k_12k # 1575222
percent_imputed_4k_12k # 52.99555
percent_imputed_8k_12k # 9.397578

print("CD3E" %in% rownames(zero_imputed_4k_12k))
print("CD3D" %in% rownames(zero_imputed_4k_12k))
print("TRDC" %in% rownames(zero_imputed_4k_12k))
print("TRGC1" %in% rownames(zero_imputed_4k_12k))
print("TRGC2" %in% rownames(zero_imputed_4k_12k))
print("CD8A" %in% rownames(zero_imputed_4k_12k))
print("CD8B" %in% rownames(zero_imputed_4k_12k))
 # ALL FALSE