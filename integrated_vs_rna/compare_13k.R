# Matrix subtraction: integrated - rna
pbmc_4k.delta_13k <- pbmc_4k.int_13k - pbmc_4k.rna_13k
pbmc_8k.delta_13k <- pbmc_8k.int_13k - pbmc_8k.rna_13k
gdt.delta <- gdt_donor1.int_13k - gdt_donor1.rna_13k
pbmc_4k.mean_delta_13k <- Matrix::rowMeans(pbmc_4k.delta_13k)
pbmc_8k.mean_delta_13k <- Matrix::rowMeans(pbmc_8k.delta_13k)
gdt.mean_delta <- Matrix::rowMeans(gdt.delta)
# Top 30 genes with highest mean shift
top_shift_genes_4k_13k <- head(sort(pbmc_4k.mean_delta_13k, decreasing = TRUE), 30)
top_shift_genes_8k_13k <- head(sort(pbmc_8k.mean_delta_13k, decreasing = TRUE), 30)
top_shift_genes_gdt <- head(sort(gdt.mean_delta, decreasing = TRUE), 30)

# Barplot
barplot(top_shift_genes_4k_13k,
        las = 2, cex.names = 0.7,
        main = "Top 30 Genes by Mean Δ (Integrated - RNA) in PBMC 4k(13k)",
        ylab = "Mean Expression Change",
        col = "steelblue")

barplot(top_shift_genes_8k_13k,
        las = 2, cex.names = 0.7,
        main = "Top 30 Genes by Mean Δ (Integrated - RNA) in PBMC 8k(13k)",
        ylab = "Mean Expression Change",
        col = "red")

barplot(top_shift_genes_gdt,
        las = 2, cex.names = 0.7,
        main = "Top 30 Genes by Mean Δ (Integrated - RNA) in GDT donor1",
        ylab = "Mean Expression Change",
        col = "green")

# Imputed zeros in integrated
zero_imputed_4k_13k <- sum(pbmc_4k.rna_13k == 0 & pbmc_4k.int_13k != 0)
zero_imputed_8k_13k <- sum(pbmc_8k.rna_13k == 0 & pbmc_8k.int_13k != 0)
zero_imputed_gdt <- sum(gdt_donor1.rna_13k == 0 & gdt_donor1.int_13k != 0)
percent_imputed_4k_13k <- zero_imputed_4k_13k / length(pbmc_4k.rna_13k) * 100
percent_imputed_8k_13k <- zero_imputed_8k_13k / length(pbmc_8k.rna_13k) * 100
percent_imputed_gdt <- zero_imputed_gdt / length(gdt_donor1.rna_13k) * 100
zero_imputed_4k_13k # 7468811
zero_imputed_8k_13k # 2033977
zero_imputed_gdt # 2324886
percent_imputed_4k_13k # 86.04621
percent_imputed_8k_13k # 12.13445
percent_imputed_gdt # 76.17582

print("CD3E" %in% rownames(zero_imputed_4k_13k)) # FALSE
print("CD3D" %in% rownames(zero_imputed_4k_13k)) # FALSE
print("TRDC" %in% rownames(zero_imputed_4k_13k)) # TRUE
print("TRGC1" %in% rownames(zero_imputed_4k_13k)) # TRUE
print("TRGC2" %in% rownames(zero_imputed_4k_13k)) # TRUE
print("CD8A" %in% rownames(zero_imputed_4k_13k)) # TRUE
print("CD8B" %in% rownames(zero_imputed_4k_13k)) # TRUE