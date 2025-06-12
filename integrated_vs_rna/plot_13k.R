library(reshape2)


# Ensure common genes and cell names are matched
common_genes_4k <- intersect(rownames(pbmc_4k.int_markersubset), rownames(pbmc_4k.rna_markersubset))
common_cells_4k <- intersect(colnames(pbmc_4k.int_markersubset), colnames(pbmc_4k.rna_markersubset))
common_genes_8k <- intersect(rownames(pbmc_8k.int_markersubset), rownames(pbmc_8k.rna_markersubset))
common_cells_8k <- intersect(colnames(pbmc_8k.int_markersubset), colnames(pbmc_8k.rna_markersubset))
common_genes_gdt <- intersect(rownames(gdt_donor1.int_markersubset), rownames(gdt_donor1.rna_markersubset))
common_cells_gdt <- intersect(colnames(gdt_donor1.int_markersubset), colnames(gdt_donor1.rna_markersubset))

# Subset and order both matrices identically
int_mat_4k <- as.matrix(pbmc_4k.int_markersubset[common_genes_4k, common_cells_4k, drop = FALSE])
rna_mat_4k <- as.matrix(pbmc_4k.rna_markersubset[common_genes_4k, common_cells_4k, drop = FALSE])
int_mat_8k <- as.matrix(pbmc_8k.int_markersubset[common_genes_8k, common_cells_8k, drop = FALSE])
rna_mat_8k <- as.matrix(pbmc_8k.rna_markersubset[common_genes_8k, common_cells_8k, drop = FALSE])
int_mat_gdt <- as.matrix(gdt_donor1.int_markersubset[common_genes_gdt, common_cells_gdt, drop = FALSE])
rna_mat_gdt <- as.matrix(gdt_donor1.rna_markersubset[common_genes_gdt, common_cells_gdt, drop = FALSE])

# Delta
delta_4k <- int_mat_4k - rna_mat_4k
delta_8k <- int_mat_8k - rna_mat_8k
delta_gdt <- int_mat_gdt - rna_mat_gdt

df_delta_4k <- reshape2::melt(delta_4k)
df_delta_4k$Assay <- "Delta"
df_delta_4k$Group <- "PBMC_4k"

df_delta_8k <- reshape2::melt(delta_8k)
df_delta_8k$Assay <- "Delta"
df_delta_8k$Group <- "PBMC_8k"

df_delta_gdt <- reshape2::melt(delta_gdt)
df_delta_gdt$Assay <- "Delta"
df_delta_gdt$Group <- "GDT"

# RNA
df_rna_4k <- reshape2::melt(as.matrix(pbmc_4k.rna_markersubset))
df_rna_4k$Assay <- "RNA"
df_rna_4k$Group <- "PBMC_4k"

df_rna_8k <- reshape2::melt(as.matrix(pbmc_8k.rna_markersubset))
df_rna_8k$Assay <- "RNA"
df_rna_8k$Group <- "PBMC_8k"

df_rna_gdt <- reshape2::melt(as.matrix(gdt_donor1.rna_markersubset))
df_rna_gdt$Assay <- "RNA"
df_rna_gdt$Group <- "GDT"

# Integrated
df_int_4k <- reshape2::melt(as.matrix(pbmc_4k.int_markersubset))
df_int_4k$Assay <- "Integrated"
df_int_4k$Group <- "PBMC_4k"

df_int_8k <- reshape2::melt(as.matrix(pbmc_8k.int_markersubset))
df_int_8k$Assay <- "Integrated"
df_int_8k$Group <- "PBMC_8k"

df_int_gdt <- reshape2::melt(as.matrix(gdt_donor1.int_markersubset))
df_int_gdt$Assay <- "Integrated"
df_int_gdt$Group <- "GDT"

# Combine
marker_plot_df <- rbind(df_rna_4k, df_int_4k,
                        df_rna_8k, df_int_8k,
                        df_rna_gdt, df_int_gdt)
colnames(marker_plot_df)[1:2] <- c("Gene", "Cell")

# Plot
ggplot(marker_plot_df, aes(x = Gene, y = value, fill = Assay)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.8)) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(title = "Marker Gene Expression in PBMC13k (8000 genes)",
       y = "Expression",
       x = "Marker Gene") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ===============================
cd_genes <- c("CD8A", "CD8B")

int_mat_4k_cd <- as.matrix(pbmc_4k.int_markersubset[cd_genes, common_cells_4k, drop = FALSE])
rna_mat_4k_cd <- as.matrix(pbmc_4k.rna_markersubset[cd_genes, common_cells_4k, drop = FALSE])
int_mat_8k_cd <- as.matrix(pbmc_8k.int_markersubset[cd_genes, common_cells_8k, drop = FALSE])
rna_mat_8k_cd <- as.matrix(pbmc_8k.rna_markersubset[cd_genes, common_cells_8k, drop = FALSE])
int_mat_gdt_cd <- as.matrix(gdt_donor1.int_markersubset[cd_genes, common_cells_gdt, drop = FALSE])
rna_mat_gdt_cd <- as.matrix(gdt_donor1.rna_markersubset[cd_genes, common_cells_gdt, drop = FALSE])

# Delta
delta_4k_cd <- int_mat_4k_cd - rna_mat_4k_cd
delta_8k_cd <- int_mat_8k_cd - rna_mat_8k_cd
delta_gdt_cd <- int_mat_gdt_cd - rna_mat_gdt_cd

df_delta_4k_cd <- reshape2::melt(delta_4k_cd)
df_delta_4k_cd$Assay <- "Delta"
df_delta_4k_cd$Group <- "PBMC_4k"

df_delta_8k_cd <- reshape2::melt(delta_8k_cd)
df_delta_8k_cd$Assay <- "Delta"
df_delta_8k_cd$Group <- "PBMC_8k"

df_delta_gdt_cd <- reshape2::melt(delta_gdt_cd)
df_delta_gdt_cd$Assay <- "Delta"
df_delta_gdt_cd$Group <- "GDT"

# RNA
df_rna_4k_cd <- reshape2::melt(as.matrix(pbmc_4k.rna_markersubset[cd_genes, , drop = FALSE]))
df_rna_4k_cd$Assay <- "RNA"
df_rna_4k_cd$Group <- "PBMC_4k"

df_rna_8k_cd <- reshape2::melt(as.matrix(pbmc_8k.rna_markersubset[cd_genes, , drop = FALSE]))
df_rna_8k_cd$Assay <- "RNA"
df_rna_8k_cd$Group <- "PBMC_8k"

df_rna_gdt_cd <- reshape2::melt(as.matrix(gdt_donor1.rna_markersubset[cd_genes, , drop = FALSE]))
df_rna_gdt_cd$Assay <- "RNA"
df_rna_gdt_cd$Group <- "GDT"

# Integrated
df_int_4k_cd <- reshape2::melt(as.matrix(pbmc_4k.int_markersubset[cd_genes, , drop = FALSE]))
df_int_4k_cd$Assay <- "Integrated"
df_int_4k_cd$Group <- "PBMC_4k"

df_int_8k_cd <- reshape2::melt(as.matrix(pbmc_8k.int_markersubset[cd_genes, , drop = FALSE]))
df_int_8k_cd$Assay <- "Integrated"
df_int_8k_cd$Group <- "PBMC_8k"

df_int_gdt_cd <- reshape2::melt(as.matrix(gdt_donor1.int_markersubset[cd_genes, , drop = FALSE]))
df_int_gdt_cd$Assay <- "Integrated"
df_int_gdt_cd$Group <- "GDT"

# Combine
marker_plot_df_cd <- rbind(df_rna_4k_cd, df_int_4k_cd,
                        df_rna_8k_cd, df_int_8k_cd,
                        df_rna_gdt_cd, df_int_gdt_cd)
colnames(marker_plot_df_cd)[1:2] <- c("Gene", "Cell")

# Plot
ggplot(marker_plot_df_cd, aes(x = Gene, y = value, fill = Assay)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.8)) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(title = "Marker Gene Expression in PBMC13k (CD8 only) using 8000 genes",
       y = "Expression",
       x = "Marker Gene") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ===============================
# Mean box plots
# 1. Compute gene-wise means
pbmc_4k.mean_rna <- Matrix::rowMeans(pbmc_4k.rna_13k)
pbmc_4k.mean_int <- Matrix::rowMeans(pbmc_4k.int_13k)
pbmc_4k.mean_delta <- pbmc_4k.mean_int - pbmc_4k.mean_rna

pbmc_8k.mean_rna <- Matrix::rowMeans(pbmc_8k.rna_13k)
pbmc_8k.mean_int <- Matrix::rowMeans(pbmc_8k.int_13k)
pbmc_8k.mean_delta <- pbmc_8k.mean_int - pbmc_8k.mean_rna

gdt.mean_rna <- Matrix::rowMeans(gdt_donor1.rna_13k)
gdt.mean_int <- Matrix::rowMeans(gdt_donor1.int_13k)
gdt.mean_delta <- gdt.mean_int - gdt.mean_rna

# 2. Combine into long-format dataframe
df_all <- rbind(
  data.frame(Gene = names(pbmc_4k.mean_rna), Mean = pbmc_4k.mean_rna, Assay = "RNA", Group = "PBMC_4k"),
  data.frame(Gene = names(pbmc_4k.mean_int), Mean = pbmc_4k.mean_int, Assay = "Integrated", Group = "PBMC_4k"),
  data.frame(Gene = names(pbmc_4k.mean_delta), Mean = pbmc_4k.mean_delta, Assay = "Delta", Group = "PBMC_4k"),
  
  data.frame(Gene = names(pbmc_8k.mean_rna), Mean = pbmc_8k.mean_rna, Assay = "RNA", Group = "PBMC_8k"),
  data.frame(Gene = names(pbmc_8k.mean_int), Mean = pbmc_8k.mean_int, Assay = "Integrated", Group = "PBMC_8k"),
  data.frame(Gene = names(pbmc_8k.mean_delta), Mean = pbmc_8k.mean_delta, Assay = "Delta", Group = "PBMC_8k"),
  
  data.frame(Gene = names(gdt.mean_rna), Mean = gdt.mean_rna, Assay = "RNA", Group = "GDT"),
  data.frame(Gene = names(gdt.mean_int), Mean = gdt.mean_int, Assay = "Integrated", Group = "GDT"),
  data.frame(Gene = names(gdt.mean_delta), Mean = gdt.mean_delta, Assay = "Delta", Group = "GDT")
)

# 3. Plot boxplot
ggplot(df_all, aes(x = Assay, y = Mean, fill = Assay)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(title = "Gene-Wise Mean Expression Across Assays", y = "Mean Expression") +
  theme_minimal(base_size = 13)
