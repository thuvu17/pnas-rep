# Box plot
pbmc_4k.rna_subset <- LayerData(part1.subset, assay = "RNA", layer = "data.subset_4k")
pbmc_8k.rna_subset <- LayerData(part1.subset, assay = "RNA", layer = "data.subset_8k")
gdt_donor1.rna_subset <- LayerData(part1.subset, assay = "RNA", layer = "data.gdt")

pbmc_4k.rna_markersubset <- pbmc_4k.rna_subset[marker_genes, , drop = FALSE]
pbmc_8k.rna_markersubset <- pbmc_8k.rna_subset[marker_genes, , drop = FALSE]
gdt_donor1.rna_markersubset <- gdt_donor1.rna_subset[marker_genes, , drop = FALSE]

# From integrated assay
cells_4k_subset <- WhichCells(part1.subset, expression = orig.ident == "pbmc4k")
cells_8k_subset <- WhichCells(part1.subset, expression = orig.ident == "pbmc8k")
cells_gdt_subset <- WhichCells(part1.subset, expression = orig.ident == "donor1")

pbmc13k.int_expr_subset <- GetAssayData(part1.subset, assay = "integrated", slot = "data")
pbmc_4k.int_subset <- pbmc13k.int_expr_subset[, cells_4k_subset]
pbmc_8k.int_subset <- pbmc13k.int_expr_subset[, cells_8k_subset]
gdt_donor1.int_subset <- pbmc13k.int_expr_subset[, cells_gdt_subset]

pbmc_4k.int_markersubset <- pbmc_4k.int_subset[marker_genes, , drop = FALSE]
pbmc_8k.int_markersubset <- pbmc_8k.int_subset[marker_genes, , drop = FALSE]
gdt_donor1.int_markersubset <- gdt_donor1.int_subset[marker_genes, , drop = FALSE]

# Define delog function
delog_matrix <- function(mat) {
  exp(mat) - 1
}

# Apply to all marker subsets (RNA)
pbmc_4k.rna_markersubset <- delog_matrix(pbmc_4k.rna_markersubset)
pbmc_8k.rna_markersubset <- delog_matrix(pbmc_8k.rna_markersubset)
gdt_donor1.rna_markersubset <- delog_matrix(gdt_donor1.rna_markersubset)

# Apply to all marker subsets (Integrated)
pbmc_4k.int_markersubset <- delog_matrix(pbmc_4k.int_markersubset)
pbmc_8k.int_markersubset <- delog_matrix(pbmc_8k.int_markersubset)
gdt_donor1.int_markersubset <- delog_matrix(gdt_donor1.int_markersubset)

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
  labs(title = "Iteration 3",
       y = "Expression",
       x = "Marker Gene") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  coord_cartesian(ylim = c(0, 20))

