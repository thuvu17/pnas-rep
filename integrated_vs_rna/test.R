# Recalculate everything cleanly
pbmc_4k.mean_rna_12k <- Matrix::rowMeans(pbmc_4k.rna_12k)
pbmc_4k.mean_int_12k <- Matrix::rowMeans(pbmc_4k.int_12k)
pbmc_4k.mean_delta_12k <- pbmc_4k.mean_int_12k - pbmc_4k.mean_rna_12k

# Get top 10 delta genes by ABSOLUTE value
top10_genes <- names(sort(abs(pbmc_4k.mean_delta_12k), decreasing = TRUE))[1:10]

# Safely align values by gene name using named indexing
rna_vals <- pbmc_4k.mean_rna_12k[top10_genes]
int_vals <- pbmc_4k.mean_int_12k[top10_genes]
delta_vals <- int_vals - rna_vals  # OR use pbmc_4k.mean_delta_12k[top10_genes] directly

# Make sure all names align
stopifnot(identical(names(rna_vals), names(int_vals)))
stopifnot(identical(names(rna_vals), names(delta_vals)))

# Construct final data frame
pbmc_4k.mean_df_top10 <- data.frame(
  Gene = rep(top10_genes, each = 3),
  Assay = rep(c("RNA", "Integrated", "Delta"), times = length(top10_genes)),
  MeanExpression = c(as.numeric(rna_vals), as.numeric(int_vals), as.numeric(delta_vals))
)
