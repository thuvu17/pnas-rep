library(ggplot2)
library(reshape2)
library(pheatmap)


# Box plots
# 1. Compute gene-wise means
pbmc_4k.mean_rna_12k <- Matrix::rowMeans(pbmc_4k.rna_12k)
pbmc_4k.mean_int_12k <- Matrix::rowMeans(pbmc_4k.int_12k)
pbmc_4k.mean_delta_12k <- pbmc_4k.mean_int_12k - pbmc_4k.mean_rna_12k

pbmc_8k.mean_rna_12k <- Matrix::rowMeans(pbmc_8k.rna_12k)
pbmc_8k.mean_int_12k <- Matrix::rowMeans(pbmc_8k.int_12k)
pbmc_8k.mean_delta_12k <- pbmc_8k.mean_int_12k - pbmc_8k.mean_rna_12k
# -----------------------------------------
# 2. Combine into a long-format data frame
pbmc_4k.mean_df <- data.frame(
  Gene = names(pbmc_4k.mean_rna_12k),
  RNA = pbmc_4k.mean_rna_12k,
  Integrated = pbmc_4k.mean_int_12k,
  Delta = pbmc_4k.mean_delta_12k
)
pbmc_4k.mean_df_melted <- reshape2::melt(pbmc_4k.mean_df, id.vars = "Gene", variable.name = "Assay", value.name = "MeanExpression")

pbmc_8k.mean_df <- data.frame(
  Gene = names(pbmc_8k.mean_rna_12k),
  RNA = pbmc_8k.mean_rna_12k,
  Integrated = pbmc_8k.mean_int_12k,
  Delta = pbmc_8k.mean_delta_12k
)
pbmc_8k.mean_df_melted <- reshape2::melt(pbmc_8k.mean_df, id.vars = "Gene", variable.name = "Assay", value.name = "MeanExpression")
# -----------------------------------------
# 3. Plot
ggplot(pbmc_4k.mean_df_melted, aes(x = Assay, y = MeanExpression, fill = Assay)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "PBMC 4k in 12k: Mean Gene Expression",
       y = "Mean Expression per Gene") +
  theme_minimal(base_size = 14)

ggplot(pbmc_8k.mean_df_melted, aes(x = Assay, y = MeanExpression, fill = Assay)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "PBMC 8k in 12k: Mean Gene Expression",
       y = "Mean Expression per Gene") +
  theme_minimal(base_size = 14)

# ===================================
# 10 largest delta
# Step 1: Get gene names by delta
top10_genes_4k_12k <- names(sort(abs(pbmc_4k.mean_delta_12k), decreasing = TRUE))[1:10]
top10_genes_8k_12k <- names(sort(abs(pbmc_4k.mean_delta_12k), decreasing = TRUE))[1:10]

# Step 2: Create clean long-format data frame using melt
pbmc_4k.mean_df <- data.frame(
  Gene = names(pbmc_4k.mean_rna_12k),
  RNA = pbmc_4k.mean_rna_12k,
  Integrated = pbmc_4k.mean_int_12k,
  Delta = pbmc_4k.mean_int_12k - pbmc_4k.mean_rna_12k
)
pbmc_4k.mean_df_top10 <- pbmc_4k.mean_df[pbmc_4k.mean_df$Gene %in% top10_genes_4k_12k, ]

pbmc_8k.mean_df <- data.frame(
  Gene = names(pbmc_8k.mean_rna_12k),
  RNA = pbmc_8k.mean_rna_12k,
  Integrated = pbmc_8k.mean_int_12k,
  Delta = pbmc_8k.mean_int_12k - pbmc_8k.mean_rna_12k
)
pbmc_8k.mean_df_top10 <- pbmc_8k.mean_df[pbmc_8k.mean_df$Gene %in% top10_genes_8k_12k, ]

# Optional: reorder Gene factor to match delta order
pbmc_4k.mean_df_top10$Gene <- factor(pbmc_4k.mean_df_top10$Gene, levels = top10_genes_4k_12k)
pbmc_4k.mean_df_melted_top10 <- reshape2::melt(pbmc_4k.mean_df_top10,
                                               id.vars = "Gene",
                                               variable.name = "Assay",
                                               value.name = "MeanExpression")

pbmc_8k.mean_df_top10$Gene <- factor(pbmc_8k.mean_df_top10$Gene, levels = top10_genes_8k_12k)
pbmc_8k.mean_df_melted_top10 <- reshape2::melt(pbmc_8k.mean_df_top10,
                                               id.vars = "Gene",
                                               variable.name = "Assay",
                                               value.name = "MeanExpression")
# Step 3: Plot
ggplot(pbmc_4k.mean_df_melted_top10, aes(x = Gene, y = MeanExpression, fill = Assay)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Top 10 Delta Genes in PBMC 4k (12k): Mean Expression by Assay",
       y = "Mean Expression", x = "Gene") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(pbmc_8k.mean_df_melted_top10, aes(x = Gene, y = MeanExpression, fill = Assay)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Top 10 Delta Genes in PBMC 8k (12k): Mean Expression by Assay",
       y = "Mean Expression", x = "Gene") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
