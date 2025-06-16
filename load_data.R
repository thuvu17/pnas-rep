library(data.table)
library(Seurat)
library(tidyverse)
library(patchwork)


# Load data
gdt_df <- fread("data/GSE128223_3donors_d1_d2_v2.tsv")
pbmc_4k <- Read10X((data.dir = "data/pbmc_4k/"))
pbmc_8k <- Read10X((data.dir = "data/pbmc_8k/"))

# Convert to matrix
rownames(gdt_df) <- gdt_df$V1  # Set gene names as rownames
gdt_df_no_gene_col <- gdt_df[,-1, with = FALSE]

# Create Seurat pbmc_4kects
pbmc_4k <- CreateSeuratObject(counts = pbmc_4k, project = "pbmc_4k", min.cells = 3, min.features = 200)
pbmc_8k <- CreateSeuratObject(counts = pbmc_8k, project = "pbmc_8k", min.cells = 3, min.features = 200)
gdt <- CreateSeuratObject(counts = as.matrix(gdt_df_no_gene_col), project = "gdt", min.cells = 3, min.features = 200)
rownames(gdt) <- rownames(gdt_df)
remove(gdt_df)
remove(gdt_df_no_gene_col)

# Adding cell types metadata
pbmc_4k$celltype <- "pbmc"
pbmc_8k$celltype <- "pbmc"
gdt$celltype <- "gdt"

# Adding subtypes
gdt.subtypes <- gsub(".*_(d[0-9]+)_.*", "\\1", colnames(gdt))
gdt <- AddMetaData(gdt, metadata = gdt.subtypes, col.name = "subtype")
pbmc_4k$subtype <- "NA"
pbmc_8k$subtype <- "NA"

gdt_donor1 <- subset(gdt, subset = orig.ident == "donor1")
# gdt_donor2 <- subset(gdt, subset = orig.ident == "donor2")
# gdt_donor3 <- subset(gdt, subset = orig.ident == "donor3")

# gdt_d1 <- subset(gdt, subset = subtype == "d1")
# gdt_d2 <- subset(gdt, subset = subtype == "d2")

# Preprocessing
pbmc_4k[["percent.mt"]] <- PercentageFeatureSet(pbmc_4k, pattern = "^MT-")
gene_counts <- LayerData(pbmc_4k[["RNA"]], layer = "counts")
gene_filter <- Matrix::rowSums(gene_counts > 0) >= 3
pbmc_4k <- pbmc_4k[gene_filter, ]
pbmc_4k <- subset(pbmc_4k, subset = nFeature_RNA > 200 & percent.mt >= 0.05)

pbmc_8k[["percent.mt"]] <- PercentageFeatureSet(pbmc_8k, pattern = "^MT-")
gene_counts <- LayerData(pbmc_8k[["RNA"]], layer = "counts")
gene_filter <- Matrix::rowSums(gene_counts > 0) >= 3
pbmc_8k <- pbmc_8k[gene_filter, ]
pbmc_8k <- subset(pbmc_8k, subset = nFeature_RNA > 200 & percent.mt >= 0.05)

gdt_donor1[["percent.mt"]] <- PercentageFeatureSet(gdt_donor1, pattern = "^MT-")
gene_counts <- LayerData(gdt_donor1[["RNA"]], layer = "counts")
gene_filter <- Matrix::rowSums(gene_counts > 0) >= 3
gdt_donor1 <- gdt_donor1[gene_filter, ]
gdt_donor1 <- subset(gdt_donor1, subset = nFeature_RNA > 200 & percent.mt >= 0.05)

pbmc_4k <- RenameCells(pbmc_4k, new.names = paste0("pbmc4k_", colnames(pbmc_4k)))
pbmc_8k <- RenameCells(pbmc_8k, new.names = paste0("pbmc8k_", colnames(pbmc_8k)))

# Add CD4 zeros to gdt_donor1
counts <- GetAssayData(gdt_donor1, assay = "RNA", slot = "counts")
# If CD4 missing, add it
if (!("CD4" %in% rownames(counts))) {
  # Create zero vector for CD4 across all cells
  zero_row <- Matrix::sparseMatrix(i = NULL, j = NULL, dims = c(1, ncol(counts)))
  rownames(zero_row) <- "CD4"
  
  # Add row to count matrix
  counts <- rbind(counts, zero_row)
  
  # Write back to object
  gdt_donor1[["RNA"]] <- CreateAssayObject(counts = counts)
}


# ======================================
# Generate tSNE plots before integration
merged_unintegrated <- merge(x=gdt_donor1, y=list(pbmc_4k, pbmc_8k))
Layers(merged_unintegrated[["RNA"]])

# Run standard analysis workflow
merged_unintegrated <- NormalizeData(merged_unintegrated)
merged_unintegrated <- FindVariableFeatures(merged_unintegrated)
merged_unintegrated <- ScaleData(merged_unintegrated)
merged_unintegrated <- RunPCA(merged_unintegrated)
merged_unintegrated <- FindNeighbors(merged_unintegrated, dims = 1:30, reduction = "pca")
merged_unintegrated <- FindClusters(merged_unintegrated, resolution = 2, cluster.name = "unintegrated_clusters")

# Visualization
merged_unintegrated <- RunTSNE(merged_unintegrated, dims = 1:30, reduction = "pca")
tsne_plot_unintegrated <- DimPlot(merged_unintegrated, reduction = "tsne.unintegrated", group.by = "orig.ident", shuffle=TRUE) + ggtitle("tSNE unintegrated")
tsne_plot_unintegrated
