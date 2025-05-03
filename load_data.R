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

# Create Seurat objects
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
gdt_donor2 <- subset(gdt, subset = orig.ident == "donor2")
gdt_donor3 <- subset(gdt, subset = orig.ident == "donor3")

gdt_d1 <- subset(gdt, subset = subtype == "d1")
gdt_d2 <- subset(gdt, subset = subtype == "d2")

# ======================================
# Generate tSNE plots before integration
merged_unintegrated <- merge(x=gdt, y=list(pbmc_4k, pbmc_8k))
Layers(merged_unintegrated[["RNA"]])

# Run standard analysis workflow
merged_unintegrated <- NormalizeData(merged_unintegrated)
merged_unintegrated <- FindVariableFeatures(merged_unintegrated)
merged_unintegrated <- ScaleData(merged_unintegrated)
merged_unintegrated <- RunPCA(merged_unintegrated)
merged_unintegrated <- FindNeighbors(merged_unintegrated, dims = 1:30, reduction = "pca")
merged_unintegrated <- FindClusters(merged_unintegrated, resolution = 2, cluster.name = "unintegrated_clusters")

# Visualization
merged_unintegrated <- RunTSNE(merged_unintegrated, dims = 1:30, reduction = "pca", reduction.name = "tsne.unintegrated")
tsne_plot_unintegrated <- DimPlot(merged_unintegrated, reduction = "tsne.unintegrated", group.by = "orig.ident", shuffle=TRUE) + ggtitle("tSNE unintegrated")
tsne_plot_unintegrated
