library(readr)


# Read cluster analysis
cluster.4k <- read_csv("data/cluster_analysis/analysis_4k/clustering/graphclust/clusters.csv")
cluster.8k <- read_csv("data/cluster_analysis/analysis_8k/clustering/graphclust/clusters.csv")
cluster.4k$Barcode <- paste0("pbmc4k_", cluster.4k$Barcode)
cluster.8k$Barcode <- paste0("pbmc8k_", cluster.8k$Barcode)

# Identify NK cells in part1.v4
part1.pbmc4k <- subset(part1.v4, cells = cluster.4k$Barcode)
part1.pbmc8k <- subset(part1.v4, cells = cluster.8k$Barcode)
Idents(pbmc_4k) <- cluster.4k$Cluster[match(Cells(pbmc_4k), cluster.4k$Barcode)]
Idents(pbmc_8k) <- cluster.8k$Cluster[match(Cells(pbmc_8k), cluster.8k$Barcode)]
# NK cluster
VlnPlot(pbmc_4k, features = c("NKG7", "CD3E"), pt.size = 0) + ggtitle("CD3E_4k")
VlnPlot(pbmc_8k, features = c("NKG7", "CD3E"), pt.size = 0) + ggtitle("CD3E_8k")
nk_cluster.4k <- 8
nk_cluster.8k <- 6 # (and 11?)
# Cytotoxic CD8 T cell cluster
VlnPlot(pbmc_4k, features = c("CD3E", "CD8A", "FCGR3A", "GNLY"), pt.size = 0) + ggtitle("GNLY_4k")
VlnPlot(pbmc_8k, features = c("CD3E", "CD8A", "FCGR3A", "GNLY"), pt.size = 0) + ggtitle("GNLY_8k")
cd8_cluster.4k <- 5
cd8_cluster.8k <- 10
# Save NK and CD8 barcodes
nk_barcodes.4k <- cluster.4k$Barcode[cluster.4k$Cluster == nk_cluster.4k]
nk_barcodes.8k <- cluster.8k$Barcode[cluster.8k$Cluster == nk_cluster.8k]
cd8_barcodes.4k <- cluster.4k$Barcode[cluster.4k$Cluster == cd8_cluster.4k]
cd8_barcodes.8k <- cluster.8k$Barcode[cluster.8k$Cluster == cd8_cluster.8k]
# Subset from part1.v4
part1.nk <- c(nk_barcodes.4k, nk_barcodes.8k) # 860 cells
part1.cd8 <- c(cd8_barcodes.4k, cd8_barcodes.8k) # 716 cells


# Add tsne embeddings
tsne_df.4k <- read_csv("data/cluster_analysis/analysis_4k/tsne/2_components/projection.csv")
tsne_df.8k <- read_csv("data/cluster_analysis/analysis_8k/tsne/2_components/projection.csv")
tsne_df.4k$Barcode <- paste0("pbmc4k_", tsne_df.4k$Barcode)
tsne_df.8k$Barcode <- paste0("pbmc8k_", tsne_df.8k$Barcode)
# 4k
tsne_mat.4k <- as.matrix(tsne_df.4k[, c("TSNE-1", "TSNE-2")])
rownames(tsne_mat.4k) <- tsne_df.4k$Barcode
colnames(tsne_mat.4k) <- c("tSNE_1", "tSNE_2")
tsne_mat.4k <- tsne_mat.4k[intersect(rownames(tsne_mat.4k), colnames(pbmc_4k)), ]
pbmc_4k[["tsne"]] <- CreateDimReducObject(
  embeddings = tsne_mat.4k,
  key = "tSNE_",
  assay = DefaultAssay(pbmc_4k)
)

DimPlot(pbmc_4k, reduction = "tsne") + ggtitle("PBMC 4k tSNE by cluster") + FeaturePlot(pbmc_4k, reduction = "tsne", features = c("CD3E", "CD8A", "NKG7", "FCGR3A", "GNLY")) + ggtitle("GNLY_4k")
# 8k
tsne_mat.8k <- as.matrix(tsne_df.8k[, c("TSNE-1", "TSNE-2")])
rownames(tsne_mat.8k) <- tsne_df.8k$Barcode
colnames(tsne_mat.8k) <- c("tSNE_1", "tSNE_2")
tsne_mat.8k <- tsne_mat.8k[intersect(rownames(tsne_mat.8k), colnames(pbmc_8k)), ]
pbmc_8k[["tsne"]] <- CreateDimReducObject(
  embeddings = tsne_mat.8k,
  key = "tSNE_",
  assay = DefaultAssay(pbmc_8k)
)
DimPlot(pbmc_8k, reduction = "tsne") + ggtitle("PBMC 8k tSNE by cluster")+ FeaturePlot(pbmc_8k, reduction = "tsne", features = c("CD3E", "CD8A", "NKG7", "FCGR3A", "GNLY")) + ggtitle("GNLY_8k")

