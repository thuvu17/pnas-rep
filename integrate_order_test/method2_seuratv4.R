# Using Seurat V4.3 integration pipeline:
# Method 2: Align as 5 separate datasets using Seurat V5.2 integration pipeline


# 3 seurat objects split by orig.ident
gdt_method2_seuratv4 <- gdt
pbmc_4k_method2_seuratv4 <- pbmc_4k
pbmc_8k_method2_seuratv4 <- pbmc_8k
datasets_method2_seuratv4 <- list(gdt_method2_seuratv4, pbmc_4k_method2_seuratv4, pbmc_8k_method2_seuratv4)

# normalize and identify variable features for each dataset independently
datasets_method2_seuratv4 <- lapply(X = datasets_method2_seuratv4, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features_method2 <- SelectIntegrationFeatures(object.list = datasets_method2_seuratv4)

# ----------------------
# Perform integration
method2_seuratv4.anchors <- FindIntegrationAnchors(object.list = datasets_method2_seuratv4, anchor.features = features_method2)
# this command creates an 'integrated' data assay
seuratv4.method2 <- IntegrateData(anchorset = method2_seuratv4.anchors)

# --------------------------------
# Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seuratv4.method2) <- "integrated"

# Run the standard workflow for visualization and clustering
seuratv4.method2 <- ScaleData(seuratv4.method2, verbose = FALSE)
seuratv4.method2 <- RunPCA(seuratv4.method2, npcs = 30, verbose = FALSE)
seuratv4.method2 <- RunTSNE(seuratv4.method2, reduction = "pca", dims = 1:30)
seuratv4.method2 <- FindNeighbors(seuratv4.method2, reduction = "pca", dims = 1:30)
seuratv4.method2 <- FindClusters(seuratv4.method2, resolution = 0.5)

# Visualization
method2_tsne_plot_integrated_v4 <- DimPlot(seuratv4.method2, reduction = "tsne", group.by = "orig.ident") + ggtitle("Method_2 tSNE integrated Seurat v4")
method2_tsne_plot_integrated_v5 + method2_tsne_plot_integrated_v4