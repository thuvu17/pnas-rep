# Using Seurat V4.3 integration pipeline

# Method 1: Align 3 donors first, then align it with pbmc
# 1.1 Align 3 donors using Seurat V5.2 integration pipeline
gdt_donor1_seuratv4 <- subset(gdt, subset = orig.ident == "donor1")
gdt_donor2_seuratv4 <- subset(gdt, subset = orig.ident == "donor2")
gdt_donor3_seuratv4 <- subset(gdt, subset = orig.ident == "donor3")
# 3 seurat objects split by orig.ident
datasets_gdt3d_seuratv4 <- list(gdt_donor1_seuratv4, gdt_donor2_seuratv4, gdt_donor3_seuratv4)

# normalize and identify variable features for each dataset independently
datasets_gdt3d_seuratv4 <- lapply(X = datasets_gdt3d_seuratv4, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features_gdt3d <- SelectIntegrationFeatures(object.list = datasets_gdt3d_seuratv4)

# ----------------------
# Perform integration
seuratv4_gdt3d.anchors <- FindIntegrationAnchors(object.list = datasets_gdt3d_seuratv4, anchor.features = features_gdt3d)
# this command creates an 'integrated' data assay
seuratv4_gdt3d.combined <- IntegrateData(anchorset = seuratv4_gdt3d.anchors)

# --------------------------------
# Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seuratv4_gdt3d.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seuratv4_gdt3d.combined <- ScaleData(seuratv4_gdt3d.combined, verbose = FALSE)
seuratv4_gdt3d.combined <- RunPCA(seuratv4_gdt3d.combined, npcs = 30, verbose = FALSE)
seuratv4_gdt3d.combined <- RunTSNE(seuratv4_gdt3d.combined, reduction = "pca", dims = 1:30)
seuratv4_gdt3d.combined <- FindNeighbors(seuratv4_gdt3d.combined, reduction = "pca", dims = 1:30)
seuratv4_gdt3d.combined <- FindClusters(seuratv4_gdt3d.combined, resolution = 0.5)

# Visualization
# tsne_plot_integrated_gdt_v4 <- DimPlot(seuratv4_gdt3d.combined, reduction = "tsne", group.by = "orig.ident") + ggtitle("tSNE gdt integrated Seurat v4")

# ------------------------------------------------------------------------
# 1.2 Align merged_gdt with pbmcs using Seurat V5.2 integration pipeline
pbmc_4k_method1_seuratv4 <- pbmc_4k
pbmc_8k_method1_seuratv4 <- pbmc_8k
pbmc_4k_method1_seuratv4 <- NormalizeData(pbmc_4k_method1_seuratv4)
pbmc_4k_method1_seuratv4 <- FindVariableFeatures(pbmc_4k_method1_seuratv4, selection.method = "vst", nfeatures = 2000)
pbmc_8k_method1_seuratv4 <- NormalizeData(pbmc_8k_method1_seuratv4)
pbmc_8k_method1_seuratv4 <- FindVariableFeatures(pbmc_8k_method1_seuratv4, selection.method = "vst", nfeatures = 2000)
# 3 seurat objects split by orig.ident
datasets_method1_seuratv4 <- list(seuratv4_gdt3d.combined, pbmc_4k_method1_seuratv4, pbmc_8k_method1_seuratv4)

# select features that are repeatedly variable across datasets for integration
features_method1 <- SelectIntegrationFeatures(object.list = datasets_method1_seuratv4)

# ----------------------
# Perform integration
seuratv4_method1.anchors <- FindIntegrationAnchors(object.list = datasets_method1_seuratv4, anchor.features = features_method1)
# this command creates an 'integrated' data assay
seuratv4.method1 <- IntegrateData(anchorset = seuratv4_method1.anchors)

# --------------------------------
# Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seuratv4.method1) <- "integrated"

# Run the standard workflow for visualization and clustering
seuratv4.method1 <- ScaleData(seuratv4.method1, verbose = FALSE)
seuratv4.method1 <- RunPCA(seuratv4.method1, npcs = 30, verbose = FALSE)
seuratv4.method1 <- RunTSNE(seuratv4.method1, reduction = "pca", dims = 1:30)
seuratv4.method1 <- FindNeighbors(seuratv4.method1, reduction = "pca", dims = 1:30)
seuratv4.method1 <- FindClusters(seuratv4.method1, resolution = 0.5)
# Visualization
# method1_tsne_plot_integrated_v4 <- DimPlot(seuratv4.method1, reduction = "tsne", group.by = "orig.ident") + ggtitle("Method_1 tSNE integrated Seurat v4")
