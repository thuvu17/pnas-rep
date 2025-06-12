# Using Seurat V5.2 data integration pipeline

# Method 1: Align 3 donors first, then align it with pbmc
# 1.1 Align 3 donors using Seurat V5.2 integration pipeline
merged_gdt <- merge(x=gdt_donor1, y=list(gdt_donor2, gdt_donor3))
# view(merged_gdt)
# class(merged_gdt[["RNA"]])

# Split the RNA measurements into layers (Already split?)
# merged_gdt[["RNA"]] <- split(merged_gdt[["RNA"]], f=merged_gdt$orig.ident)
Layers(merged_gdt[["RNA"]])

# Run standard analysis workflow
merged_gdt <- NormalizeData(merged_gdt)
merged_gdt <- FindVariableFeatures(merged_gdt)
merged_gdt <- ScaleData(merged_gdt)
merged_gdt <- RunPCA(merged_gdt)
merged_gdt <- FindNeighbors(merged_gdt, dims = 1:30, reduction = "pca")
merged_gdt <- FindClusters(merged_gdt, resolution = 2, cluster.name = "unintegrated_clusters")

# Visualization
# merged_gdt <- RunTSNE(merged_gdt, dims = 1:30, reduction = "pca", , reduction.name = "tsne.unintegrated")
# gdt_tsne_plot_unintegrated_v5 <- DimPlot(merged_gdt, reduction = "tsne.unintegrated", group.by = "orig.ident") + ggtitle("Method_1 tsne unintegrated Seurat v5")
# -------------------
# Perform integration
merged_gdt <- IntegrateLayers(object = merged_gdt, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                   verbose = FALSE)

# re-join layers after integration
merged_gdt[["RNA"]] <- JoinLayers(merged_gdt[["RNA"]])

merged_gdt <- FindNeighbors(merged_gdt, reduction = "integrated.cca", dims = 1:30)
merged_gdt <- FindClusters(merged_gdt, resolution = 1)

# Visualization
# merged_gdt <- RunTSNE(merged_gdt, dims = 1:30, reduction = "integrated.cca", reduction.name = "tsne.integrated")
# tsne_plot_integrated_gdt_v5 <- DimPlot(merged_gdt, reduction = "tsne.integrated", group.by = "orig.ident") + ggtitle("tSNE gdt integrated Seurat v5")
# tsne_plot_unintegrated_gdt_v5 + tsne_plot_integrated_gdt_v5

# ------------------------------------------------------------------------
# 1.2 Align merged_gdt with pbmcs using Seurat V5.2 integration pipeline
seuratv5.method1 <- merge(x=merged_gdt, y=list(pbmc_4k, pbmc_8k))
# view(merged_gdt)
# class(merged_gdt[["RNA"]])

# Split the RNA measurements into layers (Already split?)
# merged_gdt[["RNA"]] <- split(merged_gdt[["RNA"]], f=merged_gdt$orig.ident)
Layers(seuratv5.method1[["RNA"]])

# Run standard analysis workflow
seuratv5.method1 <- NormalizeData(seuratv5.method1)
seuratv5.method1 <- FindVariableFeatures(seuratv5.method1)
seuratv5.method1 <- ScaleData(seuratv5.method1)
seuratv5.method1 <- RunPCA(seuratv5.method1)
seuratv5.method1 <- FindNeighbors(seuratv5.method1, dims = 1:30, reduction = "pca")
seuratv5.method1 <- FindClusters(seuratv5.method1, resolution = 2, cluster.name = "unintegrated_clusters")

# -------------------
# Perform integration
seuratv5.method1 <- IntegrateLayers(object = seuratv5.method1, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                              verbose = FALSE)

# re-join layers after integration
seuratv5.method1[["RNA"]] <- JoinLayers(seuratv5.method1[["RNA"]])

seuratv5.method1 <- FindNeighbors(seuratv5.method1, reduction = "integrated.cca", dims = 1:30)
seuratv5.method1 <- FindClusters(seuratv5.method1, resolution = 1)

# Visualization
seuratv5.method1 <- RunTSNE(seuratv5.method1, dims = 1:30, reduction = "integrated.cca", reduction.name = "tsne.integrated")
# method1_tsne_plot_integrated_v5 <- DimPlot(seuratv5.method1, reduction = "tsne.integrated", group.by = "orig.ident") + ggtitle("Method_1 tsne integrated Seurat v5")
