# USing Seurat V5.2 data integration pipeline

# Method 2: Align as 5 separate datasets using Seurat V5.2 integration pipeline
seuratv5.method2 <- merge(x=gdt_donor1, y=list(gdt_donor2, gdt_donor3, pbmc_4k, pbmc_8k))
# view(merged_gdt)
# class(merged_gdt[["RNA"]])

# Split the RNA measurements into layers (Already split?)
# merged_gdt[["RNA"]] <- split(merged_gdt[["RNA"]], f=merged_gdt$orig.ident)
Layers(seuratv5.method2[["RNA"]])

# Run standard analysis workflow
seuratv5.method2 <- NormalizeData(seuratv5.method2)
seuratv5.method2 <- FindVariableFeatures(seuratv5.method2)
seuratv5.method2 <- ScaleData(seuratv5.method2)
seuratv5.method2 <- RunPCA(seuratv5.method2)
seuratv5.method2 <- FindNeighbors(seuratv5.method2, dims = 1:30, reduction = "pca")
seuratv5.method2 <- FindClusters(seuratv5.method2, resolution = 2, cluster.name = "unintegrated_clusters")

# Visualization
seuratv5.method2 <- RunTSNE(seuratv5.method2, dims = 1:30, reduction = "pca", , reduction.name = "tsne.unintegrated")
method2_tsne_plot_unintegrated_v5 <- DimPlot(seuratv5.method2, reduction = "tsne.unintegrated", group.by = "orig.ident") + ggtitle("Method_2 tsne unintegrated Seurat v5")
# -------------------
# Perform integration
seuratv5.method2 <- IntegrateLayers(object = seuratv5.method2, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                  verbose = FALSE)

# re-join layers after integration
seuratv5.method2[["RNA"]] <- JoinLayers(seuratv5.method2[["RNA"]])

seuratv5.method2 <- FindNeighbors(seuratv5.method2, reduction = "integrated.cca", dims = 1:30)
seuratv5.method2 <- FindClusters(seuratv5.method2, resolution = 1)

# Visualization
seuratv5.method2 <- RunTSNE(seuratv5.method2, dims = 1:30, reduction = "integrated.cca", reduction.name = "tsne.integrated")
method2_tsne_plot_integrated_v5 <- DimPlot(seuratv5.method2, reduction = "tsne.integrated", group.by = "orig.ident") + ggtitle("Method_2 tsne integrated Seurat v5")
method2_tsne_plot_unintegrated_v5 + method2_tsne_plot_integrated_v5
tsne_plot_integrated_v5 + method2_tsne_plot_integrated_v5