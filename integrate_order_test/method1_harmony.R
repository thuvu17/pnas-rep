# Harmony integration pipeline
# Method 1: Aligning 3 gdt first, then align that with the pbmc datasets

# !Can't really do method 1 on Harmony because it only modifies the embeddings!!
library(Seurat)
library(harmony)


# Merge gdt and pbmc datasets
harmony.method1 <- merge(x=gdt_donor1, y=list(gdt_donor2, gdt_donor3, pbmc_4k, pbmc_8k))
Layers(harmony.method1[["RNA"]])

# Run standard analysis workflow
harmony.method1 <- NormalizeData(harmony.method1)
harmony.method1 <- FindVariableFeatures(harmony.method1)
harmony.method1 <- ScaleData(harmony.method1)
harmony.method1 <- RunPCA(harmony.method1)

# Run Harmony with theta
harmony.method1 <- RunHarmony(harmony.method1, group.by.vars = c("orig.ident", "celltype"), theta = c(2, 0), verbose = FALSE)

# Run tSNE for visualization
harmony.method1 <- RunTSNE(harmony.method1, reduction = "harmony", dims = 1:30)
method1_tsne_plot_integrated_harmony <- DimPlot(harmony.method1, reduction = "tsne", group.by = "orig.ident", label = TRUE) + ggtitle("Method_1 tSNE integrated Harmony ")
# tsne_plot_unintegrated + method1_tsne_plot_integrated_harmony
