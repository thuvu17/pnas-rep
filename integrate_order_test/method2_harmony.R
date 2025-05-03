# Harmony integration pipeline
# Method 2: Aligning all 5 as separate datasets

# Prepare Seurat objects
harmony.method2 <- merge(gdt_donor1, y = c(gdt_donor2, gdt_donor3, pbmc_4k, pbmc_8k), add.cell.ids = c("donor1", "donor2", "donor3", "pbmc_4k", "pbmc_8k"))

# Find common variable features across datasets
gdt_donor1 <- NormalizeData(gdt_donor1)
gdt_donor1 <- FindVariableFeatures(gdt_donor1)
gdt_donor1 <- ScaleData(gdt_donor1)
gdt_donor1_features <- VariableFeatures(gdt_donor1)

gdt_donor2 <- NormalizeData(gdt_donor2)
gdt_donor2 <- FindVariableFeatures(gdt_donor2)
gdt_donor2 <- ScaleData(gdt_donor2)
gdt_donor2_features <- VariableFeatures(gdt_donor2)

gdt_donor3 <- NormalizeData(gdt_donor3)
gdt_donor3 <- FindVariableFeatures(gdt_donor3)
gdt_donor3 <- ScaleData(gdt_donor3)
gdt_donor3_features <- VariableFeatures(gdt_donor3)

pbmc_4k <- NormalizeData(pbmc_4k)
pbmc_4k <- FindVariableFeatures(pbmc_4k)
pbmc_4k <- ScaleData(pbmc_4k)
pbmc_4k_features <- VariableFeatures(pbmc_4k)

pbmc_8k <- NormalizeData(pbmc_8k)
pbmc_8k <- FindVariableFeatures(pbmc_8k)
pbmc_8k <- ScaleData(pbmc_8k)
pbmc_8k_features <- VariableFeatures(pbmc_8k)

common_features <- Reduce(intersect, list(gdt_donor1_features, gdt_donor2_features, gdt_donor3_features, pbmc_4k_features, pbmc_8k_features))

# Manually set the variable features
VariableFeatures(harmony.method2) <- common_features

# Standard workflow steps
harmony.method2 <- NormalizeData(harmony.method2)
harmony.method2 <- ScaleData(harmony.method2)
harmony.method2 <- RunPCA(harmony.method2)
harmony.method2 <- RunTSNE(harmony.method2, dims = 1:30)
method2_tsne_plot_unintegrated_harmony <- DimPlot(harmony.method2, reduction = "tsne", group.by = "orig.ident", label = TRUE) + ggtitle("Method 2 tSNE unintegrated Harmony")

# Run Harmony to correct batch effects (using "orig.ident" to represent batch info)
harmony.method2 <- RunHarmony(harmony.method2, group.by.vars = "orig.ident", assay = "SCT", verbose = FALSE)
harmony.method2.embeddings <- Embeddings(harmony.method2, reduction = "harmony")

# Find neighbors and clusters
harmony.method2 <- FindNeighbors(harmony.method2, reduction = "harmony", dims = 1:30)
harmony.method2 <- FindClusters(harmony.method2, resolution = 0.8)

# Run tSNE for visualization
harmony.method2 <- RunTSNE(harmony.method2, reduction = "harmony", dims = 1:30)
method2_tsne_plot_integrated_harmony <- DimPlot(harmony.method2, reduction = "tsne", group.by = "orig.ident", label = TRUE) + ggtitle("Method 2 tSNE integrated Harmony")
