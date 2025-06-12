# Integration strategy: Integrate, then only select the subset of PBMCs 
# that intersect with GDT. Then perform integration again


# Subset PBMC cells overlapping with GDT
DimPlot(part1.v4, reduction = "tsne", label = TRUE) + NoLegend() +
  ggtitle("PBMC_13k cluster analysis") +
  DimPlot(part1.v4, reduction = "tsne", 
          cells.highlight = WhichCells(part1.v4, expression = orig.ident == "donor1")) 
overlap_cluster <- c(1, 2, 3, 5, 6, 7, 8, 11, 12)
pbmc_4k.subset <- subset(part1.subset, subset = orig.ident == "pbmc4k" 
                         & seurat_clusters %in% overlap_cluster)
pbmc_8k.subset <- subset(part1.subset, subset = orig.ident == "pbmc8k" 
                         & seurat_clusters %in% overlap_cluster)
pbmc_4k.subset <- CreateSeuratObject(LayerData(part1.v4, assay = "RNA", layer = "counts.pbmc_4k")
                                     [, colnames(pbmc_4k.subset)], project = "subset_4k")
pbmc_8k.subset <- CreateSeuratObject(LayerData(part1.v4, assay = "RNA", layer = "counts.pbmc_8k")
                                     [, colnames(pbmc_8k.subset)], project = "subset_8k")

marker_genes <- c("CD3E", "CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B")
marker_genes %in% rownames(pbmc_4k.subset)
marker_genes %in% rownames(pbmc_8k.subset)

# Perform integration with subset pbmc + gdt_donor1
part1.subset.list <- c(gdt_donor1, pbmc_4k.subset, pbmc_8k.subset)
part1.subset.list <- lapply(part1.subset.list, function(obj) {
  DefaultAssay(obj) <- "RNA"
  if ("scale.data" %in% Layers(obj[["RNA"]])) {
    obj[["RNA"]]@layers[["scale.data"]] <- NULL
  }
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
})
subset.features <- SelectIntegrationFeatures(part1.subset.list, nfeatures = 2000)
subset.features <- unique(c(subset.features, marker_genes))
subset.anchors <- FindIntegrationAnchors(part1.subset.list, anchor.features = subset.features)
part1.subset <- IntegrateData(subset.anchors)

# Check if all marker genes are present
marker_genes %in% rownames(part1.subset)


DefaultAssay(part1.subset) <- "integrated"
part1.subset <- ScaleData(part1.subset)
part1.subset <- RunPCA(part1.subset, npcs = 30)
part1.subset <- RunTSNE(part1.subset, dims = 1:30)
part1.subset <- FindNeighbors(part1.subset, dims = 1:30)
part1.subset <- FindClusters(part1.subset, resolution = 1.2)

DimPlot(part1.subset, reduction = "tsne", label = TRUE) + NoLegend() +
  ggtitle("Iteration 3") + 
  DimPlot(
    part1.subset, 
    reduction = "tsne",
    cells.highlight = WhichCells(part1.subset, expression = orig.ident == "donor1"),
    cols.highlight = "red",
    cols = "pink")  # color for pbmc

