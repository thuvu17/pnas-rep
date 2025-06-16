# Integration strategy: Integrate, then only select the subset of PBMCs 
# that intersect with GDT. Then perform integration again


iter_num <- 0
# Subset PBMC cells overlapping with GDT
DimPlot(part1.v4, reduction = "tsne", label = TRUE) + NoLegend() +
  ggtitle("PBMC 13k 2000 + marker genes cluster & Highlight GDT") +
  DimPlot(part1.v4, reduction = "tsne", 
          cells.highlight = WhichCells(part1.v4, expression = orig.ident == "donor1"))

# Subset pbmc_4k and pbmc_8k cells
overlap_cluster <- c(1, 2, 4, 6, 7, 9)
overlap_4k.cells <- WhichCells(part1.v4, expression = orig.ident == "pbmc_4k" & seurat_clusters %in% overlap_cluster)
overlap_8k.cells <- WhichCells(part1.v4, expression = orig.ident == "pbmc_8k" & seurat_clusters %in% overlap_cluster)
# Extract counts layer for each
pbmc_4k.counts <- LayerData(part1.v4, assay = "RNA", layer = "counts.pbmc_4k")[, overlap_4k.cells]
pbmc_8k.counts <- LayerData(part1.v4, assay = "RNA", layer = "counts.pbmc_8k")[, overlap_8k.cells]
# Create new Seurat objects
pbmc_4k.subset <- CreateSeuratObject(pbmc_4k.counts, project = "pbmc_4k")
pbmc_8k.subset <- CreateSeuratObject(pbmc_8k.counts, project = "pbmc_8k")

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

DefaultAssay(part1.subset) <- "integrated"
part1.subset <- ScaleData(part1.subset)
part1.subset <- RunPCA(part1.subset, npcs = 30)
part1.subset <- RunTSNE(part1.subset, dims = 1:30)
part1.subset <- FindNeighbors(part1.subset, dims = 1:30)
part1.subset <- FindClusters(part1.subset, resolution = 0.6)

iter_num <- iter_num + 1
DimPlot(part1.subset, reduction = "tsne", label = TRUE) + NoLegend() +
  ggtitle(paste("PBMC 13k 2000 + marker genes & Highlight GDT: Iteration", iter_num)) + 
  DimPlot(
    part1.subset, 
    reduction = "tsne",
    cells.highlight = WhichCells(part1.subset, expression = orig.ident == "donor1"),
    cols.highlight = "red",
    cols = "pink")  # color for pbmc

