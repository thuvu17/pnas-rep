library(Seurat)
library(cluster)
library(lisi)

# ========== Utility Functions ========== #

# Compute ASW
compute_asw <- function(seurat_obj, reduction, cluster_col = "seurat_clusters") {
  embeddings <- Embeddings(seurat_obj, reduction = reduction)
  clusters <- seurat_obj[[cluster_col]][, 1]
  sil <- silhouette(as.integer(as.factor(clusters)), dist(embeddings))
  return(mean(sil[, 3]))
}

# Compute LISI
compute_lisi_group <- function(seurat_obj, reduction, label_col, cells.use) {
  embeddings <- Embeddings(seurat_obj, reduction = reduction)[cells.use, ]
  metadata <- seurat_obj@meta.data[cells.use, , drop = FALSE]
  metadata[[label_col]] <- as.character(metadata[[label_col]])
  lisi_scores <- compute_lisi(embeddings, metadata, label_col)
  return(mean(lisi_scores[[label_col]]))
}

# ========== Main Benchmark Function ========== #

evaluate_integration_metrics <- function(seurat_obj, reduction = "pca") {
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction, dims = 1:30)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  }

  # Identify cells by group
  pbmc_cells <- WhichCells(seurat_obj, idents = NULL, expression = orig.ident %in% c("pbmc_4k", "pbmc_8k"))
  gdt_cells <- WhichCells(seurat_obj, idents = NULL, expression = orig.ident %in% c("donor1", "donor2", "donor3"))

  # Compute metrics
  asw <- compute_asw(seurat_obj, reduction)
  
  ilisi_pbmc <- compute_lisi_group(seurat_obj, reduction, label_col = "orig.ident", cells.use = pbmc_cells)
  ilisi_gdt  <- compute_lisi_group(seurat_obj, reduction, label_col = "orig.ident", cells.use = gdt_cells)

  clisi <- compute_lisi_group(seurat_obj, reduction, label_col = "celltype")

  metrics <- data.frame(
    #Group = c("PBMC", "GDT"),
    ASW = c(asw)
    #iLISI = c(ilisi_pbmc, ilisi_gdt),
    #cLISI = c(clisi)
  )

  return(metrics)
}

# Run this on each of your integrated Seurat objects
metrics_v4 <- evaluate_integration_metrics(seuratv4.method1, reduction = "pca")
metrics_v5 <- evaluate_integration_metrics(seuratv5.method1, reduction = "integrated.cca")
metrics_harmony <- evaluate_integration_metrics(harmony.method1, reduction = "harmony")

print(metrics_v4)
print(metrics_v5)
print(metrics_harmony)
