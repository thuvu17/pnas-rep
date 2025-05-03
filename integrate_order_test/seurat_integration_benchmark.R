# Load required libraries
library(Seurat)
library(lisi)
library(cluster)
library(ggplot2)

# Function to compute ASW
compute_asw <- function(seurat_obj, reduction, label_col = "seurat_clusters") {
  library(cluster)
  
  # Get embeddings
  embeddings <- Embeddings(seurat_obj, reduction = reduction)
  
  # Get cluster labels as vector
  cluster_labels <- as.vector(seurat_obj[[label_col]][, 1])
  
  # Compute distance matrix
  dist_matrix <- dist(embeddings)
  
  # Compute silhouette
  sil <- silhouette(as.integer(as.factor(cluster_labels)), dist_matrix)
  mean(sil[, 3])  # Average silhouette width
}


# Function to compute iLISI and cLISI
compute_lisi_scores <- function(seurat_obj, reduction, label_col = "celltype", batch_col = "orig.ident") {
  library(lisi)
  
  # Get embeddings
  embeddings <- Embeddings(seurat_obj, reduction = reduction)
  
  # Prepare metadata
  metadata <- seurat_obj@meta.data
  metadata[[label_col]] <- as.character(metadata[[label_col]])
  metadata[[batch_col]] <- as.character(metadata[[batch_col]])
  
  # Compute LISI
  lisi_batch <- compute_lisi(embeddings, metadata, batch_col)
  lisi_celltype <- compute_lisi(embeddings, metadata, label_col)
  
  # Extract and return means
  iLISI <- mean(lisi_batch[[batch_col]])
  cLISI <- mean(lisi_celltype[[label_col]])
  
  return(list(iLISI = iLISI, cLISI = cLISI))
}



# List of Seurat objects to benchmark
method_list <- list(
  #  seuratv4_method1 = seuratv4.method1,
  # seuratv4_method2 = seuratv4.method2,
  # seuratv5_method1 = seuratv5.method1
  seuratv5_method2 = seuratv5.method2
#  harmony_method1 = harmony.method1
#    harmony_method2 = harmony.method2
)

reduction_map <- list(
  #  seuratv4_method1 = "pca",
  # seuratv4_method2 = "pca"
  # seuratv5_method1 = "integrated.cca"
  seuratv5_method2 = "integrated.cca"
#  harmony_method1 = "harmony"
#    harmony_method2 = "harmony"
)


# Create a results data frame
results <- data.frame(Method = character(), ASW = numeric(), iLISI = numeric(), cLISI = numeric(), stringsAsFactors = FALSE)

# Iterate and compute metrics
for (name in names(method_list)) {
  cat("Processing:", name, "\n")
  obj <- method_list[[name]]
  reduction <- reduction_map[[name]]
  embeddings <- Embeddings(obj, reduction = reduction)
  

  # Set default assay to integrated (assuming it's been run)
  # DefaultAssay(obj) <- "integrated"

  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    obj <- FindNeighbors(obj, reduction = reduction, dims = 1:30)
    obj <- FindClusters(obj, resolution = 1)
  }
  
  # Run PCA if not already present
  if (!"pca" %in% names(obj@reductions)) {
    obj <- RunPCA(obj)
  }

  # ASW
  asw <- compute_asw(obj, reduction)

  # LISI
  lisi <- compute_lisi_scores(obj, reduction)

  # Store results
  results <- rbind(results, data.frame(Method = name, ASW = asw, iLISI = lisi$iLISI, cLISI = lisi$cLISI))
}
  # Print results
  print(paste("ASW:", asw))
  print(paste("iLISI:", lisi$iLISI))
  print(paste("cLISI:", lisi$cLISI))
