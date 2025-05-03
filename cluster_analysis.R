# Using Seurat cluster analysis workflow: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Assigning the integrated data to "seurat"
seurat <- seuratv4.integrated
data_name <- "seuratv4"
umap_plot_cluster_name <- paste("UMAP of Integrated", data_name, "with Clusters")

# ---------------------------
# Cluster the integrated data
if (data_name == "seuratv5") {
  seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])  # Required in Seurat v5
  DefaultAssay(seurat) <- "RNA"  # Not "integrated", because integration is done in the reduction slot
} else {
  DefaultAssay(seurat) <- "integrated"
}

if (data_name == "seuratv5") {
  seurat <- FindNeighbors(seurat, reduction = "integrated.cca", dims = 1:30)
} else {
  seurat <- FindNeighbors(seurat, dims = 1:30)
}

seurat <- FindClusters(seurat, resolution = 0.3)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat), 5)
# donor1_d1_AAACCTGGTAGAGGAA_1 donor1_d1_AAACGGGCAGACACTT_1 donor1_d1_AAAGCAAAGAGTAATC_1 
# 3                            3                            5 
# donor1_d1_AAAGCAATCATGCATG_1 donor1_d1_AAAGCAATCCTCAACC_1 
# 3                            3 
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13

# ------------------------------------
# Run non-linear dimensional reductiseurat$orig.identon
if (data_name == "seuratv5") {
  seurat <- RunUMAP(seurat, reduction = "integrated.cca", dims = 1:30)
} else {
  seurat <- RunUMAP(seurat, dims = 1:10)
}

umap_plot_cluster <- DimPlot(seurat, reduction = "umap", label = TRUE, shuffle = TRUE, pt.size = 0.5) + ggtitle(umap_plot_cluster_name)
umap_plot_cluster
ggsave(paste0("plots/seuratv4_vs_v5/cluster_analysis/umap_", data_name, "_cluster.png"))
# --------------------------------------------------------------
# Finding differentially expressed features (cluster biomarkers)
if (data_name == "seuratv4") {
  seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])
}

DefaultAssay(seurat) <- "RNA"
markers_d1_vs_d2 <- FindMarkers(seurat, 
                                ident.1 = "d1", 
                                ident.2 = "d2", 
                                group.by = "subtype",
                                min.pct = 0.2, 
                                logfc.threshold = 0.5)
head(markers_d1_vs_d2, n = 5)
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# HIST1H2BB        0 -1.1238307 0.271 0.947         0
# PDGFA            0  1.6466348 0.239 0.887         0
# ZDHHC1           0 -1.7178670 0.313 0.944         0
# RP11-72M17.1     0 -0.7495296 0.927 0.314         0
# KIF18A           0  1.7952190 0.768 0.179         0

d1d2_DE_markers <- head(rownames(markers_d1_vs_d2), n = 5)
gdt_subtype_markers <- c("GZMB", "PRF1", "IL7R", "TCF7")
DC_markers <- c("FCER1A", "CST3", "LYZ", "S100A9")
platelet_markers <- c("PPBP")
paper_markers <- c("CD3E", "CD4", "CD8A", "NKG7", "CD79A", "LYZ", "CD14", "FCGR3A", "FCER1A", "PRF1", "GNLY", "CD27")

FeaturePlot(seurat, features = gdt_subtype_markers, split.by = "subtype", pt.size = 0.5) + 
  plot_annotation(title = paste(data_name, "known GDT subtype markers"))
ggsave(paste0("plots/seuratv4_vs_v5/cluster_analysis/", data_name, "_known_GDT_subtype.png"))
FeaturePlot(seurat, features = d1d2_DE_markers, split.by = "subtype", pt.size = 0.5) + 
  plot_annotation(title = paste(data_name, "d1_d2 DE markers"))
ggsave(paste0("plots/seuratv4_vs_v5/cluster_analysis/", data_name, "_d1d2_DE.png"))
FeaturePlot(seurat, features = DC_markers, pt.size = 0.5) + 
  plot_annotation(title = paste(data_name, "DC markers"))
ggsave(paste0("plots/seuratv4_vs_v5/cluster_analysis/", data_name, "_DC_markers.png"))
FeaturePlot(seurat, features = platelet_markers, pt.size = 0.5) + 
  plot_annotation(title = paste(data_name, "Platelet markers"))
ggsave(paste0("plots/seuratv4_vs_v5/cluster_analysis/", data_name, "_platelet_markers.png"))
FeaturePlot(seurat, features = paper_markers, pt.size = 0.5) + 
  plot_annotation(title = paste(data_name, "PNAS paper markers"))
ggsave(paste0("plots/seuratv4_vs_v5/cluster_analysis/", data_name, "_paper_markers.png"))
