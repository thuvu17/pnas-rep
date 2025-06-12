library(Seurat)
library(Matrix)
library(ggplot2)
library(magrittr)
library(dplyr)
library(patchwork)

# 1. Preprocessing and QC
pbmc_12k.list <- list(pbmc_4k, pbmc_8k)
pbmc_12k.list <- lapply(pbmc_12k.list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  gene_counts <- LayerData(obj[["RNA"]], layer = "counts")
  gene_filter <- Matrix::rowSums(gene_counts > 0) >= 3
  obj <- obj[gene_filter, ]
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt >= 0.05)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize")
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 8000)
  return(obj)
})

# 2. Integration using Seurat v4
pbmc_12k.features <- SelectIntegrationFeatures(pbmc_12k.list, nfeatures = 8000)
pbmc_12k.anchors <- FindIntegrationAnchors(pbmc_12k.list, anchor.features = pbmc_12k.features)
pbmc_12k.v4 <- IntegrateData(pbmc_12k.anchors)

# 3. Dimensionality reduction and clustering
DefaultAssay(pbmc_12k.v4) <- "integrated"
pbmc_12k.v4 <- ScaleData(pbmc_12k.v4)
pbmc_12k.v4 <- RunPCA(pbmc_12k.v4, npcs = 30)
pbmc_12k.v4 <- RunTSNE(pbmc_12k.v4, dims = 1:30)
pbmc_12k.v4 <- FindNeighbors(pbmc_12k.v4, dims = 1:30)
pbmc_12k.v4.lowres <- FindClusters(pbmc_12k.v4, resolution = 0.6) # low res returned 19 clusters
pbmc_12k.v4.highres <- FindClusters(pbmc_12k.v4, resolution = 1.2) # high res returned 25 clusters
DimPlot(pbmc_12k.v4.lowres, reduction = "tsne", label = TRUE) + ggtitle("0.6 res part 1 Seurat v4")
DimPlot(pbmc_12k.v4.highres, reduction = "tsne", label = TRUE) + ggtitle("1.2 res part 1 Seurat v4")
DimPlot(pbmc_12k.v4, reduction = "tsne", group.by = "orig.ident") + ggtitle("PBMC 12k Seurat v4")


# 4. Finding cluster biomarkers
# find markers for every cluster compared to all remaining cells
# Seurat v4
pbmc_12k.v4.lowres.markers <- FindAllMarkers(pbmc_12k.v4.lowres, only.pos = TRUE)
pbmc_12k.v4.lowres.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

pbmc_12k.v4.highres.markers <- FindAllMarkers(pbmc_12k.v4.highres, only.pos = TRUE)
pbmc_12k.v4.highres.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# 5. Annotate clusters using canonical markers
DefaultAssay(pbmc_12k.v4) <- "RNA"  # For marker expression
FeaturePlot(pbmc_12k.v4.lowres, features = c("CD3E", "CD8A", "TRGC1", "TRGC2", "TRDC"))
plots <- FeaturePlot(pbmc_12k.v4, features = c(
  "CD3E", "CD4", "CD8A", "CD27", "PRF1", "GNLY", "SELL","TRDC", "TRDV1", "TRDV2"
), combine = FALSE)

plots <- lapply(plots, function(p) p + NoLegend())

# Optionally combine them again for display
patchwork::wrap_plots(plots)

# effector memory CD4 and CD8 T cells (CD3E, CD4 or CD8A, CD27, PRF1, GNLY) 
# naïve and central memory CD4 and CD8 T cells (CD3E, CD4 or CD8A, SELL, CD27) 
# NK cells (CD3E-negativeNKG7 cells) 
# classic and intermediate monocytes (e.g., LYZ, S100A9, CD14)
# dendritic cells (DC) (LYZ, S100A9, FCER1A)
# intermediate monocytes (CD14−, LYZ, S100A9, FCGR3A)
# B cells (CD79A, CD79B, CD19) 

# Assign cluster labels based on expression
Idents(pbmc_12k.v4.lowres) <- "seurat_clusters"
pbmc_12k.v4.lowres <- RenameIdents(pbmc_12k.v4.lowres, 
                                `0` = "Classical Mono",
                                `1` = "Naive T Cells/GDT?",
                                `2` = "B Cells",
                                `3` = "Naive T Cells/GDT?",
                                `4` = "CD8 T Cells",
                                `5` = "CD8 T Cells",
                                `6` = "NK Cells",
                                `7` = "Naive T Cells/GDT?",
                                `8` = "Naive T Cells/GDT?",
                                `9` = "Naive T Cells/GDT?",
                                `10` = "Intermediate Mono",
                                `11` = "DC",
                                `12` = "Effector Memory CD8 T Cells",
                                `13` = "Central Memory CD4 T Cells",
                                `14` = "?",
                                `15` = "Platelets",
                                `16` = "Effector Memory CD4 T Cells"
)
DimPlot(pbmc_12k.v4.lowres, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("Part 1 low res 0.6")



# Cytotoxic CD8 T cells
pbmc_12k.cytocd8.bymarker <- WhichCells(pbmc_12k.v4, expression = CD3E > 2 & CD8A > 2 & FCGR3A > 1 & GNLY > 1)
length(pbmc_12k.cytocd8.bymarker) # 20
# NK cells
pbmc_12k.nk.bymarker <- WhichCells(pbmc_12k.v4, expression = CD3E == 0 & NKG7 > 2)
length(pbmc_12k.nk.bymarker) # 456
# GDT cells
pbmc_12k.gdt.bymarker <- WhichCells(pbmc_12k.v4, expression = CD3E > 2 & TRDC > 2 & CD4 == 0 & CD8A == 0)
length(pbmc_12k.gdt.bymarker) # 268
DimPlot(pbmc_12k.v4, reduction = "tsne",
        cells.highlight = list("Cytotoxic CD8 T" = pbmc_12k.cytocd8.bymarker,
                               "NK Cells" = pbmc_12k.nk.bymarker,
                               "GDT Cells" = pbmc_12k.gdt.bymarker),
        cols.highlight = c("blue", "green", "red"),
        cols = "lightgrey") + ggtitle("Highlighted Cell Subtypes by Gene Expression")
