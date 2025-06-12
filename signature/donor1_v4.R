library(Seurat)
library(Matrix)
library(ggplot2)
library(magrittr)
library(dplyr)
library(patchwork)

# 1. Preprocessing and QC
# Add prefix to cell names
part1.list <- list(pbmc_4k, pbmc_8k, gdt_donor1)
part1.list <- lapply(part1.list, function(obj) {
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
part1.features <- SelectIntegrationFeatures(part1.list, nfeatures = 8000)
part1.anchors <- FindIntegrationAnchors(part1.list, anchor.features = part1.features)
part1.v4 <- IntegrateData(part1.anchors)

# 3. Dimensionality reduction and clustering
DefaultAssay(part1.v4) <- "integrated"
part1.v4 <- ScaleData(part1.v4)
part1.v4 <- RunPCA(part1.v4, npcs = 30)
ElbowPlot(part1.v4) + ggtitle("Part 1 Seurat v4")

part1.v4 <- RunTSNE(part1.v4, dims = 1:30)
part1.v4 <- FindNeighbors(part1.v4, dims = 1:30)
part1.v4 <- FindClusters(part1.v4, resolution = 0.6) # low res returned 19 clusters
# part1.v4.highres <- FindClusters(part1.v4, resolution = 1.2) # high res returned 25 clusters
DimPlot(part1.v4, reduction = "tsne", label = TRUE) + ggtitle("0.6 res PBMC 13k Seurat v4")
# DimPlot(part1.v4.highres, reduction = "tsne", label = TRUE) + ggtitle("1.2 res part 1 Seurat v4")
DimPlot(part1.v4, reduction = "tsne", cells.highlight = WhichCells(part1.v4, expression = orig.ident == "donor1")) + ggtitle("PBMC 13k highlight GDT")


# 4. Finding cluster biomarkers
# find markers for every cluster compared to all remaining cells
# Seurat v4
part1.v4.lowres.markers <- FindAllMarkers(part1.v4.lowres, only.pos = TRUE)
part1.v4.lowres.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

part1.v4.highres.markers <- FindAllMarkers(part1.v4.highres, only.pos = TRUE)
part1.v4.highres.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# 5. Annotate clusters using canonical markers
DefaultAssay(part1.v4.lowres) <- "RNA"  # For marker expression
FeaturePlot(part1.v4.lowres, features = c("CD3E", "CD8A", "TRGC1", "TRGC2", "TRDC"))
plots <- FeaturePlot(part1.v4, features = c(
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
Idents(part1.v4.lowres) <- "seurat_clusters"
part1.v4.lowres <- RenameIdents(part1.v4.lowres, 
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
DimPlot(part1.v4.lowres, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("Part 1 low res 0.6")



# Cytotoxic CD8 T cells
part1.cytocd8.bymarker <- WhichCells(part1.v4.lowres, expression = CD3E > 1 & CD8A > 1 & FCGR3A > 1 & GNLY > 1)
length(part1.cytocd8.bymarker) # 248 using 8000 genes, 227 using 2000 genes in IntegrateData()
# NK cells
part1.nk.bymarker <- WhichCells(part1.v4.lowres, expression = CD3E == 0 & NKG7 > 2)
length(part1.nk.bymarker) # 588
# GDT cells
part1.gdt.bymarker <- WhichCells(part1.v4.lowres, expression = CD3E > 1 & TRDC > 1 & CD4 == 0 & CD8A == 0)
length(part1.gdt.bymarker) # 707
cd8_plot <- DimPlot(part1.v4, reduction = "tsne", cells.highlight = part1.cytocd8.bymarker) + ggtitle("part 1 Seurat v4 T CD8 cells by marker")
nk_plot <- DimPlot(part1.v4, reduction = "tsne", cells.highlight = part1.nk.bymarker) + ggtitle("part 1 Seurat v4 NK cells by marker")
gdt_plot <- DimPlot(part1.v4, reduction = "tsne", cells.highlight = part1.gdt.bymarker) + ggtitle("part 1 Seurat v4 GDT cells by marker")
markers_plot <- DimPlot(part1.v4.lowres, reduction = "tsne",
                        cells.highlight = list("Cytotoxic CD8 T" = part1.cytocd8.bymarker,
                                               "NK Cells" = part1.nk.bymarker,
                                               "GDT Cells" = part1.gdt.bymarker),
                        cols.highlight = c("blue", "green", "red"),
                        cols = "lightgrey") + ggtitle("Highlighted Cell Subtypes by Gene Expression")

cd8_plot + nk_plot + gdt_plot + markers_plot
