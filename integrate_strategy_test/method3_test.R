# Test: Identify GDT based on marker gene expression using
# 1) Raw "RNA" count
# 2) Integrated data from Seurat v4
# 3) Integrated data from the last iteration of method 3
# 4) No integration


# Set default assay
marker_genes <- c("CD3D", "CD3E", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B", "CD4")
data <- "Subset iterative"

part1.subset[["RNA"]] <- JoinLayers(part1.subset[["RNA"]])
DefaultAssay(part1.subset) <- "RNA"

# Feature plot
FeaturePlot(part1.subset, 
            features = c("CD3E", "CD8A", "CD4", "TRGC1", "TRGC2", "TRDC")) +
  plot_annotation(title = data)
ggsave("integrate_strategy_test/plots/subset/test/subset_iterative_feature_plot.png")

# Identify GDT by markers
gdt.bymarker <- WhichCells(part1.subset, 
                           expression = CD3E > 2 & TRDC > 2 & CD4 == 0 & CD8A == 0)
gdt.true <- WhichCells(part1.subset, expression = orig.ident == "donor1")
gdt.truepos <- intersect(gdt.true, gdt.bymarker)
gdt.falseneg <- setdiff(gdt.true, gdt.bymarker)

cat(length(gdt.true), "\n",
    length(gdt.bymarker), "\n",
    length(gdt.truepos), "\n",
    length(gdt.falseneg), "\n")


DimPlot(part1.subset, reduction = "tsne",
        cells.highlight = gdt.bymarker) + ggtitle(paste("Highlighted identified GDT: ", data))
ggsave("integrate_strategy_test/plots/subset/test/subset_iterative_gdt_highlight.png")
