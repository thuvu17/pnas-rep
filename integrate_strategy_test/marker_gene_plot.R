library(Seurat)
library(patchwork)


# GDT markers from this paper: https://www.pnas.org/doi/10.1073/pnas.1818488116
gdt_markers <- c('TRDC', 'TRGC1', 'TRGC2')

methods_list <- list(
  harmony_method1 = harmony.method1,
  harmony_method2 = harmony.method2
)

for (method_name in names(methods_list)) {
  obj <- methods_list[[method_name]]

  DefaultAssay(obj) <- "RNA"
  
  # Generate feature plots
  p <- FeaturePlot(obj, features = gdt_markers, pt.size = 0.3, ncol = 3) + 
    plot_annotation(title = paste(method_name, "- Marker Genes"))
  
  print(p)

  # ggsave(paste0("plots/", method_name, "_featureplots.png"), p, width = 15, height = 10)
}

harmony.method1.gdt <- subset(harmony.method1, subset = celltype == "gdt", preserve.reductions = TRUE)
harmony.method2.gdt <- subset(harmony.method2, subset = celltype == "gdt", preserve.reductions = TRUE)

harmony.method1.gdt.plot <- DimPlot(harmony.method1,
                                group.by = "celltype", 
                                reduction = "tsne",
                                cells.highlight = list("gdt" = Cells(harmony.method1.gdt)),
                                cols.highlight = "red", cols = "lightgrey") + ggtitle("Seurat v5 method 1 GDT only")
harmony.method2.gdt.plot <- DimPlot(harmony.method2,
                                     group.by = "celltype", 
                                     reduction = "tsne",
                                     cells.highlight = list("gdt" = Cells(harmony.method2.gdt)),
                                     cols.highlight = "red", cols = "lightgrey") + ggtitle("Seurat v5 method 2 GDT only")

harmony.method1.gdt.plot
harmony.method2.gdt.plot