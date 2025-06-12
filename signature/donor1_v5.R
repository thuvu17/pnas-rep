# Seurat V5.2 data integration pipeline: https://satijalab.org/seurat/articles/integration_introduction.html

part1.v5 <- merge(x=gdt_donor1, y=list(pbmc_4k, pbmc_8k))

# run standard anlaysis workflow
part1.v5 <- NormalizeData(part1.v5)
part1.v5 <- FindVariableFeatures(part1.v5)
part1.v5 <- ScaleData(part1.v5)
part1.v5 <- RunPCA(part1.v5)

part1.v5 <- IntegrateLayers(object = part1.v5, 
                            method = CCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.cca",
                            new.layer = "integrated",
                            verbose = FALSE)

part1.v5 <- FindNeighbors(part1.v5, reduction = "integrated.cca", dims = 1:30)
part1.v5 <- FindClusters(part1.v5, resolution = 1)
