# Extract only the marker genes
marker_genes <- c("CD3D", "CD3E", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B")
# --------------
# PBMC 12k
# From RNA assay
pbmc_4k.rna_12k <- LayerData(pbmc_12k.v4, assay = "RNA", layer = "data.pbmc_4k")
pbmc_8k.rna_12k <- LayerData(pbmc_12k.v4, assay = "RNA", layer = "data.pbmc_8k")
pbmc_4k.rna_markersubset <- pbmc_4k.rna_12k[marker_genes_noCD3E, , drop = FALSE]
pbmc_8k.rna_markersubset <- pbmc_8k.rna_12k[marker_genes_noCD3E, , drop = FALSE]

# From integrated assay
cells_4k_12k <- WhichCells(pbmc_12k.v4, expression = orig.ident == "pbmc_4k")
cells_8k_12k <- WhichCells(pbmc_12k.v4, expression = orig.ident == "pbmc_8k")
pbmc12k.int_expr_12k <- GetAssayData(pbmc_12k.v4, assay = "integrated", slot = "data")
pbmc_4k.int_12k <- pbmc12k.int_expr_12k[, cells_4k_12k]
pbmc_8k.int_12k <- pbmc12k.int_expr_12k[, cells_8k_12k]
# !!! CD3E is missing from the integrated assay, so we will not include it in the subset
pbmc_4k.int_markersubset <- pbmc_4k.int_12k[marker_genes_noCD3E, , drop = FALSE]
pbmc_8k.int_markersubset <- pbmc_8k.int_12k[marker_genes_noCD3E, , drop = FALSE]

# The "integrated" assay contains only the top 2,000 integration features, not the full transcriptome. 
# So we can’t compare full integrated with RNA.
# Identify shared genes
common_genes.4k_12k <- intersect(rownames(pbmc_4k.rna_12k), rownames(pbmc_4k.int_12k))
common_genes.8k_12k <- intersect(rownames(pbmc_8k.rna_12k), rownames(pbmc_8k.int_12k))
# Subset to common genes
pbmc_4k.rna_12k <- pbmc_4k.rna_12k[common_genes.4k_12k, ]
pbmc_8k.rna_12k <- pbmc_8k.rna_12k[common_genes.8k_12k, ]


# ====================================
# PBMC 13k
# From RNA assay
pbmc_4k.rna_13k <- LayerData(part1.v4, assay = "RNA", layer = "data.pbmc_4k")
pbmc_8k.rna_13k <- LayerData(part1.v4, assay = "RNA", layer = "data.pbmc_8k")
gdt_donor1.rna_13k <- LayerData(part1.v4, assay = "RNA", layer = "data.gdt")

present_markers_4k <- intersect(marker_genes, rownames(pbmc_4k.rna_13k))
present_markers_8k <- intersect(marker_genes, rownames(pbmc_8k.rna_13k))
present_markers_gdt <- intersect(marker_genes, rownames(gdt_donor1.rna_13k))
pbmc_4k.rna_markersubset <- pbmc_4k.rna_13k[present_markers_4k, , drop = FALSE]
pbmc_8k.rna_markersubset <- pbmc_8k.rna_13k[present_markers_8k, , drop = FALSE]
gdt_donor1.rna_markersubset <- gdt_donor1.rna_13k[present_markers_gdt, , drop = FALSE]

# From integrated assay
cells_4k_13k <- WhichCells(part1.v4, expression = orig.ident == "pbmc_4k")
cells_8k_13k <- WhichCells(part1.v4, expression = orig.ident == "pbmc_8k")
cells_gdt_13k <- WhichCells(part1.v4, expression = orig.ident == "donor1")

pbmc13k.int_expr_13k <- GetAssayData(part1.v4, assay = "integrated", slot = "data")
pbmc_4k.int_13k <- pbmc13k.int_expr_13k[, cells_4k_13k]
pbmc_8k.int_13k <- pbmc13k.int_expr_13k[, cells_8k_13k]
gdt_donor1.int_13k <- pbmc13k.int_expr_13k[, cells_gdt_13k]

# The "integrated" assay contains only the top 8,000 integration features, not the full transcriptome. 
# So we can’t compare full integrated with RNA.
# Identify shared genes
common_genes.4k_13k <- intersect(rownames(pbmc_4k.rna_13k), rownames(pbmc_4k.int_13k))
common_genes.8k_13k <- intersect(rownames(pbmc_8k.rna_13k), rownames(pbmc_8k.int_13k))
common_genes.gdt <- intersect(rownames(gdt_donor1.rna_13k), rownames(gdt_donor1.int_13k))
# Subset to common genes
pbmc_4k.rna_13k <- pbmc_4k.rna_13k[common_genes.4k_13k, ]
pbmc_8k.rna_13k <- pbmc_8k.rna_13k[common_genes.8k_13k, ]
gdt_donor1.rna_13k <- gdt_donor1.rna_13k[common_genes.gdt, ]
gdt_donor1.int_13k <- gdt_donor1.int_13k[common_genes.gdt, ]

present_markers_4k <- intersect(marker_genes, rownames(pbmc_4k.int_13k))
present_markers_8k <- intersect(marker_genes, rownames(pbmc_8k.int_13k))
present_markers_gdt <- intersect(marker_genes, rownames(gdt_donor1.int_13k))
pbmc_4k.int_markersubset <- pbmc_4k.int_13k[marker_genes, , drop = FALSE]
pbmc_8k.int_markersubset <- pbmc_8k.int_13k[marker_genes, , drop = FALSE]
gdt_donor1.int_markersubset <- gdt_donor1.int_13k[marker_genes, , drop = FALSE]

# Define delog function
delog_matrix <- function(mat) {
  exp(mat) - 1
}

# Apply to all marker subsets (RNA)
pbmc_4k.rna_markersubset <- delog_matrix(pbmc_4k.rna_markersubset)
pbmc_8k.rna_markersubset <- delog_matrix(pbmc_8k.rna_markersubset)
gdt_donor1.rna_markersubset <- delog_matrix(gdt_donor1.rna_markersubset)

# Apply to all marker subsets (Integrated)
pbmc_4k.int_markersubset <- delog_matrix(pbmc_4k.int_markersubset)
pbmc_8k.int_markersubset <- delog_matrix(pbmc_8k.int_markersubset)
gdt_donor1.int_markersubset <- delog_matrix(gdt_donor1.int_markersubset)

