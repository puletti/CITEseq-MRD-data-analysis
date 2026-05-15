if (!dir.exists(gene_cell_filtering_outdir)) {
  dir.create(gene_cell_filtering_outdir, recursive = TRUE)
}

message("proceding with normalization: ", method_neighbors)
if (method_norm == "classic"){
  DefaultAssay(seu) <- "RNA"
} else {
  DefaultAssay(seu) <- "SCT"
}

if (method_neighbors == "multimodal") {
  seu <- FindMultiModalNeighbors(
    seu, reduction.list = list("pca", "apca"), 
    dims.list = list(1:n_pcs_RNA, 1:n_pcs_ADT), modality.weight.name = "RNA.weight"
  )
  seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", n.components = 2L)
  
} else if (method_neighbors == "pca") {
  seu = RunUMAP(seu, dims = 1:50, n.neighbors = 30, min.dist = 0.6, spread = 3, negative.sample.rate = 4, reduction = "harmony.RNA")
}




seu = RunUMAP(seu, dims = 1:n_pcs_RNA, reduction.name = "RNA.umap", reduction.key = "rnaUMAP_", n.components = 3)


seu = RunUMAP(seu, dims = 1:n_pcs_RNA, reduction.name = "SCT.umap",reduction.key = "sctUMAP_", n.components = 3, assay = 'SCT')