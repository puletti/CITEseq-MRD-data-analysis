if (!dir.exists(dim_reduction_outdir)) {
  dir.create(dim_reduction_outdir, recursive = TRUE)
}


message("proceding with method: ", method_neighbors)

if (method_norm == "classic"){
  DefaultAssay(seu) <- "RNA"
} else {
  DefaultAssay(seu) <- "SCT"
}
message("default assay: ", DefaultAssay(seu))

if (harmony == TRUE) {
  message("HARMONIZATION...")
  reduction_prefix <- "harmony"
  seu <- harmony::RunHarmony(object = seu,
                      group.by.vars = c(batchvar, batchvar2),
                      use.assay = DefaultAssay(seu),
                      reduction.save = paste0(reduction_prefix, as.character(DefaultAssay(seu))),
                      theta = c(4, 4))
  #
  message("HARMONIZATION DONE!")
  print(seu@reductions)
} else {
  reduction_prefix <- ""
}

if (method_neighbors == "multimodal") {
  seu <- FindMultiModalNeighbors(
    seu, reduction.list = list("pca", "apca"), 
    dims.list = list(1:n_pcs_RNA, 1:n_pcs_ADT), modality.weight.name = "RNA.weight"
  )
  seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", n.components = 2L)
  
} else if (method_neighbors == "pca") {
  seu <- RunUMAP(seu,
                dims = 1:n_pcs_RNA,
                reduction = paste0(reduction_prefix, as.character(DefaultAssay(seu))),
                reduction.name = paste0(reduction_prefix, as.character(DefaultAssay(seu)), "_umap"),
                n.neighbors = 30,
                min.dist = 0.6,
                spread = 3,
                negative.sample.rate = 4)
  p1 <- DimPlot(seu, 
                reduction = paste0(reduction_prefix, DefaultAssay(seu), "_umap"),
                group.by = "DonorID",
                label = FALSE,
                label.size = 4,
                repel = TRUE,
                cols = palette_DonorID,
                pt.size = 0.5,
                shuffle = TRUE,
                label.box = TRUE,
                label.color = "white") + theme_bw()
  ggsave(
    filename = paste0(dim_reduction_outdir, reduction_prefix, DefaultAssay(seu), "_umap.png"),
    plot = p1,
    width = 14,
    height = 7,  
    dpi = 300
  )
} 




#seu = RunUMAP(seu, dims = 1:n_pcs_RNA, reduction.name = "RNA.umap", reduction.key = "rnaUMAP_", n.components = 3)
#seu = RunUMAP(seu, dims = 1:n_pcs_RNA, reduction.name = "SCT.umap",reduction.key = "sctUMAP_", n.components = 3, assay = 'SCT')