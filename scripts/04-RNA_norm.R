if (!dir.exists(RNA_norm_outdir)) {
  dir.create(RNA_norm_outdir, recursive = TRUE)
}

message("proceding with normalization: ", method_norm)

seu$raw_counts_per_cell <- Matrix::colSums(GetAssayData(seu, slot = "counts"))
meds <- seu@meta.data %>%
  group_by(Description) %>%
  summarise(median_value = median(raw_counts_per_cell, na.rm = TRUE))

if (method_norm == "classic") {
  DefaultAssay(seu) <- 'RNA'
  
  seu <- NormalizeData(seu, scale.factor = median(seu$nCount_RNA))
  seu = FindVariableFeatures(seu, selection.method_norm = "vst", nfeatures = n_var_features)
  seu = ScaleData(seu, features = VariableFeatures(seu), vars.to.regress = vars_to_regress)
  seu = RunPCA(seu, features = VariableFeatures(object = seu), npcs = 50, assay = 'RNA')
  #seu = RunUMAP(seu, dims = 1:n_pcs, reduction.name = "RNA.umap", n.components = 3)
  
  ### PLOTS #######
  
  #checking correlations between principale components and sequencing depth
  
  message("checking correlation between Principal Components and sequencing depth")
  
  png(file=paste0(RNA_norm_outdir, "corr_pcs_seqDepth.png"),
      width=1200, height=700)
  plot(cor(Embeddings(seu, "pca")[,1:30], seu$nCount_RNA))
  abline(h=c(0.1, -0.1), col='orange', lty = 2)
  abline(h=c(0.3, -0.3), col='red', lty = 2)
  dev.off()
  
  
  p1= ggplot(seu@meta.data, aes(x = Description, y = raw_counts_per_cell, fill = Description)) +
    scale_fill_manual(values = palette_description) +
    geom_boxplot(outlier.size = 0.5) +
    geom_text(data = meds,
              aes(x = Description, y = median_value, label = round(median_value, 2)),
              vjust = -0.5, color = "black", fontface = "bold") +
    labs(y = "Raw sequencing depth (counts per cell)",
         x = "Sample") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = 'none') 
  
  seu$norm_counts_per_cell <- Matrix::colSums(GetAssayData(seu, slot = "data"))
  meds_norm <- seu@meta.data %>%
    group_by(Description) %>%
    summarise(median_value = median(norm_counts_per_cell, na.rm = TRUE))
  
  p2= ggplot(seu@meta.data, aes(x = Description, y = norm_counts_per_cell, fill = Description)) +
    scale_fill_manual(values = palette_description) +
    geom_boxplot(outlier.size = 0.5) +
    geom_text(data = meds_norm,
              aes(x = Description, y = median_value, label = round(median_value, 2)),
              vjust = -0.5, color = "black", fontface = "bold") +
    labs(y = "Normalized sequencing depth (log-normalized counts per cell)",
         x = "Sample") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = 'none') 
  pcombined <- p1 | p2
  
  ggsave(
    filename = paste0(RNA_norm_outdir, "seqDepht_rawVSnorm_RNA.png"),
    plot = pcombined,   
    height = 9,  
    dpi = 300
  )
} else if (method_norm == "SCTransform split") {
  DefaultAssay(seu) <- 'RNA'
  
  seu_list <- SplitObject(seu, split.by = samples)
  seu_list <- lapply(names(seu_list), function(sample_name) {
    message(paste0("Processing ", sample_name, " (", ncol(seu_list[[sample_name]]), " cells)"))
    SCTransform(seu_list[[sample_name]],
                method_norm = "glmGamPoi",
                vst.flavor = "v2", assay='RNA', n_cells = 2000, exclude_poisson = T,
                return.only.var.genes = FALSE,
                verbose =T)
  })
  # Restore names
  names(seu_list) <- levels(seu$DonorID)
  # Select features that are variable across multiple samples
  features <- SelectIntegrationFeatures(
    object.list = seu_list,
    nfeatures = 3000
  )
  message(paste0("Selected ", length(features), " integration features"))
  # Prepare SCT data for integration
  message("Preparing SCT integration...")
  seu_list <- PrepSCTIntegration(
    object.list = seu_list,
    anchor.features = features
  )
  message("PrepSCTIntegration complete!")
  # Merge all samples back into single object
  message("Merging samples...")
  seu <- merge(
    seu_list[[1]],
    y = seu_list[2:length(seu_list)],
    merge.data = TRUE
  )
  
  DefaultAssay(seu) <- 'SCT'
  seu = RunPCA(seu, dims = 1:30, features = features, assay = 'SCT')
  
  #seurat_rna = RunUMAP(seurat_rna, dims = 1:30, reduction.name = "SCT.umap", n.components = 2, assay = 'SCT')
  
  ### PLOTS #######
  
  #checking correlations between principale components and sequencing depth
  
  message("checking correlation between Principal Components and sequencing depth")
  
  png(file=paste0(RNA_norm_outdir, "corr_pcs_seqDepth.png"),
      width=1200, height=700)
  plot(cor(Embeddings(seu, "pca")[,1:30], seu$nCount_RNA))
  abline(h=c(0.1, -0.1), col='orange', lty = 2)
  abline(h=c(0.3, -0.3), col='red', lty = 2)
  dev.off()
  
  p1= ggplot(seu@meta.data, aes(x = Description, y = raw_counts_per_cell, fill = Description)) +
    scale_fill_manual(values = palette_description) +
    geom_boxplot(outlier.size = 0.5) +
    geom_text(data = meds,
              aes(x = Description, y = median_value, label = round(median_value, 2)),
              vjust = -0.5, color = "black", fontface = "bold") +
    labs(y = "Raw sequencing depth (counts per cell)",
         x = "Sample") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = 'none') 
  
  seu$norm_counts_per_cell <- Matrix::colSums(GetAssayData(seu, slot = "data", assay = "SCT"))
  meds_norm <- seu@meta.data %>%
    group_by(Description) %>%
    summarise(median_value = median(norm_counts_per_cell, na.rm = TRUE))
  
  p2= ggplot(seu@meta.data, aes(x = Description, y = norm_counts_per_cell, fill = Description)) +
    scale_fill_manual(values = palette_description) +
    geom_boxplot(outlier.size = 0.5) +
    geom_text(data = meds_norm,
              aes(x = Description, y = median_value, label = round(median_value, 2)),
              vjust = -0.5, color = "black", fontface = "bold") +
    labs(y = "Normalized sequencing depth (log-normalized counts per cell)",
         x = "Sample") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = 'none') 
  pcombined <- p1 | p2
  
  ggsave(
    filename = paste0(RNA_norm_outdir, "seqDepht_rawVSnorm_RNA.png"),
    plot = pcombined,   
    height = 9,  
    dpi = 300
  )
} else if (method_norm == "SCTransform_classic") {
  DefaultAssay(seu) <- 'RNA'
  seu <- SCTransform(seu,
                     assay = 'RNA',
                     verbose = T)
  DefaultAssay(seu) <- 'SCT'
  seu = RunPCA(seu, features = VariableFeatures(object = seu), npcs = 30)
  
  ### PLOTS #######
  
  #checking correlations between principale components and sequencing depth
  
  message("checking correlation between Principal Components and sequencing depth")
  
  png(file=paste0(RNA_norm_outdir, "corr_pcs_seqDepth.png"),
      width=1200, height=700)
  plot(cor(Embeddings(seu, "pca")[,1:30], seu$nCount_RNA))
  abline(h=c(0.1, -0.1), col='orange', lty = 2)
  abline(h=c(0.3, -0.3), col='red', lty = 2)
  dev.off()
  
  p1= ggplot(seu@meta.data, aes(x = Description, y = raw_counts_per_cell, fill = Description)) +
    scale_fill_manual(values = palette_description) +
    geom_boxplot(outlier.size = 0.5) +
    geom_text(data = meds,
              aes(x = Description, y = median_value, label = round(median_value, 2)),
              vjust = -0.5, color = "black", fontface = "bold") +
    labs(y = "Raw sequencing depth (counts per cell)",
         x = "Sample") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = 'none') 
  
  seu$norm_counts_per_cell <- Matrix::colSums(GetAssayData(seu, slot = "data", assay = "SCT"))
  meds_norm <- seu@meta.data %>%
    group_by(Description) %>%
    summarise(median_value = median(norm_counts_per_cell, na.rm = TRUE))
  
  p2= ggplot(seu@meta.data, aes(x = Description, y = norm_counts_per_cell, fill = Description)) +
    scale_fill_manual(values = palette_description) +
    geom_boxplot(outlier.size = 0.5) +
    geom_text(data = meds_norm,
              aes(x = Description, y = median_value, label = round(median_value, 2)),
              vjust = -0.5, color = "black", fontface = "bold") +
    labs(y = "Normalized sequencing depth (log-normalized counts per cell)",
         x = "Sample") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = 'none') 
  pcombined <- p1 | p2
  
  ggsave(
    filename = paste0(RNA_norm_outdir, "seqDepht_rawVSnorm_RNA.png"),
    plot = pcombined,   
    height = 9,  
    dpi = 300
  )
} else {
  errorCondition("invalid normalization method_norm provided! (classic or SCTransform)")
}

