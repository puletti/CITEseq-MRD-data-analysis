if (!dir.exists(RNA_norm_outdir)) {
  dir.create(RNA_norm_outdir, recursive = TRUE)
}

DefaultAssay(seu) <- 'RNA'
seu$raw_counts_per_cell <- Matrix::colSums(GetAssayData(seu, slot = "counts"))

meds <- seu@meta.data %>%
  group_by(Description) %>%
  summarise(median_value = median(raw_counts_per_cell, na.rm = TRUE))

seu <- NormalizeData(seu, scale.factor = median(seu$nCount_RNA))

seu$norm_counts_per_cell <- Matrix::colSums(GetAssayData(seu, slot = "data"))


seu = FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu = ScaleData(seu)
seu = RunPCA(seu, features = VariableFeatures(object = seu))


### PLOTS #######

meds <- seu@meta.data %>%
  group_by(Description) %>%
  summarise(median_value = median(raw_counts_per_cell, na.rm = TRUE))

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
  filename = paste0(RNA_norm_outdir, "seqDepht_rawVSnorm.png"),
  plot = pcombined,   
  height = 9,  
  dpi = 300
)