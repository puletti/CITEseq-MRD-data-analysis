if (!dir.exists(ADT_norm_outdir)) {
  dir.create(ADT_norm_outdir, recursive = TRUE)
}

DefaultAssay(seu) <- 'CITE'
VariableFeatures(seu) <- rownames(seu[["CITE"]])
seu <- NormalizeData(seu, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

norm_mat_ADT <- GetAssayData(seu, slot = "data")

# Sum across genes per cell
seu$norm_counts_per_cell_ADT <- Matrix::colSums(norm_mat_ADT)

meds_ADT <- seu@meta.data %>%
  group_by(Description) %>%
  summarise(median_value = median(norm_counts_per_cell_ADT, na.rm = TRUE))

p1 <- ggplot(seu@meta.data, aes(x = Description, y = norm_counts_per_cell_ADT, fill = Description)) +
  scale_fill_manual(values = palette_description) +
  geom_boxplot(outlier.size = 0.5) +
  geom_text(data = meds_ADT,
            aes(x = Description, y = median_value, label = round(median_value, 2)),
            vjust = -0.5, color = "black", fontface = "bold") +
  labs(y = "Normalized sequencing depth (log-normalized counts per cell)",
       x = "Sample") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = 'none') 

ggsave(
  filename = paste0(ADT_norm_outdir, "seqDepht_norm_ADT.png"),
  plot = p1,   
  height = 9,  
  dpi = 300
)