library(Seurat)
library("tidyverse")
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(plotly)
library(ggsci)
library(viridis)

if (!dir.exists(QC_ADT_outdir)) {
  dir.create(QC_ADT_outdir, recursive = TRUE)
}

counts_ADT <- GetAssayData(seu, assay = 'CITE', layer = 'counts')

colors = viridis(12, option = "turbo")

#protein counts per sample ---------------

meta <- seu@meta.data %>%
  mutate(total_counts = Matrix::colSums(GetAssayData(seu, assay = "CITE", slot = "counts")))
# Reorder Description by median total_counts
meta$Description <- reorder(meta$Description, meta$total_counts, FUN = median)
# Compute medians for labeling
meds <- meta %>%
  group_by(Description) %>%
  summarize(median_value = median(total_counts))

qc3 <- ggplot(meta, aes(x = Description, y = total_counts, fill = Description)) + 
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = palette_description) +
  geom_text(data = meds,
            aes(x = Description, y = median_value, label = round(median_value, 0)),
            color = "white",           # white text contrasts with colored box
            fontface = "bold",
            vjust = 0.5,            
            size = 4.5) +
  labs(y = "RAW ADT sequencing depth (counts per cell)",
       x = "Sample") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = 'none') +
  ylim(0, 3000)

#seu$nCount_ADT

qc4 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nCount_ADT", group.by = "DonorID", shuffle = FALSE, pt.size = 0.02, smooth = TRUE, cols = palette_DonorID) + geom_smooth(method = "lm", se = TRUE) + ylim(1,5000)
#qc5 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nCount_ADT", group.by = "RunID.x", shuffle = FALSE, pt.size = 0.02) + geom_smooth(method = "lm", se = TRUE) + ylim(1,4000)
qc6 <- VlnPlot(seu, features = "nCount_ADT", group.by = "DonorID", cols = palette_DonorID, pt.size = 0.3, alpha = 0.05)+ ylim(1,10000)

adt_df = as.data.frame(t(as.matrix(counts_ADT)))
adt_df$sample = meta$Description

adt_mean = adt_df %>%
  group_by(sample) %>%
  summarise(across(starts_with("CITE-"), mean))

adt_long <- adt_mean %>%
  pivot_longer(cols = starts_with("CITE-"), names_to = "marker", values_to = "expr") %>%
  group_by(sample) %>%
  mutate(percent = expr / sum(expr) * 100)
adt_mat_norm <- adt_long %>%
  select(sample, marker, percent) %>%
  pivot_wider(names_from = marker, values_from = percent) %>%
  column_to_rownames("sample") %>%
  as.matrix()

colors = viridis(11, option = "H") 
qc7 <- ggplot(adt_long, aes(x = sample, y = percent, fill = marker)) +
  scale_fill_manual(values = colors) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Relative ADT composition per sample",
    x = "Sample",
    y = "Relative ADT expression (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45),
    legend.position = "right"
  ) 


qc8 <- pheatmap(adt_mat_norm, scale = "row",
         main = "ADT composition heatmap (row-scaled per marker)",
         cluster_rows = TRUE, cluster_cols = TRUE, breaks = c(seq(-2,3, length.out =100)))


####################### PLOTS
#total log-counts per protein
png(file=paste0(QC_ADT_outdir, "total_log-count_per_protein.png"),
    width=1200, height=700)
boxplot(as.matrix(t(log(counts_ADT+1))), col = colors, horizontal = FALSE, main = "total log-counts per feature")
dev.off()

#total RAW-counts per protein
png(file=paste0(QC_ADT_outdir, "total_RAW-count_per_protein.png"),
    width=1200, height=700)
boxplot(as.matrix(t(counts_ADT+1)), col = colors, horizontal = FALSE, main = "total log-counts per feature", ylim = c(0,1000))
dev.off()

ggsave(
  filename = paste0(QC_ADT_outdir, "RAW_ADT_sequencing_depth.png"),
  plot = qc3,   # set a suitable width
  height = 9,   # depends on your layout
  dpi = 300
)
ggsave(
  filename = paste0(QC_ADT_outdir, "scatter_RNA_vs_ADT.png"),
  plot = qc4,   # set a suitable width
  height = 9,   # depends on your layout
  dpi = 300
)
#ggsave(
#  filename = paste0(QC_ADT_outdir, ""),
#  plot = qc5,   # set a suitable width
#  height = 9,   # depends on your layout
#  dpi = 300
#)
ggsave(
  filename = paste0(QC_ADT_outdir, "vln_nCountADT_perDonor.png"),
  plot = qc6,   # set a suitable width
  height = 9,   # depends on your layout
  dpi = 300
)
ggsave(
  filename = paste0(QC_ADT_outdir, "ADT_composition_perSample.png"),
  plot = qc7,   # set a suitable width
  height = 9,   # depends on your layout
  dpi = 300
)
ggsave(
  filename = paste0(QC_ADT_outdir, "ADT_composition_heatmap.png"),
  plot = qc8,   # set a suitable width
  height = 9,   # depends on your layout
  dpi = 300
)