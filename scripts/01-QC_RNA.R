library(Seurat, quietly = T, verbose = F, warn.conflicts = F)
library("tidyverse", quietly = T, verbose = F, warn.conflicts = F)
library(RColorBrewer, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(pheatmap, quietly = T, verbose = F, warn.conflicts = F)
library(plotly, quietly = T, verbose = F, warn.conflicts = F)
library(ggsci, quietly = T, verbose = F, warn.conflicts = F)
library(viridis, quietly = T, verbose = F, warn.conflicts = F)

if (!dir.exists(QC_RNA_outdir)) {
  dir.create(QC_RNA_outdir, recursive = TRUE)
}

print(paste0("INITIAL OBJECT DIMENSIONs: ", dim(seu)))

#REMOVE MRD2792?? YES
seu = subset(seu, subset = DonorID != "MRD2792")

### COUNTS check
counts <- GetAssayData(seu, assay = 'RNA', layer = 'counts')
genes_per_cell <- Matrix::colSums(counts>0) # count a gene only if it has non-zero reads mapped.
counts_per_cell <- Matrix::colSums(counts)

png(filename=paste0(QC_RNA_outdir, "genes_per_cell_ordered.png"),
    width=1200, height=700)
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
dev.off()

# now replot with the thresholds being shown:
png(file=paste0(QC_RNA_outdir, "genes_per_cell_ordered_thr.png"),
    width=1200, height=700)
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
abline(h=max_nFeature_per_cell, col='blue')
abline(h=min_nFeature_per_cell, col='magenta')  # lower threshold
dev.off()



### MITO check
mito_genes <- grep("^mt-", rownames(counts) , ignore.case=T, value=T)
mito_gene_read_counts = Matrix::colSums(counts[mito_genes,])
pct_mito = mito_gene_read_counts / counts_per_cell * 100

png(file=paste0(QC_RNA_outdir, "mito_perc_ordered.png"),
    width=1200, height=700)
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab = 
       "percentage mitochondrial counts")
dev.off()


png(file=paste0(QC_RNA_outdir, "mito_perc_ordered_thr.png"),
    width=1200, height=700)
plot(sort(pct_mito))
abline(h=max_mito, col='red')
dev.off()



### RIBOSOMAL check
seu[["percent.ribo"]] = PercentageFeatureSet(seu, pattern = "^RPS|^RPL")
ribo_genes = grep("^RPS|^RPL", rownames(counts), ignore.case = T, value = T)
ribo_gene_read_counts = Matrix::colSums(counts[ribo_genes,])
pct_ribo = ribo_gene_read_counts / counts_per_cell * 100

png(file=paste0(QC_RNA_outdir, "ribo_perc_ordered_thr.png"),
    width=1200, height=700)
plot(sort(pct_ribo), xlab = "cells sorted by percentage ribosomal counts", ylab = 
       "percentage ribosomal counts")
abline(h=max_ribo, col='red')
dev.off()


#####


qc_by_sample = seu@meta.data %>%
  group_by(Description) %>%
  summarise(n_cells = n(),
            median_nCount = median(nCount_RNA),
            median_nFeature = median(nFeature_RNA),
            median_pctMT = median(percent_MT))
#print(qc_by_sample)
qc_long <- qc_by_sample %>%
  select(Description, median_nCount, median_nFeature, median_pctMT) %>%
  pivot_longer(cols = -Description, names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    median_nCount = "Median nCount",
                    median_nFeature = "Median nFeature",
                    median_pctMT = "Median %MT"),
    # scale values so they align with secondary axis
    value_scaled = case_when(
      metric == "Median %MT" ~ value * 100,  # boost %MT
      TRUE ~ value
    )
  )

qc0=ggplot() +
  geom_col(
    data = qc_by_sample,
    aes(x = Description, y = n_cells, fill = Description),
    alpha = 0.7
  ) +
  geom_line(
    data = qc_long,
    aes(x = Description, y = value_scaled,
        color = metric, group = metric),
    size = 1.2) +
  geom_point(
    data = qc_long,
    aes(x = Description, y = value_scaled, color = metric),
    size = 2
  ) +
  scale_fill_manual(values = palette_description) + # <- bar colors
  scale_y_continuous(
    name = "Number of cells",
    sec.axis = sec_axis(
      ~./1000,
      name = "Median metrics (nCount, nFeature, %MT×1000)"
    )
  ) +
  theme_minimal() +
  labs(
    title = "QC Overview per Sample",
    x = "Sample",
    color = "Metric"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

colors = viridis(12, option = "turbo")
qc1= FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE) + NoLegend()
qc2= FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent_MT", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE)
qc2b= FeatureScatter(seu, feature1 = "percent_MT", feature2 = "nFeature_RNA", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE)

qc_combined <- qc0 | qc1  |qc2

ggsave(
  filename = paste0(QC_RNA_outdir, "qc_overview.png"),
  plot = qc_combined,
  width = 17,   # set a suitable width
  height = 7,   # depends on your layout
  dpi = 300
)

qc3 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3, group.by = "DonorID", cols = palette_DonorID)
ggsave(
  filename = paste0(QC_RNA_outdir, "qc_overview_violin.png"),
  plot = qc3,
  width = 17,   # set a suitable width
  height = 7,   # depends on your layout
  dpi = 300
)