library(scater)


DefaultAssay(seu) <- 'RNA'
### SUBSETTING
####### option 1: adaptive for each sample (based on median and Median Absolute Deviation)
dim(seurat_rna)

meta <- seurat_rna@meta.data
sample_thresholds <- meta %>% group_by(Description) %>%
  summarise(med_genes = median(nFeature_RNA),
            mad_genes = mad(nFeature_RNA),
            med_counts = median(nCount_RNA),
            mad_counts = mad(nCount_RNA),
            med_pctmt = median(percent_MT),
            mad_pctmt = mad(percent_MT)) %>%
  mutate(min_genes = pmax(200, med_genes - 3*mad_genes),
         max_counts = med_counts + 3*mad_counts,
         max_pctmt = med_pctmt + 3*mad_pctmt)
sample_thresholds


seurat_rna@meta.data <- left_join(seurat_rna@meta.data, sample_thresholds, by = "Description")
to_keep <- with(seurat_rna@meta.data,
                nFeature_RNA >= min_genes &
                  nCount_RNA <= max_counts &
                  percent_MT <= max_pctmt)

seurat_rna@graphs = list()

DefaultAssay(seurat_rna) = "RNA"
rownames(seurat_rna@meta.data) = seurat_rna@meta.data$CellID
seurat_rna = subset(seurat_rna, cells = rownames(seurat_rna@meta.data)[to_keep])

seurat_rna = subset(seurat_rna, subset = !is.na(Description))
dim(seurat_rna)

####### option 2: hardcode: same threshold for all the samples
seu = subset(seu,
                    subset = nFeature_RNA > min_nFeature_per_cell &
                      nFeature_RNA < max_nFeature_per_cell &
                      percent_MT < max_mito &
                      percent.ribo < max_ribo &
                      nCount_RNA < max_nCount)



dim(seu)



### PLOTS
qc3= FeatureScatter(seurat_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE) + NoLegend()
qc4= FeatureScatter(seurat_rna, feature1 = "nCount_RNA", feature2 = "percent_MT", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE) + NoLegend()
qc5= FeatureScatter(seurat_rna, feature1 = "percent_MT", feature2 = "nFeature_RNA", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE)

qc6 <- VlnPlot(seurat_rna, features = c("nFeature_RNA", "nCount_RNA", "percent_MT"),pt.size = 0.2, alpha = 0.02, group.by = "Description", ncol = 3, cols = palette_description)

qc_combined <- qc3 |qc4 | qc5

ggsave(
  filename = paste0(gene_cell_filtering_outdir, "qc_filtered_overview.png"),
  plot = qc_combined,   
  height = 8,  
  dpi = 300
)

ggsave(
  filename = paste0(gene_cell_filtering_outdir, "qc_filtered_Vln.png"),
  plot = qc6,   
  height = 8,  
  dpi = 300
)


# some cells are labeled as NA. remove them.
seu = subset(seu, subset = !is.na(Description))
dim(seu)




## GENE FILTERING

sce = as.SingleCellExperiment(seu)

sce = addPerCellQC(sce)

#colnames(colData(sce))
plotHighestExprs(sce, exprs_values = "counts", colour_cells_by = "Description", n= 35) +theme(legend.key.size = unit(0.5, "cm"))+ guides(color = guide_legend(override.aes = list(size = 10)))

FeaturePlot(seurat_rna, features = "MALAT1", slot = "counts", max.cutoff = 1000)
VlnPlot(seurat_rna, features = "MALAT1", group.by = "Description", pt.size = 0.01, alpha = 0.02, cols = palette_description, slot = "counts") + ylim(1,1250)

FeaturePlot(seurat_rna, features = "IGKC", slot = "counts", max.cutoff = 20)
Plot_Density_Custom(seurat_rna, features = c('NEAT1'),reduction = 'SCT.umap', method='wkde')
Plot_Density_Custom(seurat_rna, features = c('RPL13', 'RPLP1', 'RPS18', "RPS6", "RPS2", "RPS27", "RPL41", "RPL10", "RPL32", "RPS3A", "EEF1A1", "RPL37"),reduction = 'SCT.umap', method='wkde')

## B2M    0.019416325  52587
## TMSB4X 0.016583278  44914
## EEF1A1 0.012944217  35058
## RPL21  0.010354856  28045
## RPS27  0.010097877  27349
## RPL13  0.009351309  25327
## RPL13A 0.008643877  23411
## RPL10  0.008105920  21954
## RPLP1  0.007909124  21421
## RPL34  0.007891401  21373
## RPS12  0.007672083  20779
## RPS18  0.007560208  20476
## RPS3A  0.007439842  20150
## RPL32  0.007433934  20134
## RPS6   0.007406612  20060
## RPS27A 0.007243415  19618
## RPS2   0.007188770  19470
## RPL41  0.006997143  18951
## RPL11  0.006277897  17003

genes_to_remove <- c("MALAT1", "IGKC", "MT-CO3","MT-CO2", "MT-CO1", "FTH1")
blacklist_patterns <- c("^RPS", "^RPL", "^MT-")
blacklist_genes <- grep(paste(blacklist_patterns, collapse="|"), rownames(seurat_rna), value=TRUE)
blacklist_genes <- unique(c(genes_to_remove, blacklist_genes))
rna_assay <- subset(seurat_rna, features = setdiff(rownames(seurat_rna), blacklist_genes))

sce = as.SingleCellExperiment(rna_assay)

sce = addPerCellQC(sce)

#colnames(colData(sce))
plotHighestExprs(sce, exprs_values = "counts", colour_cells_by = "Description", n= 35) +theme(legend.key.size = unit(0.5, "cm"))+ guides(color = guide_legend(override.aes = list(size = 10)))

seurat_rna[["RNA"]] <- rna_assay[["RNA"]]