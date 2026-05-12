library(scater, quietly = T, verbose = F, warn.conflicts = F)

if (!dir.exists(gene_cell_filtering_outdir)) {
  dir.create(gene_cell_filtering_outdir, recursive = TRUE)
}

DefaultAssay(seu) <- 'RNA'
### SUBSETTING
####### option 2: hardcode: same threshold for all the samples
print(dim(seu))
print(max_nCount)
seu = subset(seu,
             subset = nFeature_RNA > min_nFeature_per_cell &
                      nFeature_RNA < max_nFeature_per_cell &
                      percent_MT < max_mito &
                      percent.ribo < max_ribo &
                      nCount_RNA < max_nCount
             )

dim(seu)

### PLOTS
qc3= FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE) + NoLegend()
qc4= FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent_MT", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE) + NoLegend()
qc5= FeatureScatter(seu, feature1 = "percent_MT", feature2 = "nFeature_RNA", group.by = "Description", jitter = TRUE,cols = palette_description, shuffle = TRUE)

qc6 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent_MT"),pt.size = 0, group.by = "Description", ncol = 3, cols = palette_description)

qc_combined <- qc3 |qc4 | qc5

ggsave(
  filename = paste0(gene_cell_filtering_outdir, "qc_filtered_overview.png"),
  plot = qc_combined,
  height=7,
  width=17,  
  dpi = 300
)

ggsave(
  filename = paste0(gene_cell_filtering_outdir, "qc_filtered_Vln.png"),
  plot = qc6,   
  height = 8,  
  dpi = 300
)

# some cells are labeled as NA. remove them.
#seu = subset(seu, subset = !is.na(Description))
dim(seu)

## GENE FILTERING

sce = as.SingleCellExperiment(seu)

sce = addPerCellQC(sce)

#colnames(colData(sce))
#plotHighestExprs(sce, exprs_values = "counts", colour_cells_by = "Description", n= 35) +theme(legend.key.size = unit(0.5, "cm"))+ guides(color = guide_legend(override.aes = list(size = 10)))

#FeaturePlot(seu, features = "MALAT1", slot = "counts", max.cutoff = 1000)
#VlnPlot(seu, features = "MALAT1", group.by = "Description", pt.size = 0.01, alpha = 0.02, cols = palette_description, slot = "counts") + ylim(1,1250)

#FeaturePlot(seu, features = "IGKC", slot = "counts", max.cutoff = 20)
#Plot_Density_Custom(seu, features = c('NEAT1'),reduction = 'SCT.umap', method='wkde')
#Plot_Density_Custom(seu, features = c('RPL13', 'RPLP1', 'RPS18', "RPS6", "RPS2", "RPS27", "RPL41", "RPL10", "RPL32", "RPS3A", "EEF1A1", "RPL37"),reduction = 'SCT.umap', method='wkde')

genes_to_remove <- c("MALAT1", "IGKC", "MT-CO3","MT-CO2", "MT-CO1", "FTH1")
blacklist_patterns <- c("^RPS", "^RPL", "^MT-")
blacklist_genes <- grep(paste(blacklist_patterns, collapse="|"), rownames(seu), value=TRUE)
blacklist_genes <- unique(c(genes_to_remove, blacklist_genes))
rna_assay <- subset(seu, features = setdiff(rownames(seu), blacklist_genes))

sce = as.SingleCellExperiment(rna_assay)

sce = addPerCellQC(sce)

#colnames(colData(sce))
p1 <- plotHighestExprs(sce, exprs_values = "counts", colour_cells_by = "Description", n= 35) +theme(legend.key.size = unit(0.5, "cm"))+ guides(color = guide_legend(override.aes = list(size = 10)))
ggsave(
  filename = paste0(gene_cell_filtering_outdir, "HighExpr_genes.png"),
  plot = p1,  
  height = 8,   # depends on your layout
  dpi = 300
)
seu[["RNA"]] <- rna_assay[["RNA"]]