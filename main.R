library(yaml)

#library(Seurat)
#library(harmony)
#library("tidyverse")
#library(scCustomize)

#library(RColorBrewer)
#library(ggplot2)
#library(pheatmap)
#library(plotly)
#library(ggsci)
#library(viridis)
#library(ComplexHeatmap)

#library(dittoSeq)



config <- yaml::read_yaml("./config.yml")
setwd(config$options$home_dir)
.libPaths(config$options$libPath)

# DIRECTORIES
resources_dir <- config$paths$resources
scripts_dir <- config$paths$scripts_dir
input_dir <- config$paths$input_dir
output_dir <- config$paths$output_dir

#STEPS
steps <- config$steps
#1.   qc on RNA 
assign(x = paste0(steps[1], "_outdir"), value = paste0(output_dir, "01-", steps[1], "/"))
#2.   qc on ADT
assign(x = paste0(steps[2], "_outdir"), value = paste0(output_dir, "02-", steps[2], "/"))
#3.   cell and gene filtering based on QC and gene expression
assign(x = paste0(steps[3], "_outdir"), value = paste0(output_dir, "03-", steps[3], "/"))
#4.   RNA normalization
assign(x = paste0(steps[4], "_outdir"), value = paste0(output_dir, "04-", steps[4], "/"))
#5.   ADT normamlization
assign(x = paste0(steps[5], "_outdir"), value = paste0(output_dir, "05-", steps[5], "/"))


# VARIABLES
## QC
min_nFeature_per_cell <- as.numeric(config$analysis$qc$min_nFeature)  #750
max_nFeature_per_cell <- as.numeric(config$analysis$qc$max_nFeature) #6000
max_mito <- as.numeric(config$analysis$qc$max_mito) #15
max_ribo <- as.numeric(config$analysis$qc$max_ribo) #30
max_nCount <- as.numeric(config$analysis$qc$max_ncount) #30000

## Norm
method <- as.character(config$analysis$norm$method)
n_var_features <- as.numeric(config$analysis$norm$n_var_features) #3000
vars_to_regress <- as.character(config$analysis$norm$vars_to_regress)
samples <- as.character(config$analysis$norm$samples)
n_pcs <- as.numeric(config$analysis$norm$n_pcs)

# FILES
seu <- config$files$input
seu <- readRDS(paste0(input_dir, seu))

#PALETTES
palette_description <- c(
  # Group 1 (blue shades)
  "#1f77b4", "#6baed6", "#c6dbef",
  # Group 2 (green shades)
  "#2ca02c", "#74c476", "#c7e9c0",
  # Group 3 (orange shades)
  "#ff7f0e", "#fdae6b", "#fee6ce",
  # Group 4 (purple shades)
  "#8c2d99", "#af8dc3", "#e7c4e8"
)
palette_DonorID <- c(
  # Group 1 (blue)
  "#6baed6",
  # Group 2 (green)
  "#74c476",
  # Group 3 (orange)
  "#fdae6b",
  # Group 4 (purple)
  "#af8dc3"
)
print("executing RNA QC step...")
source(paste0(scripts_dir, "QC_RNA.R"))
print("RNA QC step complete!")
print("executing ADT QC step...")
source(paste0(scripts_dir, "QC_ADT.R"))
print("ADT QC step complete!")
print("executing Gene and Cell filtering step...")
source(paste0(scripts_dir, "Gene_Cell_filtering.R"))
print("Gene and Cell filtering step complete!")
print("executing RNA normalization step...")
source(paste0(scripts_dir, "RNA_norm.R"))
print("RNA normalization step complete!")
print("executing ADT normalization step...")
source(paste0(scripts_dir, "ADT_norm.R"))
print("ADT normalization step complete!")
