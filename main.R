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
setwd("/home/P100068/lindas")


config <- yaml::read_yaml("./config.yml")

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

# VARIABLES
## QC
min_nFeature_per_cell <- as.numeric(config$analysis$qc$min_nFeature)  #750
max_nFeature_per_cell <- as.numeric(config$analysis$qc$max_nFeature) #6000
max_mito <- as.numeric(config$analysis$qc$max_mito) #15
max_ribo <- as.numeric(config$analysis$qc$max_ribo) #30

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

source(paste0(scripts_dir, "QC_RNA.R"))
source(paste0(scripts_dir, "QC_ADT.R"))
