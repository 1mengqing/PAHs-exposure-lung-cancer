# Load necessary packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(celldex)
library(BiocParallel)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
library(ggsci)
library(plyr)
library(dplyr)

# Function to create and preprocess Seurat objects
# @param data_dir Directory containing the data
# @param project_name Project name for the Seurat object
# @param group_name Group name for the metadata
# @return A preprocessed Seurat object
create_and_preprocess_seurat <- function(data_dir, project_name, group_name) {
  data <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = data, project = project_name, min.cells = 3, min.features = 200)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj@meta.data$group <- group_name
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
  return(seurat_obj)
}

# Function to merge multiple Seurat objects
# @param seurat_objs List of Seurat objects to merge
# @return A merged Seurat object
merge_seurat_objects <- function(seurat_objs) {
  scRNA <- do.call(merge, c(list(x = seurat_objs[[1]]), seurat_objs[-1]))
  return(scRNA)
}

# Function for initial quality control and normalization
# @param scRNA A Seurat object
# @return A Seurat object with normalized data
initial_qc_normalization <- function(scRNA) {
  scRNA <- NormalizeData(object = scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
  scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(scRNA)
  scRNA <- ScaleData(scRNA, features = all.genes)
  return(scRNA)
}

# Function to run PCA and Harmony
# @param scRNA A Seurat object
# @return A Seurat object with PCA and Harmony results
run_pca_harmony <- function(scRNA) {
  scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
  system.time({scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")})
  return(scRNA)
}

# Function to find neighbors and clusters, and run UMAP and TSNE
# @param scRNA A Seurat object
# @param dims Dimensions to use for clustering
# @param resolution Resolution for clustering
# @return A Seurat object with clustering and dimensionality reduction results
run_clustering <- function(scRNA, dims = 1:20, resolution = 0.5) {
  scRNA <- FindNeighbors(scRNA, dims = dims, reduction = "harmony")
  scRNA <- FindClusters(scRNA, resolution = resolution)
  scRNA <- RunUMAP(scRNA, dims = dims, reduction = "harmony")
  scRNA <- RunTSNE(scRNA, dims = dims, reduction = "harmony")
  return(scRNA)
}

# Function to annotate cells using SingleR
# @param scRNA A Seurat object
# @param ref_data List of reference datasets for annotation
# @return A Seurat object with annotated cell types
annotate_cells <- function(scRNA, ref_data) {
  data <- GetAssayData(scRNA, slot = "data")
  pred.cluster <- SingleR(test = data, ref = ref_data, clusters = Idents(scRNA))
  celltype <- data.frame(ClusterID = rownames(pred.cluster), celltype_SR = pred.cluster$labels, stringsAsFactors = FALSE)
  scRNA[['celltype_R']] <- celltype$celltype_SR[match(Idents(scRNA), celltype$ClusterID)]
  return(scRNA)
}

# Function to generate plots
# @param scRNA A Seurat object
# @param file_prefix Prefix for the output files
generate_plots <- function(scRNA, file_prefix) {
  pdf(paste0(file_prefix, "_nFeature_RNA.pdf"), width = 10, height = 6)
  VlnPlot(object = scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
  
  pdf(paste0(file_prefix, "_pca.pdf"), width = 10, height = 8)
  ElbowPlot(scRNA)
  dev.off()
  
  pdf(paste0(file_prefix, "_tsne_umap.pdf"), width = 10, height = 8)
  DimPlot(scRNA, reduction = "umap", raster = FALSE) + DimPlot(scRNA, reduction = "tsne", raster = FALSE)
  dev.off()
}

# Function to find and save differential genes
# @param scRNA A Seurat object
# @param output_dir Directory to save the output files
find_save_diff_genes <- function(scRNA, output_dir) {
  subc <- levels(x = scRNA@active.ident)
  for (l in subc) {
    cluster.markers <- FindMarkers(object = scRNA, ident.1 = l, min.pct = 0.25, logfc.threshold = 0.25)
    write.table(cluster.markers, file = file.path(output_dir, paste0(l, '_diffgenes.xls')), sep = '\t', quote = FALSE, row.names = TRUE)
  }
}
