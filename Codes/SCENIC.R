# Load necessary libraries
load_required_libraries <- function() {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SCENIC)
    library(AUCell)
    library(RcisTarget)
    library(GENIE3)
    library(KernSmooth)
    library(BiocParallel)
    library(ggplot2)
    library(data.table)
    library(grid)
    library(ComplexHeatmap)
    library(tidyverse)
    library(doRNG)
    library(doMC)
  })
}

# Step 1: Data Preprocessing
preprocess_data <- function(data_file, sample_size = 3255, output_dir = "int") {
  setwd(output_dir) # Set working directory
  
  # Load Seurat object
  SeuratObject <- readRDS(data_file)
  Idents(SeuratObject) <- 'orig.ident'
  
  # Sample a subset of cells to reduce computation
  subcell <- sample(colnames(SeuratObject), sample_size)
  scRNAsub <- SeuratObject[, subcell]
  
  # Extract expression matrix and cell information
  exprMat <- scRNAsub@assays$RNA@counts
  cellInfo <- SeuratObject@meta.data %>%
    select(orig.ident, nCount_RNA, nFeature_RNA, seurat_clusters, group) %>%
    set_names("Sample", "nUMI", "nGene", "Clusters", "Group")
  
  # Save the processed data
  dir.create(output_dir, showWarnings = FALSE)
  saveRDS(exprMat, file = file.path(output_dir, "exprMat.Rds"))
  saveRDS(cellInfo, file = file.path(output_dir, "cellInfo.Rds"))
  
  # Define sample colors
  colVars <- list(Sample = c("B1-Lung" = "forestgreen", 
                             "B2-Lung" = "darkorange",
                             "C1-Lung" = "#FF97FF",
                             "C2-Lung" = "#FFA15A"))
  colVars$Sample <- colVars$Sample[intersect(names(colVars$Sample), cellInfo$Sample)]
  saveRDS(colVars, file = file.path(output_dir, "colVars.Rds"))
}

# Step 2: SCENIC Setup
setup_scenic <- function(dbDir, output_dir = "int", org = "mgi", nCores = 50) {
  # Load motif annotations
  data(list = "motifAnnotations_mgi_v9", package = "RcisTarget")
  motifAnnotations_mgi <- motifAnnotations_mgi_v9
  
  # Initialize SCENIC options
  scenicOptions <- initializeScenic(org = org, 
                                    dbDir = dbDir, 
                                    dbs = c("mm9-500bp-upstream-7species.mc9nr.feather", 
                                            "mm9-tss-centered-10kb-7species.mc9nr.feather"), 
                                    datasetTitle = "SCENIC example on mice blood", 
                                    nCores = nCores)
  
  # Modify input dataset info
  scenicOptions@inputDatasetInfo$cellInfo <- file.path(output_dir, "cellInfo.Rds")
  scenicOptions@inputDatasetInfo$colVars <- file.path(output_dir, "colVars.Rds")
  
  # Save SCENIC options
  saveRDS(scenicOptions, file = file.path(output_dir, "scenicOptions.Rds"))
}

# Step 3: Coexpression Network Inference
run_coexpression_analysis <- function(exprMat_file, scenicOptions_file, nCores = 200) {
  # Load expression matrix and SCENIC options
  exprMat <- as.matrix(readRDS(exprMat_file))
  scenicOptions <- readRDS(scenicOptions_file)
  
  # Set the number of cores
  scenicOptions@settings$nCores <- nCores
  
  # Filter genes
  genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions)
  exprMat_filtered <- exprMat[genesKept, ]
  
  # Compute correlation and run GENIE3 or GRNBoost2
  runCorrelation(exprMat_filtered, scenicOptions)
  exprMat_filtered_log <- log2(exprMat_filtered + 1)
  
  # Run GENIE3
  runGenie3(exprMat_filtered_log, scenicOptions)
  
  # Save updated SCENIC options
  saveRDS(scenicOptions, file = scenicOptions_file)
}

# Step 4: Run SCENIC
run_scenic_analysis <- function(exprMat_file, scenicOptions_file, nCores = 9) {
  # Load expression matrix and SCENIC options
  exprMat <- as.matrix(readRDS(exprMat_file))
  exprMat_log <- log2(exprMat + 1)
  scenicOptions <- readRDS(scenicOptions_file)
  
  # Set the number of cores
  scenicOptions@settings$nCores <- nCores
  
  # Run SCENIC pipeline
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions, nTopTfs = c(5, 10, 20, 30, 50), corrThr = 0.02)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, minGenes = 5)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
  scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
  
  # Save final SCENIC options
  saveRDS(scenicOptions, file = scenicOptions_file)
}

# Main function to run the entire pipeline
run_pipeline <- function(data_file, dbDir, output_dir = "int", sample_size = 3255, nCores_coexpression = 200, nCores_scenic = 9) {
  # Load required libraries
  load_required_libraries()
  
  # Step 1: Preprocess data
  preprocess_data(data_file = data_file, sample_size = sample_size, output_dir = output_dir)
  
  # Step 2: Setup SCENIC
  setup_scenic(dbDir = dbDir, output_dir = output_dir)
  
  # Step 3: Run coexpression network inference
  run_coexpression_analysis(exprMat_file = file.path(output_dir, "exprMat.Rds"), 
                            scenicOptions_file = file.path(output_dir, "scenicOptions.Rds"), 
                            nCores = nCores_coexpression)
  
  # Step 4: Run SCENIC analysis
  run_scenic_analysis(exprMat_file = file.path(output_dir, "exprMat.Rds"), 
                      scenicOptions_file = file.path(output_dir, "scenicOptions.Rds"), 
                      nCores = nCores_scenic)
}

# Execute the pipeline
run_pipeline(data_file = "/work/user/myan/Mice/mice-lung-pbmc/Neutrophils/Lung/Neutrophils_Lung.rds",
             dbDir = "/home/myan/work/TF/cisTarget_databases/mm9/",
             output_dir = "/work/user/myan/Mice/mice-lung-pbmc/Neutrophils/Lung/tf")
