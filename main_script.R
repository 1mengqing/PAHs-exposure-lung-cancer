# Load the functions
source("seurat_functions.R")

# Main script to process data
data_dirs <- c("/mnt/data/userdata/svip017/biny/data/raw/P/PBMC01_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/P/PBMC07_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/P/PBMC18_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/P/EX14PBMC_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/P/EX16PBMC_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/P/EX28PBMC_filtered_matrix",
               "/mnt/data/userdata/svip017/biny/data/raw/control_P/Y1_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/control_B/Y2_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/control_B/Y3_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/control_P/Y4_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/control_P/Y5_filtered_matrix", 
               "/mnt/data/userdata/svip017/biny/data/raw/control_P/Y6_filtered_matrix")

project_names <- c("P1", "P2", "P3", "P4", "P5", "P6", "C1", "C2", "C3", "C4", "C5", "C6")
group_names <- c(rep("PAHs", 6), rep("Control", 6))

# Create and process Seurat objects
seurat_objs <- mapply(create_and_preprocess_seurat, data_dirs, project_names, group_names, SIMPLIFY = FALSE)
scRNA <- merge_seurat_objects(seurat_objs)
saveRDS(scRNA, file = 'scRNA12.rds')

# Quality control and normalization
scRNA <- initial_qc_normalization(scRNA)

# PCA and Harmony
scRNA <- run_pca_harmony(scRNA)

# Clustering
scRNA <- run_clustering(scRNA)

# Load reference data for annotation
Mouse.se <- MouseRNAseqData()
ImmGen.se <- ImmGenData()
ref_data <- list(Mouse.se = Mouse.se, ImmGen.se = ImmGen.se)

# Annotate cells
scRNA <- annotate_cells(scRNA, ref_data)

# Generate plots
generate_plots(scRNA, "scRNA")

# Find and save differential genes
find_save_diff_genes(scRNA, ".")

# Save final object
saveRDS(scRNA, file = "scRNA_harmony_annotated.rds")
