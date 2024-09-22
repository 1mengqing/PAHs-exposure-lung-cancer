# Load required libraries
load_required_libraries <- function() {
  library(Seurat)
  library(infercnv)
  library(dplyr)
}

# Function to run infercnv analysis
run_infercnv_analysis <- function(seurat_file, gene_order_file, output_dir, ref_groups = c("B cells", "T cells"), cutoff = 0.1, num_threads = 20) {
  # Load required libraries
  load_required_libraries()
  
  # Step 1: Load Seurat object
  scRNA <- readRDS(seurat_file)
  scRNA <- subset(scRNA, subset = orig.ident %in% c("B2-Lung"))
  print(dim(scRNA))  # Check dimensions of the subsetted data
  
  # Step 2: Extract counts matrix and cell annotations
  counts_matrix <- GetAssayData(object = scRNA, assay = "RNA", slot = "counts")
  annotation_cell <- scRNA@meta.data %>% select(Celltype_TOTAL)
  
  # Step 3: Load gene order information
  gene_order <- read.table(gene_order_file, row.names = 1)
  
  # Step 4: Create infercnv object
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                       annotations_file = annotation_cell,
                                       delim = "\t",
                                       gene_order_file = gene_order,
                                       ref_group_names = ref_groups)
  
  # Step 5: Run infercnv analysis
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff = cutoff,  # Default: 0.1 for 10x Genomics
                                out_dir = output_dir,
                                denoise = TRUE,
                                HMM = FALSE,
                                num_threads = num_threads,
                                output_format = "pdf",
                                cluster_by_groups = TRUE,
                                write_expr_matrix = TRUE)
  
  return(infercnv_obj)
}

# Example usage:
infercnv_result <- run_infercnv_analysis(
  seurat_file = "/home/myan/work/Mice/mice-lung-pbmc/scRNA8R-96443.rds",
  gene_order_file = "/home/myan/work/Mice/Tissue/total/infercnv/mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_name.infercnv_positions.txt",
  output_dir = "/home/myan/work/Mice/mice-lung-pbmc/infercnv"
)
