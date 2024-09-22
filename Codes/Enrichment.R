# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(GSEABase)
library(GSVA)
library(org.Mm.eg.db)
library(msigdbr)
library(pheatmap)
library(ggplot2)

# Function to load data
load_data <- function(file_path) {
  data <- readRDS(file_path)
  return(data)
}

# Function to subset data for Monocytes
subset_monocytes <- function(data) {
  Mo <- subset(data, subset = Celltype == "Monocytes")
  saveRDS(Mo, "Monocytes.rds")
  return(Mo)
}

# Function to perform differential expression analysis
perform_deg <- function(pro1) {
  cells1 <- subset(pro1@meta.data, group %in% c("Bap")) %>% rownames()
  cells2 <- subset(pro1@meta.data, group %in% c("Control")) %>% rownames()
  deg <- FindMarkers(pro1, ident.1 = cells1, ident.2 = cells2)
  deg <- data.frame(gene = rownames(deg), deg)
  return(deg)
}

# Function to classify genes as up or down-regulated
classify_genes <- function(deg) {
  k1 <- (deg$p_val_adj < 0.05) & (deg$avg_log2FC < -0.5)
  k2 <- (deg$p_val_adj < 0.05) & (deg$avg_log2FC > 0.5)
  change <- ifelse(k1, "down", ifelse(k2, "up", "stable"))
  deg$change <- change
  return(deg)
}

# Function to convert gene names to Entrez IDs
convert_gene_ids <- function(deg) {
  s2e <- bitr(deg$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  deg <- inner_join(deg, s2e, by = c("gene" = "SYMBOL"))
  return(deg)
}

# Function to perform GO enrichment analysis
perform_go_analysis <- function(gene_up, gene_down) {
  go_up <- enrichGO(gene = gene_up, OrgDb = "org.Mm.eg.db", keyType = 'SYMBOL', ont = "all")
  go_down <- enrichGO(gene = gene_down, OrgDb = "org.Mm.eg.db", keyType = 'SYMBOL', ont = "all")
  write.table(go_up, file = "go_up_N.csv", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(go_down, file = "go_down_N.csv", sep = "\t", quote = FALSE, row.names = FALSE)
  return(list(go_up = go_up, go_down = go_down))
}

# Function to plot GO enrichment results
plot_go_results <- function(go_up, go_down) {
  pdf(file = "1GO-up.pdf", width = 12, height = 10)
  dotplot(go_up, split = "ONTOLOGY", label_format = 70) + facet_grid(ONTOLOGY ~ ., scale = "free")
  dev.off()
  
  pdf(file = "2GO-down.pdf", width = 12, height = 10)
  dotplot(go_down, split = "ONTOLOGY", label_format = 70) + facet_grid(ONTOLOGY ~ ., scale = "free")
  dev.off()
}

# Function to perform KEGG enrichment analysis
perform_kegg_analysis <- function(deg) {
  gene_up_KEGG <- deg[deg$change == 'up', 'ENTREZID']
  gene_down_KEGG <- deg[deg$change == 'down', 'ENTREZID']
  gene_all <- deg[, 'ENTREZID']
  
  kk_up <- enrichKEGG(gene = gene_up_KEGG, keyType = 'kegg', organism = 'mmu', pAdjustMethod = 'BH', universe = gene_all, qvalueCutoff = 0.25)
  kk_down <- enrichKEGG(gene = gene_down_KEGG, keyType = 'kegg', organism = 'mmu', pAdjustMethod = 'BH', universe = gene_all, qvalueCutoff = 0.25)
  
  write.table(kk_up, file = "KEGG_up_B1C.N.csv", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(kk_down, file = "KEGG_down_N.csv", sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(list(kk_up = kk_up, kk_down = kk_down))
}

# Function to plot KEGG results
plot_kegg_results <- function(kk_up, kk_down) {
  pdf(file = "3barplot_up.pdf", width = 8, height = 12)
  barplot(kk_up, drop = TRUE, showCategory = 20)
  dev.off()
  
  pdf(file = "4barplot_down.pdf", width = 8, height = 12)
  barplot(kk_down, drop = TRUE, showCategory = 20)
  dev.off()
}

# Function to perform GSEA analysis
perform_gsea_analysis <- function(deg) {
  geneList <- deg$avg_log2FC
  names(geneList) <- deg$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  genesets <- msigdbr(species = "Mus musculus", category = "H") %>%
    dplyr::select("gs_name", "gene_symbol") %>%
    as.data.frame()
  
  egmt <- GSEA(geneList, TERM2GENE = genesets, verbose = FALSE, pvalueCutoff = 0.05)
  saveRDS(egmt, "egmt.rds")
  
  pdf(file = "gseapot_H.pdf", width = 12, height = 10)
  gseaplot2(egmt, geneSetID = c(1, 2, 3, 4, 5))
  dev.off()
  
  return(egmt)
}

# Main function to run the analysis pipeline
run_analysis_pipeline <- function(data_file) {
  # Load data
  data <- load_data(data_file)
  
  # Subset for Monocytes
  pro1 <- subset_monocytes(data)
  
  # Perform differential expression analysis
  deg <- perform_deg(pro1)
  
  # Classify genes as up or down-regulated
  deg <- classify_genes(deg)
  
  # Convert gene names to Entrez IDs
  deg <- convert_gene_ids(deg)
  
  # Perform GO analysis
  gene_up <- deg[deg$change == 'up', 'gene']
  gene_down <- deg[deg$change == 'down', 'gene']
  go_results <- perform_go_analysis(gene_up, gene_down)
  
  # Plot GO results
  plot_go_results(go_results$go_up, go_results$go_down)
  
  # Perform KEGG analysis
  kegg_results <- perform_kegg_analysis(deg)
  
  # Plot KEGG results
  plot_kegg_results(kegg_results$kk_up, kegg_results$kk_down)
  
  # Perform GSEA analysis
  perform_gsea_analysis(deg)
}

# Run the analysis pipeline
run_analysis_pipeline("/Users/mq/Downloads/Mice:Blood/mergedata.rds")
