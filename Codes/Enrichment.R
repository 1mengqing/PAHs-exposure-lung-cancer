rm(list = ls())

setwd("/Users/mq/Downloads/Mice:Blood/Monocytes")


library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler) #GO,KEGG,GSEA
library(enrichplot) #GO,KEGG,GSEA
library(GSEABase) #GSEA
library(GSVA) #GSVS
library(DOSE)
library(topGO)
#library(pathview) 
library(msigdbr) #数据库
library(pheatmap)
library(Rgraphviz)
library(ggplot2)




data <- readRDS("/Users/mq/Downloads/Mice:Blood/mergedata.rds")
dim(data)
head(data)
Mo <- subset(data, subset = Celltype == "Monocytes")
#Ncells <- subset(pro1,  subset = Celltype =="Neutrophils")
#Tcells <- subset(pro1,  subset = Celltype =="T cells")
#dim(Tcells)
saveRDS(Mo,"Monocytes.rds")
pro1 <- Mo
dim(pro1)
# group --------------------------------------------------------------------
cells1 <- subset(pro1@meta.data, group %in% c("Bap"))  %>% rownames()
cells2 <- subset(pro1@meta.data, group %in% c("Control"))  %>% rownames()
deg <- FindMarkers(pro1, ident.1 = cells1, ident.2 = cells2)
deg <- data.frame(gene = rownames(deg), deg)#

head(deg)
dim(deg)
#Distinguish between up- and down-regulated genes
colnames(deg)
k1 = (deg$p_val_adj < 0.05)&(deg$avg_log2FC < -0.5) 
k2 = (deg$p_val_adj < 0.05)&(deg$avg_log2FC > 0.5) 
table(k1)
table(k2)

change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg$change <- change
head(deg)
table(deg$change)

#change name to gene id
s2e <- bitr(deg$gene, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Mm.eg.db)#

head(s2e)
deg <- inner_join(deg,s2e,by=c("gene"="SYMBOL")) #
head(deg)

#select gene 

gene_up = deg[deg$change == 'up','gene'] #
gene_down = deg[deg$change == 'down','gene'] #
gene_diff = c(gene_up,gene_down)

pvalueFilter=0.05       
qvalueFilter=0.05       

#color
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}



go_up<- enrichGO(gene = gene_up, OrgDb="org.Mm.eg.db", keyType  = 'SYMBOL',ont="all")
go_down <- enrichGO(gene = gene_down, OrgDb="org.Mm.eg.db", keyType  = 'SYMBOL',ont="all")

write.table(go_up, file="go_up_N.csv", sep="\t", quote=F, row.names = F)
write.table(go_down, file="go_down_N.csv", sep="\t", quote=F, row.names = F)

dim(go_up)
dim(go_down)
#plot

pdf(file = "1GO-up.pdf",width = 12,height = 10)
p <- dotplot(go_up, split="ONTOLOGY",label_format = 70) +facet_grid(ONTOLOGY~., scale="free")
p
dev.off()
pdf(file = "2GO-down.pdf",width = 12,height = 10)
p <- dotplot(go_down, split="ONTOLOGY",label_format = 70) +facet_grid(ONTOLOGY~., scale="free")
p
dev.off()


#KEGG enrichment
#KEGG

gene_all = deg[,'ENTREZID']
gene_up_KEGG = deg[deg$change == 'up','ENTREZID']
gene_down_KEGG = deg[deg$change == 'down','ENTREZID']
gene_diff_KEGG = c(gene_up_KEGG,gene_down_KEGG)

#remove.packages('clusterProfiler')
#devtools::install_github("YuLab-SMU/clusterProfiler")
#install.packages("clusterProfiler", dependencies = TRUE)
#library(clusterProfiler)

kk.up <- enrichKEGG(
  gene = gene_up_KEGG, 
  keyType = 'kegg',  
  organism = 'mmu',  
  pAdjustMethod = 'BH', 
  universe= gene_all,
  #pvalueCutoff = 1,  
  qvalueCutoff = 0.25,  
)

kk.down <- enrichKEGG(
  gene = gene_down_KEGG,        
  keyType = 'kegg',           
  organism = 'mmu',           
  pAdjustMethod = 'BH',       
  universe     = gene_all,
  pvalueCutoff = 1,         
  qvalueCutoff = 0.25,      
)

dim(kk.up)
dim(kk.down)

ekegg_up <- setReadable(kk.up, OrgDb = org.Mm.eg.db, keyType="ENTREZID") #
ekegg_down <- setReadable(kk.down, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
#

KEGG_up=as.data.frame(ekegg_up)
KEGG_down=as.data.frame(ekegg_down)

write.table(KEGG_up, file="KEGG_up_B1C.N.csv", sep="\t", quote=F, row.names = F)
write.table(KEGG_down, file="KEGG_down_N.csv", sep="\t", quote=F, row.names = F)
# showNum=15
# if(nrow(KEGG)<showNum){
#   showNum=nrow(KEGG)
# }

#柱状图
pdf(file="3barplot_up.pdf", width=8, height=12)
barplot(kk.up, drop=TRUE, showCategory=20)
dev.off()

pdf(file="4barplot_down.pdf", width=8, height=12)
barplot(kk.down, drop=TRUE, showCategory=20)
dev.off()


#GSEA

geneList = deg$avg_log2FC
names(geneList) = deg$ENTREZID
geneList = sort(geneList,decreasing = T)
#genesets <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select("gs_name","gene_symbol" )%>% as.data.frame()
genesets <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select("gs_name","gene_symbol" )%>% as.data.frame()
#msigdbr_species()

names(geneList) = deg$gene
egmt <- GSEA(geneList, TERM2GENE=genesets,verbose=F,pvalueCutoff = 0.05) #
y=data.frame(egmt) #
head(y)
y
saveRDS(egmt,"egmt.rds")

pdf(file = "gseapot_H.pdf",width = 12,height = 10)
gseaplot2(egmt, geneSetID = c(1,2,3,4,5))
gseaplot2
dev.off()
