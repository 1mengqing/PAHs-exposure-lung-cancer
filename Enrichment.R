#load library
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
library(msigdbr) #database
library(pheatmap)
library(Rgraphviz)
#load data
pb <- readRDS("/public/home/yangbin/Analysis1/Ben/Seurat/filter.pb.rds")
pb <- readRDS("pb_new.rds")
head(pb@meta.data)
# names(pbmc@meta.data)[4] <- "celltype"
Idents(pb) <- "group"
#different genes
cells1 <- subset(pb@meta.data, group %in% c("B"))  %>% rownames()
cells2 <- subset(pb@meta.data, group %in%  c("C"))  %>% rownames()
deg <- FindMarkers(pb, ident.1 = cells1, ident.2 = cells2)
deg <- data.frame(gene = rownames(deg), deg)
head(deg)


deg <- read.table("deg.txt", header=T, sep="\t", check.names=F)     

#Filter differential genes
deg1 <- deg
#logFC_t=0.5
#P.Value_t = 0.05
k1 = (deg1$p_val_adj < 0.05)&(deg1$avg_log2FC < -0.5) #down gene
k2 = (deg1$p_val_adj < 0.05)&(deg1$avg_log2FC > 0.5) #up gene
table(k1)
table(k2)

#Grouping up/down
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg1$change <- change
head(deg1)
table(deg1$change)


write.table(deg1, file="deg1.txt", sep="\t", quote=F, row.names = F)
write.table(deg1, file="degg1.txt", sep="\t", quote=F, row.names = T)

#DEGs
#NK_cell"     "Monocyte"    "Neutrophils" "HSC_-G-CSF"  "T_cells"    "Myelocyte"   "B_cell"
T_cell <- subset(pb, celltype == "T_cells")

T_cells1 <- subset(T_cell@meta.data, group %in% c("B"))  %>% rownames()
T_cells2 <- subset(T_cell@meta.data, group %in%  c("C"))  %>% rownames()
T_deg <- FindMarkers(T_cell, ident.1 = T_cells1, ident.2 = T_cells2)
T_deg <- data.frame(gene = rownames(T_deg), T_deg)
head(T_deg)
write.table(T_deg, file="T_deg.txt", sep="\t", quote=F, row.names = F)

B_cell <- subset(pb, celltype == "B_cell")
diff_B <- FindMarkers(B_cell, 
                      group.by = "group",
                      ident.1 ="B",
                      ident.2="C")  
write.table(diff_B, file="diff_B.txt", sep="\t", quote=F, row.names = F)

N_cell <- subset(pb, celltype == "Neutrophils")
diff_N <- FindMarkers(N_cell, 
                      group.by = "group",
                      ident.1 ="B",
                      ident.2="C")  
write.table(diff_N, file="diff_N.txt", sep="\t", quote=F, row.names = F)

Mono_cell <- subset(pb, celltype == "Monocyte")
diff_Mono <- FindMarkers(Mono_cell, 
                         group.by = "group",
                         ident.1 ="B",
                         ident.2="C")  
write.table(diff_Mono, file="diff_Mono.txt", sep="\t", quote=F, row.names = F)

NK_cell <- subset(pb, celltype == "NK_cell")
diff_NK <- FindMarkers(NK_cell, 
                       group.by = "group",
                       ident.1 ="B",
                       ident.2="C")  
write.table(diff_NK, file="diff_NK.txt", sep="\t", quote=F, row.names = F)


My_cell <- subset(pb, celltype == "Myelocyte")
diff_My <- FindMarkers(My_cell, 
                       group.by = "group",
                       ident.1 ="B",
                       ident.2="C")  
write.table(diff_My, file="diff_My.txt", sep="\t", quote=F, row.names = F)

saveRDS(pb, file = "pb_new.rds")
saveRDS(T_cell, file = "T_cell.rds")
saveRDS(B_cell, file = "B_cell.rds")
saveRDS(Mono_cell, file = "Mono_cell.rds")
saveRDS(N_cell, file = "N_cell.rds")
saveRDS(NK_cell, file = "NK_cell.rds")
saveRDS(My_cell, file = "My_cell.rds")
#Gene name conversion, ID conversion, gene_symbol into ENTREZID.
#keytypes(org.Hs.eg.db) #View genes
s2e <- bitr(deg1$gene, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
head(s2e)
deg1 <- inner_join(deg1,s2e,by=c("gene"="SYMBOL")) 
head(deg1)

#GO
gene_up = deg1[deg1$change == 'up','gene'] 
gene_down = deg1[deg1$change == 'down','gene'] 
gene_diff = c(gene_up,gene_down)

#KEGG
gene_all = deg1[,'ENTREZID']
gene_up_KEGG = deg1[deg1$change == 'up','ENTREZID']
gene_down_KEGG = deg1[deg1$change == 'down','ENTREZID']
gene_diff_KEGG = c(gene_up_KEGG,gene_down_KEGG)

#plot
#go
#cellular components
ego_CC <- enrichGO(gene          = gene_up,
                   keyType       = 'SYMBOL', 
                   OrgDb         = org.Hs.eg.db,  
                   ont           = "CC",
                   pAdjustMethod = "BH", 
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

#biological process
ego_BP <- enrichGO(gene          = gene_up,
                   OrgDb          = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

#molecular function
ego_MF <- enrichGO(gene          = gene_up,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

save(ego_CC,ego_BP,ego_MF,file = "GO.Rdata")


pdf(file = "GO_BP_plotGO.pdf",width = 12,height = 10)
plotGOgraph(ego_BP) #topGO包
dev.off()

pdf(file = "GO_BP_goplot.pdf",width = 12,height = 10)
goplot(ego_BP)
dev.off()


#Combining cellular components, molecular functions, and biological processes
go <- enrichGO(gene = gene_up, OrgDb = "org.Hs.eg.db", keyType  = 'SYMBOL',ont="all")

head(go)
#plot
dotplot(ego_CC, showCategory=30)#showCategory
barplot(ego_CC) #top8
p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p

p=dotplot(ALL,x="count",
          showCategory = 14,colorBy="pvalue") #showCategory
#4 KEGG
# KEGG
#Enrichment of up-regulated genes
kk.up <- enrichKEGG(gene         = gene_up_KEGG, #entrzeid
                    organism     = 'hsa',
                    universe     = gene_all, ##
                    pvalueCutoff = 0.9, ##
                    qvalueCutoff = 0.9)
#Enrichment of down-regulated genes
kk.down <- enrichKEGG(gene         =  gene_down_KEGG,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
kk.diff <- enrichKEGG(gene         = gene_diff_KEGG,
                      organism     = 'hsa',
                      pvalueCutoff = 0.9)
save(kk.diff,kk.down,kk.up,file = "kegg.Rdata")

#
ekegg <- setReadable(kk.up, OrgDb = org.Hs.eg.db, keyType="ENTREZID") 
head(ekegg)

#plot
p1 <- barplot(ekegg, showCategory=10)
p2 <- dotplot(ekegg, showCategory=10)
plotc = p1/p2
plotc
ggsave(plotc, file='plotc.pdf', width=12, height=10)
#cnetplot 
#circluar
plot2 <- enrichplot::cnetplot(ekegg,circular=FALSE,colorEdge = TRUE)
ggsave(plot2, file='plot2.pdf', width=12, height=10) 

plot3 <- enrichplot::cnetplot(ekegg,circular=TRUE,colorEdge = TRUE)
ggsave(plot3, file='plot3.pdf', width=12, height=10) 

#Gene-pathway correlation heat map
enrichplot::heatplot(ekegg,showCategory = 5)

#5 GSEA
#5.1 GSEA-clusterProfiler
#
geneList = deg1$avg_log2FC 
names(geneList) = deg1$ENTREZID
geneList = sort(geneList,decreasing = T)
geneList[1:10]

## gsea-KEGG
kk_gse <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  keyType = 'kegg',
  minGSSize = 10,
  pvalueCutoff = 0.9,
  verbose = FALSE  
)

kk_gse=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')  #基因ID转换
sortkk<-kk_gse[order(kk_gse$enrichmentScore, decreasing = T),] #排序
head(sortkk)
#plot
pdf(file = "gseapot.pdf",width = 12,height = 10)
gseaplot2(kk_gse, 
          "hsa03010", 
          color = "firebrick")
dev.off()
#top4 pathway
pdf(file = "gseapot2.pdf",width = 12,height = 10)
gseaplot2(kk_gse, row.names(sortkk)[1:4]) 
dev.off()

#5.2 GSEA-GSEAbase
#database
genesets <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select("gs_name","gene_symbol" )%>% as.data.frame()
#https://www.gsea-msigdb.org/gsea/msigdb c2 database


#
egmt <- GSEA(geneLists, TERM2GENE=genesets,verbose=F,pvalueCutoff = 0.5) 
y=data.frame(egmt) 
head(y)
#gseaplot
pdf(file = "gseapot.pdf",width = 12,height = 10)
gseaplot2(egmt, geneSetID = 1, title = egmt$Description[1])

