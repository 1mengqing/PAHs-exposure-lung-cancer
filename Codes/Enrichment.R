
#Load package
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
library(msigdbr) #
library(pheatmap)
library(Rgraphviz)
#Load data
scRNA <- readRDS("~/biny/ymq/12sample/2024/117757/scRNA_T-cells.rds")

head(scRNA@meta.data)
# names(pbmc@meta.data)[4] <- "celltype"
Idents(scRNA) <- "group"
#Obtain differential genes between groups
cells1 <- subset(scRNA@meta.data, group %in% c("PAHs"))  %>% rownames()
cells2 <- subset(scRNA@meta.data, group %in%  c("Control"))  %>% rownames()
deg <- FindMarkers(scRNA, ident.1 = cells1, ident.2 = cells2)
deg <- data.frame(gene = rownames(deg), deg)
head(deg)

write.table(deg, file="deg.txt", sep="\t", quote=F, row.names = F)
deg <- read.table("deg.txt", header=T, sep="\t", check.names=F)     #读

#筛选基因
deg1 <- deg

deg1$change=ifelse(deg1$p_val_adj>0.05,"stable",
                   ifelse(deg1$avg_log2FC>0.5 ,"up",
                          ifelse(deg1$avg_log2FC<(-0.5),"down","stable")))
table(deg1$change)



write.table(deg1, file="deg1.txt", sep="\t", quote=F, row.names = F)



#gene_symbol is converted to ENTREZID
#keytypes(org.Hs.eg.db) 
s2e <- bitr(deg1$gene, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
head(s2e)

deg1 <- inner_join(deg1,s2e,by=c("gene"="SYMBOL")) 
head(deg1)

#GO enrichment[Symbol]
gene_up = deg1[deg1$change == 'up','gene'] #
gene_down = deg1[deg1$change == 'down','gene'] 
gene_diff  = c(gene_up,gene_down)

#KEGG enrichment[ENTREZID]

gene_all = deg1[,'ENTREZID']
gene_up_KEGG = deg1[deg1$change == 'up','ENTREZID']
gene_down_KEGG = deg1[deg1$change == 'down','ENTREZID']
gene_diff_KEGG = c(gene_up_KEGG,gene_down_KEGG)

#GO enrichment 
#CC
ego_CC <- enrichGO(gene          = gene_down,
                   keyType       = 'SYMBOL', 
                   OrgDb         = org.Hs.eg.db,  #database
                   ont           = "CC",
                   pAdjustMethod = "BH", #
                   pvalueCutoff  = 0.05, #
                   qvalueCutoff  = 1) #

#BP
ego_BP <- enrichGO(gene          = gene_down,
                   OrgDb          = org.Hs.eg.db,
                   #universe = 
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 1)

#分子功能
ego_MF <- enrichGO(gene          = gene_down,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 1)

save(ego_CC,ego_BP,ego_MF,file = "GO_down.Rdata")


pdf(file = "GO_BP_plotGO.pdf",width = 12,height = 10)
plotGOgraph(ego_BP) #topGO包
dev.off()

pdf(file = "GO_BP_goplot.pdf",width = 12,height = 10)
goplot(ego_BP)
dev.off()


#细胞组分、分子功能、生物学过程
go_up<- enrichGO(gene = gene_up, OrgDb = "org.Hs.eg.db", keyType  = 'SYMBOL',ont="all")

go_down<- enrichGO(gene = gene_down, OrgDb = "org.Hs.eg.db", keyType  = 'SYMBOL',ont="all")


#plot
dotplot(ego_CC, showCategory=30)
barplot(ego_CC) #top8
pdf(file = "1.2GO-up.pdf",width = 12,height = 10)
p <- dotplot(go_up, split="ONTOLOGY",label_format = 70) +facet_grid(ONTOLOGY~., scale="free")
p
dev.off()
pdf(file = "1.3GO-down.pdf",width = 12,height = 10)
p <- dotplot(go_down, split="ONTOLOGY",label_format = 70) +facet_grid(ONTOLOGY~., scale="free")
p
dev.off()

p=dotplot(go,x="count", showCategory = 14,color="pvalue") #
p

#4 KEGG enrichment

#Enrichment of up-regulated genes
kk.up <- enrichKEGG(gene         = gene_up_KEGG, #注意这里只能用 entrzeid
                    organism     = 'hsa',
                    #universe     = gene_all, ##背景基因集，可省
                    pvalueCutoff = 0.05, ##指定 p 值阈值，不显著的值将不显示在结果中
                    qvalueCutoff = 1)
#Enrichment of down-regulated genes
kk.down <- enrichKEGG(gene         =  gene_down_KEGG,
                      organism     = 'hsa',
                      # universe     = gene_all,
                      pvalueCutoff = 0.05,
                      qvalueCutoff =1)
kk.diff <- enrichKEGG(gene         = gene_diff_KEGG,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
save(kk.diff,kk.down,kk.up,file = "kegg.Rdata")

ekegg1 <- setReadable(kk.up, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #将基因的entrzeid转换成SYMBOL
head(ekegg1)
ekegg2 <- setReadable(kk.down, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #将基因的entrzeid转换成SYMBOL
head(ekegg2)
#plot
p1 <- barplot(ekegg1, showCategory=10,label_format = 70)
p2 <- barplot(ekegg2, showCategory=10,label_format = 70)
p1
p2
ggsave(p1, file='plotc1.pdf',width = 12,height = 10)
ggsave(p2, file='plotc2.pdf',width = 12,height = 10)


#3 GSEA
#3.1 GSEA富集-使用clusterProfiler
deg1
geneList = deg1$avg_log2FC 
names(geneList) = deg1$ENTREZID
deg1$ENTREZID
geneList = sort(geneList,decreasing = T)
geneList[1:10]
geneList

## gsea-KEGG
kk_gse <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  keyType = 'kegg',
  minGSSize = 3,
  pvalueCutoff = 0.05,
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
#the top4 pathway
pdf(file = "gseapot2.pdf",width = 12,height = 10)
gseaplot2(kk_gse, row.names(sortkk)[1:4]) 
dev.off()

#5.2 GSEA-GSEAbase
#database
genesets <- msigdbr(species = "Homo sapiens", category = "C3") %>% dplyr::select("gs_name","gene_symbol" )%>% as.data.frame()
#https://www.gsea-msigdb.org/gsea/msigdb 不同数据集包含的内容，在这里我们选择c2数据集：（专家）校验基因集合，基于通路、文献等：
#可以使用msigdbr_species() 和 msigdbr_collections()查看支持的物种和基因集类别。
msigdbr_species()
#plot
names(geneList) = deg1$gene
egmt <- GSEA(geneList, TERM2GENE=genesets,verbose=F,pvalueCutoff = 0.05)
y=data.frame(egmt) 
head(y)
#gseaplot
pdf(file = "gseapot_C3.pdf",width = 12,height = 10)
gseaplot2(egmt, geneSetID = 1, title = egmt$Description[1])
dev.off()

