rm(list = ls())
setwd()
getwd()
#load package----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(tidyverse)
library(reshape2)
library(patchwork)
library(Seurat)
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

#Load data--------------------------------------------------------------------
scRNA<- readRDS(file = "/Volumes/YMQ/mac/scRNA12R.rds")
view(scRNA@meta.data)
dim(scRNA)

#Cell filter------------------------------------------------------------------
#scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 15)
pdf("1.1.nFeature_RNA.pdf",width = 10,height = 6)
VlnPlot(object = scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
head(scRNA@meta.data)
scRNA


#Remove mitochondrial and ribosomal genes
counts <- GetAssayData(scRNA, assay = "RNA")
mt.genes <- rownames(scRNA)[grep("^MT-",rownames(scRNA))]
rb.genes <- rownames(scRNA)[grep("^RPL|^RPS|^MRPL|^MRPS",rownames(scRNA))]
dim(scRNA)
view(scRNA@meta.data)
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes))),]
scRNA <- subset(scRNA, features = rownames(counts))
dim(scRNA)

#Sequencing depth--------------------------------------------------------------
plot1 <- FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
pdf("2.0 featureCor.pdf",width = 10,height = 6)
plot1+plot2
dev.off()

#gene expression---------------------------------------------------------------
pdf(file=" 3.0 colSums.pdf",width=10,height=6)
hist(colSums(scRNA$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")
dev.off()

#standardization

scRNA <- NormalizeData(object = scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#Average gene expression after normalization
pdf(file=" 3.1 colSums.pdf",width=10,height=6)
hist(colSums(scRNA$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")
dev.off()
#Looking for hypervariable genes-----------------------------------------------

scRNA <- FindVariableFeatures(object = scRNA,selection.method = "vst", nfeatures = 2000)

#The 10 most mutated genes
top10 <- head(x = VariableFeatures(object = scRNA), 10)
pdf(file=" 4.0.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

CombinePlots(plots = list(plot1, plot2))
dev.off()

pdf("4.1.LabelPoints.pdf",width = 10,height = 6)
plot2
dev.off()

#Data normalization------------------------------------------------------------
all.gene <- rownames(scRNA)
scRNA <- ScaleData(scRNA,features = all.gene)


#Data dimensionality reduction-------------------------------------------------
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))

#Select dimensions
# method 1 ElbowPlot
ElbowPlot(scRNA)
pdf(file=" 5.0.pca.pdf",width=10,height=8)
ElbowPlot(scRNA)
dev.off()
plot
pdf(file=" 5.1.1pcaGene.pdf",width=20,height=20)
VizDimLoadings(object = scRNA, dims = 1:20, reduction = "pca",nfeatures = 20)
dev.off()
#method 2 Heatmap
DimHeatmap(scRNA, dims = 1:20,  cells = 500, balanced = TRUE)
#method 3 JackStraw
scRNA <- JackStraw(scRNA, num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA,dims = 1:20)
pdf(file=" 5.2.pca.pdf",width=10,height=8)
JackStrawPlot(scRNA, dims = 1:20)
dev.off()

#Remove batch effect-----------------------------------------------------------
system.time({scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")})
names(scRNA@reductions)

#reduction---------------------------------------------------------------------
#Calculate the nearest neighbor graph based on the Euclidean distance in PCA space and optimize the distance weight between any two cells 
#(enter the PC dimension obtained in the previous step)
scRNA <- FindNeighbors(scRNA, dims = 1:20,reduction = "harmony")

#Optimize the model. 
#The resolution parameter determines the number of clusters obtained by downstream cluster analysis. For about 3k cells, setting it to 0.4-1.2 can get better results (official description); 
#if the amount of data increases, this parameter should also be increased appropriately. .
scRNA <- FindClusters(scRNA, resolution = 0.5)
##Nonlinear dimensionality reduction
scRNA <- RunUMAP(scRNA,dims = 1:20,reduction = "harmony")
scRNA <- RunTSNE(scRNA,dims = 1:20,reduction = "harmony")

#Find grouping parameters
res.used <- seq(0.5,1.5,by=0.2)
res.used
# Loop over and perform clustering of different resolutions
for(i in res.used){
  sce <- FindClusters(object = scRNA, verbose = T, resolution = res.used)
}
# Make plot
install.packages("clustree")
library(clustree)
pdf("Epithelium Cells.tree.pdf")
clustree(sce)
dev.off()
#theme(legend.position = “bottom”)
# scale_color_brewer(palette = “Set1") +
# scale_edge_color_continuous(low = “grey80”, high = “red”)

#plot
pdf("6.0.seurat_clusters_umap.pdf",width = 10,height = 8)
DimPlot(scRNA, reduction = "umap",raster=FALSE)
dev.off()

pdf("6.1.seurat_clusters_tsne.pdf",width = 10,height = 8)
DimPlot(scRNA, reduction = "tsne",raster=FALSE)
dev.off()

#group
pdf("6.2.umap_group.pdf",width = 8,height = 6)
DimPlot(scRNA, reduction = "umap", group.by = 'group')
dev.off()

pdf("6.3.tsne_group.pdf",width = 8,height = 6)
DimPlot(pb, reduction = "tsne", group.by = 'group')
dev.off()
head(scRNA@meta.data)


#Find and save differential genes----------------------------------------------
subc<-levels(x= scRNA@active.ident)
for (l in subc){
  cluster.markers <- FindMarkers(object= scRNA, ident.1=l, min.pct=0.25,logfc.threshold = 0.25)
  write.table(cluster.markers,file=paste(l,'_diffgenes.xls',sep=''),sep='\t',quote=F,row.names=T)}
dev.off()
head(scRNA)

#top10 differential genes
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pb, features = top20$gene) + NoLegend()#展示前10个标记基因的热图
saveRDS(scRNA, file = "scRNA_harmony.rds")

# devtools::install_github("eddelbuettel/harmony",force = TRUE)
# devtools::install_version("Rcpp", version = "1.0.7", repos = "http://cran.us.r-project.org/")
# packageVersion("harmony")
# packageVersion("Rcpp")

#cell annotation------------------------------------------------------------
#singleR
#load package
library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)
#database
#Mouse.se=MouseRNAseqData()##mouse
hpca.se=HumanPrimaryCellAtlasData() ##human
Blue.se=BlueprintEncodeData()##human
Immune.se=DatabaseImmuneCellExpressionData()##human
Nover.se=NovershternHematopoieticData()##human
MonacoIm.se=MonacoImmuneData()##human
#ImmGen.se=ImmGenData() #(mouse)

saveRDS(hpca.se,'hpca.se.rds')
saveRDS(Blue.se,'Blue.se.rds')
saveRDS(Immune.se,'Immune.se.rds')
saveRDS(Nover.se,'Nover.se.rds')
saveRDS(MonacoIm.se,'MonacoIm.se.rds')

#load data
scRNA_harmony<- readRDS("scRNA_harmony.rds")
data <- GetAssayData(scRNA_harmony, slot="data") ##Get normalized matrix
#annotate by cluster :resolution=0.5:RNA_snn_res.0.5
Idents(scRNA_harmony) <- 'RNA_snn_res.0.5'
pred.cluster <- SingleR(test = data, ref = list(Bu=Blue.se, hp=hpca.se, Im=Immune.se, No=Nover.se, Mo=MonacoIm.se), clusters=Idents(scRNA_harmony), labels = list(Blue.se$label.main,hpca.se$label.main,Immune.se$label.main,Nover.se$label.main,MonacoIm.se$label.main)) #按照细胞群进行注释
#pred.cluster <- SingleR(test = data, ref = list(HPCA=hpca.se), clusters=Idents(scRNA_harmony),labels = list(hpca.se$label.main)) #按照细胞群进行注释
celltype <-data.frame(ClusterID=rownames(pred.cluster),celltype_SR=pred.cluster$labels,stringsAsFactors = F)
head(celltype)
celltype
scRNA_harmony[['celltype_R']]<- celltype$celltype_SR[match(Idents(scRNA_harmony), celltype$ClusterID)]
head(scRNA_harmony@meta.data)
#plot
pdf(file="7.1.umap.pdf",width=12,height = 6)
DimPlot(scRNA_harmony, group.by = c("celltype_R"),reduction = "umap",raster=FALSE,label=TRUE) #画图
dev.off()
pdf(file="7.2.tsne.pdf",width=12,height = 6)
DimPlot(scRNA_harmony, group.by = c("celltype_R"),reduction = "tsne",raster=FALSE,label=TRUE) #画图
dev.off()

#group
view(scRNA_harmony@meta.data)
Idents(scRNA_harmony) <- 'group'
CON <- subset(scRNA_harmony, idents="Control")  
PAH <- subset(scRNA_harmony, idents="PAHs")  

plot1 <- DimPlot(CON, group.by = c("celltype_R"),reduction = "umap",raster=FALSE,label=TRUE)
plot2 <- DimPlot(PAH, group.by = c("celltype_R"),reduction = "umap",raster=FALSE,label=TRUE)
pdf(file="7.3.umap.pdf",width=12,height = 6)
plot1+plot2
dev.off()

plot1 <- DimPlot(CON, group.by = c("celltype_R"),reduction = "tsne",raster=FALSE,label=TRUE)
plot2 <- DimPlot(PAH, group.by = c("celltype_R"),reduction = "tsne",raster=FALSE,label=TRUE)
pdf(file="7.4.tsne.pdf",width=12,height = 6)
plot3+plot4
dev.off()

Idents(scRNA_harmony) <- 'celltype_R'
pdf(file="7.5 DOTPLOT.pdf",width=12,height = 10)
DotPlot(object = scRNA_harmony, feature = c("FCGR3B","CSF3R","LYZ","CD68","CD163","CD14","FGFBP2","FCG3RA","CX3CR1","CD8A","CD3D","CD3E","CD79A","MS4A1","II6","Gata2","Cpa3","FCER1A","CST3"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

Idents(scRNA_harmony) <- 'RNA_snn_res.0.5'

pdf(file="7.6 featureplot.pdf",width=12,height = 10)
FeaturePlot(object = scRNA_harmony, feature = c("FCGR3B","CSF3R","LYZ","CD68","CD163","CD14","FGFBP2","FCG3RA","CX3CR1","CD8A","CD3D","CD3E","CD79A","MS4A1","II6","Gata2","Cpa3","FCER1A","CST3"),raster=FALSE)
dev.off()


#Rating of cell annotations-----------------------------------------------------

library(SummarizedExperiment)
library(MatrixGenerics)
library(matrixStats)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(Biobase)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scuttle)

print(plotScoreHeatmap(pred.cluster))

pdf("8.0 singleR-plotScoreHeatmap.pdf",width=20,height = 10)
print(plotScoreHeatmap(pred.cluster))
dev.off()

#Manual annotation-cellmarker
Idents(scRNA)<-"RNA_snn_res.1.1"
Celltype<-c("Naive B cell","Naive B cell","Memory B cell","Naive B cell","Naive B cell", "Memory B cell","Activated T cell","Neutrophils","CD8+ Tem","Erythrocyte","B cell progenitor")
names(Celltype) <- levels(scRNA)
scRNA <- RenameIdents(scRNA,Celltype) 
scRNA[["Celltype"]] <- Idents(object = scRNA)
view(scRNA)

#plot--------------------------------------------------------------------------
alldata=scRNA_harmony
display.brewer.all()
colourCount = length(table(scRNA_harmony@meta.data$seurat_clusters))
pdf("9.1.seurat_clusters_tsne.pdf")
DimPlot(alldata, 
        reduction = "tsne", raster=FALSE,
        cols = colorRampPalette(brewer.pal(20, "Paired"))(colourCount),#当颜色参考版颜色不够
        group.by  = "seurat_clusters",
        label = T) +  
  NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())
dev.off()

pdf("9.2.seurat_clusters_umap.pdf")
DimPlot(alldata, 
        reduction = "umap", raster=FALSE,
        cols = colorRampPalette(brewer.pal(20, "Paired"))(colourCount),#当颜色参考版颜色不够
        group.by  = "seurat_clusters",
        label = T) +  
  NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())
dev.off()


alldata=scRNA_harmony
display.brewer.all()
colourCount = length(table(scRNA_harmony@meta.data$celltype))
pdf("9.3.cell_type_tsne.pdf")
DimPlot(alldata, 
        reduction = "tsne", raster=FALSE,
        cols = colorRampPalette(brewer.pal(20, "Paired"))(colourCount),#当颜色参考版颜色不够
        group.by  = "celltype_R",
        label = T) +  
  NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())
dev.off()

pdf("9.4.cell_type_umap.pdf")
DimPlot(alldata, 
        reduction = "umap", raster=FALSE,
        cols = colorRampPalette(brewer.pal(20, "Paired"))(colourCount),#当颜色参考版颜色不够
        group.by  = "celltype_R",
        label = T) +  
  NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())
dev.off()

alldata=scRNA_harmony
display.brewer.all()
colourCount = length(table(scRNA_harmony@meta.data$orig.ident))
pdf("9.5.orig.ident_tsne.pdf")
DimPlot(alldata, raster=FALSE,
        reduction = "tsne", 
        cols = colorRampPalette(brewer.pal(6, "Paired"))(colourCount),#当颜色参考版颜色不够
        group.by  = "orig.ident",
        label = T) +  
  NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())
dev.off()

pdf("9.6.orig.ident_umap.pdf")
DimPlot(alldata, 
        reduction = "umap", raster=FALSE,
        
        group.by  = "orig.ident",
        label = T) +  
  NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())
dev.off()

a=scRNA_harmony@meta.data
alldata=scRNA_harmony
display.brewer.all()
colourCount = length(table(scRNA_harmony@meta.data$group))
pdf("9.7.group_tsne.pdf")
DimPlot(alldata, 
        reduction = "tsne", raster=FALSE,
        cols = colorRampPalette(brewer.pal(3, "Paired"))(colourCount),#当颜色参考版颜色不够
        group.by  = "group",
        label = T) +  
  NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())
dev.off()

pdf("9.8.group_umap.pdf")
DimPlot(alldata, 
        reduction = "umap", raster=FALSE,
        cols = colorRampPalette(brewer.pal(3, "Paired"))(colourCount),#当颜色参考版颜色不够
        group.by  = "group",
        label = T) +  
  NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),   
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())

dev.off()

#Cell type proportion
table(scRNA_harmony$orig.ident)
prop.table(table(Idents(scRNA_harmony)))
Idents(scRNA_harmony) <- 'celltype_R'
#Count Cell numbers of different celltype in each group
table(Idents(scRNA_harmony), scRNA_harmony$orig.ident)
#Calculate the proportion of different cell populations in each group of samples
Cellratio <- prop.table(table(Idents(scRNA_harmony), scRNA_harmony$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(table(scRNA_harmony@meta.data$celltype))
allcolour=c(colorRampPalette(brewer.pal(20, "Paired"))(colourCount))
library(ggplot2)
pdf (file="9.9.bar.pdf",width=12,height=6)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()
scRNA_harmony@meta.data

#group
table(scRNA_harmony$group)
prop.table(table(Idents(scRNA_harmony)))
Idents(scRNA_harmony) <- 'celltype_R'
table(Idents(scRNA_harmony), scRNA_harmony$group)
Cellratio <- prop.table(table(Idents(scRNA_harmony), scRNA_harmony$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(table(scRNA_harmony@meta.data$celltype))
allcolour=c(colorRampPalette(brewer.pal(20, "Paired"))(colourCount))
library(ggplot2)
pdf (file="9.10.group.bar.pdf",width=12,height=6)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()

#save---------------------------------------------------------------------------
rate<-table(Idents(scRNA_harmony), scRNA_harmony$orig.ident)
write.table(rate,file="01.cluster-rate.txt",sep="\t",row.names=T,quote=F)
write.table(scRNA_harmony@meta.data,file="02.cluster-cell.txt",sep="\t",row.names=F,quote=F)

saveRDS(scRNA_harmony, file ='scRNA12R.rds')

scRNA12R <- readRDS("scRNA12R.rds")
view(scRNA12R@meta.data)

#marker gene--------------------------------------------------------------------

table(Idents(scRNA12R))
Idents(scRNA12R) <- "celltype_R"
table(Idents(scRNA12R))    

cluster.markers <- FindMarkers(object=scRNA12R,ident.1="Neutrophils",ident.2 = "Monocytes",ident.3 = "NK cells",ident.4 = "CD8+ T-cells",ident.5 = " B-cells",ident.6 = " Basophils ",ident.7 = " Dendritic cells  ", logfc.threshold= 0.25,min.pct = 0.25)
top10 <- cluster.markers %>% group_by(ident.1="Neutrophils",ident.2 = "Monocytes",ident.3 = "NK cells",ident.4 = "CD8+ T-cells",ident.5 = " B-cells",ident.6 = " Basophils ",ident.7 = " Dendritic cells") %>% top_n(n = 10, wt = avg_log2FC)

top10
pdf("10.0.cluster.markers.pdf",width = 10,height = 8)
DoHeatmap(object = scRNA12R, features = rownames(top10), group.by = (ident.1="Neutrophils",ident.2 = "Monocytes",ident.3 = "NK cells",ident.4 = "CD8+ T-cells",ident.5 = " B-cells",ident.6 = " Basophils ",ident.7 = " Dendritic cells")
          dev.off()
      
          VlnPlot(scRNA12R, features = c("MS4A1", "CD79A"),raster=FALSE) 
          FeaturePlot(scRNA12R, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),raster=FALSE)
          

