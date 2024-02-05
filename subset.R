#subset
#Neutrophils 
#Monocytes
#NK cells
#CD8+ T cells
#B cells
#Dendritic cells
#Basophils
#View data and set save path
setwd("/mnt/data/userdata/svip017/biny/ymq/12sample/2024")
head(scRNA_harmony@meta.data)
dim(scRNA_harmony)
Idents(scRNA_harmony) <- "Celltype_R"
table(scRNA_harmony$celltype_R)
#Extract cell subsets
colnames(scRNA_harmony@meta.data)
Neu <- subset(scRNA_harmony, subset = celltype_R == "Neutrophils")
dim(Neu)
saveRDS(Neu, file = "scRNA_Neutrophils.rds")




#seurat
pro1 <- Neu
dim(pro1)
pro1 <- FindVariableFeatures(pro1)
pro1 <- RunPCA(object=pro1,features = VariableFeatures(object = pro1))
pro1 <- FindNeighbors(pro1, reduction = "pca", dims = 1:20)
pro1 <- FindClusters(pro1,resolution = 0.8, algorithm = 1)
pro1 <- RunTSNE(object=pro1,dims.use=1:20,do.fast=TRUE,check_duplicates = FALSE)
pro1 <- RunUMAP(pro1, reduction = "pca", dims = 1:20)

#singleR
data <- GetAssayData(pro1, slot="data")
view(pro1@meta.data)
Idents(pro1) <- 'RNA_snn_res.0.8'
pred.cluster <- SingleR(test = data, ref = list(Bu=Blue.se, hp=hpca.se, Mo=MonacoIm.se), clusters=Idents(pro1), labels = list(Blue.se$label.fine,hpca.se$label.fine,MonacoIm.se$label.fine)) #
celltype <-data.frame(ClusterID=rownames(pred.cluster),celltype_SR=pred.cluster$labels,stringsAsFactors = F)
celltype
pro1[['celltype_NR']]<- celltype$celltype_SR[match(Idents(pro1), celltype$ClusterID)]
head(pro1)

#Subset annotations are added to the original meta
head(scRNA_harmony)
Idents(scRNA_harmony) <- "new_celltype"
Idents(pro1) <- "celltype_BR"
dim(scRNA_harmony)
dim(pro1)

Idents(scRNA_harmony,cells= colnames(pro1)) <- Idents(pro1)
head(scRNA_harmony)
scRNA_harmony$new_celltype <- Idents(scRNA_harmony)

view(scRNA_harmony@meta.data)

saveRDS(scRNA_harmony,"12NKR.rds")