#转录因子分析 ----------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c("zoo", "mixtools","DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"),ask = F,update = F) 
remotes::install_github("bokeh/rbokeh")
install.packages("doMC", repos="http://R-Forge.R-project.org")

if (!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

#arrow,doRNG,RcisTarget包需要把软件包下载下来后本地安装。

library(Seurat)
library(tidyverse)
library(foreach)
library(doParallel)
library(SCopeLoomR)
library(zoo)
library(mixtools)
library(rbokeh)
library(DT)
library(NMF)
library(ComplexHeatmap)
library(rbokeh)
library(R2HTML)
library(Rtsne)
library(doMC)
library(doRNG)
library(RcisTarget)
library(AUCell)
library(GENIE3)
library(arrow)

if (!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")
library(SCENIC)

# 5.1 #准备表达矩阵 -------------------------------------------------------------
#准备表达矩阵
pbmc <- readRDS(file = "/Volumes/YMQ/PHD/Neutrophils/scRNA_Neucell.rds")
pbmc <- pro1
head(pbmc@meta.data)
Idents(pbmc) <- "RNA_snn_res.0.5"
levels(pbmc)

#为了节省计算资源，可随机抽取1000个细胞的数据子集，也可用所以细胞
subcell <- sample(colnames(pbmc),1000)
scRNAsub <- pbmc[,subcell]

exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
dim(exprMat)

cellInfo <-  scRNAsub@meta.data[,c(9,3,2)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
saveRDS(cellInfo, file="int/cellInfo.Rds")

#设置颜色
colVars <- list(CellType=c("B-cells"="#FF97FF",
                           "CD4+ T-cells"="#FFA15A",
                           "CD8+ T-cells"="#FF6692",
                           "Monocytes"="#AB63FA",
                           "NK cells"="#2E91E5"))

colVars$CellType <-colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
head(colVars)
saveRDS(colVars, file="int/colVars.Rds")

# 5.2 设置分析环境 --------------------------------------------------------------
#下载参考基因数据库
# dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
#              "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather")
# # mc9nr: Motif collection version 9: 24k motifs
# 
# dir.create("cisTarget_databases");   #创建一个文件夹保存数据库
# setwd("cisTarget_databases")
# #如果3个参考数据库都想下载，每次设置变量dbFiles后，都要运行以下代码
# for(featherURL in dbFiles)
# {
#   download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
# }

##手动下载参考基因数据库（网址如上），放在cisTarget_databases文件夹中
mydbDIR <- "/Volumes/YMQ/cisTarget_databases"
list.files(path = mydbDIR)
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather",
           "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")


data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=4,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "SCENIC")

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds" #添加注释信息
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds" #添加细胞颜色

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


# 5.3 共表达网络计算 -------------------------------------------------------------
##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)

exprMat_filtered <- exprMat[genesKept, ]
head(exprMat_filtered)

##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)

##TF-Targets相关性回归分析,#推断TF靶点
exprMat_filtered_log <- log2(exprMat_filtered+1)

set.seed(1234)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 4)
#这一步消耗的计算资源非常大，个人电脑需要几个小时的运行时间
#查看
Gen <- readRDS(file = "F:\\Single cell train\\SCENIC\\int\\1.4_GENIE3_linkList.Rds")
#TF是转录因子名称，Target是潜在靶基因的名字，weight是TF与Target之间的相关性权重。
head(Gen)


# 5.4 推断共表达模块 -------------------------------------------------------------
# 推断共表达模块
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
# scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
Gen2 <- readRDS(file = "F:\\Single cell train\\SCENIC\\int\\1.6_tfModules_asDF.Rds")
#method是上面提到的6种方法，corr是runCorrelation(exprMat_filtered, scenicOptions)命令得到的。
#1代表激活，-1代表抑制，0代表中性，SCENIC只会采用corr值为1的数据用于后续分析
head(Gen2)
saveRDS(scenicOptions,'scenicOptions.rds') 


#Motif验证共表达模块，推断转录调控网络（regulon）
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量
#可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
# #后续消耗内存很大，关闭后重启
# scenicOptions <- readRDS(file = "F:\\Single cell train\\SCENIC\\scenicOptions.rds")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget"))

saveRDS(scenicOptions,'scenicOptions2.rds') 

MotifEnrich <- read.table(file = "F:\\Single cell train\\SCENIC\\output\\Step2_MotifEnrichment.tsv",header=T, sep="\t")
regulonTarget <- read.table(file = "F:\\Single cell train\\SCENIC\\output\\Step2_regulonTargetsInfo.tsv",header=T, sep="\t")

# 5.5 regulon活性评分与可视化 -----------------------------------------------------
##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
exprMat_all <- as.matrix(pbmc@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat)
scenicOptions_all <- runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

#查看调整阈值
#使用shiny互动调整阈值
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
savedSelections <- shiny::runApp(aucellApp)
#保存调整后的阈值
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#生成二分类的regulon活性矩阵
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)

saveRDS(scenicOptions,'scenicOptions.rds') 


# 5.6 可视化分析 ---------------------------------------------------------------
#寻找核心基因和转录因子
cellInfo <- readRDS("int/cellinfo.Rds")
head(cellInfo)

Celltype = subset(cellInfo,select = 'CellType')
Celltype <- Celltype %>% arrange(Celltype$CellType)
cell_dd <- rownames(Celltype)

head(Celltype)

Regulon <-readRDS("F:\\Single cell train\\SCENIC\\int\\3.4_regulonAUC.Rds")
Regulon <- Regulon@assays@data@listData$AUC
head(Regulon)

Regulon_all <- Regulon[,cell_dd] #选择目标细胞
head(Regulon_all)

#画图
library(pheatmap)
library(eoffice)

pheatmap(Regulon_all, show_colnames=F, 
         annotation_col=Celltype)
topptx(filename = 'myAUCmatrix_heatmap.pptx',height = 8,width = 6)#保存为PPT

#整体看一下，然后挑选你觉得有意义且差异大的转录因子继续做热图
cols <- readRDS('F:\\Single cell train\\SCENIC\\int\\colVars.Rds') #自定义的颜色
cols

##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(pbmc, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')
head(scRNAauc@meta.data)
##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(pbmc, BINmatrix)
scRNAbin@assays$integrated <- NULL
head(scRNAbin@meta.data)
saveRDS(scRNAbin, 'scRNAbin.rds')

##目标转录因子，利用Seurat可视化AUC,
dir.create('scenic_seurat')
#FeaturePlot
p1 = FeaturePlot(scRNAauc, features='ATF1_11g', label=T, reduction = 'tsne')
p2 = FeaturePlot(scRNAbin, features='ATF1_11g', label=T, reduction = 'tsne')
p3 = DimPlot(pbmc, reduction = 'tsne',label=T)
plotc = p1|p2|p3
ggsave('scenic_seurat/CEBPB_extended_2290g.png', plotc, width=14 ,height=4)

#RidgePlot&VlnPlot
p1 = RidgePlot(scRNAauc, features = "CEBPB_extended_71g", group.by="celltype") + 
  theme(legend.position='none')
p2 = VlnPlot(scRNAauc, features = "CEBPB_extended_71g", pt.size = 0, group.by="celltype") + 
  theme(legend.position='none')
plotc = p1 + p2
ggsave('scenic_seurat/Ridge-Vln_CEBPB_extended_2290g.png', plotc, width=10, height=8)


#热图
library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'CellType')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#挑选部分感兴趣的regulons
my.regulons <- c('ETS1_2372g','ETV7_981g','IRF7_239g','XBP1_854g','ATF4_37g',
                 'KLF13_78g','ATF6_129g','CREB3L2_619g','TAGLN2_13g',
                 'STAT1_extended_1808g','CEBPB_extended_2290g','IRF5_extended_422g',
                 'SPI1_1606g','HMGA1_14g','SPIB_1866g','IRF8_348g','BCL11A_136g',
                 'EBF1_40g','MAF_45g','BATF_131g','FOXP3_55g','TBX21_388g',
                 'EOMES_extended_101g','TCF7_extended_31g','LEF1_extended_49g')
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
#使用regulon原始AUC值绘制热图
pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype,
         filename = 'scenic_seurat/myAUCmatrix_heatmap.png',
         width = 6, height = 5)
#使用regulon二进制AUC值绘制热图
pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         filename = 'scenic_seurat/myBINmatrix_heatmap.png',
         color = colorRampPalette(colors = c("white","black"))(100),
         width = 6, height = 5)