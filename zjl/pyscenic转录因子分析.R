getwd()
setwd("~/pyscenic")

####1.下载pyscenic分析的数据集####
## SCENIC需要一些依赖包，先安装好
R.Version()
R.Version()$version.string
#[1] "R version 4.3.1 (2023-06-16)"
#查看版本
BiocManager::version()
#[1] '3.18'
BiocManager::install(c("AUCell", "RcisTarget","GENIE3","zoo", "mixtools", "rbokeh","DT", "NMF", "pheatmap", "R2HTML", "Rtsne","doMC", "doRNG","scRNAseq"))
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")
#BiocManager::install("SeuratObject")
#BiocManager::install("Seurat")
#check
library(SCENIC)
packageVersion("SCENIC")

library(Seurat)
library(SeuratObject)
#remotes::install_github("satijalab/seurat-data")
AvailableData()
InstallData("pbmc3k")

#上面不行  直接参考：https://blog.csdn.net/wmm131333/article/details/134316164
#wget http://seurat.nygenome.org/src/contrib/pbmc3k.SeuratData_3.1.4.tar.gz
# 手动安装该数据集
install.packages('~/pyscenic/pbmc3k.SeuratData_3.1.4.tar.gz', repos = NULL, type = "source")
# 仅首次需要手动导入，后续可以直接加载
library(pbmc3k.SeuratData)	
# 加载该数据集
data("pbmc3k")
pbmc3k
pbmc3k = UpdateSeuratObject(pbmc3k)
pbmc3k

#或者读入测序数据
#pbmc3k = readRDS("./pbmc3k.test.seurat.Rds")
#pbmc3k
#注意矩阵一定要转置，不然会报错
write.csv(t(as.matrix(pbmc3k@assays$RNA@counts)),file = "for.scenic.data.csv")

####2.下面的处置转入python：jupyter->pyscenic####

####3.接 jupyter:out_SCENIC.loom导入R里面进行可视化####
##可视化
rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
install.packages("pacman")
library(pacman)
p_load(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)

#提取 out_SCENIC.loom 信息:
###3.1.提取 out_SCENIC.loom 信息
loom <- open_loom('out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)

###3.2 加载SeuratData
library(pbmc3k.SeuratData)	
# 加载该数据集
data("pbmc3k")
seurat.data = pbmc3k
#接下来直接跑seurat流程
seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = F)

n.pcs = 30
seurat.data <- seurat.data %>% 
  RunUMAP(reduction = "pca", dims = 1:n.pcs, verbose = F) %>% 
  FindNeighbors(reduction = "pca", k.param = 10, dims = 1:n.pcs)

# 这里有自带的注释
seurat.data$seurat_annotations[is.na(seurat.data$seurat_annotations)] = "B"
Idents(seurat.data) <- "seurat_annotations"
DimPlot(seurat.data,reduction = "umap",label=T ) 

###3.3 可视化
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$seurat_annotations))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$seurat_annotations)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

#保存一下
save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')

##B细胞有两个非常出名的转录因子，TCF4(+) 以及NR2C1(+)，接下来就可以对这两个进行简单的可视化。
##首先我们需要把这两个转录因子活性信息 添加到降维聚类分群后的的seurat对象里面。
regulonsToPlot = c('TCF4(+)','NR2C1(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))

# Vis
p1 = DotPlot(seurat.data, features = unique(regulonsToPlot)) + RotatedAxis()
p2 = RidgePlot(seurat.data, features = regulonsToPlot , ncol = 2) 
p3 = VlnPlot(seurat.data, features = regulonsToPlot,pt.size = 0)
p4 = FeaturePlot(seurat.data,features = regulonsToPlot)

wrap_plots(p1,p2,p3,p4)
#可选
#scenic_res = assay(sub_regulonAUC) %>% as.matrix()
#seurat.data[["scenic"]] <- SeuratObject::CreateAssayObject(counts = scenic_res)
#seurat.data <- SeuratObject::SetAssayData(seurat.data, slot = "scale.data",
#                                  new.data = scenic_res, assay = "scenic")

###3.4.亚群特异性转录因子分析及可视化
###3.4.1. TF活性均值
# 看看不同单细胞亚群的转录因子活性平均值
# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
#[1] 209   9
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
#热图查看TF分布：
Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
#可以看到，确实每个单细胞亚群都是有 自己的特异性的激活的转录因子。

###3.4.2. rss查看特异TF
#不过，SCENIC包自己提供了一个 calcRSS函数，帮助我们来挑选各个单细胞亚群特异性的转录因子，全称是：Calculates the regulon specificity score
#参考文章：The RSS was first used by Suo et al. in: Revealing the Critical Regulators of Cell Identity in the Mouse Cell Atlas. Cell Reports (2018). doi: 10.1016/j.celrep.2018.10.045 运行超级简单。
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

###3.4.3 其他查看TF方式
rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)
#这个图好看一点