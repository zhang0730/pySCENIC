#安装pyscenic，并创建分析环境
conda create -n pyscenic python=3.7
conda activate pyscenic #激活pyscenic 环境

pip list
#安装依赖包
conda install -y numpy
conda install -y -c anaconda cytoolz
conda install -y scanpy
#安装pyscenic
pip install pyscenic

#==================================================================================================
######数据文件准备
#1.TF注释
#鼠的
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
#人的
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
 
#2.转录组因子列表
#下载地址
#https://github.com/aertslab/pySCENIC/tree/master/resources
#人的文件名：hs_hgnc_tfs.txt，复制为txt文件即可
#鼠的文件名：mm_mgi_tfs.txt，复制为txt文件即可
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt

 
#3.reference数据库，之前一些网上教程的链接文件已经不行了，因为做了更新，跑的时候会出错，我是根据报错选择了下面的文件
#鼠的
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_clustered.genes_vs_motifs.rankings.feather
#人的
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather



#==================================================================================================
######分析文件准备

#文件准备
#pyscenic的输入文件是行为基因名，列为细胞ID的矩阵，所以在seurat对象中导出矩阵的时候需要转置一下，可以用标准化矩阵，也可以用counts矩阵，影响不大！
#表达矩阵、meta----R中进行
write.csv(t(as.matrix(scRNA@assays$RNA@counts)),file = "sce_exp.csv")
#cellInfo <- sce@meta.data[,c("celltype","nCount_RNA","nFeature_RNA")]
#colnames(cellInfo) <- c('CellType', 'nGene' ,'nUMI')
#head(cellInfo)
#write.csv(cellInfo, file = "cellInfo.csv")
 
#转化为loom文件，Linux下的python脚本(用中号服务器)
##注意在响应的虚拟环境下要安装需要的模块
#编辑脚本
vim trans.py
#输入以下内容
import os, sys
os.getcwd()
os.listdir(os.getcwd())
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("sce_exp.csv");#R中导出的表达矩阵
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sce.loom",x.X.transpose(),row_attrs,col_attrs)
 
#保存并退出
#运行trans.py
python trans.py
ls
#这样在文件夹中会出现sce.loom文件，就是接下来输入pyscenic的文件。

#====================================================================================
###分析
##分析第一步：GRN---运行完得到sce.adj.csv文件

pyscenic grn --num_workers 10 \
  --sparse \
  --method grnboost2 \
  --output sce.adj.csv \
  sce.loom \
  allTFs_hg38.txt
  #这一步的目的
  #推断转录因子与提供的表达矩阵基因的共表达模块，基于grnboost2，R中时GENIE3
  
###分析第二步：RcisTarget---运行完得到sce.regulons.csv文件
pyscenic ctx --num_workers 10 \
  --output sce.regulons.csv \
  --expression_mtx_fname sce.loom \
  --all_modules \
  --mask_dropouts \
  --mode "dask_multiprocessing" \
  --min_genes 10 \
  --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
  sce.adj.csv \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
  #这一步的目的
  #进行TF-motif富集分析，识别直接靶标
  #得到转录因子(TF)与其对应的直接作用的靶点,称为regulon(每一个regulon是1个TF和其调控的靶基因)
  
#####分析第三步：AUCell---运行完得到sce_SCENIC.loom文件，即分析结果
pyscenic aucell --num_workers 3 \
  --output sce_SCENIC.loom \
  sce.loom \
  sce.regulons.csv
  #这一步的目的
  #使用AUCell对每个细胞的每个regulon活性进行评分。


