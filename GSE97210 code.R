#一、下载数据
rm(list = ls()) 
options(stringsAsFactors = F)
setwd("E:/Desktop/新")
library(GEOmirror)
library(GEOquery)
library(limma)
library(extrafont)
library(ggpubr)
gset <- getGEO("GSE97210", destdir=".", AnnotGPL = F, getGPL = F)
gset   
dat=exprs(gset[[1]])
pd=pData(gset[[1]])
dat=dat[apply(dat,1,sd)>0,]
dat[dat<0]=1
dat=normalizeBetweenArrays(dat)

#二、分组&注释
table(pd$source_name_ch1)
##分组
group_list=ifelse(grepl('Normal arterial intima',pd$title),'Normal','AS')
table(group_list)
#设置参考水平，对照在前，处理在后
group_list = factor(group_list,
                    levels = c('Normal','AS'))
table(group_list)

# GPL16956 注释
library(openxlsx)
ann=read.xlsx("2.xlsx",sheet=1,colNames = TRUE)##seqmap比对得到2.xlsx
colnames(ann)
dat=as.data.frame(dat)
dat$ID=rownames(dat)
id_symbol_expr<-na.omit(merge(x=dat,y=ann[c('ID','lncRNA_ID')],by='ID',all.x=T))
#这里处理一个探针对应不同gene名的情况，随便取一个就行
symbol<-lapply(id_symbol_expr$lncRNA_ID,FUN = function(x){strsplit(x,'///')[[1]][1]})
id_symbol_expr$lncRNA_ID<-as.character(id_symbol_expr$lncRNA_ID)
#去除掉NA值，就是说有些探针对应不到已知的基因上，至少在GPL文件中没有对应关系
ids=id_symbol_expr[id_symbol_expr$lncRNA_ID != 'NA',]
ids=ids[ids$ID %in%  rownames(dat),]
dat=dat[ids$ID,] 
#处理一个gene名对应多个探针情况，中位数排序取最大
#去掉dat最后一列的ID
b=dat[,-7]
ids$median=apply(b,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$lncRNA_ID,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$lncRNA_ID),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat=dat[ids$ID,] #新的ids取出ID这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$lncRNA_ID#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
dat=dat[,-7]
dim(dat) 
dat[1:4,1:4] 
table(group_list) #table函数，查看group_list中的分组个数

#三、非匹配分析
library(limma)
design=model.matrix(~group_list)
fit=lmFit(dat,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t=2
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

#Figure S1A PCA
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,group_list)
library("FactoMineR")
library("factoextra")
dat.pca<-PCA(dat[,-ncol(dat)],graph=FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(dat.pca,
             geom.ind="point",
             col.ind=dat$group_list,
             palette=c("#00AFBB","#E7B800"),
             addEllipses=TRUE,
             legend.title="Groups"
)

#Figure S1B volcano plot
nrDEG=deg
attach(nrDEG)
library(ggpubr)
library(ggthemes)
library(ggplot2)
df=nrDEG
df$g=ifelse(df$P.Value>0.05,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
            ifelse( df$logFC >2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                    ifelse( df$logFC < -2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
)
table(df$g)

ggplot(df, aes(x = logFC, y = -log10(P.Value), colour=g)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 坐标轴
  labs(x="log2FoldChange",
       y="-log10(P.Value)")+
  theme_bw()+
  ggtitle("GSE97210")+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

#Figure 1A heatmap
#ctrl+C复制GSE97210 TOP10下调热图.xlsx中数据
dat=read.table("clipboard",header = T)
library(readxl)
library(openxlsx)
library(RColorBrewer) 
library(pheatmap)
annotation_col = data.frame(Sample=factor(c(rep("Normal",3),rep("AS",3))))#创建分组列
row.names(annotation_col) = colnames(dat) #这一行必须有，否则会报错：Error in check.length("fill") :  'gpar' element 'fill' must not be length 0

ann_colors = list(Sample = c('Normal'="blue", 'AS'="red")) #定义分组颜色
pheatmap(dat,
         scale = "row", # 按行归一化，查看因子在不同样本中的分布情况
         cluster_cols = FALSE, clustering_distance_rows = "correlation", #取消列聚类，表示行聚类使用皮尔森相关系数聚类
         treeheight_row = 30, # 设置行聚类树高
         cutree_rows =2, #根据样品列聚类情况将热图的行方向隔开为3份
         main="GSE97210", # 设置图形标题
         show_colnames = T, # 设置行列标签的显示
         show_rownames = T,
         cellwidth = 28,cellheight = 20,
         border="white", # 设置边框为白色
         legend = T, # FALSE去除图例; T显示图例
         fontsize_row = 10, # 分别设置行列标签字体大小
         fontsize_col = 10,
         angle_col = 90, # 设置标签显示角度
         annotation_col = annotation_col, #显示样品列的分组信息及图例
         annotation_colors = ann_colors, #使用annotation_colors参数设定样品列分组的颜色
)
