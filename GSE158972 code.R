#一、下载数据
rm(list = ls()) 
options(stringsAsFactors = F)
setwd("E:/Desktop/新")
library(GEOmirror)
library(GEOquery)
library(limma)
library(openxlsx)
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE158972", "file=GSE158972_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID
write.table(tbl,"12345.txt",sep='\t')####矩阵文件输出

dat=read.xlsx("GSE158972矩阵.xlsx",sheet=1,colNames = TRUE)
rownames(dat)=dat[,1]
dat=dat[,-1]
range(dat)
dat=log2(dat+1)
dat=normalizeBetweenArrays(dat)
#二、分组&注释
##分组
pd=read.csv('GSE158972.csv',header=T,check.names = F)
group_list=ifelse(grepl('control',pd$title),'control','si-CARMN')
table(group_list)
#设置参考水平，对照在前，处理在后
group_list = factor(group_list,
                    levels = c('control','si-CARMN'))
table(group_list)

# GPL20301 注释
library(openxlsx)
ann=read.xlsx("GSE158972注释文件.xlsx",sheet=1,colNames = TRUE)
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
table(group_list)

#三、非匹配分析
library(limma)
design=model.matrix(~group_list)
fit=lmFit(dat,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg,"GSE158972.logFC.前.txt",sep='\t')
## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t=1
deg$g=ifelse(deg$adj.P.Val>0.01,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

#Figure 2A volcano plot
nrDEG=deg
attach(nrDEG)
library(ggpubr)
library(ggthemes)
library(ggplot2)
df=nrDEG
df$g=ifelse(df$adj.P.Val>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
            ifelse( df$logFC >1,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                    ifelse( df$logFC < -1,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
)
table(df$g)

ggplot(df, aes(x = logFC, y = -log10(P.Value), colour=g)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 坐标轴
  labs(x="log2FoldChange",
       y="-log10(P.Value)")+
  theme_bw()+
  ggtitle("GSE158972")+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

#Figure 2B KEGG functional enrichment analysis
library(openxlsx)
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)
#自噬结合前富集分析
info <- read.xlsx( "GO10.xlsx", rowNames = F,colNames = T)
#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
dotplot(KEGG,label_format=70)

#Figure 5B-D GO and KEGG functional enrichment analyses
info <- read.xlsx( "GO8.xlsx", rowNames = F,colNames = T)

#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

GO<-enrichGO( gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T)

KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)
library(ggplot2)
#气泡图好用
dotplot(GO, split="ONTOLOGY",showCategory = 10,label_format=200)+ 
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(panel.grid = element_blank())+#修改主题
  theme(axis.title =element_text(size = 12, color = 'black'),
        axis.text.y =element_text(size = 12),
        legend.title=element_text(size=12))+
  scale_color_gradient(high="#FC8D62",low="#4b5cc4")#设置颜色

#KEGG弦图+弦表图
genedata<-data.frame(ID=info$gene_symbol,logFC=info$log2FoldChange)
KEGG<-read.xlsx("sc2.xlsx",sheet=1,colNames = TRUE)
KEGG$geneID <-str_replace_all(KEGG$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
names(KEGG)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
KEGG$Category = "BP"#分类信息

circ_BP<-GOplot::circle_dat(KEGG,genedata) #GOplot导入数据格式整理

chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框

GOChord(data = chord_BP,#弦图
        title = 'KEGG',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小
GOCircle(circ_BP) #弦表图
browseKEGG(KEGG,'hsa04140')