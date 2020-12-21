library(dplyr)
library(DESeq2)
library(patchwork)
library(ggplot2)
library(metap)
library(ape)
library(AnnotationHub)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Mmu.eg.db)
library(org.Hs.eg.db)
library(ggfortify)
library(stringr)
library(biomaRt)
library(RColorBrewer)
library(ClusterR)
library(Seurat)
library(monocle3)
library(pheatmap)
library(UpSetR)
rm(list = ls())

human52d<-read.table("/data2/bulk/T1human 52d/Quant.genes.results",header=T, row.names=1, com='', quote='',
                     check.names=F, sep="\t")
human63d<-read.table("/data2/bulk/T2human 63d/Quant.genes.results",header=T, row.names=1, com='', quote='',
                     check.names=F, sep="\t")
human77d<-read.table("/data2/bulk/T3human 77d/Quant.genes.results",header=T, row.names=1, com='', quote='',
                     check.names=F, sep="\t")
moca40d<-read.table("/data2/bulk/T1monkey 40d/Quant.genes.results",header=T, row.names=1, com='', quote='',
                     check.names=F, sep="\t")
moue15pc<-read.csv("/data2/bulk/all_compare.csv")
head(moue15pc)
a<-as.data.frame(matrix(NA,dim(human52d)[1],2))
#a[,1]<-"52d"
a[,1]<-rownames(a)<-rownames(human52d)
b<-as.data.frame(matrix(NA,dim(human52d)[1],2))
d<-as.data.frame(matrix(NA,dim(moca40d)[1],2))
c<-b
####a[,2]<-human52d$TPM
a[,2]<-human52d$expected_count
#b[,1]<-"63d"
b[,2]<-human63d$expected_count
b[,1]<-rownames(b)<-rownames(human63d)
c[,2]<-human77d
#c[,1]<-"77d"
v<-d[,1]<-rownames(moca40d)
c[,2]<-human77d$expected_count
d[,2]<-moca40d$expected_count
c[,1]<-s<-rownames(c)<-rownames(human77d)
n<-cbind(a[2],b[2],c[2])
#m<-rbind(a,b,c)
#colnames(m)<-c("batch","tpm")
#e<-colnames(n)<-c("52dTPM","63dTPM","77dTPM")
########################################
###human
#head(listDatasets(ensembl)
ensembl = useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
hsa_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values =s, mart = mart)
rownames(hsa_symbols)<-hsa_symbols$ensembl_gene_id
#####monkey
ensembl = useMart("ensembl")
mart<-useDataset("mfascicularis_gene_ensembl",useMart("ensembl"))
moca_symbol<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                   filters = 'ensembl_gene_id', values =v, mart = mart)
?getBM
#####mouse
##ensembl = useMart("ensembl")
###mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
###moca_symbol<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values =v, mart = mart)

##############counts to tpm
#countToTpm <- function(counts, effLen)
#{
#  rate <- log(counts) - log(effLen)
#  denom <- log(sum(exp(rate)))
#  exp(rate - denom + log(1e6))
#}

#countToFpkm <- function(counts, effLen)
#{
#  N <- sum(counts)
# exp( log(counts) + log(1e9) - log(effLen) - log(N) )
#}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  
}

#countToEffCounts <- function(counts, len, effLen)
#{
#  counts * (len / effLen)
#}
#s1<-moue15pc[,c(17:21)]
s1<-moue15pc[,c(2:6)]
#wt1tpm <- apply(s1,2,fpkmToTpm)
#head(wt1tpm)
#dim(wt1tpm)
#s1$en<-rownames(s1)<-moue15pc$gene_id
#colnames(s1)<-gsub("fpkm","mou15dtpm",colnames(s1))
###########
result<-c()
for (i in moue15pc$gene_name){
  i<-toupper(i)
  result<-c(result,i)
}
moue15pc$re<-result
s1$external_gene_name<-result
#############
############combine
m1<-unique(intersect(moca_symbol$external_gene_name,hsa_symbols$external_gene_name))
m1<-intersect(result,m1)
b1<-subset(hsa_symbols,external_gene_name%in%m1)
b1<-b1[!duplicated(b1$external_gene_name),]

n<-n[rownames(b1),]
n$en<-rownames(n)
n$external_gene_name<-rownames(n)<-b1$external_gene_name
rownames(b1)<-b1$external_gene_name

colnames(n)<-c("hsa52d","hsa63d","hsa77d","en","external_gene_name")
head(n)

#############################3
c1<-subset(moca_symbol,external_gene_name%in%m1)
c1<-c1[!duplicated(c1$external_gene_name),]
rownames(c1)<-c1$ensembl_gene_id
moca_n<-moca40d[rownames(c1),]
rownames(moca_n)<-c1$external_gene_name
moca_n$external_gene_name<-c1$external_gene_name
moca_n<-moca_n[,c(5,7)]
colnames(moca_n)<-c("moca40d","external_gene_name")
###########################################

s1<-subset(s1,external_gene_name%in%m1)
s1<-s1[!duplicated(s1$external_gene_name),]

rownames(s1)<-s1$external_gene_name
##############combine
combine<-merge(n,moca_n,by="external_gene_name")
combine<-merge(combine,s1,by="external_gene_name")
write.csv(combine,file="combine1.csv",quote = F)
#s<-stringr::str_to_title(s)
###eg <- bitr(s, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db"); head(eg)
###eh <- bitr(r, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mmu.eg.db"); head(eh)
###unique(eg$SYMBOL)
#############################
combine<-read.csv("combine1.csv")
################################二代分析
####假设有批次，去批次之后硬算
#人鼠，人猴，猴鼠，就硬算
#人鼠
library(DESeq2)
library(limma)
d<-combine[,c(3:5,7:10)]

head(d)
d$external_gene_name<-combine$external_gene_name
d<-d[!duplicated(d$external_gene_name),]
rownames(d)<-d$external_gene_name
d<-d[,-10]
# d<-as.matrix(d[,-10])

hsa_mou_data<-d[,c(1:3,5:7)]
# hsa_mou_data<-as.matrix(hsa_mou_data)


hsa_mou_data<-ceiling(hsa_mou_data)
head(hsa_mou_data)
sample1<-as.data.frame(matrix(NA,6,1))
rownames(sample1)<-colnames(hsa_mou_data)

colnames(sample1)<-"condition"
sample1$condition<-c()

sample1$condition<-c(c(sample1$conditions,rep("hsa",3)),rep("mou",3))
#sample1$batch<-as.factor(c("1","2","3",rep("4",5)))
#sample1$individual<-factor(gsub("^W.*","mou",rownames(sample1)))
#sample1<-sample1[c(4,1,5,2,6,3,7)]
dds <- DESeqDataSetFromMatrix(countData = hsa_mou_data,
                              colData = sample1,
                              design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
#assay(dds) <- limma::removeBatchEffect(assay(dds), vsd$batch)
?removeBatchEffect
dds <- DESeq(dds)
###################norm
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
# 标准化后的数据输出
write.table(normalized_counts, file="ehbio_trans.Count_matrix.xls.DESeq2.normalized.xls",
            quote=F, sep="\t", row.names=T, col.names=T)

# log转换后的结果
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]


write.table(rlogMat, file="ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls",
            quote=F, sep="\t", row.names=T, col.names=T)
# 生成颜色
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# 计算相关性pearson correlation
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))

# 层级聚类
hc <- hcluster(t(rlogMat), method="pearson")

# 热图绘制
## 在命令行下运行时，需要去除下面开头的#号，把输出的图保存到文件中
## 输出到文件时，dev.off()命令是关闭输出，从而完成图片的写入。如果不做这一步，则图片则不能写入，从而不能打开
## 如果在Rstudio或其它可视化界面时，可以直接把图输出到屏幕
#pdf("ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)
pheatmap(pearson_cor,scale = "none",show_rownames = T,cluster_cols = F)
r<-as.list(rlogMat)
M<-upset(r,nsets = 6)
?upset
res <- results(dds)
summary(res)
res=res[order(res$padj),]
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
####################
boxplot(as.data.frame(normalized_counts),main="Original")
###################################################


expmatrix <- as.matrix(rlogMat)

celldata <- as.data.frame(colnames(rlogMat))
rownames(celldata) <- celldata[, 1]

genedata <- as.data.frame(rownames(rlogMat))
rownames(genedata) <- genedata[, 1]
cds <- CreateSeuratObject(expression_data = expmatrix,
                         cell_metadata = celldata,
                         gene_metadata = genedata)
cds<-CreateSeuratObject(counts = expmatrix,project = 'cds',min.cells = 3,min.features = 200)
cds<-new_cell_data_set(expmatrix)
str(expmatrix)
cds<-new_cell_data_set()