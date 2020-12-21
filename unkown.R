library(dplyr)
library(Seurat)
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
library(dplyr)
library(patchwork)
library(ggplot2)
library(AnnotationHub)
library(VennDiagram)
library(UpSetR)
library(pheatmap)
library(clusterProfiler)

human52d<-read.table("/data2/bulk/T1human 52d/Quant.genes.results",header=T, com='', quote='',
                     check.names=F, sep="\t")
human63d<-read.table("/data2/bulk/T2human 63d/Quant.genes.results",header=T, com='', quote='',
                     check.names=F, sep="\t")
human77d<-read.table("/data2/bulk/T3human 77d/Quant.genes.results",header=T, com='', quote='',
                     check.names=F, sep="\t")
moca40d<-read.table("/data2/bulk/T1monkey 40d/Quant.genes.results",header=T, com='', quote='',
                    check.names=F, sep="\t")
moue15pc<-read.csv("/data2/bulk/all_compare.csv",header=T,row.names = 1)
hsa52d<-human52d[,c(1,6)]
hsa63d<-human63d[,c(1,6)]
hsa77d<-human77d[,c(1,6)]
moc40d<-moca40d[,c(1,6)]
rownames(hsa52d)<-hsa52d$gene_id
rownames(hsa63d)<-hsa52d$gene_id
rownames(hsa77d)<-hsa52d$gene_id
rownames(moc40d)<-moc40d$gene_id
mou15pc<-moue15pc[,c(16:20,46)]


#m<-rbind(a,b,c)
#colnames(m)<-c("batch","tpm")
#e<-colnames(n)<-c("52dTPM","63dTPM","77dTPM")
########################################
###human
#head(listDatasets(ensembl)
ensembl = useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
hsa_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values =hsa52d$gene_id, mart = mart)
rownames(hsa_symbols)<-hsa_symbols$ensembl_gene_id
#####monkey
ensembl = useMart("ensembl")
mart<-useDataset("mfascicularis_gene_ensembl",useMart("ensembl"))
moca_symbol<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                   filters = 'ensembl_gene_id', values =moc40d$gene_id, mart = mart)
#####mouse
##ensembl = useMart("ensembl")
###mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
###moca_symbol<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
#

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
wttpm <- apply(mou15pc[,1:5],2,fpkmToTpm)
mou15pc[,1:5]<-wttpm
#head(wt1tpm)
#dim(wt1tpm)
#s1$en<-rownames(s1)<-moue15pc$gene_id
colnames(mou15pc)<-gsub("fpkm","TPM",colnames(mou15pc))
###########
result<-c()
for (i in mou15pc$gene_name){
  i<-toupper(i)
  result<-c(result,i)
}
mou15pc$gene_id<-result
#############
hsa_symbols<-subset(hsa_symbols,ensembl_gene_id%in%hsa52d$gene_id)
hsa52d<-hsa52d[hsa_symbols$ensembl_gene_id,]
hsa63d<-hsa52d[hsa_symbols$ensembl_gene_id,]
hsa77d<-hsa77d[hsa_symbols$ensembl_gene_id,]
a<-hsa_symbols$external_gene_name
hsa52d$external_gene_id<-a
hsa63d$external_gene_id<-a
hsa77d$external_gene_id<-a
hsa52d<-hsa52d[!duplicated(hsa52d$external_gene_name),]
hsa63d<-hsa63d[!duplicated(hsa63d$external_gene_name),]
hsa77d<-hsa77d[!duplicated(hsa77d$external_gene_name),]
rownames(hsa63d)<-rownames(hsa77d)<-rownames(hsa52d)<-hsa52d$external_gene_id
moca_symbol<-subset(moca_symbol,ensembl_gene_id%in%moc40d$gene_id)
moc40d<-moc40d[moca_symbol$ensembl_gene_id,]
moc40d$external_gene_id<-moca_symbol$external_gene_name
moc40d<-moc40d[!duplicated(moca40d$external_gene_name),]
rownames(moc40d)<-moc40d$external_gene_id
############combine

# # m1<-unique(intersect(moca_symbol$external_gene_name,hsa_symbols$external_gene_name))
# # m1<-intersect(result,m1)
# # b1<-subset(hsa_symbols,external_gene_name%in%m1)
# # b1<-b1[!duplicated(b1$external_gene_name),]
# # 
# # n<-n[rownames(b1),]
# # n$en<-rownames(n)
# # n$external_gene_name<-rownames(n)<-b1$external_gene_name
# # rownames(b1)<-b1$external_gene_name
# # 
# # colnames(n)<-c("hsa52d","hsa63d","hsa77d","en","external_gene_name")
# # head(n)
# 
# #############################3
# c1<-subset(moca_symbol,external_gene_name%in%m1)
# c1<-c1[!duplicated(c1$external_gene_name),]
# rownames(c1)<-c1$ensembl_gene_id
# moca_n<-moca40d[rownames(c1),]
# rownames(moca_n)<-c1$external_gene_name
# moca_n$external_gene_name<-c1$external_gene_name
# moca_n<-moca_n[,c(5,7)]
# colnames(moca_n)<-c("moca40d","external_gene_name")
# ###########################################
# 
# s1<-subset(s1,external_gene_name%in%m1)
# s1<-s1[!duplicated(s1$external_gene_name),]
# 
# rownames(s1)<-s1$external_gene_name
# ##############combine
# combine<-merge(n,moca_n,by="external_gene_name")
# combine<-merge(combine,s1,by="external_gene_name")
# write.csv(combine,file="combine1.csv",quote = F)
