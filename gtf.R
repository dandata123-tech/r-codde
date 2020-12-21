library(rtracklayer)
library(dplyr)
rm(list = ls())
gtf <- rtracklayer::import("/data3/cellranger/cellranger3.1/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
hsamcfmmu <- read.csv(file = "/home/user/hsamcfmmu.csv", 
                      header = T, encoding = 'utf-8')

dim(hsamcfmmu)
gtf<-as.data.frame(gtf)
vc  hsamcfmmu <- hsamcfmmu[!duplicated(hsamcfmmu[, 1]), ]
hsamcfmmu <- hsamcfmmu[!duplicated(hsamcfmmu[, 2]), ]
rownames(hsamcfmmu)<-hsamcfmmu$Gene.name
gtf<-gtf[!duplicated(gtf$gene_id),]
rownames(gtf)<-gtf$gene_name
a<-intersect(gtf$gene_name,hsamcfmmu$Gene.name)
hsamcfmmu<-subset(hsamcfmmu,Gene.name%in%a)
gtf<-gtf[,c(10,12,14)]
gtf<-subset(gtf,gene_name%in%a)
colnames(hsamcfmmu)[2]<-"gene_name"
hsamcfmmu<-merge(hsamcfmmu,gtf,by="gene_name")

for (i in 1:length(c)) {
  c$all[i]<-paste(c$gene_name[i],c$gene_biotype[i],sep = "|")
}

# #gtf=dplyr::select(gtf,c(gene_name,gene_id,gene_biotype))
# 
# for (i in 1:length(hsamcfmmu)) {
#   current<-hsamcfmmu[i,]$Gene.stable.ID
#   cu<-hsamcfmmu[i,]$Gene.name
#   a<-gtf[which(current == gtf$gene_id),]$gene_biotype
#   cur<-paste(cu,a,sep = "|")
#   hsamcfmmu[i,]$Gene.name<-cur
# }
# hsamcfmmu=dplyr::inner_join(a,hsamcfmmu,by ="gene_id")
# c[1:5,1:5]
# current<-current[a,]
# current<-current[!duplicated(rownames(current)),]
# dim(current)
# one<-one2oneorth[a,][,c(2,6)]
# c=dplyr::inner_join(b,LIHCdata1,by ="gene_id")c[1:5,1:5]
# #one<-one2oneorth[!duplicated(subset(one2oneorth,Mouse.gene.name%in%a)),]
# rownames(current)<-one$Gene.name