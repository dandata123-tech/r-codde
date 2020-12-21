library(Biobase)
library(BisqueRNA)
library(dplyr)
library(Seurat)
library(patchwork)
library(metap)
options(stringsAsFactors = F)
library(DoubletFinder)
library(randomForest)
rm(list = ls())
moue15pc<-read.csv("/data2/bulk/all_compare.csv")
moue15pc<-moue15pc[!duplicated(moue15pc$gene_name),]
#rownames(moue15pc)<-moue15pc$gene_id
rownames(moue15pc)<-moue15pc$gene_name
moue15pc<-moue15pc[,2:6]
moue15pc<-as.matrix(moue15pc)
#moue15pc<-t(moue15pc)
mue14annotation.data<-read.table(file = '/home/user/Downloads/GSE123335_E14_combined_matrix_ClusterAnnotations.txt.gz',sep = '\t')
mue14.data<-read.table("/home/user/Downloads/GSE123335_E14_combined_matrix.txt.gz")
colnames(mue14annotation.data)<-mue14annotation.data[1,]

mue14annotation.data<-mue14annotation.data[-1,]
barcode.name<-mue14.data[1,-1]
feature.name<-mue14.data[-1,1]
matrix.name<-mue14.data[-1,-1]
rownames(matrix.name)<-as.vector(feature.name)
colnames(matrix.name)<-as.vector(barcode.name)
rownames(matrix.name)<-gsub("_","-",rownames(matrix.name))
mue14<-CreateSeuratObject(counts = matrix.name,project = 'mue14',min.cells = 3,min.features = 200)
head(mue14$nCount_RNA)
mue14[["cellid"]]<-gsub("-",".",rownames(as.data.frame(mue14$nCount_RNA)))
head(mue14$cellid)
class(mue14annotation.data$CellID)
m<-intersect(mue14$cellid,mue14annotation.data$CellID)
mue14<-subset(mue14,cellid%in%m)

mue14annotation.data<-subset(mue14annotation.data,CellID%in%m)
mue14[["cluster"]]<-mue14annotation.data$Cluster
mue14[["cellid"]]<-rownames(as.data.frame(mue14$nCount_RNA))
feature<-rownames(mue14@assays$RNA)
mue14 <- NormalizeData(mue14, normalization.method = "LogNormalize", scale.factor = 10000)
mue14<- FindVariableFeatures(mue14, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mue14)
mue14<- ScaleData(mue14, features = all.genes)
mue14 <- RunPCA(mue14, features = VariableFeatures(object = mue14))

#mue14<-mue14[intersect(feature,rownames(moue15pc)),]
#moue15pc<-moue15pc[intersect(feature,rownames(moue15pc)),]
#################################3
bulk.eset <- Biobase::ExpressionSet(assayData = moue15pc)
bulk.eset[["sampleNames"]]<-colnames(moue15pc)
bulk.eset[["SubjectName"]]<-colnames(moue15pc)
# bulk.sample<-Biobase:sampleNames(bulk.est)
###############################################################3
# sample.ids <- mue14$cellid

# sc.meta <- data.frame(labelDescription=c("SubjectName",
#                                          "cellType"),
#                       row.names=c("SubjectName",
#                                   "cellType"))
# sc.pdata <- new("AnnotatedDataFrame",
#                 data=sc.pheno,
#                 varMetadata=sc.meta)
# l<-mue14[['RNA']]@counts
# colnames(l)<-mue14$cellid
# rownames(l)<-mue14
# head(l)
#sc.eset <- Biobase::ExpressionSet(assayData=as.matrix(mue14[['RNA']]@counts))
# sc.eset@phenoData<-sc.pdata

sc.eset <- BisqueRNA::SeuratToExpressionSet(mue14, delimiter="-", position=2, version="v3")
sc.eset$cellType<-mue14$cluster<-gsub("\\[.*\\]","",mue14$cluster)
sc.eset$cellType<-mue14$cluster<-gsub("^SVZ.*","SVZ",mue14$cluster)
sc.eset$cellType<-mue14$cluster<-gsub("^RG.*","RG",mue14$cluster)
sc.eset$cellType<-mue14$cluster<-gsub("^I.*","Int",mue14$cluster)
#GenerateSCReference(sc.eset,cell.types="cellType")
#sc.sample<-base::levels(base::factor(sc.eset[['SubjectName']]))
sc.eset[["sampleNames"]]<-colnames(mue14)

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)
ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)
g<-round(as.data.frame(ref.based.estimates),2)
x<-res$transformed.bulk
sc<-specc(x,centers=12)
  ?kmeas
r <- cor(as.vector(ref.based.estimates), 
         as.vector(true.props[row.names(ref.based.estimates),colnames(ref.based.estimates)]))
reknitr::knit_print(sprintf("R: %f", r))
#res <- BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers, weighted=F)
saveRDS(res,"res.rds")
########################################################################
# library(xbioc)
# library(MuSiC)
# 
# # Estimate cell type proportions
# 
# # Estimate cell type proportions
# Est.prop.GSE50244 = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'cellType',
#                                samples = 'SubjectName', select.ct = NULL, verbose = F)
# saveRDS(sc.eset,"sc.eset")
# saveRDS(bulk.eset,"bulk.eset")

