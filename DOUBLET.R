library(dplyr)
library(Seurat)
library(patchwork)
library(metap)
options(stringsAsFactors = F)
library(DoubletFinder)
library(randomForest)
rm(list = ls())

mue14annotation.data<-read.table(file = '/home/user/Downloads/GSE123335_E14_combined_matrix_ClusterAnnotations.txt.gz',sep = '\t')
mue14.data<-read.table("/home/user/Downloads/GSE123335_E14_combined_matrix.txt.gz")

barcode.name<-mue14.data[1,-1]
feature.name<-mue14.data[-1,1]
matrix.name<-mue14.data[-1,-1]
rownames(matrix.name)<-as.vector(feature.name)
colnames(matrix.name)<-as.vector(barcode.name)

rownames(matrix.name)<-gsub("_","-",rownames(matrix.name))
mue14<-CreateSeuratObject(counts = matrix.name,project = 'mue14',min.cells = 3,min.features = 200)
colnames(mue14$nCount_RNA)
mue14[["cellid"]]<-gsub("-",".",rownames(as.data.frame(mue14$nCount_RNA)))
dim(mue14)#21955  7376
#QC
mue14[["percent.mt"]]<-PercentageFeatureSet(object = mue14,pattern = "^MT-")
mue14[["percent.ribo"]]<-PercentageFeatureSet(object = mue14,pattern = "^RP[SL]")
plot1 <- FeatureScatter(mue14, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mue14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
VlnPlot(mue14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mue14 <- subset(mue14, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
mue14 <- NormalizeData(mue14, normalization.method = "LogNormalize", scale.factor = 10000)
mue14<- FindVariableFeatures(mue14, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mue14), 10)
plot3 <- VariableFeaturePlot(mue14)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot4
all.genes <- rownames(mue14)
mue14<- ScaleData(mue14, features = all.genes)
####################################
mue14 <- RunPCA(mue14, features = VariableFeatures(object = mue14))
print(mue14[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mue14, dims = 1:2, reduction = "pca")
DimPlot(mue14, reduction = "pca")
DimHeatmap(mue14, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mue14, dims = 1:15, cells = 500, balanced = TRUE)
mue14 <- JackStraw(mue14, num.replicate = 100)
mue14 <- ScoreJackStraw(mue14, dims = 1:20)
JackStrawPlot(mue14, dims = 1:15)
ElbowPlot(mue14,ndims = 40, reduction = "pca")
mue14 <- FindNeighbors(object = mue14, dims = 1:10)
mue14 <- FindClusters(object = mue14, resolution = 0.2)
head(Idents(mue14),5)
# non-linear reduction
# t-SNE
mue14 <- RunTSNE(object = mue14, dims.use = 1:10, do.fast = TRUE)
DimPlot(mue14, reduction = "tsne", label = T)
# UMAP
mue14 <- RunUMAP(mue14, dims = 1:10)
f1<-DimPlot(mue14, reduction = "umap")
FeaturePlot(mue14, features = c('PTPRC', 'PAX6', 'PDGFRA', 'NEUROD2', 'GAD1', 'AQP4'))
png(file="umap.png")
print(f1)

dev.off()
######
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.mue14 <- paramSweep_v3(mue14, PCs = 1:30, sct = FALSE)
sweep.stats_pfc77 <- summarizeSweep(sweep.mue14 , GT = FALSE)
bcmvn_pfc77 <- find.pK(sweep.stats_pfc77)
## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
##sweep.mue14 <- paramSweep_v3(mue14, PCs = 1:30, sct = FALSE)
##gt.calls <- mue14@meta.data[rownames(sweep.mue14[[1]]), "GT"]  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
##sweep.stats_pfc77 <- summarizeSweep(sweep.mue14, GT = TRUE, GT.calls = gt.calls)
##bcmvn_pfc77 <- find.pK(sweep.stats_pfc77)
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(mue14@meta.data$ClusteringResults)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(mue14@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
mue14 <- doubletFinder_v3(mue14, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mue14 <- doubletFinder_v3(mue14, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
## Plot results --------------------------------------------------------------------------------------------------------------
mue14@meta.data[,"DF_hi.lo"] <- mue14@meta.data$DF.classifications_0.25_0.09_643
mue14@meta.data$DF_hi.lo[which(mue14@meta.data$DF_hi.lo == "Doublet" & mue14@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mue14, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mou<-subset(mue14,subset=DF_hi.lo=="Singlet")

