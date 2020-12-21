library(dplyr)
library(reticulate)
use_python("/home/user/anaconda3/envs/R_env/bin/python3",required = TRUE)
library(sceasy)
library(reticulate)
use_condaenv('/home/user/anaconda3/envs/R_env')
library(Seurat)
library(patchwork)
library(ggplot2)
library(metap)
library(UpSetR)
library(clusterProfiler)
library(AnnotationHub)
library(org.Hs.eg.db)
library(pheatmap)
library(rlist)
library(DoubletFinder)
library(SeuratWrappers)
library(harmony)
library(liger)
library(conos)
rm(list = ls())


ctxge141.data <- Read10X("~/Downloads/GSM4610591_E14_CtxGE_Emx1_N1/E14_CtxGE_Emx1_N1")
ctxge142.data <- Read10X(data.dir = "/home/user/Downloads/GSM4610592_E14_CtxGE_Emx1_N2/E14_CtxGE_Emx1_N2")
ctxge143.data <- Read10X(data.dir = "~/Downloads/GSM4610591_E14_CtxGE_Emx1_N1/E14_CtxGE_Emx1_N1")
hsamcfmmu <- read.csv(file = "/home/user/hsamcfmmu.csv", 
                      header = T, encoding = 'utf-8')

orth <- hsamcfmmu %>% filter(Crab.eating.macaque.gene.stable.ID != "" & 
                               Crab.eating.macaque.gene.name != "" & 
                               Mouse.gene.stable.ID != "" & 
                               Mouse.gene.name != "")

orthh <- orth[!duplicated(orth[, 1]), ]
orthhf <- orthh[!duplicated(orthh[, 3]), ]
orthhfu <- orthhf[!duplicated(orthhf[, 5]), ]
# 把三个物种中重复的编号找出，默认保留第一个，其余全部去掉
# （？？？）存在id一样，但是gene name 不一样的条目，也是神奇
# （？？？）神他妈有基因不同id，标注了orth，但用同一个基因名
orthhfuh <- orthhfu[!duplicated(orthhfu[, 2]), ]
orthhfuhf <- orthhfuh[!duplicated(orthhfuh[, 4]), ]
orthhfuhfu <- orthhfuhf[!duplicated(orthhfuhf[, 6]), ]
one2oneorth <- orthhfuhfu
dim(ctxge141.data)
# [1] 27999  2432
dim(ctxge142.data)
# [1] 27999  4204
dim(ctxge143.data)
# [1] 27999  5597
# colnames(ctxge141.data) <- paste('N1', seq(1, ncol(ctxge141.data)), sep = '_')
# colnames(ctxge142.data) <- paste('N2', seq(1, ncol(ctxge142.data)), sep = '_')
# colnames(ctxge143.data) <- paste('N3', seq(1, ncol(ctxge143.data)), sep = '_')
# # 重命名列名，既是避免barcode 重复，也为了后面方便添加identity
# ctxge14.data <- cbind(ctxge141.data, ctxge142.data, ctxge143.data)
# dim(ctxge14.data)
# [1] 27999 12233
x<-list(mmu141=ctxge141.data,mmu142=ctxge142.data,mmu143=ctxge143.data)
for (i in 1:length(x)) {
  current<- x[[i]]
  dim(current)
  rownames(one2oneorth) <- one2oneorth[, 6]
  a<-unique(intersect(rownames(current), one2oneorth[, 6]))
  current<-current[a,]
  current<-current[!duplicated(rownames(current)),]
  dim(current)
  one<-one2oneorth[a,][,c(2,6)]
  #one<-one2oneorth[!duplicated(subset(one2oneorth,Mouse.gene.name%in%a)),]
  rownames(current)<-one$Gene.name
  current <- CreateSeuratObject(counts = current, 
                                project = "mmu14_pfc14", 
                                min.cells = 12, 
                                min.features = 800)
  current[["percent.mt"]] <- PercentageFeatureSet(object = current, pattern = '^MT-')
  current<- subset(current, subset = percent.mt < 10 & 
                     percent.mt < quantile(current@meta.data$percent.mt, probs=seq(0, 1, 0.05))[20] & 
                     nCount_RNA < 10*nFeature_RNA)
  #current <- subset(current,!duplicated(current[unique(intersect(rownames(current), one2oneorth[, 6])), ]))
  x[[i]]<-current
}

#b <- CreateSeuratObject(counts = mmu14afqc, project = 'mmu14')
# 15410 features across 11051 samples within 1 assay 

# # 15410 features across 10471 samples within 1 assay 
# mmu14 <- ctxge14[intersect(rownames(ctxge14), one2oneorth[, 6]), ]
# #mmu14afqc <- ctxge14@assays$RNA@counts
# mmu14afqc <- mmu14@assays$RNA@counts
# rownames(one2oneorth) <- one2oneorth[, 6]
# rownames(mmu14afqc) <- one2oneorth[rownames(mmu14afqc), 2]
# b <- CreateSeuratObject(counts = mmu14afqc, project = 'mmu14')
# b@meta.data$percent.mt <- mmu14@meta.data$percent.mt
# mmu14 <- b
# # 11657 features across 10471 samples within 1 assay 
# mmu14@meta.data$authorident <- mmu14@meta.data$orig.ident
# mmu14@meta.data$orig.ident <- factor(rep('mmu14', ncol(mmu14)))
# levels(mmu14@active.ident)[levels(mmu14@active.ident) %in% c('N1', 'N2', 'N3')] <- 'mmu14'
# rm(ctxge141.data, ctxge142.data, ctxge143.data, ctxge14.data, ctxge14, mmu14afqc, b)


# import mmu pfc18
pfc18.data <- read.csv(file = "~/Downloads/GSM3494230_E18_cerebral_cortex.csv",
                       header = T,
                       row.names = 1,
                       encoding = 'utf-8')
head(pfc18.data)[1:10,1:10]
rownames(one2oneorth) <- one2oneorth[, 6]
a<-unique(intersect(rownames(pfc18.data), one2oneorth[, 6]))
pfc18.data<-pfc18.data[a,]
pfc18.data<-pfc18.data[!duplicated(rownames(pfc18.data)),]
dim(pfc18.data)
one<-one2oneorth[a,][,c(2,6)]
#one<-one2oneorth[!duplicated(subset(one2oneorth,Mouse.gene.name%in%a)),]
rownames(pfc18.data)<-one$Gene.name
# dim(pfc18.data)
# # [1] 16179  6753
# dgcge <- as(as.matrix(pfc18.data), 'dgCMatrix')
mmu <- CreateSeuratObject(counts = pfc18.data,
                          project = 'mmu14_pfc18',
                          min.cells = 6,
                          min.features = 800)
# 14889 features across 6709 samples within 1 assay
mmu[["percent.mt"]] <- PercentageFeatureSet(object = mmu, pattern = '^MT-')
mmu <- subset(mmu, subset = percent.mt < 7 & nCount_RNA < 5*nFeature_RNA)
# # 14889 features across 6698 samples within 1 assay 
# mmu <- mmu[intersect(rownames(mmu), one2oneorth[, 6]), ]
# mmuafqc <- mmu@assays$RNA@counts
# rownames(one2oneorth) <- one2oneorth[, 6]
# rownames(mmuafqc) <- one2oneorth[rownames(mmuafqc), 2]
# a <- CreateSeuratObject(counts = mmuafqc, project = 'mmu18')
# # 由于intersect之后，每个细胞中表达的基因会变少，因此不能添加min feature 语句
# # （？？？）但为什么min cell 语句也会导致基因数目减少，讲道理没有动细胞啊
# a@meta.data$percent.mt <- mmu@meta.data$percent.mt
# mmu18 <- a
# # 11420 features across 6698 samples within 1 assay 
# rm(a, mmu, dgcge, mmuafqc, pfc18.data)


### import mcf pfc45
pfc45.data <- Read10X(data.dir = "/data2/scRNA_data/mcf45PFC/filtered_feature_bc_matrix")
current<-pfc45.data
rownames(one2oneorth) <- one2oneorth[, 4]
a<-unique(intersect(rownames(current), one2oneorth[, 4]))
current<-current[a,]
current<-current[!duplicated(rownames(current)),]
dim(current)
one<-one2oneorth[a,][,c(2,4)]
rownames(current)<-one$Gene.name
# 
# dim(current)
# rownames(one2oneorth) <- one2oneorth[, 4]<-as.character(one2oneorth$Crab.eating.macaque.gene.name)
# current<-current[unique(intersect(rownames(current), one2oneorth[, 4])),]
# current<-current[!duplicated(rownames(current)),]
# one2oneorth<-one2oneorth[!duplicated(subset(one2oneorth,Crab.eating.macaque.gene.name%in%unique(intersect(rownames(current), one2oneorth[, 4])))),]
# rownames(current)<-one2oneorth$Gene.name
# colnames(pfc45.data) <- paste0('mcfpfc45_', seq(1, ncol(pfc45.data)))
# dim(pfc45.data)
# # [1] 29324 12704
pfc45 <- CreateSeuratObject(counts = current, 
                            project = 'mcf_pfc45', 
                            min.cells = 12, 
                            min.features = 800)
# mtgene <- c('ENSMFAG00000000001', 'ENSMFAG00000000002', 'ENSMFAG00000000003', 
#             'ENSMFAG00000000004', 'ENSMFAG00000000005', 'ENSMFAG00000000006', 
#             'ENSMFAG00000000007', 'ENSMFAG00000000008', 'ENSMFAG00000000009', 
#             'ENSMFAG00000000010', 'ENSMFAG00000000011', 'ENSMFAG00000000012', 
#             'ENSMFAG00000000013', 'ENSMFAG00000000014', 'ENSMFAG00000000015', 
#             'ENSMFAG00000000016', 'ENSMFAG00000000017', 'ENSMFAG00000000018', 
#             'ENSMFAG00000000019', 'ENSMFAG00000000020', 'ENSMFAG00000000021', 
#             'ENSMFAG00000000022', 'ENSMFAG00000000023', 'ENSMFAG00000000024', 
#             'ND1', 'ND2', 'ND3', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 
#             'ATP6', 'ATP8', 'ND4L', 'CYTB')
# remainmt <- intersect(mtgene, rownames(pfc45))
pfc45[["percent.mt"]] <- PercentageFeatureSet(object = pfc45, pattern = '^MT-')
pfc45 <- subset(pfc45, subset = percent.mt < 15 & 
                  percent.mt < quantile(pfc45@meta.data$percent.mt, probs=seq(0, 1, 0.05))[20] & 
                  nCount_RNA < 5*nFeature_RNA)
pfc45[["group"]]<-"pfc45"
# 13071 features across 9643 samples within 1 assay 
# pfc45 <- pfc45[intersect(rownames(pfc45), one2oneorth[, 4]), ]
# pfc45afqc <- pfc45@assays$RNA@counts
# rownames(one2oneorth) <- one2oneorth[, 4]
# rownames(pfc45afqc) <- one2oneorth[rownames(pfc45afqc), 2]
# a <- CreateSeuratObject(counts = pfc45afqc, project = 'mcf45')
# a@meta.data$percent.mt <- pfc45@meta.data$percent.mt
# mcf45 <- a
# # 10479 features across 9923 samples within 1 assay
# rm(a, pfc45, pfc45afqc,pfc45.data)


### import mcf pfc53
pfc53.data <- Read10X(data.dir = "/data2/scRNA_data/mcf53PFC/filtered_feature_bc_matrix")
# colnames(pfc53.data) <- paste0('mcfpfc53_', seq(1, ncol(pfc53.data)))
dim(pfc53.data)
# [1] 29324 10021
current<-pfc53.data
rownames(one2oneorth) <- one2oneorth[, 4]
a<-unique(intersect(rownames(current), one2oneorth[, 4]))
current<-current[a,]
current<-current[!duplicated(rownames(current)),]
dim(current)
one<-one2oneorth[a,][,c(2,4)]
rownames(current)<-one$Gene.name
pfc53 <- CreateSeuratObject(counts = current, 
                            project = 'mcf_pfc53', 
                            min.cells = 10, 
                            min.features = 800)
pfc53[["percent.mt"]] <- PercentageFeatureSet(object = pfc53, pattern = '^MT-')
pfc53 <- subset(pfc53, subset = percent.mt < 15 & 
                  percent.mt < quantile(pfc53@meta.data$percent.mt, probs=seq(0, 1, 0.05))[20] &  
                  nCount_RNA < 5*nFeature_RNA)
pfc53[["group"]]<-"pfc53"
# 13432 features across 8003 samples within 1 assay 
# pfc53 <- pfc53[intersect(rownames(pfc53), one2oneorth[, 4]), ]
# pfc53afqc <- pfc53@assays$RNA@counts
# rownames(one2oneorth) <- one2oneorth[, 4]
# rownames(pfc53afqc) <- one2oneorth[rownames(pfc53afqc), 2]
# a <- CreateSeuratObject(counts = pfc53afqc, project = 'mcf53')
# a@meta.data$percent.mt <- pfc53@meta.data$percent.mt
# mcf53 <- a
# # 11052 features across 8003 samples within 1 assay  
# rm(a, pfc53, pfc53afqc, pfc53.data)


### import hsa pfc77
pfc77.data <- Read10X(data.dir = "/data2/scRNA_data/hsa77pfc/filtered_feature_bc_matrix")
dim(pfc77.data)
# [1] 36601  7376
current<-pfc77.data
rownames(one2oneorth) <- one2oneorth[, 2]
a<-unique(intersect(rownames(current), one2oneorth[, 2]))
current<-current[a,]
current<-current[!duplicated(rownames(current)),]
dim(current)
one<-one2oneorth[a,][,c(2,4)]
rownames(current)<-one$Gene.name
pfc77 <- CreateSeuratObject(counts = pfc77.data, 
                            project = 'hsa_pfc77', 
                            min.cells = 7, 
                            min.features = 800)
pfc77[["percent.mt"]] <- PercentageFeatureSet(object = pfc77, pattern = '^MT-')
pfc77 <- subset(pfc77, subset = percent.mt < 15 & 
                  percent.mt < quantile(pfc77@meta.data$percent.mt, probs=seq(0, 1, 0.05))[20] &
                  nCount_RNA < 7*nFeature_RNA)
pfc77[["group"]]<-"pfc77"
# 19555 features across 6555 samples within 1 assay 
# pfc77 <- pfc77[intersect(rownames(pfc77), one2oneorth[, 2]), ]
# hsa <- pfc77
# # 12584 features across 6555 samples within 1 assay  
# rm(pfc77, pfc77.data)
############
x<-list.append(x,pattern = mmu,pfc45,pfc53,pfc77)

# 
# for (i in x) {
#   current<-i
#   current<-NormalizeData(current,normalization.method = "LogNormalize")
#   current<-FindVariableFeatures(current,selection.method = "vst", nfeatures = 5000)
#   all.genes <- rownames(current)
#   current<- ScaleData(current, features = all.genes)
#   current <- RunPCA(current, npcs = 50)
#   current <- FindNeighbors(current, 
#                            reduction = "pca", 
#                            dims = 1:19)
#   current <- FindClusters(current, 
#                           resolution = 0.7)
#   current <- RunUMAP(current, 
#                      reduction = "pca", 
#                      dims = 1:19)
#   i<-current
# }
# for (i in 1:length(x)) {
#   current<-x[[i]]
#   current<-NormalizeData(current,normalization.method = "LogNormalize")
#   current<-FindVariableFeatures(current,selection.method = "vst", nfeatures = 5000)
#  
# }
for (i in 1:length(x)) {
  current<-x[[i]]
  current<-NormalizeData(current,normalization.method = "LogNormalize")
  current<-FindVariableFeatures(current,selection.method = "vst", nfeatures = 5000)
  all.genes <- rownames(current)
  current<- ScaleData(current, features = all.genes)
  current <- RunPCA(current, npcs = 50)
  current <- FindNeighbors(current, 
                           reduction = "pca", 
                           dims = 1:19)
  current <- FindClusters(current, 
                          resolution = 0.7)
  current <- RunUMAP(current, 
                     reduction = "pca", 
                     dims = 1:19)
  x[[i]]<-current
}

#########################################33
for (i in 1:length(x)) {
  current<-x[[i]]
  sweep.mue14 <- paramSweep_v3(current, PCs = 1:30, sct = FALSE)
  sweep.stats_pfc77 <- summarizeSweep(sweep.mue14 , GT = FALSE)
  bcmvn_pfc77 <- find.pK(sweep.stats_pfc77)
  homotypic.prop <- modelHomotypic(current@meta.data$ClusteringResults)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(current@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  current <- doubletFinder_v3(current, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  current <- doubletFinder(current, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = grep('^pANN',names(current@meta.data[grep('^pANN',names(current@meta.data))])))
  x[[i]]<-current
}
change<-list()
doubletplot<-list()
for (i in 1:length(x)) {
  current<-x[[i]]
  current@meta.data[,"DF_hi.lo"] <- current@meta.data[[grep('^DF.classifications',names(current@meta.data))[1]]]
  current@meta.data$DF_hi.lo[which(current@meta.data$DF_hi.lo=="Doublet"& current@meta.data[[grep('^DF.classifications',names(current@meta.data))[2]]]=="Singlet")] <- "Doublet_lo"
  current@meta.data$DF_hi.lo[which(current@meta.data$DF_hi.lo=="Doublet")]<-"Doublet_hi"
  current@meta.data$DF_hi.lo[which(current@meta.data$DF_hi.lo=="Doublet_lo"|current@meta.data$DF_hi.lo=="Singlet")]<-"single"
  current1<-subset(current,subset=DF_hi.lo=="single")
  
  doubletplot[[i]]<-UMAPPlot(current, group.by="DF_hi.lo",label.color=c("black","gold","red"))
  x[[i]]<-current
  change[[i]]<-current1
}

m1<-do.call(gridExtra::grid.arrange,c(doubletplot))
ggsave("doublet.pdf",m1)
#########################
# 
# for (i in 1:length(x)) {
#   current<-change[[i]]
#   current<-NormalizeData(current,normalization.method = "LogNormalize")
#   current<-FindVariableFeatures(current,selection.method = "vst", nfeatures = 5000)
#   all.genes <- rownames(current)
#   current<- ScaleData(current, features = all.genes)
#   current <- RunPCA(current, npcs = 50)
#   current <- FindNeighbors(current, 
#                            reduction = "pca", 
#                            dims = 1:19)
#   current <- FindClusters(current, 
#                           resolution = 0.7)
#   current <- RunUMAP(current, 
#                      reduction = "pca", 
#                      dims = 1:19)
#   change[[i]]<-current
# }
# for (i in 1:length(change)) {
#   current<-change[[i]]
#   assign(paste(current@project.name, "loom", sep = "."),as.loom(current, filename = paste(names(current@meta.data), "loom", sep = "."), verbose = FALSE))
#   write.csv(current@meta.data,paste(current@project.name, "csv", sep = "."))
# }


#######################
# for (i in 1:length(x)) {
#   current<-change[[i]]
#   current<-current[one2oneorth,]
#   change[[i]]<-current
# }
comeb<-merge(change[[1]],y=c(change[[2]],change[[3]],change[[4]],change[[5]],change[[6]],change[[7]]),project="combine",merge.data=TRUE)
# combine<-merge(x[[1]],y=c(x[[2]],x[[3]],x[[4]],x[[5]],x[[6]]),project="combine",merge.data=TRUE)
combine<-comeb
#################harmony
harmonydata <- NormalizeData(combine) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
harmonydata <- RunHarmony(harmonydata, group.by.vars = "orig.ident")
harmonydata <- RunUMAP(harmonydata, reduction = "harmony", dims = 1:30)
harmonydata <- FindNeighbors(harmonydata, reduction = "harmony", dims = 1:30) %>% FindClusters()
p1<-DimPlot(harmonydata, group.by = "orig.ident")+labs(title = "harmony")
p1
######################liger
# Please update your `liger` version to 0.5.0 or above before following this tutorial.
ligerdata <- NormalizeData(combine)
ligerdata <- FindVariableFeatures(ligerdata)
ligerdata <- ScaleData(ligerdata , split.by = "orig.ident", do.center = FALSE)
ligerdata <- RunOptimizeALS(ligerdata , k = 20, lambda = 5, split.by = "orig.ident")
ligerdata <- RunQuantileNorm(ligerdata , split.by = "orig.ident")
ligerdata <- RunQuantileNorm(ligerdata, split.by = "orig.ident")
# You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
# `RunQuantileNorm` according to your needs
ligerdata <- FindNeighbors(ligerdata, reduction = "iNMF", dims = 1:20)
ligerdata <- FindClusters(ligerdata, resolution = 0.4)
# Dimensional reduction and plotting
ligerdata <- RunUMAP(ligerdata, dims = 1:ncol(ligerdata[["iNMF"]]), reduction = "iNMF")
p2<-DimPlot(ligerdata, group.by = "orig.ident")+labs(title = "liger")
p2

###############################fastmnn
fastmnndata <- NormalizeData(combine)
fastmnndata <- FindVariableFeatures(fastmnndata)
fastmnndata <- RunFastMNN(object.list = SplitObject(fastmnndata, split.by = "orig.ident"))
fastmnndata <- RunUMAP(fastmnndata, reduction = "mnn", dims = 1:19)
fastmnndata <- FindNeighbors(fastmnndata, reduction = "mnn", dims = 1:19)
fastmnndata <- FindClusters(fastmnndata)
p3<-DimPlot(fastmnndata, group.by ="orig.ident")+labs(title = "fastmnn")
##############################SCTransform
# combine.panel <- SplitObject(combine, split.by = "orig.ident")
# combine.panel<-combine.panel[names(combine.panel)]
# for (i in 1:length(combine.panel)) {
#   combine.panel[[i]] <- SCTransform(combine.panel[[i]], verbose = FALSE)
# }
# combine.features <- SelectIntegrationFeatures(object.list = combine.panel, nfeatures = 3000)
# combine.list <- PrepSCTIntegration(object.list = combine.panel, anchor.features = combine.features)
# combine.list <- PrepSCTIntegration(object.list = combine.panel, anchor.features = all.genes)
# combine.anchors <- FindIntegrationAnchors(object.list = combine.panel, normalization.method = "SCT", 
#                                        anchor.features = combine.features)
# combine.panel <- FindIntegrationAnchors(combine.panel, normalization.method = "SCT", 
#                                            anchor.features = combine.features, verbose = FALSE)
# combine.panel <- IntegrateData(anchorset = combine.anchors, normalization.method = "SCT", 
#                                      verbose = FALSE)
#####################################3
combine.panel <- SplitObject(combine, split.by = "orig.ident")
for (i in 1:length(combine.panel)) {
  combine.panel[[i]] <- NormalizeData(combine.panel[[i]], verbose = FALSE)
  combine.panel[[i]] <- FindVariableFeatures(combine.panel[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}
bothpfc.anchors <- FindIntegrationAnchors(object.list = combine.panel, 
                                          anchor.features = 2000, 
                                          scale = T, 
                                          reduction = 'cca', 
                                          dims = 1:30, 
                                          k.filter = 200)
bothpfc.combined <- IntegrateData(anchorset = bothpfc.anchors, dims = 1:10)
DefaultAssay(bothpfc.combined) <- "integrated"
bothpfc.combined  <- ScaleData(bothpfc.combined , verbose = FALSE)
bothpfc.combined  <- RunPCA(bothpfc.combined , npcs = 19, verbose = FALSE)
bothpfc.combined <- RunUMAP(bothpfc.combined, reduction = "pca", dims = 1:19)
p4 <- DimPlot(bothpfc.combined , reduction = "umap", group.by = "orig.ident")+labs(title = "seurat3.0")
p4
#############################simspec
library(simspec)
# Import packages and data and data preprocessing
simspecdata <- NormalizeData(combine)
simspecdata <- FindVariableFeatures(simspecdata, nfeatures = 3000)
simspecdata <- ScaleData(simspecdata)
simspecdata <- RunPCA(simspecdata, npcs = 20)
simspecdata <- RunUMAP(simspecdata, dims = 1:20)

p5<-UMAPPlot(simspecdata, group.by = "orig.ident") + labs(title = "merge")
p5
#CSS calculation with kernel probability transformation
simspecdata <- cluster_sim_spectrum(object =simspecdata, label_tag = "orig.ident",
                                cluster_resolution = 0.4,
                                corr_method = "pearson",
                                spectrum_type = "corr_kernel")
simspecdata <- RunUMAP(simspecdata,reduction = "css", dims = 1:ncol(Embeddings(simspecdata, "css")))
#The dimensionality of CSS increases when more samples/batches are to integrate. 
#The high dimensionality can slow down the later analysis. Optionally, 
#we can further reduce the dimensionality by applying an extra PCA
simspecdata <- run_PCA(simspecdata, reduction = "css", npcs = 10)
simspecdata <- RunUMAP(simspecdata, reduction = "css_pca", dims = 1:10)
simspecdata <- FindNeighbors(simspecdata, reduction = "css_pca", dims = 1:10)
simspecdata <- FindClusters(simspecdata, resolution = 0.3)
p6<-UMAPPlot(simspecdata, group.by = "orig.ident") +labs(title = "simspec")
p6
plot_grid(p1+p2+p3+p4+p5+p6)
?plot_grid
################################3

# mou<-subset(mue14,subset=DF_hi.lo=="Singlet")
# main.loom <- as.loom(bothpfc.anchors, filename = "/data1/main.loom", verbose = FALSE)
# write.csv(both@meta.data,'/data1/mian.csv')
# 
# main1.loom <- as.loom(bothpfc.anchors, filename = "/data1/main1.loom", verbose = FALSE)
# write.csv(bothpfc.anchor@meta.data,'/data1/mian1.csv')
# pfc45@meta.data$
#   ################################


loompy <- reticulate::import('loompy')
repl_python()
