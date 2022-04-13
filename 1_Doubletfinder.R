################## remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') ##
# install doubletfinder
rm(list=ls())
library(dplyr)
library(cowplot)
library(Seurat)
library(harmony)
library(DoubletFinder)
library(ggplot2) #  visualization
doublet_rate<-0.016
##################### 
Stroke.data <- Read10X(data.dir = "./Stroke/")
Stroke.data <- CreateSeuratObject(counts = Stroke.data, project = "Stroke", min.cells = 3)  ##18377  6796
load_number<-dim(Stroke.data)
Stroke.data[["percent.mt"]] <- PercentageFeatureSet(Stroke.data, pattern = "^mt-")
Stroke.data <- subset(Stroke.data, subset = nFeature_RNA > 200  & percent.mt < 20)  ##18377  6559
MT_number<-dim(Stroke.data)
Stroke.data <- NormalizeData(Stroke.data)
Stroke.data <- FindVariableFeatures(Stroke.data, selection.method = "vst", nfeatures = 2000)
Stroke.data <- ScaleData(Stroke.data)
Stroke.data <- RunPCA(Stroke.data)
Stroke.data <- RunUMAP(Stroke.data, dims = 1:10)
Stroke.data <- RunTSNE(Stroke.data, dims = 1:10)
#####pre-process####################
sweep.res.list <- paramSweep_v3(Stroke.data, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
######## Homotypic Doublet Proportion Estimate    ##ref:  A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart
nExp_poi <- round(doublet_rate*ncol(Stroke.data))        ## Assuming 1.6% doublet 
Stroke.data <- doubletFinder_v3(Stroke.data, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# table(Stroke.data$DF.classifications_0.25_0.005_328)
pdf("1_doubletFinder_Stroke.tsne.pdf")
DimPlot(Stroke.data,pt.size = 1,label=TRUE, label.size = 5,reduction = "tsne",group.by = "DF.classifications_0.25_0.3_105" )+theme(aspect.ratio = 1)
dev.off()
Stroke.data<- subset(Stroke.data, DF.classifications_0.25_0.005_105=="Singlet" )  ##18377  6231
DoubletFinder_number<-dim(Stroke.data)

##############
write.table(as.matrix(GetAssayData(object = Stroke.data, slot = "counts")), 
            'Stroke.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)
######################################################################################
######################################################################################
Sham.data <- Read10X(data.dir = "./Sham/")
Sham.data <- CreateSeuratObject(counts = Sham.data, project = "Sham", min.cells = 3) 
load_number<-dim(Sham.data)
Sham.data[["percent.mt"]] <- PercentageFeatureSet(Sham.data, pattern = "^mt-")
Sham.data <- subset(Sham.data, subset = nFeature_RNA > 200  & percent.mt < 20)  

MT_number<-dim(Sham.data)
Sham.data <- NormalizeData(Sham.data)
Sham.data <- FindVariableFeatures(Sham.data, selection.method = "vst", nfeatures = 2000)
Sham.data <- ScaleData(Sham.data)
Sham.data <- RunPCA(Sham.data)
Sham.data <- RunUMAP(Sham.data, dims = 1:10)
Sham.data <- RunTSNE(Sham.data, dims = 1:10)
#####pre-process####################
sweep.res.list <- paramSweep_v3(Sham.data, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
######## Homotypic Doublet Proportion Estimate    ##ref:  A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart
nExp_poi <- round(doublet_rate*ncol(Sham.data))        ## Assuming 1.6% doublet 
Sham.data <- doubletFinder_v3(Sham.data, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# table(Sham.data$DF.classifications_0.25_0.005_328)
pdf("1_doubletFinder_Sham.tsne.pdf")
DimPlot(Sham.data,pt.size = 1,label=TRUE, label.size = 5,reduction = "tsne",group.by = "DF.classifications_0.25_0.07_137" )+theme(aspect.ratio = 1)
dev.off()
Sham.data<- subset(Sham.data, DF.classifications_0.25_0.15_137=="Singlet" )  ##18377  6231
DoubletFinder_number<-dim(Sham.data)

##############
write.table(as.matrix(GetAssayData(object = Sham.data, slot = "counts")), 
            'Sham.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)
############
