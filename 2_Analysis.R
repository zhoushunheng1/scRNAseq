################https://portals.broadinstitute.org/harmony/SeuratV3.html
rm(list=ls())
library(dplyr)
library(cowplot)
library(Seurat)
library(ggplot2)
resolution<-2  ## set the resolution
####load the expression matrix of the doubletfinder################ 
Stroke.sparse<-read.csv("./Stroke/Stroke.csv")  ## read CSV
Sham.sparse<-read.csv("./Sham/Sham.csv")

common_genes<-intersect(rownames(Stroke.sparse),rownames(Sham.sparse))
Stroke_specific_genes<-setdiff(rownames(Stroke.sparse),rownames(Sham.sparse))
Sham_specific_genes<-setdiff(rownames(Sham.sparse),rownames(Stroke.sparse))
####merge the two expression matrix########################
Stroke.sparse_common<-Stroke.sparse[common_genes,]
Stroke.sparse_specific<-Stroke.sparse[Stroke_specific_genes,]
Stroke.sparse_zero<-matrix(0,ncol=ncol(Stroke.sparse),nrow=length(Sham_specific_genes))
colnames(Stroke.sparse_zero)<-colnames(Stroke.sparse_common)
rownames(Stroke.sparse_zero)<-Sham_specific_genes
Stroke.sparse<-rbind(Stroke.sparse_common,Stroke.sparse_specific,Stroke.sparse_zero)
###############################
Sham.sparse_common<-Sham.sparse[common_genes,]
Sham.sparse_specific<-Sham.sparse[Sham_specific_genes,]
Sham.sparse_zero<-matrix(0,ncol=ncol(Sham.sparse),nrow=length(Stroke_specific_genes))
colnames(Sham.sparse_zero)<-colnames(Sham.sparse_common)
rownames(Sham.sparse_zero)<-Stroke_specific_genes
Sham.sparse<-rbind(Sham.sparse_common,Sham.sparse_zero,Sham.sparse_specific)
##############################
colnames(Sham.sparse)<- sub(pattern = "1", replacement = "2", x = colnames(Sham.sparse))##unique colname
####create the seurat object####################################
Stroke <- CreateSeuratObject(counts = cbind(Stroke.sparse, Sham.sparse), project = "Stroke", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = Stroke@var.genes, npcs = 50, verbose = FALSE)
	
Stroke@meta.data$stim <- c(rep("Stroke", ncol(Stroke.sparse)), rep("Sham", ncol(Sham.sparse)))#

#######################dimension reduction##############
Stroke <- Stroke %>% 
	RunTSNE( dims = 1:20) %>% 
	FindNeighbors( dims = 1:20) %>% 
	FindClusters(resolution = resolution) %>%   identity()
###################################################################################
################################################################
##########################################################
# find all markers of clusters###
Stroke.markers <- FindAllMarkers(Stroke,  min.pct = 0.25, logfc.threshold = 0.25)

top20markers <- Stroke.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  
write.table(top20markers,"6_top20_markers.resolution2.txt",sep="\t",quote=F)
save(Stroke,file="Stroke.resolution2.rds")
#########################remove cluster and reclustering ############################

Stroke$previous<-Stroke$RNA_snn_res.2
idents= c(0,2,3,19,22)  ##remove  cells##
cells.use = rownames(Stroke@meta.data)[!Idents(Stroke) %in% idents]
Stroke<- subset(x = Stroke, cells =cells.use , invert = FALSE)  ##get count mat##

resolution<-2
Stroke <- Stroke %>% 
    RunTSNE(dims = 1:20) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = resolution) %>%     identity()
save(Stroke,file="Stroke.recluster.rds")
Stroke.markers <- FindAllMarkers(Stroke, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20_markers<-Stroke.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  ##top 1 of each cluster 

new.cluster.ids <-c("Microglia","BAMs","Microglia","Microglia","Microglia","Microglia","Astrocytes","Astrocytes","Monocytes","Astrocytes","Oligodendrocytes","Microglia","DCs","Microglia","Granulocytes","Astrocytes","OPCs","NK&T cells","OPCs","Endothelial cells","Astrocytes","Endothelial cells","Fibroblast-like cells","Astrocytes","Oligodendrocytes","Endothelial cells","DCs","Mural cells","Mural cells","OPCs","Mural cells","Mural cells","Oligodendrocytes","Endothelial cells","Microglia","NK&T cells")
names(new.cluster.ids) <- levels(Stroke)
Stroke <- RenameIdents(Stroke, new.cluster.ids) #
save(Stroke,file="Stroke.rename.rds") ## 
