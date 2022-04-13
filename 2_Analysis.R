################https://portals.broadinstitute.org/harmony/SeuratV3.html
rm(list=ls())
library(dplyr)
library(cowplot)
library(Seurat)
library(ggplot2)
library(harmony)
resolution<-2  ## set the resolution
####load the expression matrix of the doubletfinder################ 
Stroke.sparse<-read.csv("/mnt/zhoushunheng/stroke/2seurat/Result3_doubletfinder/Stroke.csv")  ## read CSV
Sham.sparse<-read.csv("/mnt/zhoushunheng/stroke/2seurat/Result3_doubletfinder/Sham.csv")
file<-"/mnt/zhoushunheng/stroke/2seurat/result5_doubletfinder_harmony_refinedClusterR3/no_harmony/"
dir.create(file)
setwd(file)
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

#######################umap和tsne降维##############
Stroke <- Stroke %>% 
    RunUMAP( dims = 1:20) %>% 
    RunTSNE( dims = 1:20) %>% 
    FindNeighbors( dims = 1:20) %>% 
    FindClusters(resolution = resolution) %>%   identity()

options(repr.plot.height = 4, repr.plot.width = 10)
pdf("2tsne_split_by_group.pdf") # pre harmony  correction
DimPlot(Stroke, reduction = "tsne", group.by = "stim", pt.size = .1)
dev.off()
pdf("2umap_split_by_cluster.pdf") # pre harmony  correction
DimPlot(Stroke, reduction = "umap", label = TRUE, pt.size = .1)
dev.off()
pdf("2tSNE_split_by_cluster.pdf") # pre harmony  correction
DimPlot(Stroke, reduction = "tsne", label = TRUE, pt.size = .1)
dev.off()
###################################################################################
Stroke <- BuildClusterTree(object = Stroke)
pdf("3_cluster_tree.pdf")
PlotClusterTree(object = Stroke)
dev.off()
cells<-table(Idents(Stroke))   #num of cells per cluster
write.table(cells,"4_Num_of_cell_of_each_cluster.txt",sep="\t",quote=F)
################################################################
markergene_list<-c("Ptgds","Plp1","Hexb","Siglech","Mfge8","Gfap","H2-Aa","Ms4a7","Pf4","Pdgfra","Nnat","S100a9","S100a8","Ly6c1","Igfbp7","Cd3e","Cd3g","Col1a1","Col3a1","Tagln","Acta2")
pdf("4_FeaturePlot_of_selected_gene.pdf")
FeaturePlot(Stroke, reduction="umap",features = markergene_list,pt.size = .0001 )
dev.off()  
pdf("4_vlnPlot_of_selected_gene.pdf")
VlnPlot(Stroke, features = markergene_list,ncol=2,pt.size = 0)
dev.off() 
###########对选择的基因进行表达可视化##############
p1 <- FeaturePlot(Stroke, features = c("S100a8","S100a9","Meg3"), combine = FALSE )
library(ggplot2)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'blue'),  limits = c(0, 3))
p2 <- lapply(p1, function (x) x + fix.sc)
pdf("4_FeaturePlot_of_given_gene_zero_start.pdf")
CombinePlots(p2)
dev.off()###
##########看MT基因比例，ncount，nFeature在不同亚群的分布#########################################################
Stroke[["percent.mt"]] <- PercentageFeatureSet(Stroke, pattern = "^mt-")  #percent of mt
pdf("5_QC_MT_Feature_Plot.pdf")
FeaturePlot(Stroke, reduction="umap",features = "percent.mt",pt.size = 0.1)
dev.off()##
pdf("5-QC_NCount_Plot.pdf")
FeaturePlot(Stroke, reduction="umap",features = c("nCount_RNA"),pt.size = 0.1)
dev.off()
pdf("5-QC_NFeature_Plot.pdf")
FeaturePlot(Stroke, reduction="umap",features = c("nFeature_RNA"),pt.size = 0.1)
dev.off()
################################################################
# find all markers of clusters###
Stroke.markers <- FindAllMarkers(Stroke,  min.pct = 0.25, logfc.threshold = 0.25)
top1markers<-Stroke.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)  ##top 1 of each cluster 
write.table(Stroke.markers,"6_allmarkers.resolution2.txt",sep="\t",quote=F)  ##save the findmarker result
######### https://satijalab.org/seurat/seurat_clustering_tutorial_part2.html   

top20markers <- Stroke.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  
write.table(top20markers,"6_top20_markers.resolution2.txt",sep="\t",quote=F)
save(Stroke,file="Stroke.resolution2.rds")
#########################删除cluster  recluster ############################
library(ggplot2)
library(Seurat)
library(dplyr)
setwd("/mnt/zhoushunheng/stroke/2seurat/result5_doubletfinder_harmony_refinedClusterR3/no_harmony/")
load("/mnt/zhoushunheng/stroke/2seurat/result5_doubletfinder_harmony_refinedClusterR3/no_harmony/Stroke.resolution2.rds")
Stroke$old<-Stroke$RNA_snn_res.2
# idents= c(0,2,3,19)  ##remove the unknown cells##
idents= c(0,2,3,19,22)  ##remove the unknown cells##
# idents= c(0,1,2)  ##remove the unknown cells##
cells.use = rownames(Stroke@meta.data)[!Idents(Stroke) %in% idents]
Stroke<- subset(x = Stroke, cells =cells.use , invert = FALSE)  ##get count mat##

resolution<-2
Stroke <- Stroke %>% 
    RunUMAP(dims = 1:20) %>% 
    RunTSNE(dims = 1:20) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = resolution) %>%     identity()
save(Stroke,file="Stroke.recluster.rds")
# write.table(table(Stroke$old,Idents(Stroke)),"remove0231922.txt",sep="\t",quote=FALSE)
# DimPlot(Stroke,reduction="tsne")
table(Stroke$old,Idents(Stroke))
# find all markers of clusters###
Stroke.markers <- FindAllMarkers(Stroke, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20_markers<-Stroke.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  ##top 1 of each cluster 
write.table(top20_markers,"top20_markers.recluster.txt",sep="\t",quote=F)  ##save the findmarker result


# Stroke<-RunTSNE(Stroke,seed.use=3)
DimPlot(Stroke, reduction = "tsne", label = TRUE)
DimPlot(Stroke, reduction = "tsne",group.by="stim", label = TRUE)


load("Stroke.recluster.rds")

# new.cluster.ids <-c("Microglia","BAMs","Microglia","Microglia","Microglia","Astrocytes","Astrocytes","Microglia","Astrocytes","Monocytes","Oligodendrocytes","DCs","Microglia","Microglia","Astrocytes","Granulocytes","OPCs","NK&T cells","OPCs","BAMs","Endothelial cells","Macrophages","Astrocytes","Endothelial cells","Fibroblast-like cells","Astrocytes","Oligodendrocytes","Endothelial cells","DCs","Mural cells","Mural cells","OPCs","Mural cells","Mural cells","Oligodendrocytes","Microglia","Endothelial cells","NK&T cells")
new.cluster.ids <-c("Microglia","BAMs","Microglia","Microglia","Microglia","Microglia","Astrocytes","Astrocytes","Monocytes","Astrocytes","Oligodendrocytes","Microglia","DCs","Microglia","Granulocytes","Astrocytes","OPCs","NK&T cells","OPCs","Endothelial cells","Astrocytes","Endothelial cells","Fibroblast-like cells","Astrocytes","Oligodendrocytes","Endothelial cells","DCs","Mural cells","Mural cells","OPCs","Mural cells","Mural cells","Oligodendrocytes","Endothelial cells","Microglia","NK&T cells")
names(new.cluster.ids) <- levels(Stroke)
Stroke <- RenameIdents(Stroke, new.cluster.ids)  #combine and rename cell clusters 
save(Stroke,file="Stroke.rename.rds") ## ##  as.matrix(GetAssayData(object = Stroke, slot = "counts"))[1,1]




