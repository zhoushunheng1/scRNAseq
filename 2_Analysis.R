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
# find all markers of clusters###
Stroke.markers <- FindAllMarkers(Stroke,  min.pct = 0.25, logfc.threshold = 0.25)

top20markers <- Stroke.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  
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
##rename each cluster based on marker genes#
new.cluster.ids <-c("Microglia","BAMs","Microglia","Microglia","Microglia","Microglia","Astrocytes","Astrocytes","Monocytes","Astrocytes","Oligodendrocytes","Microglia","DCs","Microglia","Granulocytes","Astrocytes","OPCs","NK&T cells","OPCs","Endothelial cells","Astrocytes","Endothelial cells","Fibroblast-like cells","Astrocytes","Oligodendrocytes","Endothelial cells","DCs","Mural cells","Mural cells","OPCs","Mural cells","Mural cells","Oligodendrocytes","Endothelial cells","Microglia","NK&T cells")
names(new.cluster.ids) <- levels(Stroke)
Stroke <- RenameIdents(Stroke, new.cluster.ids) #
save(Stroke,file="Stroke.rename.rds") ##
################Figure 3A ######################################################################################
cols<-c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0","#999999")
x<-as.data.frame(Stroke@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(Stroke)
cell.type<-as.data.frame(cell.type)
cell.type$orig.ident<-rownames(cell.type)  
c<-merge(x,cell.type,by='orig.ident')
c$cell.type<- factor(c$cell.type,  
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))

xlim<-c(min(c$tSNE_1)-3,max(c$tSNE_1)+2)
ylim<-c(min(c$tSNE_2)-2,max(c$tSNE_2)+2)
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))

stat_density <- ggplot(c,aes(tSNE_1,tSNE_2,colour=cell.type)) + xlim(xlim)+ ylim(ylim)+
	stat_density_2d(size = 0.25, colour = "#B5BABE",h=bwidth ,contour_var = "ndensity") 
	
gg =stat_density+geom_point(data = c, aes(x = tSNE_1, y = tSNE_2),size=0.4,stroke=0.8) +  ##node size, stroke size ##
	scale_color_manual(values=cols)+
	       
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+2,mean(ylim)),legend.justification = c(0, 1))+
	guides(colour = guide_legend(override.aes = list(size=2)))+ 
	theme(legend.text=element_text(size=8),legend.title=element_blank())+
	 theme_classic()	

ggsave(gg,file="Figure3A-tSNE1.pdf",width=7.5,height=6)
#################
Stroke_Stroke<- subset(Stroke, stim == "Stroke")  ##
x<-as.data.frame(Stroke_Stroke@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(Stroke_Stroke)
cell.type<-as.data.frame(cell.type)

cell.type$orig.ident<-rownames(cell.type)  

c<-merge(x,cell.type,by='orig.ident')
c$cell.type<- factor(c$cell.type,  
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))

bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
##
gg =stat_density+geom_point(data = c, aes(x = tSNE_1, y = tSNE_2),size=0.4,stroke=0.8) +  ##node size, stroke size ##
	scale_color_manual(values=cols)+
	       
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+2,mean(ylim)),legend.justification = c(0, 1))+
	guides(colour = guide_legend(override.aes = list(size=2)))+ 
	theme(legend.text=element_text(size=8),legend.title=element_blank())+
	 theme_classic()	
	
ggsave(gg,file="Figure3A-PT_tSNE.pdf",width=7.5,height=6)
###############
Stroke_Sham<- subset(Stroke, stim == "Sham")  ##
x<-as.data.frame(Stroke_Sham@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(Stroke_Sham)
cell.type<-as.data.frame(cell.type)

cell.type$orig.ident<-rownames(cell.type)  

c<-merge(x,cell.type,by='orig.ident')
c$cell.type<- factor(c$cell.type,     #######sham  remove                                                           Granulocytes######
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))
cols<-                              c("#FED023",         "#5492CD", "#E12F8B", "#833C64",       "#23936A","#FF7F05",             "#33BEB6",   "#4B4E4F",          "#5B5FAA",    "#9D73B0")
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
##
gg =stat_density+geom_point(data = c, aes(x = tSNE_1, y = tSNE_2),size=0.4,stroke=0.8) +  ##node size, stroke size ##
	scale_color_manual(values=cols)+
	       
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+2,mean(ylim)),legend.justification = c(0, 1))+
	guides(colour = guide_legend(override.aes = list(size=2)))+ 
	theme(legend.text=element_text(size=8),legend.title=element_blank())+
	 theme_classic()
ggsave(gg,file="Figure3A-Sham_tSNE.pdf",width=7.5,height=6)
############
cell_order<-c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells")

Stroke@active.ident <- factor(Stroke@active.ident, 
                            levels=cell_order)
Num<-table(Idents(Stroke),Stroke$stim)
Group <- c(rep("Stroke" , length(cell_order)) , rep("Sham" ,  length(cell_order))  )
cell.type <- c(cell_order,cell_order)
percent <- c(Num[,2],Num[,1])
data <- data.frame(Group,cell.type ,percent)
data$cell.type <- factor(data$cell.type,levels =cell_order)
# Stacked + percent
cols<-c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0","#999999")
gg<-ggplot(data, aes(fill=cell.type, y=percent, x=Group)) + 
    scale_fill_manual(values=cols)+
	geom_bar( position="fill",stat="identity",width = 0.7)+
	theme_classic()
		
ggsave(gg,file="Figure3A_proportion.pdf")	
#################################################################

