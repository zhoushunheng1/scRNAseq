################## remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') ### install doubletfinder
#####remove potential doublets#############################################################
rm(list=ls())
library(dplyr)
library(cowplot)
library(Seurat)
library(DoubletFinder)
library(ggplot2) # 
doublet_rate<-0.016
##################### 
PT.data <- Read10X(data.dir = "./PT/")
PT.data <- CreateSeuratObject(counts = PT.data, project = "PT", min.cells = 3)  ##18377  6796
load_number<-dim(PT.data)
PT.data[["percent.mt"]] <- PercentageFeatureSet(PT.data, pattern = "^mt-")
PT.data <- subset(PT.data, subset = nFeature_RNA > 200  & percent.mt < 20)  ##18377  6559
MT_number<-dim(PT.data)
PT.data <- NormalizeData(PT.data)
PT.data <- FindVariableFeatures(PT.data, selection.method = "vst", nfeatures = 2000)
PT.data <- ScaleData(PT.data)
PT.data <- RunPCA(PT.data)
PT.data <- RunUMAP(PT.data, dims = 1:10)
PT.data <- RunTSNE(PT.data, dims = 1:10)
#####pre-process####################
sweep.res.list <- paramSweep_v3(PT.data, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
######## Homotypic Doublet Proportion Estimate    ##ref:  A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart
nExp_poi <- round(doublet_rate*ncol(PT.data))        ## Assuming 1.6% doublet 
PT.data <- doubletFinder_v3(PT.data, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pdf("Supplementary.Fig.S4A.PT.pdf")
DimPlot(PT.data,pt.size = 1,label=TRUE, label.size = 5,reduction = "tsne",group.by = "DF.classifications_0.25_0.3_105" )+theme(aspect.ratio = 1)
dev.off()
PT.data<- subset(PT.data, DF.classifications_0.25_0.005_105=="Singlet" )  ##18377  6231
DoubletFinder_number<-dim(PT.data)

write.table(as.matrix(GetAssayData(object = PT.data, slot = "counts")), 
            'PT.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)
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
#####pre-process#############################################################################
sweep.res.list <- paramSweep_v3(Sham.data, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
######## Homotypic Doublet Proportion Estimate    ##ref:  A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart
nExp_poi <- round(doublet_rate*ncol(Sham.data))        ## Assuming 1.6% doublet 
Sham.data <- doubletFinder_v3(Sham.data, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# table(Sham.data$DF.classifications_0.25_0.005_328)
pdf("Supplementary.Fig.S4A.Sham.pdf")
DimPlot(Sham.data,pt.size = 1,label=TRUE, label.size = 5,reduction = "tsne",group.by = "DF.classifications_0.25_0.07_137" )+theme(aspect.ratio = 1)
dev.off()
Sham.data<- subset(Sham.data, DF.classifications_0.25_0.15_137=="Singlet" )  ##18377  6231
DoubletFinder_number<-dim(Sham.data)

##############
write.table(as.matrix(GetAssayData(object = Sham.data, slot = "counts")), 
            'Sham.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)
##########process data with seurat###################################################################################
rm(list=ls())
library(dplyr)
library(cowplot)
library(Seurat)
resolution<-2  ## set the resolution
####load the expression matrix of the doubletfinder################ 
PT.sparse<-read.csv("./PT/PT.csv")  ## read CSV
Sham.sparse<-read.csv("./Sham/Sham.csv")

common_genes<-intersect(rownames(PT.sparse),rownames(Sham.sparse))
PT_specific_genes<-setdiff(rownames(PT.sparse),rownames(Sham.sparse))
Sham_specific_genes<-setdiff(rownames(Sham.sparse),rownames(PT.sparse))
####merge the two expression matrix########################
PT.sparse_common<-PT.sparse[common_genes,]
PT.sparse_specific<-PT.sparse[PT_specific_genes,]
PT.sparse_zero<-matrix(0,ncol=ncol(PT.sparse),nrow=length(Sham_specific_genes))
colnames(PT.sparse_zero)<-colnames(PT.sparse_common)
rownames(PT.sparse_zero)<-Sham_specific_genes
PT.sparse<-rbind(PT.sparse_common,PT.sparse_specific,PT.sparse_zero)
###############################
Sham.sparse_common<-Sham.sparse[common_genes,]
Sham.sparse_specific<-Sham.sparse[Sham_specific_genes,]
Sham.sparse_zero<-matrix(0,ncol=ncol(Sham.sparse),nrow=length(PT_specific_genes))
colnames(Sham.sparse_zero)<-colnames(Sham.sparse_common)
rownames(Sham.sparse_zero)<-PT_specific_genes
Sham.sparse<-rbind(Sham.sparse_common,Sham.sparse_zero,Sham.sparse_specific)
##############################
colnames(Sham.sparse)<- sub(pattern = "1", replacement = "2", x = colnames(Sham.sparse))##unique colname
####create the seurat object####################################
PT <- CreateSeuratObject(counts = cbind(PT.sparse, Sham.sparse), project = "PT", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = PT@var.genes, npcs = 50, verbose = FALSE)

PT@meta.data$stim <- c(rep("PT", ncol(PT.sparse)), rep("Sham", ncol(Sham.sparse)))#

#######################dimension reduction##############
PT <- PT %>% 
	RunTSNE( dims = 1:20) %>% 
	FindNeighbors( dims = 1:20) %>% 
	FindClusters(resolution = resolution) %>%   identity()
# find all markers of clusters###
PT.markers <- FindAllMarkers(PT,  min.pct = 0.25, logfc.threshold = 0.25)

top20markers <- PT.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  
#########################remove cluster and reclustering ############################

PT$previous<-PT$RNA_snn_res.2
idents= c(0,2,3,19,22)  ##remove  cells##
cells.use = rownames(PT@meta.data)[!Idents(PT) %in% idents]
PT<- subset(x = PT, cells =cells.use , invert = FALSE)  ##get count mat##

resolution<-2
PT <- PT %>% 
    RunTSNE(dims = 1:20) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = resolution) %>%     identity()
save(PT,file="PT.recluster.rds")
PT.markers <- FindAllMarkers(PT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20_markers<-PT.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  ##top 1 of each cluster 
##rename each cluster based on marker genes#
new.cluster.ids <-c("Microglia","BAMs","Microglia","Microglia","Microglia","Microglia","Astrocytes","Astrocytes","Monocytes","Astrocytes","Oligodendrocytes","Microglia","DCs","Microglia","Granulocytes","Astrocytes","OPCs","NK&T cells","OPCs","Endothelial cells","Astrocytes","Endothelial cells","Fibroblast-like cells","Astrocytes","Oligodendrocytes","Endothelial cells","DCs","Mural cells","Mural cells","OPCs","Mural cells","Mural cells","Oligodendrocytes","Endothelial cells","Microglia","NK&T cells")
names(new.cluster.ids) <- levels(PT)
PT <- RenameIdents(PT, new.cluster.ids) #
save(PT,file="PT.rename.rds") ##
################Figure 3A ######################################################################################
cols<-c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0","#999999")
x<-as.data.frame(PT@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(PT)
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
PT_PT<- subset(PT, stim == "PT")  ##
x<-as.data.frame(PT_PT@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(PT_PT)
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
PT_Sham<- subset(PT, stim == "Sham")  ##
x<-as.data.frame(PT_Sham@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(PT_Sham)
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

PT@active.ident <- factor(PT@active.ident, 
                            levels=cell_order)
Num<-table(Idents(PT),PT$stim)
Group <- c(rep("PT" , length(cell_order)) , rep("Sham" ,  length(cell_order))  )
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
#########Figure 3B Dotplot######################################################################################
PT@active.ident <- factor(PT@active.ident,  
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))
dotplot_markergene<-c("Tmem119","Siglech","Ms4a7","Pf4","Mfge8","Aqp4","Pdgfra","Nnat","Ptgds","Plp1","Ccr2","Ly6c2","H2-Aa","Cd74","S100a9","S100a8","Cd3e","Cd3g","Ly6c1","Igfbp7","Tagln","Acta2","Col1a1","Col3a1")
pdf("Figure3B.Dotplot.pdf",width=11,height=5) ######dotplot of marker gene
DotPlot(object = PT, features = dotplot_markergene,cols =c("#4F7936","#FECA0A"))+theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.4,
 colour=c("#FED023","#FED023","#CE5227","#CE5227","#5492CD","#5492CD","#E12F8B","#E12F8B","#833C64","#833C64","#23936A","#23936A","#FF7F05","#FF7F05","#CD9201","#CD9201","#33BEB6","#33BEB6","#4B4E4F","#4B4E4F","#5B5FAA","#5B5FAA","#9D73B0","#9D73B0","#999999","#999999")),   
  axis.text.y = element_text(  colour=c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0","#999999")))
dev.off()

####### draw the barplot###################################
Count_No_Feature<-function(cell_type){
	new_data <- subset(x = PT, idents = c(cell_type), invert = FALSE)
	return(c(mean(new_data$nFeature_RNA),##gene##
			mean(new_data$nCount_RNA))  )##UMI ##
}
cell_list<-as.character(unique(Idents(PT)))
cell_type_nFeature_result<-c()
for(i in 1:length(cell_list)){
	cell_type_nFeature_result<-rbind(cell_type_nFeature_result,c(Count_No_Feature(cell_list[i])))
}
colnames(cell_type_nFeature_result)<-c("mean.detected.gene","mean.UMI"); rownames(cell_type_nFeature_result)<-cell_list;
barplot_data <- data.frame(cell.type=as.character(rownames(cell_type_nFeature_result)),mean.UMI=cell_type_nFeature_result[,2])
barplot_data$cell.type <- factor(barplot_data$cell.type,levels = rev(c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","Macrophages","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells")))
gg<-ggplot(barplot_data, aes( x=cell.type,y=mean.UMI)) + 
    geom_bar(stat = "identity",width = 0.7,fill=c("#999999"))+
	theme_classic()
ggsave(gg,file="Figure3B_barplot.pdf",width=12,height=5)
markergene_list<-c("Siglech","Ms4a7","Aqp4","Pdgfra","Ptgds","Ly6c2","H2-Aa","S100a9","Cd3e","Ly6c1","Acta2","Col1a1")
pdf("Figure3C.FeaturePlot.pdf")
FeaturePlot(PT, reduction="tsne",features = markergene_list,pt.size = .001, cols = c("#a1a3a6", "#840228") )
dev.off()  
###################################################################
pltdata<-Stroke$nFeature_RNA
sd<-sd(pltdata)
pdf("FigS4B.nFeature.sd.pdf")

hist(pltdata, breaks = 60, main = "Distribution of unique genes per cell", col = "#E5FBB3", border="#74c69d",cex.main=2,
        xlab = expression("Number of unique genes per cell"), ylab = "Number of cells", cex.lab = 1.7, cex.axis = 1.7)
abline(v=mean(pltdata), col="black",lty=5,lwd=3)	
abline(v=c(mean(pltdata)-sd,mean(pltdata)+sd),lty=5, col="red",lwd=3)	
legend("topright", c(paste("mean =",round(mean(pltdata)),sep=" "), paste("SD =",round(sd),sep=" ")), col=c("black", "red"), bty="n",lwd=3)

dev.off()
########################
pltdata<-Stroke$nCount_RNA
sd<-sd(pltdata)
pdf("FigS3_histplot_nCount_RNA.sd.pdf")

hist(pltdata, breaks = 60, main = "Distribution of unique transcripts per cell", col = "#f0ead2", border="#dde5b6",cex.main=2,
        xlab = expression("Number of unique transcripts per cell"), ylab = "Number of cells", cex.lab = 1.7, cex.axis = 1.7)
abline(v=mean(pltdata), col="black",lty=5,lwd=3)	
abline(v=c(mean(pltdata)-sd,mean(pltdata)+sd),lty=5, col="red",lwd=3)	
legend("topright", c(paste("mean =",round(mean(pltdata)),sep=" "), paste("SD =",round(sd),sep=" ")), col=c("black", "red"), bty="n",lwd=3)

dev.off()
###


##################################################################################################################################
##########cell cell communication analysis############################################################
##update the cellchatDB###########













