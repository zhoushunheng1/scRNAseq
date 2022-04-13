##3.6.0  https://github.com/sqjin/CellChat
# https://mp.weixin.qq.com/s/Pz_pEv4RooGKKzNBphWnRw
# devtools::install_github("sqjin/CellChat") 
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html ##
####首先需要更新cellchatDB数据库，添加Lgals1 Lgasl3 Lgasl9
rm(list=ls())
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE)
###更新cellchatDB###
setwd("/mnt/zhoushunheng/stroke/2seurat/result6_cell_cell_communication/cellchatDBupdate")
options(stringsAsFactors = FALSE)
interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = 'geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
setwd("/mnt/zhoushunheng/package/CellChat-master/") # This is the folder of CellChat package downloaded from Github
CellChatDB.mouse <- CellChatDB
usethis::use_data(CellChatDB.mouse, overwrite = TRUE)



options(stringsAsFactors = FALSE)
setwd("/mnt/zhoushunheng/stroke/2seurat/result6_cell_cell_communication/")

load("/mnt/zhoushunheng/stroke/2seurat/result5_doubletfinder_harmony_refinedClusterR3/no_harmony/Stroke.rename.rds")

Stroke@active.ident <- factor(Stroke@active.ident,  
                             # levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","Macrophages","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))


result_file<-paste("/mnt/zhoushunheng/stroke/2seurat/result6_cell_cell_communication/","joint_analysis",sep="")
unlink(result_file, recursive = TRUE)
dir.create(result_file)
setwd(result_file)

Stroke_PT<-subset(x = Stroke, subset= stim =="Stroke" , invert = FALSE)  ##get count mat##
data_input_PT <- GetAssayData(Stroke_PT, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Stroke_PT)
meta_PT <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
##create CellChat object#
cellchat_PT <- createCellChat(object = data_input_PT, meta = meta_PT, group.by = "labels")
######add meta data##
cellchat_PT <- addMeta(cellchat_PT, meta = meta_PT, meta.name = "labels")
cellchat_PT <- setIdent(cellchat_PT, ident.use = "labels") # set "labels" as default cell identity

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling  "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 
cellchat_PT@DB <- CellChatDB.use
cellchat_PT <- subsetData(cellchat_PT) # subset the expression data of signaling genes for saving computation cost
######
future::plan("multiprocess", workers = 4) # do parallel

cellchat_PT <- identifyOverExpressedGenes(cellchat_PT)
cellchat_PT <- identifyOverExpressedInteractions(cellchat_PT)
cellchat_PT <- projectData(cellchat_PT, PPI.mouse)
###2 Inference of cell-cell communication network ###########################################################################################################
cellchat_PT <- computeCommunProb(cellchat_PT, raw.use = TRUE)
# Filter
cellchat_PT <- filterCommunication(cellchat_PT, min.cells = 10)
#######
cellchat_PT <- computeCommunProbPathway(cellchat_PT)
cellchat_PT <- netAnalysis_computeCentrality(cellchat_PT, slot.name = "netP")

cellchat_PT <- aggregateNet(cellchat_PT)
###############################################################################################
Stroke_Sham<-subset(x = Stroke, subset= stim =="Sham" , invert = FALSE)  ##get count mat##
data_input_Sham <- GetAssayData(Stroke_Sham, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Stroke_Sham)
meta_Sham <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
##create CellChat object#
cellchat_Sham <- createCellChat(object = data_input_Sham, meta = meta_Sham, group.by = "labels")
######add meta data##
cellchat_Sham <- addMeta(cellchat_Sham, meta = meta_Sham, meta.name = "labels")
cellchat_Sham <- setIdent(cellchat_Sham, ident.use = "labels") # set "labels" as default cell identity

###############
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling  "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 
cellchat_Sham@DB <- CellChatDB.use
cellchat_Sham <- subsetData(cellchat_Sham) # subset the expression data of signaling genes for saving computation cost
######
future::plan("multiprocess", workers = 4) # do parallel

cellchat_Sham <- identifyOverExpressedGenes(cellchat_Sham)
cellchat_Sham <- identifyOverExpressedInteractions(cellchat_Sham)
cellchat_Sham <- projectData(cellchat_Sham, PPI.mouse)
###2 Inference of cell-cell communication network ###########################################################################################################
cellchat_Sham <- computeCommunProb(cellchat_Sham, raw.use = TRUE)
# Filter
cellchat_Sham <- filterCommunication(cellchat_Sham, min.cells = 10)
#######
cellchat_Sham <- computeCommunProbPathway(cellchat_Sham)
cellchat_Sham <- netAnalysis_computeCentrality(cellchat_Sham, slot.name = "netP")
cellchat_Sham <- aggregateNet(cellchat_Sham)




######################### Lift up CellChat object and merge together
group.new = levels(cellchat_PT@idents)  ###全的那个细胞类型####
cellchat_Sham <- liftCellChat(cellchat_Sham, group.new)
#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...
object.list <- list(PT = cellchat_PT, Sham = cellchat_Sham)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)





###由于细胞类不一致 只能计算结构相似性   要用root  否则需要更新anaconda
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("structural-similarity.pdf")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()
#> 2D visualization of signaling networks from datasets 1 2




weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Number of interactions.pdf")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Interaction weights.pdf")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]))
}
dev.off()




library(ComplexHeatmap)

i=1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i],font.size = 6, width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1],font.size = 6, width = 5, height = 6)
pdf("outgoing.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
##########
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i],font.size = 6, width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1],font.size = 6, width = 5, height = 6, color.heatmap = "GnBu")
pdf("incoming.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i],font.size = 6, width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1],font.size = 6, width = 5, height = 6, color.heatmap = "OrRd")
pdf("all.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()




