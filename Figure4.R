########### Fig4a
Tcell2@commands
Tcell2 = SCTransform(Tcell2)
Tcell2 = RunPCA(Tcell2, features = VariableFeatures(Tcell2))
Tcell2 = RunUMAP(Tcell2, reduction = "harmony", dims = 1:20)
Tcell2 = RunTSNE(Tcell2, reduction = "harmony", dims = 1:20)
Tcell2 = FindNeighbors(Tcell2, reduction = "harmony", dims = 1:20)
Tcell2 = FindClusters(Tcell2, resolution = 0.8)

Bcell = SCTransform(Bcell)
Bcell = RunPCA(Bcell, features = VariableFeatures(Tcell2))
Bcell = RunUMAP(Bcell, reduction = "harmony", dims = 1:20)
Bcell = RunTSNE(Bcell, reduction = "harmony", dims = 1:20)
Bcell = FindNeighbors(Bcell, reduction = "harmony", dims = 1:20)
Bcell = FindClusters(Bcell, resolution = 0.8)

Mono = SCTransform(Mono)
Mono = RunPCA(Mono, features = VariableFeatures(Tcell2))
Mono = RunUMAP(Mono, reduction = "harmony", dims = 1:20)
Mono = RunTSNE(Mono, reduction = "harmony", dims = 1:20)
Mono = FindNeighbors(Mono, reduction = "harmony", dims = 1:20)
Mono = FindClusters(Mono, resolution = 0.8)

DimPlot(Tcell2,group.by = 'celltype',label = T)
DimPlot(Bcell,group.by = 'celltype',label = T)
DimPlot(Mono,group.by = 'celltype',label = T)



########### Fig4e
library(Seurat)
library(dplyr)
library(AUCell)
library(SCENIC)

scenicOptions <- initializeScenic(org="hgnc", dbDir="/share/pub/like/script/UVM/SCENIC/cisTarget_databases")

exprMat = data@assays$RNA@counts
exprMat = as.matrix(exprMat)
cellinfo = data.frame(data@meta.data)
cellinfo = cellinfo[,c("celltype")]
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir="/share/pub/like/script/UVM/SCENIC/cisTarget_databases" , nCores=1) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

pheatmap::pheatmap(regulonActivity_byCellType,show_colnames = F,cluster_cols = F,border = F)




########### Fig4F-G
library(Seurat)
library(dplyr)
library(CellChat)
data4 = subset(data3@meta.data, celltype %in% c("tumor cell","Mono","Tcell","Bcell"))
data4 = subset(data3, cells=row.names(data4))

data_list = SplitObject(data4, split.by = "ITH") 

data.input  <- data_list$'IT_high'@assays$SCT@data
identity = data.frame(group =data_list$'IT_high'$celltype2, row.names = names(data_list$'IT_high'$celltype2)) 
unique(identity$group)
cellchat <- createCellChat(data.input, meta = identity, group.by = "group")
cellchat <- addMeta(cellchat, meta = identity)
cellchat <- setIdent(cellchat, ident.use = "group")
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.human 
dplyr::glimpse(CellChatDB$interaction)   # Show the structure of the database
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_IT_high = cellchat

data.input  <- data_list$'IT_low'@assays$SCT@data
identity = data.frame(group =data_list$'IT_low'$celltype2, row.names = names(data_list$'IT_low'$celltype2)) 
unique(identity$group)
cellchat <- createCellChat(data.input, meta = identity, group.by = "group")
cellchat <- addMeta(cellchat, meta = identity)
cellchat <- setIdent(cellchat, ident.use = "group")
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.human 
dplyr::glimpse(CellChatDB$interaction)   # Show the structure of the database
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_IT_low = cellchat

object.list <- list(IT_high = cellchat_IT_high, IT_low = cellchat_IT_low)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netAnalysis_signalingRole_scatter(object.list[[1]], title = names(object.list)[1])
netAnalysis_signalingRole_scatter(object.list[[2]], title = names(object.list)[2])

netVisual_bubble(cellchat, sources.use = 'tumor cell',
                 targets.use = c('Bcell','CD4T','CD8TEX','Prol TEX','DC',
                                 'TAM','Monocyte','NK','PlasmaCell','Treg'),comparison = c(1, 2))




































