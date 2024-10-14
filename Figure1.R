library(Seurat)
library(dplyr)
library(ggpubr)
library(harmony)

LH16_3814 = Read10X_h5("GSM4107899_LH16.3814_raw_gene_bc_matrices_h5.h5")
LH17_364 = Read10X_h5("GSM4107900_LH17.364_raw_gene_bc_matrices_h5.h5")
LH17_530 = Read10X_h5("GSM4107901_LH17.530_raw_gene_bc_matrices_h5.h5")
LH17_3222 = Read10X_h5("GSM4107902_LH17.3222_raw_gene_bc_matrices_h5.h5")
LH17_3554 = Read10X_h5("GSM4107903_LH17.3554_raw_gene_bc_matrices_h5.h5")
LH18_277 = Read10X_h5("GSM4107904_LH18.277_raw_gene_bc_matrices_h5.h5")

LH16_3814 = CreateSeuratObject(counts = LH16_3814, project = "LH16_3814", min.features = 100)
LH17_364 = CreateSeuratObject(counts = LH17_364, project = "LH17_364", min.features = 100)
LH17_530 = CreateSeuratObject(counts = LH17_530, project = "LH17_530", min.features = 100)
LH17_3222 = CreateSeuratObject(counts = LH17_3222, project = "LH17_3222", min.features = 100)
LH17_3554 = CreateSeuratObject(counts = LH17_3554, project = "LH17_3554", min.features = 100)
LH18_277 = CreateSeuratObject(counts = LH18_277, project = "LH18_277", min.features =100)

UMM069 = Read10X("UMM069")
UMM067L = Read10X("UMM067L")
UMM066 = Read10X("UMM066")
UMM065 = Read10X("UMM065")
UMM064 = Read10X("UMM064")
UMM063 = Read10X("UMM063")
UMM062 = Read10X("UMM062")
UMM061 = Read10X("UMM061")
UMM059 = Read10X("UMM059")
UMM041L = Read10X("UMM041L")
BSSR0022 = Read10X("BSSR0022")

UMM069 = CreateSeuratObject(counts = UMM069, project = "UMM069", min.features = 200)
UMM067L = CreateSeuratObject(counts = UMM067L, project = "UMM067L", min.features = 200)
UMM066 = CreateSeuratObject(counts = UMM066, project = "UMM066", min.features = 200)
UMM065 = CreateSeuratObject(counts = UMM065, project = "UMM065", min.features = 200)
UMM064 = CreateSeuratObject(counts = UMM064, project = "UMM064", min.features = 200)
UMM063 = CreateSeuratObject(counts = UMM063, project = "UMM063", min.features = 200)
UMM062 = CreateSeuratObject(counts = UMM062, project = "UMM062", min.features = 200)
UMM061 = CreateSeuratObject(counts = UMM061, project = "UMM061", min.features = 200)
UMM059 = CreateSeuratObject(counts = UMM059, project = "UMM059", min.features = 200)
UMM041L = CreateSeuratObject(counts = UMM041L, project = "UMM041L", min.features = 200)
BSSR0022 = CreateSeuratObject(counts = BSSR0022, project = "BSSR0022", min.features = 200)

GSE138665  = merge(LH16_3814,y=c(LH17_364,LH17_530,LH17_3222,LH17_3554,LH18_277))
GSE139829 = merge(UMM069,y=c(UMM067L,UMM066,UMM065,UMM064,UMM063,UMM062,UMM061,UMM059,UMM041L,BSSR0022))


data  = merge(GSE139829,y= GSE138665)
data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-") 
data[["percent.rb"]] = PercentageFeatureSet(data, pattern = "^RP[SL]")
data2 = subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10) 
data2$cb = substring(rownames(data2@meta.data),1,20)
load('data_3.Rdata')

data3@meta.data$celltype='Malignant cell'
data3@meta.data[data3$seurat_clusters %in% c(28),'celltype'] = 'Photoreceptor cell'
data3@meta.data[data3$seurat_clusters %in% c(26),'celltype'] = 'Endothelial cell'
data3@meta.data[data3$seurat_clusters %in% c(2,6,19,24,17),'celltype'] = 'T/ NK cell'
data3@meta.data[data3$seurat_clusters %in% c(15,25,20),'celltype'] = 'B/ Plasma cell'
data3@meta.data[data3$seurat_clusters %in% c(27),'celltype'] = 'Fibroblast'
data3@meta.data[data3$seurat_clusters %in% c(7,23),'celltype'] = 'Myeloid cell'


########### Fig1a-c
library(ggpubr)
DimPlot(data3,reduction = 'tsne',group.by = 'seurat_clusters',label = T)+NoLegend()
DimPlot(data3,reduction = 'tsne',group.by = 'celltype',label = T,
       cols= c("#5CB7BC","#CF9043","#53B400","pink","#A084B7","#D7709F","#54AFDE")) +NoLegend()

data3$seurat_clusters = factor(data3$seurat_clusters,levels = c(28,26,27,23,7,20,15,25,19,6,2,24,17,13,3,18,5,
                                                                8,14,11,10,4,30,31,22,16,29,1,0,21,12,32,9) %>% rev())
Idents(data3) = 'seurat_clusters'
DotPlot(data3,
        features =  c("MLANA","MITF","CD3D","CD3E","MS4A1","CD38","SDC1",
                      "C1QA","CD14","CD68","DCN","C1R","C1S",
                      "PECAM1","CLDN5","CD34","RCVRN") %>% rev(),              
        dot.scale = 8,cols = c( "white","#CF292C"),
        cluster.idents=F) + RotatedAxis()+xlab(NULL)+ylab(NULL)+coord_flip()


########### Supplementary Fig1a
DimPlot(data3,reduction = 'tsne',group.by = 'orig.ident',label = F)+ggtitle('Patient origin')


########### Fig1d
Idents(data3) = 'seurat_clusters'
marker = FindAllMarkers(data3)
top10 = marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
AverageExp = AverageExpression(data3,features=  top10$gene,slot = "data")         
AverageExp = AverageExp$SCT
coorda = psych::corr.test(AverageExp,AverageExp,method="spearman") 
pheatmap::pheatmap(coorda$r,clustering_method = "ward")




########### Fig1e-f
library(harmony)
load('tumor.Rdata')

DimPlot(tumor,group.by = 'SCT_snn_res.0.5',label = T,reduction = 'tsne')
DimPlot(tumor,group.by = 'orig.ident',label = F,reduction = 'tsne')



########### Fig1g
library(ggalluvial)
vaccinations
tumor$orig.ident
plot_data = tumor@meta.data

ggplot(data = plot_data,
       aes(axis1 = SCT_snn_res.0.5, axis2 = orig.ident)) +
  geom_alluvium() +
  geom_stratum()



########### Fig1f
library(monocle)
library(Seurat)
library(dplyr)
load('cds_1.Rdata')
plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.0001,show_branch_points=F)+ 
  theme(legend.position = "right") + scale_colour_gradient(low = "#FAF5A9", high = "#C52F2D")
plot_cell_trajectory(cds, color_by = "pure",cell_size = 0.0001,show_branch_points=F)+ 
  theme(legend.position = "right") + scale_color_manual(values=c("#3FB8EA","#7DC241"))


########### Supplementary Fig1b
ITH = data.frame(row.names = as.character(unique(tumor$orig.ident)))
for (sample in as.character(unique(tumor$orig.ident))) {
  
  tumor_sample = subset(tumor@meta.data, orig.ident == sample) 
  tumor_sample = subset(tumor, cells=row.names(tumor_sample))
  
  pc_data = data.frame(row.names = rownames(tumor_sample@reductions$pca@cell.embeddings))
  for (j in rownames(tumor_sample@reductions$pca@cell.embeddings)) {
    for (i in 1:50) {
      pc = (tumor_sample@reductions$pca@cell.embeddings[j,i] - mean(tumor_sample@reductions$pca@cell.embeddings[,i]))^2
      pc_data[j,i] = pc
    }  
    pc_data[j,51] = sqrt(rowSums(pc_data[j,1:30]))
  }  
  
  ITH[sample,1] = sum(pc_data[,51])/dim(pc_data)[1]
}


########### Supplementary Fig1c
ITH$class = ifelse(ITH$V1 > median(ITH$V1),'ITMHhi','ITMHlo')
library(grafify)
library(ggpubr)
library(ggthemes)

plot_scatterbox(ITH,   class,    V1,       
                symsize = 2, fontsize = 25,
                jitter = 0)+stat_compare_means()


########### Supplementary Fig1d
tumor@meta.data[tumor@meta.data$orig.ident== "LH16_3814","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_364","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_530","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_3222","GEP_class"] = "1b"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_3554","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "LH18_277","GEP_class"] = "1b"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM059","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM061","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM062","GEP_class"] = "1a"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM063","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM064","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM065","GEP_class"] = "1a"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM066","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM069","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "BSSR0022","GEP_class"] = "1b"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM067L","GEP_class"] = "2"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM041L","GEP_class"] = "2"


tumor$Sex = 'M'
tumor@meta.data[tumor@meta.data$orig.ident== "LH16_3814","Sex"] = "F"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM059","Sex"] = "F"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM061","Sex"] = "F"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM065","Sex"] = "F"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM069","Sex"] = "F"
tumor@meta.data[tumor@meta.data$orig.ident== "BSSR0022","Sex"] = "F"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM041L","Sex"] = "F"

tumor@meta.data[tumor@meta.data$orig.ident== "LH16_3814","Tumor_diameter"] = "14"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_364","Tumor_diameter"] = "18"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_530","Tumor_diameter"] = "19"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_3222","Tumor_diameter"] = "10"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_3554","Tumor_diameter"] = "15"
tumor@meta.data[tumor@meta.data$orig.ident== "LH18_277","Tumor_diameter"] = "17"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM059","Tumor_diameter"] = "15"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM061","Tumor_diameter"] = "19"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM062","Tumor_diameter"] = "9"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM063","Tumor_diameter"] = "20"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM064","Tumor_diameter"] = "20"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM065","Tumor_diameter"] = "10"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM066","Tumor_diameter"] = "18.4"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM069","Tumor_diameter"] = "19.8"
tumor@meta.data[tumor@meta.data$orig.ident== "BSSR0022","Tumor_diameter"] = "NA"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM067L","Tumor_diameter"] = "NA"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM041L","Tumor_diameter"] = "NA"
tumor$Histological = 'Mix'
tumor@meta.data[tumor@meta.data$orig.ident== "UMM066","Histological"] = "Epithelioid"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM066","Histological"] = "Epithelioid"
tumor@meta.data[tumor@meta.data$orig.ident== "LH16_3814","Histological"] = "Spindle"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_364","Histological"] = "Epithelioid"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_3222","Histological"] = "Spindle"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_3554","Histological"] = "Spindle"
tumor@meta.data[tumor@meta.data$orig.ident== "LH18_277","Histological"] = "Spindle"

tumor$metastasis = 'No'
tumor@meta.data[tumor@meta.data$orig.ident== "BSSR0022","metastasis"] = "Yes"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM067L","metastasis"] = "Yes"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM041L","metastasis"] = "Yes"
tumor@meta.data[tumor@meta.data$orig.ident== "LH16_3814","Age"] = "84"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_364","Age"] = "69"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_530","Age"] = "84"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_3222","Age"] = "65"
tumor@meta.data[tumor@meta.data$orig.ident== "LH17_3554","Age"] = "85"
tumor@meta.data[tumor@meta.data$orig.ident== "LH18_277","Age"] = "31"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM059","Age"] = "86"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM061","Age"] = "71"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM062","Age"] = "69"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM063","Age"] = "66"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM064","Age"] = "77"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM065","Age"] = "44"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM066","Age"] = "53"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM069","Age"] = "80"
tumor@meta.data[tumor@meta.data$orig.ident== "BSSR0022","Age"] = "68"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM067L","Age"] = "73"
tumor@meta.data[tumor@meta.data$orig.ident== "UMM041L","Age"] = "63"


metadata =group_by(tumor@meta.data[,c('orig.ident','GEP_class','Sex','Tumor_diameter','Histological','metastasis','Age')], orig.ident) %>% unique()
metadata = data.frame(metadata)
rownames(metadata) = metadata$orig.ident
metadata = metadata[rownames(ITH),]
metadata$ITH = ITH$V1

metadata$ITH2 = metadata$ITH
pheatmap::pheatmap(metadata[,c('ITH','ITH2')] %>% t(),annotation_col = metadata[,c(2:7)],cluster_cols = F)

metadata$orig.ident = factor(metadata$orig.ident,levels = metadata$orig.ident)
ggplot(metadata, aes(x=orig.ident, y= ITH)) +
  geom_bar(position = "dodge",stat = "identity",width=0.6)+ ylim(0,40)


########### Supplementary Fig1e
tumor$ITH = 'ITHMlo'
tumor@meta.data[tumor$orig.ident %in% rownames(filter(ITH,class == 'ITMHhi')),'ITH'] = 'ITMHhi'

DimPlot(tumor,reduction = 'tsne',group.by = 'ITH',cols = c('#5098D3','#D45F3A'))




########### Supplementary Fig1f
library(ggthemes)
tsne = data.frame(row.names = as.character(unique(tumor$orig.ident)))
for (i in as.character(unique(tumor$orig.ident))) {
  
  a = tumor@reductions$tsne@cell.embeddings[rownames(filter(tumor@meta.data,orig.ident == i)),1:2]
  tsne[i,'tSNE_1'] = mean(a[,1] %>% as.numeric())
  tsne[i,'tSNE_2'] = mean(a[,2] %>% as.numeric())
  
}

ggplot(data=tsne,aes(x= tSNE_1,y= tSNE_2))+geom_point(shape=21,size= 6) +theme_few()+
  geom_text(aes(tSNE_1, tSNE_2, label = rownames(tsne)))


