########### Fig3a
library(Seurat)
library(dplyr)
library(infercnv)
library(tidyverse)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=groupinfo,
                                    delim="\t",
                                    gene_order_file= geneInfor,
                                    ref_group_names=c("T cell")) 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,              
                             out_dir= "/share/pub/like/script/Infercnv/out", 
                             cluster_by_groups=TRUE)


########### Fig3c
library(maftools)
library(dplyr)
library(SummarizedExperiment)
library(TCGAbiolinks) 
library(maftools)
library(dplyr)
maf = read.maf("bcgsc.ca_UVM.IlluminaHiSeq_DNASeq.1.somatic.maf")
maf@data$Tumor_Sample_Barcode = str_sub(maf@data$Tumor_Sample_Barcode, start = 1, end = 16)

UVM_maf_1 = subsetMaf(maf,tsb = rownames(filter(TCGA_sur,class == "IT_low")))
UVM_maf_2 = subsetMaf(maf,tsb = rownames(filter(TCGA_sur,class == "IT_high")))

plotmafSummary(maf = UVM_maf_1,  
               rmOutlier = TRUE,
               showBarcodes = T,
               textSize = 0.4,
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE )
getSampleSummary(UVM_maf_1)
oncoplot(maf = UVM_maf_1,                
         top = 10,                                
         fontSize = 0.5,                       
         showTumorSampleBarcodes = F)     

plotmafSummary(maf = UVM_maf_2,  
               rmOutlier = TRUE,
               showBarcodes = T,
               textSize = 0.4,
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE )
getSampleSummary(UVM_maf_2)
oncoplot(maf = UVM_maf_2,                
         top = 10,                                
         fontSize = 0.5,                       
         showTumorSampleBarcodes = F)  


########### Fig3f
ggplot(data,aes_string(x = data$ITH_score,y = data$P2)) +
  geom_point(size = 2,color = 'black',alpha = 0.5) +
  theme_bw() +theme(axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    axis.ticks.length = unit(0.25,'cm'),
                    axis.ticks = element_line(size = 1),
                    panel.border = element_rect(size = 1.5),
                    panel.grid = element_blank()) +  
  geom_smooth(method = 'lm',se = T,size = 1.5) + 
  stat_cor(method = "spearman",digits = 3,size=6)

ggplot(data,aes_string(x = data$ITH_score,y = data$P3)) +
  geom_point(size = 2,color = 'black',alpha = 0.5) +
  theme_bw() +theme(axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    axis.ticks.length = unit(0.25,'cm'),
                    axis.ticks = element_line(size = 1),
                    panel.border = element_rect(size = 1.5),
                    panel.grid = element_blank()) +  
  geom_smooth(method = 'lm',se = T,size = 1.5) + 
  stat_cor(method = "spearman",digits = 3,size=6)

ggplot(data,aes_string(x = data$ITH_score,y = data$P4)) +
  geom_point(size = 2,color = 'black',alpha = 0.5) +
  theme_bw() +theme(axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    axis.ticks.length = unit(0.25,'cm'),
                    axis.ticks = element_line(size = 1),
                    panel.border = element_rect(size = 1.5),
                    panel.grid = element_blank()) +  
  geom_smooth(method = 'lm',se = T,size = 1.5) + 
  stat_cor(method = "spearman",digits = 3,size=6)




























