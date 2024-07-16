rm(list=ls())
library(Seurat)
library(S4Vectors)
library(ggExtra)
library(dplyr) 
library(ggplot2)
setwd("/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q3_merge/")
load("/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q3_merge/Merge_pre_D1.Rdata")



##QC step:
# Discard cells: in each sample
# (1) with fewer than 1,000 transcripts
# (2) with counts (transcripts/ cell) greater than 3 median absolute deviation (MAD) 
#     away from the median (removing potential doublets and debris) 
# (3) with genes (genes/cell) greater than 3 MAD away from the median 
#     (to remove low-abundance genes or genes with high dropout rate)
# (4) mt % > 10%
# (5) with genes fewer than 500
# (6) doublet

Remove_Cells<-NULL
for (i in unique(Merge$sample)){
  meta_data<-Merge@meta.data
  meta_data<-meta_data[meta_data$sample==i,]
  
  genes_cut<-median(meta_data$nFeature_RNA)+3*mad(meta_data$nFeature_RNA)
  counts_cut<-median(meta_data$nCount_RNA)+3*mad(meta_data$nCount_RNA)
  print(paste0(i," Genes_cut: ",genes_cut))
  print(paste0(i," counts_cut: ",counts_cut))
  
  Cells<-rownames(meta_data)[meta_data$nCount_RNA>counts_cut | meta_data$nFeature_RNA>genes_cut]
  Remove_Cells<-c(Cells,Remove_Cells)
  
}
  
length(Remove_Cells)  
Merge<-Merge[,!colnames(Merge)%in%Remove_Cells]
filtered_data <- subset(x = Merge, subset= nFeature_RNA > 500  & percent.mt < 10 & nCount_RNA > 1000 & doublet_prediction==0)


meta_data<-filtered_data@meta.data
p<-meta_data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() 
  
p
p<-ggMarginal(p, 
           xparams = list(fill = 4),
           yparams = list(fill = 3))
p
p;ggsave("/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q4_gene_cell_mito/nFeature_RNA_nCount_RNA.pdf",width = 9,height = 5)

filtered_data[["percent.ribo"]] <- PercentageFeatureSet(filtered_data, pattern = "^Rps|^Rpl")

##Gene-level filtering
##keep only genes which are expressed in 30 or more cells.
#counts <- GetAssayData(object = filtered_data, slot = "counts")
#counts[1:10,1:10]
##hist(as.matrix(counts))

#####################################
#### Suerat version 5 of above ######
############# testing ###############
filtered_data <- JoinLayers(filtered_data)
#counts <- GetAssayData(object = filtered_data, slot = "counts")
#counts[1:10,1:10]

#counts <- GetAssayData(object = filtered_data, slot = "counts")
counts <- LayerData(filtered_data, assay = "RNA", layer = "counts")






nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 30
filtered_counts <- counts[keep_genes, ]

dim(filtered_data)
filtered_data<-filtered_data[keep_genes,]

# #Remove ribosomal and mitochondrial genes from the analysis
# remove <- rownames(filtered_data)[grep("^Rps|^Rpl|^Mrps|^Mrpl|^mt-", rownames(filtered_data))]
# filtered_data<-filtered_data[!rownames(filtered_data)%in%remove,]
# 
# dim(filtered_data)
# 
# 
# filtered_data


save(filtered_data,file="/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/QCed_data_D1.Rdata")

