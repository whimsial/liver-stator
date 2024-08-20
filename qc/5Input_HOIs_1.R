rm(list = ls())
setwd("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/")
set.seed(1234)
library(Seurat)
library(sctransform)
library(randomcoloR)
library(ggplot2)
library(loomR)
library(xlsx)
library(rjson)
library(xlsx)
require(stringr)

############################
########## Subset1 #########
############################

load("./run_subset1/QCed_input_D1_subset1.Rdata")
QCed_input_D1_subset1

##Color
table(QCed_input_D1_subset1$sample)

QCed_input_D1_subset1 <- NormalizeData(QCed_input_D1_subset1, normalization.method = "LogNormalize", scale.factor = 10000)
QCed_input_D1_subset1 <- FindVariableFeatures(QCed_input_D1_subset1, selection.method = "vst", nfeatures = 1000)
table(QCed_input_D1_subset1$sample)
# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(QCed_input_D1_subset1), 20)
top20
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(QCed_input_D1_subset1)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2


colnames(QCed_input_D1_subset1@meta.data)

meta_data<-QCed_input_D1_subset1@meta.data[c("sample","Group")]
meta_data
write.csv(meta_data,file="/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/run_subset1/meta.csv")


gene<-VariableFeatures(QCed_input_D1_subset1)
#count<-QCed_input_D1_subset1@assays[["RNA"]]@counts #older version of seurat
# seurat V5: access the counts matrix from the RNA assay
count<-QCed_input_D1_subset1[["RNA"]]$counts
count<-as.matrix(count)
count<-count[gene,]
count<-t(count)

dim(count)

sub_matrix<-count
result <- fromJSON(file = "/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/ReAnalysis/jsonData.json")
expressed_count<-sort(colSums(sub_matrix),decreasing = T)
print(table(expressed_count>0))
genes<-names(expressed_count)
genes
n="run_subset1"

#dir.create(paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n))
write.table(rbind(genes),file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/genes.csv"),append = F,quote = F,row.names = F,col.names = F,sep = ",")
sub_matrix<-sub_matrix[,genes]
print(dim(sub_matrix))
expressed_count
write.csv(expressed_count,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/expressed_count.csv"),append = F,quote = F,row.names = T,col.names = T)
cluster<-rownames(sub_matrix)
cluster<-data.frame(row=cluster,cluster=1)
result_sub<-result
result_sub$rawDataPath<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/count_matrix.csv")
result_sub$nGenes<-length(genes)
result_sub$nCells<-dim(sub_matrix)[1]
result_sub$userGenes<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/genes.csv")
result_sub$clusterFile<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/clusters.csv")
jsonData <- toJSON(result_sub)
write(jsonData, paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/jsonData.json"),ncolumns = 1)
write.csv(cluster,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/clusters.csv"),append = F,quote = F,row.names = F,col.names = T)
write.csv(sub_matrix,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/count_matrix.csv"),append = F,quote = F,row.names = T,col.names = T)

## Get the full count matrix ##

#srt<-get(load("/Users/yaoyuelin/Documents/GoF_DNMT3A/ReAnalysed/filtered_data_normalised.Rdata"))
srt <- QCed_input_D1_subset1
count<- LayerData(srt, assay = "RNA", layer = "counts")
count<-as.matrix(count)
count<-t(count)
write.csv(count,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/shiny_app_count_martix_all_run_subset1.csv"),append = F,quote = F,row.names = T,col.names = T)



############################
########## Subset2 #########
############################

load("./run_subset2/QCed_input_D1_subset2.Rdata")
QCed_input_D1_subset2

##Color
table(QCed_input_D1_subset2$sample)

QCed_input_D1_subset2 <- NormalizeData(QCed_input_D1_subset2, normalization.method = "LogNormalize", scale.factor = 10000)
QCed_input_D1_subset2 <- FindVariableFeatures(QCed_input_D1_subset2, selection.method = "vst", nfeatures = 1000)
table(QCed_input_D1_subset2$sample)
# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(QCed_input_D1_subset2), 20)
top20
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(QCed_input_D1_subset2)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2


colnames(QCed_input_D1_subset2@meta.data)

meta_data<-QCed_input_D1_subset2@meta.data[c("sample","Group")]
meta_data
write.csv(meta_data,file="/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/run_subset2/meta.csv")


gene<-VariableFeatures(QCed_input_D1_subset2)
#count<-QCed_input_D1_subset2@assays[["RNA"]]@counts #older version of seurat
# seurat V5: access the counts matrix from the RNA assay
count<-QCed_input_D1_subset2[["RNA"]]$counts
count<-as.matrix(count)
count<-count[gene,]
count<-t(count)

dim(count)

sub_matrix<-count
result <- fromJSON(file = "/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/ReAnalysis/jsonData.json")
expressed_count<-sort(colSums(sub_matrix),decreasing = T)
print(table(expressed_count>0))
genes<-names(expressed_count)
genes
n="run_subset2"

#dir.create(paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n))
write.table(rbind(genes),file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/genes.csv"),append = F,quote = F,row.names = F,col.names = F,sep = ",")
sub_matrix<-sub_matrix[,genes]
print(dim(sub_matrix))
expressed_count
write.csv(expressed_count,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/expressed_count.csv"),append = F,quote = F,row.names = T,col.names = T)
cluster<-rownames(sub_matrix)
cluster<-data.frame(row=cluster,cluster=1)
result_sub<-result
result_sub$rawDataPath<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/count_matrix.csv")
result_sub$nGenes<-length(genes)
result_sub$nCells<-dim(sub_matrix)[1]
result_sub$userGenes<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/genes.csv")
result_sub$clusterFile<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/clusters.csv")
jsonData <- toJSON(result_sub)
write(jsonData, paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/jsonData.json"),ncolumns = 1)
write.csv(cluster,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/clusters.csv"),append = F,quote = F,row.names = F,col.names = T)
write.csv(sub_matrix,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/count_matrix.csv"),append = F,quote = F,row.names = T,col.names = T)

## Get the full count matrix ##

#srt<-get(load("/Users/yaoyuelin/Documents/GoF_DNMT3A/ReAnalysed/filtered_data_normalised.Rdata"))
srt <- QCed_input_D1_subset2
count<- LayerData(srt, assay = "RNA", layer = "counts")
count<-as.matrix(count)
count<-t(count)
write.csv(count,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/shiny_app_count_martix_all_run_subset2.csv"),append = F,quote = F,row.names = T,col.names = T)






############################
########## Run1 ##########
############################

load("QCed_input_D1_Run1.Rdata")
QCed_input_D1_Run1

##Color
table(QCed_input_D1_Run1$sample)

QCed_input_D1_Run1 <- NormalizeData(QCed_input_D1_Run1, normalization.method = "LogNormalize", scale.factor = 10000)
QCed_input_D1_Run1 <- FindVariableFeatures(QCed_input_D1_Run1, selection.method = "vst", nfeatures = 1000)
table(QCed_input_D1_Run1$sample)
# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(QCed_input_D1_Run1), 20)
top20
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(QCed_input_D1_Run1)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2


colnames(QCed_input_D1_Run1@meta.data)

meta_data<-QCed_input_D1_Run1@meta.data[c("sample","Group")]
meta_data
write.csv(meta_data,file="/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/meta.csv")


gene<-VariableFeatures(QCed_input_D1_Run1)
#count<-QCed_input_D1_Run1@assays[["RNA"]]@counts #older version of seurat
# seurat V5: access the counts matrix from the RNA assay
count<-QCed_input_D1_Run1[["RNA"]]$counts
count<-as.matrix(count)
count<-count[gene,]
count<-t(count)

dim(count)

sub_matrix<-count
result <- fromJSON(file = "/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/ReAnalysis/jsonData.json")
expressed_count<-sort(colSums(sub_matrix),decreasing = T)
print(table(expressed_count>0))
genes<-names(expressed_count)
genes
n="run1"

dir.create(paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n))
write.table(rbind(genes),file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/genes.csv"),append = F,quote = F,row.names = F,col.names = F,sep = ",")
sub_matrix<-sub_matrix[,genes]
print(dim(sub_matrix))
expressed_count
write.csv(expressed_count,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/expressed_count.csv"),append = F,quote = F,row.names = T,col.names = T)
cluster<-rownames(sub_matrix)
cluster<-data.frame(row=cluster,cluster=1)
result_sub<-result
result_sub$rawDataPath<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/count_matrix.csv")
result_sub$nGenes<-length(genes)
result_sub$nCells<-dim(sub_matrix)[1]
result_sub$userGenes<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/genes.csv")
result_sub$clusterFile<-paste0("/exports/igmm/eddie/khamseh-lab/ava/ME_scRNA/analysis/",n,"/clusters.csv")
jsonData <- toJSON(result_sub)
write(jsonData, paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/jsonData.json"),ncolumns = 1)
write.csv(cluster,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/clusters.csv"),append = F,quote = F,row.names = F,col.names = T)
write.csv(sub_matrix,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/",n,"/count_matrix.csv"),append = F,quote = F,row.names = T,col.names = T)

## Get the full count matrix ##

#srt<-get(load("/Users/yaoyuelin/Documents/GoF_DNMT3A/ReAnalysed/filtered_data_normalised.Rdata"))
srt <- QCed_input_D1_Run1
count<- LayerData(srt, assay = "RNA", layer = "counts")
count<-as.matrix(count)
count<-t(count)
write.csv(count,file=paste0("/Volumes/khamseh-lab/ava/ME_scRNA/data_analysis/run1/shiny_app_count_martix_all_run1.csv"),append = F,quote = F,row.names = T,col.names = T)






# ############################
# ########## Yuelin ##########
# ############################
#
#
# #### analysed data: /Users/yaoyuelin/Desktop/GoF_DNMT3A/Mouse_WT_and_GOF_DNMT3A_Progenitors/Analysis_03/04_WT_and_Mutant_Merged_All_Lineages/05_Chosen_Parameters/Mouse_WT_and_Mutant_All_Lineage_Cell_Cycle_Regresssed_02_dim_20_res_0.5.RDS
# rm(list = ls())
# srt<-readRDS("/Users/yaoyuelin/Desktop/GoF_DNMT3A/Mouse_WT_and_GOF_DNMT3A_Progenitors/Analysis_03/04_WT_and_Mutant_Merged_All_Lineages/05_Chosen_Parameters/Non_Cell_Cycle_Regressed/Mouse_WT_and_Mutant_All_Lineage_02_dim_20_res_0.5.RDS")
#
# table(srt$condition)
#
# srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
# srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 1000)
# #"Pf4"    "Hbb-bs" "Hba-a1" "Hbb-bt" "Prss34" "Prg2"
# # Identify the 20 most highly variable genes
# top20 <- head(VariableFeatures(srt), 20)
# top20
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(srt)
# plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
# plot1 + plot2
#
# colnames(srt@meta.data)
# srt$condition
# srt$names<-paste0(srt$names_WT,srt$names_Mutant)
# srt$names<-gsub("NA","",srt$names)
# srt$Names<-paste0(srt$condition,"-",srt$names)
# meta_data<-srt@meta.data[c("Names","condition","names")]
# meta_data
# write.csv(meta_data,file="/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/Analysed/meta.csv")
# write.csv(meta_data,file="/Users/yaoyuelin/Desktop/GoF_DNMT3A/Analysed/Meta.csv")
# write.csv(srt@reductions[["umap"]]@cell.embeddings,file="/Users/yaoyuelin/Desktop/GoF_DNMT3A/Analysed/umap.csv")
#
# gene<-VariableFeatures(srt)
# count<-srt@assays[["RNA"]]@counts
# count<-as.matrix(count)
# count<-count[gene,]
# count<-t(count)
#
# dim(count)
#
# sub_matrix<-count
# result <- fromJSON(file = "/Volumes/khamseh-lab/Yuelin/pancancer/SKSC/jsonData.json")
# expressed_count<-sort(colSums(sub_matrix),decreasing = T)
# print(table(expressed_count>0))
# genes<-names(expressed_count)
# genes
# n="Analysed"
#
# dir.create(paste0("/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/",n))
# write.table(rbind(genes),file=paste0("/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/",n,"/genes.csv"),append = F,quote = F,row.names = F,col.names = F,sep = ",")
# sub_matrix<-sub_matrix[,genes]
# print(dim(sub_matrix))
# expressed_count
# write.csv(expressed_count,file=paste0("/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/",n,"/expressed_count.csv"),append = F,quote = F,row.names = T,col.names = T)
# cluster<-rownames(sub_matrix)
# cluster<-data.frame(row=cluster,cluster=1)
# result_sub<-result
# result_sub$rawDataPath<-paste0("/exports/eddie/scratch/s1914230/",n,"/count_matrix.csv")
# result_sub$nGenes<-length(genes)
# result_sub$nCells<-dim(sub_matrix)[1]
# result_sub$userGenes<-paste0("/exports/eddie/scratch/s1914230/",n,"/genes.csv")
# result_sub$clusterFile<-paste0("/exports/eddie/scratch/s1914230/",n,"/clusters.csv")
# jsonData <- toJSON(result_sub)
# write(jsonData, paste0("/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/",n,"/jsonData.json"),ncolumns = 1)
# write.csv(cluster,file=paste0("/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/",n,"/clusters.csv"),append = F,quote = F,row.names = F,col.names = T)
# write.csv(sub_matrix,file=paste0("/Volumes/khamseh-lab/Yuelin/GoF_DNMT3A_Single_Cell_Jordan_Analysis/",n,"/count_matrix.csv"),append = F,quote = F,row.names = T,col.names = T)
#
# #####################################################################
#
#
#
# srt<-get(load("/Users/yaoyuelin/Documents/GoF_DNMT3A/ReAnalysed/filtered_data_normalised.Rdata"))
# count<-srt@assays[["RNA"]]@counts
# count<-as.matrix(count)
# count<-t(count)
# write.csv(count,file=paste0("/Users/yaoyuelin/Documents/GoF_DNMT3A/Yuelin_GoF_DNMT3A/ReAnalysed/Stator/ShinyApp_Input/count_martix_all.csv"),append = F,quote = F,row.names = T,col.names = T)
#
# srt<-readRDS("/Users/yaoyuelin/Documents/GoF_DNMT3A/Mouse_WT_and_GOF_DNMT3A_Progenitors/Analysis_03/04_WT_and_Mutant_Merged_All_Lineages/05_Chosen_Parameters/Non_Cell_Cycle_Regressed/Mouse_WT_and_Mutant_All_Lineage_02_dim_20_res_0.5.RDS")
# count<-srt@assays[["RNA"]]@counts
# count<-as.matrix(count)
# count<-t(count)
# write.csv(count,file=paste0("/Users/yaoyuelin/Documents/GoF_DNMT3A/Yuelin_GoF_DNMT3A/Analysed/Stator/ShinyApp_Input/count_martix_all.csv"),append = F,quote = F,row.names = T,col.names = T)
#
#
#
