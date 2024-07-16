library(Seurat)
library(xlsx) 
require(stringr)

setwd("/Volumes/khamseh-lab/ava/ME_scRNA/raw_data/")
# all samples
files<-c("GSM6603110_COR-7349-D1","GSM6603113_COR-6215-D1","GSM6603115_COR-3717-D1","GSM6603116_COR-1126-D1","GSM6603118_COR-8330-D1","GSM6603121_COR-8763-D1","GSM6603123_COR-8843-D1","GSM6603124_COR-9546-D1","GSM6603126_COR-5723-D1","GSM6603129_COR-1481-D1","GSM6603131_COR-5218-D1","GSM6603133_COR-8468-D1","GSM6603134_COR-6302-D1","GSM6603137_COR-9198-D1","GSM6603138_COR-7256-D1","GSM6603140_COR-3195-D1","GSM6603142_COR-6176-D1","GSM6603145_COR-9818-D1","GSM6603146_COR-2465-D1","GSM6603149_COR-7414-D1","GSM6603150_COR-2674-D1","GSM6603152_COR-8066-D1","GSM6603154_COR-6685-D1","GSM6603156_COR-4856-D1","GSM6603158_COR-1203-D1","GSM6603160_COR-7101-D1","GSM6603162_COR-5623-D1","GSM6603164_COR-4181-D1","GSM6603166_COR-3427-D1","GSM6603168_COR-3315-D1","GSM6603170_COR-5878-D1","GSM6603172_COR-5554-D1","GSM6603175_COR-5136-D1","GSM6603176_COR-3985-D1","GSM6603178_COR-1002-D1","GSM6603180_COR-6913-D1","GSM6603182_COR-5087-D1","GSM6603184_COR-4616-D1","GSM6603186_COR-5429-D1","GSM6603188_COR-6308-D1","GSM6603190_COR-9348-D1","GSM6603192_COR-4994-D1","GSM6603194_COR-8852-D1","GSM6603196_COR-3309-D1","GSM6603198_COR-8111-D1","GSM6603200_COR-9643-D1","GSM6603202_COR-8625-D1","GSM6603204_COR-3722-D1","GSM6603206_COR-4271-D1","GSM6603208_COR-5797-D1","GSM6603210_COR-7717-D1","GSM6603212_COR-6711-D1","GSM6603214_COR-1141-D1","GSM6603216_COR-4617-D1","GSM6603218_COR-9868-D1","GSM6603220_COR-4209-D1","GSM6603222_COR-9633-D1","GSM6603224_COR-5449-D1")
# first sample on the list
sample="GSM6603110_COR-7349-D1"
{
  
  doublet_prediction<-read.csv(file = paste0("/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q2_doublet/",sample,"_scrublet_EDR0.008_PredictedDoublets.csv"),header = F)
  doublet_score<-read.csv(file = paste0("/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q2_doublet/",sample,"_scrublet_EDR0.008_DoubletScores.csv"),header = F)
  
  ##barcodes<-read.table(paste0("Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/barcodes.tsv"),header = F)
  ##barcodes
  
  ##rownames(doublet_prediction)<-barcodes$V1
  ##system(paste0("gzip ","Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/*"))
  #pre<- Read10X(data.dir = paste0("Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/"))
  pre<- Read10X_h5(paste0('/Volumes/khamseh-lab/ava/ME_scRNA/raw_data/',sample,'_filtered_feature_bc_matrix.h5'), use.names = TRUE, unique.features = TRUE)
  pre <- CreateSeuratObject(counts = pre, project = sample)
  
  #table(colnames(pre)==barcodes$V1)
  pre$barcodes<- rownames(pre@meta.data)
  pre@meta.data
  
  pre$sample<-sample
  scvelo_index<-paste0(pre$sample,":",pre$barcodes)
  scvelo_index
  scvelo_index<-gsub(pattern = '-1',replacement = "x",scvelo_index)
  
  # different samples may have the same barcode, so add the name of the sample to the barcode
  pre<-RenameCells(object = pre,new.names =scvelo_index )
  
  #doublet_prediction$V1<-1
  # 0 means to doublet. Change threshold if not satisfied with the original threshold
  #doublet_prediction$V1[doublet_score$V1<0.35]<-0
  
  pre[["doublet_prediction"]]<-doublet_prediction$V1
  pre[["doublet_score"]]<-doublet_score$V1
  pre[["percent.mt"]] <- PercentageFeatureSet(pre, pattern = "^mt-")
  pre$log10GenesPerUMI <- log10(pre$nFeature_RNA) / log10(pre$nCount_RNA)
  
  #pdf(paste0("/Users/yaoyuelin/Desktop/QC/",sample,"_vlnPlot.pdf"),width = 9,height = 5)
  #p<-VlnPlot(pre, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
  #print(p)
  
  #dev.off()
  
  
  #plot1 <- FeatureScatter(pre, feature1 = "nCount_RNA", feature2 = "percent.mt")
  #plot2 <- FeatureScatter(pre, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  #pdf(paste0("/Users/yaoyuelin/Desktop/QC/",sample,"_feature_feature_relationships.pdf"),width = 8,height = 4)
  #p<-plot1 + plot2
  #print(p)
  #dev.off()
  
  
}


Merge<-pre
files
#sample="MAB04-2-Mutant"
for (sample in files[-1]){
  
  print(sample)
  doublet_prediction<-read.csv(file = paste0("/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q2_doublet/",sample,"_scrublet_EDR0.008_PredictedDoublets.csv"),header = F)
  doublet_score<-read.csv(file = paste0("/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q2_doublet/",sample,"_scrublet_EDR0.008_DoubletScores.csv"),header = F)
  
  ##barcodes<-read.table(paste0("Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/barcodes.tsv"),header = F)
  ##barcodes
  
  ##rownames(doublet_prediction)<-barcodes$V1
  ##system(paste0("gzip ","Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/*"))
  #pre<- Read10X(data.dir = paste0("Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/"))
  pre<- Read10X_h5(paste0('/Volumes/khamseh-lab/ava/ME_scRNA/raw_data/',sample,'_filtered_feature_bc_matrix.h5'), use.names = TRUE, unique.features = TRUE)
  pre <- CreateSeuratObject(counts = pre, project = sample)
  
  #table(colnames(pre)==barcodes$V1)
  pre$barcodes<- rownames(pre@meta.data)
  pre@meta.data
  
  pre$sample<-sample
  scvelo_index<-paste0(pre$sample,":",pre$barcodes)
  scvelo_index
  scvelo_index<-gsub(pattern = '-1',replacement = "x",scvelo_index)
  
  # different samples may have the same barcode, so add the name of the sample to the barcode
  pre<-RenameCells(object = pre,new.names =scvelo_index )
  
  #doublet_prediction$V1<-1
  # 0 means to doublet. Change threshold if not satisfied with the original threshold
  #doublet_prediction$V1[doublet_score$V1<0.35]<-0
  
  pre[["doublet_prediction"]]<-doublet_prediction$V1
  pre[["doublet_score"]]<-doublet_score$V1
  pre[["percent.mt"]] <- PercentageFeatureSet(pre, pattern = "^mt-")
  pre$log10GenesPerUMI <- log10(pre$nFeature_RNA) / log10(pre$nCount_RNA)
  
  
  Merge<-merge(x=pre,y=Merge) 
  
  
}

#save(Merge,file="/Users/akhamseh/Dropbox/interactions/ME_CFS/ME_scRNA_seq/QC3/Merge_pre_D1.Rdata")
save(Merge,file="/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q3_merge/Merge_pre_D1.Rdata")


######################################
######### Merge with metadata ########
######################################

# load meta data
meta_data <- read.xlsx("/Volumes/khamseh-lab/ava/ME_scRNA/meta_data/ME_metadata.xlsx",sheetIndex=1, colNames = TRUE)

setwd("/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q3_merge/")
load("Merge_pre_D1.Rdata")

#setwd("~/Dropbox/interactions/ME_CFS/ME_scRNA_seq/QC3/")
#load("Merge_pre_D1.Rdata")

Merge

# The last two string in sample name gives the time point (D1 or D2)
Merge$time <- str_sub(Merge$sample, start= -2)

# Strings to extract individual ID (characters 12-19 inclusive)
Merge$cor_id <- str_sub(Merge$sample, start= 12, end=19)

# add information on phenotype, under column name "Group"
Merge$Group <- meta_data$phenotype[match(Merge$cor_id, meta_data$cor_id)]

# add information on age, under column name "age_binned"
Merge$age_binned <- meta_data$age_binned[match(Merge$cor_id, meta_data$cor_id)]

# add information on bmi, under column name "bmi_binned"
Merge$bmi_binned <- meta_data$bmi_binned[match(Merge$cor_id, meta_data$cor_id)]

# add information on mecfs_duration_binned, under column name "mecfs_duration_binned"
Merge$mecfs_duration_binned <- meta_data$mecfs_duration_binned[match(Merge$cor_id, meta_data$cor_id)]

# add information on general health, under column name "gh_binned"
Merge$gh_binned <- meta_data$gh_binned[match(Merge$cor_id, meta_data$cor_id)]

# add information on post exercise, under column name "pcs_binned"
Merge$pcs_binned <- meta_data$pcs_binned[match(Merge$cor_id, meta_data$cor_id)]

# add information on Multidimensional Fatigue Inventory (MFI-20), under column name "mfi20_total_binned"
Merge$mfi20_total_binned <- meta_data$mfi20_total_binned[match(Merge$cor_id, meta_data$cor_id)]

# add information on post-exertional malaise, under column name "pem_max_delta_binned"
Merge$pem_max_delta_binned <- meta_data$pem_max_delta_binned[match(Merge$cor_id, meta_data$cor_id)]


#save(Merge,file="/Users/akhamseh/Dropbox/interactions/ME_CFS/ME_scRNA_seq/QC3/Merge_pre_D1.Rdata")
save(Merge,file="/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/Q3_merge/Merge_pre_D1.Rdata")





##############################
######### From Yuelin ########
##############################


Merge<-pre
files
sample="MAB04-2-Mutant"
for (sample in files[-1]){
  
  print(sample)
  doublet_prediction<-read.csv(file = paste0("/Users/yaoyuelin/Desktop/GoF_DNMT3A/",sample,"_scrublet_EDR0.08_PredictedDoublets.csv"),header = F)
  doublet_score<-read.csv(file = paste0("/Users/yaoyuelin/Desktop/GoF_DNMT3A/",sample,"_scrublet_EDR0.08_DoubletScores.csv"),header = F)
  
  barcodes<-read.table(paste0("Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/barcodes.tsv"),header = F)
  barcodes
  rownames(doublet_prediction)<-barcodes$V1
  
  system(paste0("gzip ","Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/*"))
  pre<- Read10X(data.dir = paste0("Raw_Data/",sample,"/outs/filtered_feature_bc_matrix/"))
  pre <- CreateSeuratObject(counts = pre, project = sample)
  
  table(colnames(pre)==barcodes$V1)
  pre$barcodes<-barcodes$V1
  pre@meta.data
  
  pre$sample<-sample
  scvelo_index<-paste0(pre$sample,":",pre$barcodes)
  scvelo_index
  scvelo_index<-gsub(pattern = '-1',replacement = "x",scvelo_index)
  
  
  pre<-RenameCells(object = pre,new.names =scvelo_index )
  
  
  doublet_prediction$V1<-1
  doublet_prediction$V1[doublet_score$V1<0.26]<-0
  
  pre[["doublet_prediction"]]<-doublet_prediction$V1
  pre[["doublet_score"]]<-doublet_score$V1
  pre[["percent.mt"]] <- PercentageFeatureSet(pre, pattern = "^mt-")
  pre$log10GenesPerUMI <- log10(pre$nFeature_RNA) / log10(pre$nCount_RNA)
  
  #pdf(paste0("/Users/yaoyuelin/Desktop/QC/",sample,"_vlnPlot.pdf"),width = 9,height = 5)
  #p<-VlnPlot(pre, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
  #print(p)
  
  #dev.off()
  
  
  #plot1 <- FeatureScatter(pre, feature1 = "nCount_RNA", feature2 = "percent.mt")
  #plot2 <- FeatureScatter(pre, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  #pdf(paste0("/Users/yaoyuelin/Desktop/QC/",sample,"_feature_feature_relationships.pdf"),width = 8,height = 4)
  #p<-plot1 + plot2
  #print(p)
  #dev.off()
  
  Merge<-merge(x=pre,y=Merge) 
  
  
  
}


Merge$Group<-Merge$sample
load("/Users/yaoyuelin/Desktop/GoF_DNMT3A/emptydrop.Rdata")
Merge<-Merge[,colnames(Merge)%in%Cell]
save(Merge,file="/Users/yaoyuelin/Desktop/GoF_DNMT3A/Merge_pre.Rdata")

