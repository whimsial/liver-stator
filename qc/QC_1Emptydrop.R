## BiocManager::install("DropletUtils")
## load 
library(DropletUtils)

## Set working directory
setwd("/Volumes/khamseh-lab/ava/ME_scRNA/")
meta_data <- read.csv("me_metadata.txt")
files<-c()

sample="MAB04-1-WT"
Cell=NULL
for (sample in files){
  
  print(sample)
  dir.name<-paste0("./Raw_Data/",sample,"/outs/raw_feature_bc_matrix/")
  dir.name
  #system(paste0("gunzip ",dir.name,"/*"))
  barcodes<-read.table(paste0(dir.name,"/barcodes.tsv"),header = F)
  barcodes
  sce <- read10xCounts(dir.name)
  set.seed(100)
  my.count<-counts(sce)
  colnames(my.count)<-paste0(sample,":",barcodes$V1)
  colnames(my.count)<-gsub(pattern = '-1',replacement = "x",colnames(my.count))
  print("starting emptyDrops")
  e.out <- emptyDrops(my.count)
  head(e.out)
  #https://support.bioconductor.org/p/123554/#123562
  is.cell <- e.out$FDR <= 0.01
  is.cell[is.na(is.cell)]<-FALSE
  True_Cell<-colnames(my.count)[is.cell]
  Cell<-c(Cell,True_Cell)
  print("done")
  
}

length(Cell)
save(Cell,file="/Users/yaoyuelin/Desktop/GoF_DNMT3A/emptydrop.Rdata")

