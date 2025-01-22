## This script is to run QC and pre-processing in Stator with already downloaded
## data. Here we use example of ME/CFS and LongCovid dataset provided by HiFiBio
## Install required dependencies on HPC.
#' Minimum required R version is 4.3.1 (4.4 recommended) and gcc version 9.4.0
#' On EDDIE run: `module load roslin/R/4.4.0`
#' Some packages can be installed from CRAN, others require BiocManager, or
#' installation from GitHub. Many packages are built from source so this step
#' may take some time (and sometimes break) but it has to be executed only once.
#'
#' Installing R packages with complex dependnecies on Eddie is a tedius and
#' complicated task and many things may throw errors. The script
#' `dependencies.R` contains code which I used to install the packages. The
#' order of installation and specific package versions are what worked for me.
#'
#' @param install.dependencies Logical specifying whether to install required
#'        packages.
#' @param my.lib Full path to the writable R library. The first user library is
#'        used by default (change this accordingly to where your R library)
## -----------------------------------------------------------------------------
install.dependencies <- FALSE ## set to TRUE to install
## To ensure installation on the store rather in the personal folder
my.lib <- "/exports/igmm/eddie/ponting-lab/sbraich2/Rlibrary"
.libPaths(c(my.lib, .libPaths()))

if (install.dependencies) source("dependencies.R")

## Load R packages.
library(data.table)
library(Seurat)
library(stringr)
library(DropletUtils)
library(HDF5Array)
library(biomaRt)
library(ggplot2)
library(cowplot)
library (scDblFinder)

read.file <- function(file.path) {
    # Read the file using fread
    dt <- fread(file.path)  # You can add more parameters to fread as needed
    dt[,file.path := file.path]
    return(dt)
}

select.clonotypes <- function(selected.clonotypes, this.cell.map, this.pattern) {
    this.clonotypes <- all.data.clonotypes[grepl(this.pattern, file.path), .(barcode, chain, raw_clonotype_id, umis, file.path)]
    this.clonotypes[,sample:=this.cell.map[grepl(this.pattern, sample)]$sample]
    selected.clonotypes <- rbind(selected.clonotypes, this.clonotypes)
    return(selected.clonotypes)
}

#' @param root.dir Full path to large storage where the project data is to be
#'        stored.
#' @param working.dir Full path to the working directory (typically the cloned
#'        repository).
#' @param map Selected samples from the dataset that will be included
#'        in the analysis mapped to the paths.
#' @param sample.map Size of map X 3 includes paths to all .mtx triplets.
#' @param metadata Samples of interest mapped to the directory 
#' @param cell.map Each cell mapped to condition: "healthy", "desease", etc. 
## -----------------------------------------------------------------------------

root.dir <- "/exports/igmm/eddie/ponting-lab/ava/ME_CSF_Hifibio/raw_data/"
working.dir <- "/exports/igmm/eddie/ponting-lab/sbraich2/liver-stator/"
# ## Load helper functions
source(file.path(working.dir, "rnaseq.functions.R"))

## Reading and unzipping step. Select samples, extract and process files
## -----------------------------------------------------------------------------
## Selected samples(ME samples: MECFS001 - no CD8CD4, HD015v3 - broken data for CD8CD4)
map <- data.table(sample=c("HD0069V1", "HD066V1", "D044", "D085", 
                                  "HD055V2", "HD053V1", "HD010V2", "HD015v3", 
                                  "HD034V2", "HD033V1", "LCOVID001", "FSDD817V2"), 
                         condition=c("healthy", "healthy", "healthy", "healthy", 
                                     "me", "me", "me", "me", 
                                     "lc", "lc", "lc", "lc"))

## Map the exsisting archives to the selected samples
all.archives <- data.table()
for (idx in seq_len(nrow(map))) {
    this.sample <- map[idx, sample]
    this.pattern <- sprintf("_%s_.*GEX.zip", this.sample)
    this.archives <- data.table(sample=this.sample, 
                                path=list.files(root.dir, pattern=this.pattern, 
                                full.names=TRUE))
    all.archives <- rbind(all.archives, this.archives)
}
## Join map with path to archives
map <- map[all.archives, on="sample"]
## Add paths to data directories
map[, sample.dir := gsub(".zip", "", path)]
## Add actual sample names based on names of directories
## because we have data with different cell types for a single sample
map[, sample.full := paste0(basename(sample.dir), "_", sample)]

## Remove HD015v3 CD8CD4 cells because the data is broken
this.pattern <- "_HD015v3_.*_CD8CD4_"
map <- map[!grepl(this.pattern, sample.dir),]

## extract archives
for (idx in seq_len(nrow(map))) {
    cmd <- sprintf("unzip %s -d %s", map[idx, path], root.dir)
    ## Uncomment this for running unzip folders
    #system(cmd)
}
## We have nested archives with .mtx triplet files
## extracted nested archives
nested.map <- data.table()
for (idx in seq_len(nrow(map))) {
    this.sample <- map[idx, sample]
    this.sample.full <- map[idx, sample.full]
    this.sample.dir <- map[idx, sample.dir]
    this.dir <- file.path(this.sample.dir, "outs/raw_feature_bc_matrix")
    this.pattern <- "*.gz"
    this.map <- data.table(sample=this.sample, sample.full=this.sample.full,
                           nested.archives=list.files(this.dir, pattern=this.pattern, full.names=TRUE))
    this.map[, nested.extracts := gsub(".gz", "", nested.archives)]
    for (idx2 in seq_len(nrow(this.map))) {
        cmd <- sprintf("gzip -c %s > %s", this.map[idx2, nested.archives], this.map[idx2, nested.extracts])
        ## Uncomment this for running unzip .mtx triplets
        #system(cmd)
    }
    nested.map <- rbind(nested.map, this.map)
}

## `Merge dirnames and mtx file paths
unique.map <- unique(map, by="sample.full")
sample.map <- unique.map[nested.map, on="sample.full"]
sample.map[, sample.dir := dirname(nested.archives)]

mtx.files <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
## scrublet threshold
doublet.threshold <- 0.15
this.file.type <- ".mtx"

## STEP 1. Read 10X data from sample files.
#'
#' @param mtx.files list of a triplet of file names usually
#'        `barcodes.tsv`, `genes.tsv` (here `features.tsv`) and `matrix.mtx`
#' @param this.extracts A list of full paths to mtx.files  which 
#'         contain counts for a given sample.
#' @param this.sample.dir Full path to the directory containing counts data for
#'         a given sample.
## -----------------------------------------------------------------------------
## STEP 2. Detect empty drops in scRNA-seq data.
#'
#' @param data 10x dataset for a single sample loaded using read.10x.data
#'        function.
#' @param this.sample Sample ID whose counts data is in the `data` variable.
## -----------------------------------------------------------------------------
## STEP 3. Detect doublets.
## Here we have two options, the second one is recommended.
## OPTION 1: scDblFinder
## Use package scDblFinder https://f1000research.com/articles/10-979/v2 to filter
## out doublets in the dataset
## OPTION 2: scrublet
## First we run Jupyter notebook (`doublets.ipynb` in this repository) manually
## to work out suitable threshold and then run the python script for all samples

#' @param data$sce input data for scDblFinder is in SingleCellExperiment format
#' @param doublet.threshold Numeric specifying the doublet cutoff. Default 0.15
#'        seems to work quite well for several tested datasets.
#' @param this.sample.dir Full path to the directory containing counts data for
#'        a given sample.
## -----------------------------------------------------------------------------
for (this.sample.dir in sample.map[, unique(sample.dir)]) {
    this.map <- sample.map[sample.dir==eval(this.sample.dir)]
    this.sample <- this.map[1, sample]
    data <- read.10x.data(this.map$nested.extracts, this.sample.dir, mtx.files)
    true.cells <- remove.emptydrops(sce=data$sce, sample=this.sample,
                                    sample.barcodes=data$sample.barcodes$V1)
    true.cells.ss <- remove.sample.from.barcode(true.cells)
    ## Save barcodes of the true cells for Scrublet
    true.cells.dt <- data.table(true.cells=true.cells.ss)
    true.cells.file <- file.path(this.sample.dir,
                                 "true_cells.csv")
    fwrite(true.cells.dt, file=true.cells.file, col.names=FALSE)
    
    ## run python script `doublet.py` from R using system command
    msg(bold, "Detecting doublets")
    cmd <- sprintf("python3 doublet.py --sample_dir %s --data_type %s \\
                    --doublet_threshold %f", this.sample.dir, this.file.type,
                    doublet.threshold)
    system(cmd)
}

## From all processd data now we want to select CD8 and CD4CD8 cells only
## Create metadata table with CD8 and CD4CD8 cells only
this.pattern <- "_All_"
subsample.map <- sample.map[!grepl(this.pattern, sample.dir), .(sample, sample.full, sample.dir)]
metadata <- subsample.map[, unique(.SD), .SDcols=c("sample", "sample.full", "sample.dir")]

## For all sells we run only selected samples
sel.samples <- c("HD066V1", "HD053V1", "HD033V1")
add.metadata <- data.table()
for (this.sample in sel.samples){
    this.pattern <- sprintf("_%s_.*_All_", this.sample)
    subsample.map <- sample.map[grepl(this.pattern, sample.dir), .(sample, sample.full, sample.dir)]
    single.subsample <- subsample.map[, unique(.SD), .SDcols=c("sample", "sample.full", "sample.dir")]
    add.metadata <- rbind(add.metadata, single.subsample)
}
metadata <- rbind(metadata, add.metadata)

## Add samples for all ME, LC and heathy controls, HD066V1, HD053V1 and HD033V1, respectively

msg(bold, "Reading single cell data to Seurat and merging")
seurat.all <- process.samples.and.merge(metadata, output.dir=root.dir, doublet.method.scDblFinder=FALSE)

## FIXME: Sort out barcodes issues in a general way
seurat.all@meta.data$barcodes <- rownames(seurat.all@meta.data)
## rewrite indents with new extended sample identifiers
new.idents <- gsub("(.*)_.*", "\\1", rownames(seurat.all@meta.data))
seurat.all@meta.data$orig.ident <- new.idents

msg(bold, "Performing QC of a merged Seurat object")
seurat.qc <- qc.seurat(seurat.all, output.dir=root.dir)

## Select highly variable genes
seurat.for.stator <- process.variable.genes(seurat.qc,
                                            output.dir=root.dir)
hvg <- seurat.for.stator@assays$RNA@var.features

## Cell map with all samples and conditions and runs
cell.map <- metadata[sample.map, on="sample.dir", condition:=condition]
cell.map[, run := c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,1,1,2,2,3,3,4,4,5,5,5)]
# cell.map[, sample := paste0(basename(dirname(dirname(sample.dir))), "_", sample)]
cell.map[, sample.simple := sample]
cell.map[, sample := sample.full]
cell.map[, sample.full := NULL]


## Number of cells to keep in each run
n.cells.keep <- 20000

## Sample cells for each run
all.runs <- list()
for (this.run in cell.map[, unique(run)]) {
    this.cell.map <- cell.map[run==eval(this.run)]
    all.runs[this.run] <- sample.seurat(seurat.for.stator, output.dir=root.dir,
                                        cell.map=this.cell.map, 
                                        n.cells.keep=n.cells.keep)
    this.run.file <- file.path(root.dir, 
                               sprintf("%s.seurat.for.stator.RDS", this.run))
    saveRDS(all.runs[this.run], file=this.run.file)
}

genes.file <- file.path(root.dir, "genes.csv")
write.table(hvg, file=genes.file, sep=",", quote=FALSE,
            col.names=FALSE, row.names=FALSE)
            
## Process the data for Stator
for (this.run in cell.map[, unique(run)]){
    seurat.file <- file.path(root.dir, 
                               sprintf("%s.seurat.for.stator.RDS", this.run))
    if (file.exists(seurat.file)) {
        msg.txt <- sprintf("Read the detected %s.seurat.for.stator.RDS file", this.run)
        msg(warn, msg.txt)
        this.seurat <- readRDS(seurat.file)[[1]]
    } else {
        this.seurat <- all.runs[[this.run]]
    }
    ## extract sparse count matrix, convert to dense matrix
    counts <- GetAssayData(object=this.seurat, slot="counts")
    summary(rowSums(counts))
    summary(colSums(counts))

    counts <- t(counts)
    counts <- as.matrix(counts)

    ## write counts and selected genes to files
    counts.file <- file.path(root.dir, sprintf("%s.counts.csv", this.run))
    write.csv(counts, file=counts.file, append=FALSE, quote=FALSE,
              row.names=TRUE, col.names=TRUE)
    gc()
}

## Read unfiltered clontypes data
clonotypes.raw <- fread(file.path(working.dir, "MECFS001_PBMC_TCRs_all_contig_annotations.csv"), 
                    header=TRUE)
clonotypes.cell <- clonotypes.raw[clonotypes.raw[,is_cell==TRUE]]
clonotypes.filtered <- clonotypes.cell[clonotypes.cell[,raw_clonotype_id!=""]]

add.clonotypes <- TRUE
if(add.clonotypes){
    ## Read clonotypes data filemanes
    all.clonotypes <- list.files(path=file.path(root.dir, "clonotypes/"), pattern=".csv", all.files=TRUE, 
                      full.names=TRUE)
    ## Read all data
    all.list.clonotypes <- lapply(all.clonotypes, read.file)
    all.data.clonotypes <- rbindlist(all.list.clonotypes)
}
## Write metadata into csv
for (this.run in cell.map[, unique(run)]){
    seurat.file <- file.path(root.dir, 
                               sprintf("%s.seurat.for.stator.RDS", this.run))
    if (file.exists(seurat.file)) {
        msg.txt <- sprintf("Read the detected %s.seurat.for.stator.RDS file", this.run)
        msg(warn, msg.txt)
        this.seurat <- readRDS(seurat.file)[[1]]
    } else {
        this.seurat <- all.runs[[this.run]]
    }
    this.dt <- data.table(sample=this.seurat$orig.ident, barcode=colnames(this.seurat))
    this.dt[, celltype := gsub(".*(_All_|_CD8_|_CD8CD4_|_CD4CD8_).*", "\\1", barcode)]
    this.dt[celltype=="_CD8CD4_", celltype := "_CD4CD8_"]
    this.merge <- this.dt[cell.map, on="sample", nomatch=NULL]
    this.export <- this.merge[,.(barcode, condition, celltype)]
    if(add.clonotypes){
        ## Read clontypes data
        metadata.file <- file.path(root.dir, sprintf("%s.metadata.clonotypes.csv", this.run))
        ## Cellmap for each run
        this.cell.map <- cell.map[run==eval(this.run)]
        ## Selected samples
        this.samples <- unique(this.cell.map$sample.simple)
        ## Select corresponding clonotypes and assign sample ids
        selected.clonotypes <- data.table()
        for (sample in this.samples){
            if(this.run != 5){
                ## CD8 cells
                this.pattern <- sprintf("%s_.*_CD8_.*", sample)
                selected.clonotypes <- select.clonotypes(selected.clonotypes, this.cell.map, this.pattern)
                ## CD4CD8 cells
                this.pattern <- sprintf("%s_.*(_CD8CD4_|_CD4CD8_).*", sample)
                selected.clonotypes <- select.clonotypes(selected.clonotypes, this.cell.map, this.pattern)
            } else {
                ## All cells
                this.pattern <- sprintf("%s_.*_All_.*", sample)
                selected.clonotypes <- select.clonotypes(selected.clonotypes, this.cell.map, this.pattern)
            }
        }
    selected.clonotypes <- selected.clonotypes[full.barcode:=paste0()]
    ## Filter out clonotypes with less than 30 umis
    filtered.colontypes <- selected.clonotypes[selected.clonotypes[,umis>30]]
    data.clonotypes[" ":= ]
    } else {
        metadata.file <- file.path(root.dir, sprintf("%s.metadata.csv", this.run))
    }
    setnames(this.export, c(1, 2, 3), c(" ", "Cell.State", "Cell.Types"))
    write.table(this.export, file=metadata.file, sep = ",", quote=FALSE,
              row.names=FALSE, col.names=TRUE)
    }
