## Example R script to download, process and QC scRNA-seq data from GEO.
## This simplified script assumes that the data comes from a single study
## and there are no major inconsistencies between sample in terms of
## experimental protocols. If you need to perform data integration from studies
## with different protocols, refer to `pipeline.R`.
## -----------------------------------------------------------------------------

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


## Preparation step.
## Set GEO ID, URLs, and study directories and file names.
#'
#' @param download.studies set downloading from GEO to TRUE/FALSE 
#' @param root.dir Full path to large storage where the project data is to be
#'        stored.
#' @param working.dir Full path to the working directory (typically the cloned
#'        repository).
#' @param study.name String specifying a descriptive study name. This will be
#'        used to create study subdirectory in root.dir.
#' @param geo.id Unique ID of the study in GEO database.
#' @param data.url URL for the raw study data on GEO FTP server.
#' @param filelist.url URL pointing to the file on GEO FTP server which contains
#'        list of raw data files.
#' @param softfile.url URL pointing to the file on GEO FTP server which contains
#'        sample metadata for the study.
## -----------------------------------------------------------------------------
download.studies <- FALSE
root.dir <- "/exports/igmm/eddie/ponting-lab/ava/ME_CSF_Hifibio/raw_data/"
working.dir <- "/exports/igmm/eddie/ponting-lab/sbraich2/liver-stator/"
# ## Load helper functions
source(file.path(working.dir, "rnaseq.functions.R"))

## (Skip if you have data already available by setting download.studies <- FALSE)
if(download.studies){
    ## Here we keep example from liver Ramachandran et al. 2019
    ## https://www.nature.com/articles/s41586-019-1631-3
    study.name <- "Ramachandran"
    geo.id <- "GSE136103"
    data.url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/GSE136103_RAW.tar"
    filelist.url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/filelist.txt"
    softfile.url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/soft/GSE136103_family.soft.gz"
    study.dir <- file.path(root.dir, study.name)
    study.file <- file.path(study.dir, basename(data.url))
    filelist.file <- file.path(study.dir, basename(filelist.url))
    soft.file <- file.path(study.dir, basename(softfile.url))
    
    download.study(study.name, study.file, data.url, filelist.file, filelist.url)
    
    ## GEO step. Download GEO dataset, extract and process files
    ## -----------------------------------------------------------------------------
    ## download specified GEO dataset
    ## read and populate metadata
    meta.dt <- fread(filelist.file)
    meta.dt[, study.name := study.name]
    meta.dt <- meta.dt[, .(study=eval(study.name),
                           geo.id=eval(geo.id),
                           data.url=eval(data.url),
                           filelist.url=eval(filelist.url),
                           softfile.url=eval(softfile.url),
                           study.dir=eval(study.dir),
                           study.file=eval(study.file),
                           filelist.file=eval(filelist.file),
                           soft.file=eval(soft.file),
                           sample.file=Name)]
    meta.dt[, sample := gsub("(GSM\\d+).*", "\\1", sample.file)]
    meta.dt[, sample.dir := file.path(study.dir, sample)]
    meta.dt <- meta.dt[!grep(".tar", sample)]

    ## assign file types (mtx and h5 files are supported), if `sample.file` column
    ## of `meta.dt` contains files with .mtx extension then use "mtx". Alternatively
    ## use "h5"
    meta.dt[, file.type := "mtx"]

    ## for this example, we simplify things by using only 1 sample from the
    ## liver atlas
    meta.dt <- meta.dt[sample=="GSM4041150"]

    ## save metadata to a file
    meta.file <- file.path(root.dir, "metadata.txt")
    fwrite(meta.dt, file=meta.file, sep="\t")

    ## Extract files into respective subdirectories
    msg(bold, "Extracting downloaded archive")
    cmd <- sprintf("tar -xvf %s -C %s", study.file, study.dir)
    system(cmd)

    msg(bold, "Extracting sample files into subdirectories")
    for (idx in seq_len(nrow(meta.dt))) {
        this.sample <- meta.dt[idx, sample]
        msg.txt <- sprintf("Extracting sample %s", this.sample)
        msg(info, msg.txt)
        this.sample.dir <- meta.dt[sample==eval(this.sample), unique(sample.dir)]
        this.archives <- meta.dt[sample==eval(this.sample),
                                 file.path(unique(study.dir), sample.file)]
                                 this.extracts <- extract.samples(this.archives,
                                                                  this.sample.dir)
        meta.dt[sample==eval(this.sample), extract := this.extracts]
    }
    
} else {
    ## Select a sample of interest manually (.mtx,.tsv)
    mtx.files <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
    this.extracts <- file.path(root.dir, 
        "PLRC_HD055V2_fresh_PBMC_CD8_10xscGEX/outs/raw_feature_bc_matrix",
        mtx.files)
    this.sample.dir <- file.path(root.dir, 
        "PLRC_HD055V2_fresh_PBMC_CD8_10xscGEX/outs/raw_feature_bc_matrix")
        # ## h5 raw_feature_bc_matrix.h5
        # this.sample.dir <- file.path(root.dir,
        #     "PLRC_HD055V2_fresh_PBMC_CD8_10xscGEX/outs/")
        # this.extracts <- file.path(root.dir,
        #     "PLRC_HD055V2_fresh_PBMC_CD8_10xscGEX/outs/filtered_feature_bc_matrix.h5")
}

## STEP 1. Read 10X data from sample files.
#'
#' @param mtx.files list of a triplet of file names usually
#'        `barcodes.tsv`, `genes.tsv` (here `features.tsv`) and `matrix.mtx`
#' @param this.extracts A list of full paths to mtx.files  which 
#'         contain counts for a given sample.
#' @param this.sample.dir Full path to the directory containing counts data for
#'         a given sample.
## -----------------------------------------------------------------------------
data <- read.10x.data(this.extracts, this.sample.dir, mtx.files)

## STEP 2. Detect empty drops in scRNA-seq data.
#'
#' @param data 10x dataset for a single sample loaded using read.10x.data
#'        function.
#' @param this.sample Sample ID whose counts data is in the `data` variable.
## -----------------------------------------------------------------------------
msg(bold, "Removing empty drops")
this.sample <- "HD055V2"
## Save barcodes of true cells
true.cells <- remove.emptydrops(sce=data$sce, sample=this.sample,
                                sample.barcodes=data$sample.barcodes$V1)
                                
## Remove sample names from the barcodes
true.cells.ss <- remove.sample.from.barcode(true.cells)

## Filter cells
sce_filtered <- data$sce[, data$sce$Barcode %in% true.cells.ss]

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
msg(bold, "Detecting doublets")
## The flag to select the method, scDblFinder is default
doublet.method.scDblFinder <- FALSE
if (doublet.method.scDblFinder){
    ## scDblFinder impelentation (DEFAULT)
    find.doublets.scDblFinder(sce_filtered, this.sample.dir)
} else {
    ## scrublet impelentation
    doublet.threshold <- 0.15
    this.file.type <- meta.dt[sample==eval(this.sample), unique(file.type)]
    
    ## TODO: run srcublet on true cells only, now it runs on all cells
    ## run python script `doublet.py` from R using system command
    msg(bold, "Detecting doublets")
    cmd <- sprintf("python3 doublet.py --sample_dir %s --data_type %s \\
                    --doublet_threshold %f", this.sample.dir, this.file.type,
                    doublet.threshold)
    system(cmd)
}

## STEP 4. Convert counts data for all sample to Seurat objects and then merge
## those objects together.
#'
#' @param metadata Data table with columns `sample` and `sample.dir` which are
#'        contain sample ID and full path to the sample directory. If the
#'        preparation step from this script was run, then `meta.dt` variable
#'        can be supplied. Alternatively, create this variable manually as shown
#'        below.
#' @param root.dir Full path to large storage where the project data is to be
#'        stored.
## -----------------------------------------------------------------------------
metadata <- meta.dt[, .(sample, sample.dir)]
## TEST: create metadata with
## metadata <- data.table(sample=this.sample, sample.dir=this.sample.dir)
msg(bold, "Reading single cell data to Seurat and merging")
seurat.all <- process.samples.and.merge(metadata, true.cells, output.dir=root.dir)

# Remove sample id from barcodes
seurat.all@meta.data$barcodes <- remove.sample.from.barcode(seurat.all@meta.data$barcodes)

## STEP 5. Filter out genes/cells that do not pass QC thresholds
#'
#' @param seurat.all Seurat object from the previous step.
#' @param root.dir Full path to large storage where the project data is to be
#'        stored.
## -----------------------------------------------------------------------------
msg(bold, "Performing QC of a merged Seurat object")
seurat.qc <- qc.seurat(seurat.all, output.dir=root.dir)

## STEP 6: Sample XXX cells from the Seurat object after the QC
#'
#' @param seurat.qc Seurat object from the previous step.
#' @param cell.map Data table containing a mapping between samples and
#'        conditions. It should have the following columns: `sample` (sample
#'        identifier), `condition` (e.g., "healthy", "disease")
#' @param n.cells.keep Integer specifying the number of cells to be sampled.
## -----------------------------------------------------------------------------
cell.map <- data.table(sample="GSM4041150", condition="healthy")
n.cells.keep <- 5000

msg(bold, "Sampling cells")
seurat.for.stator <- sample.seurat(seurat.qc, output.dir=root.dir,
                                   cell.map=cell.map, n.cells.keep=n.cells.keep)

## repeat the QC steps on a randomly selected subset of cells
seurat.for.stator <- qc.seurat(seurat.for.stator, output.dir=root.dir)

## STEP 7: Prepare data for Stator
#'
#' @param seurat.for.stator Seurat object from the previous step.
## -----------------------------------------------------------------------------
seurat.for.stator <- process.variable.genes(seurat.for.stator,
                                            output.dir=root.dir)
hvg <- seurat.for.stator@assays$RNA@var.features

## extract sparse count matrix, convert to dense matrix
counts <- GetAssayData(object=seurat.for.stator, slot="counts")
summary(rowSums(counts))
summary(colSums(counts))

counts <- t(counts)
counts <- as.matrix(counts)

## write counts and selected genes to files
counts.file <- file.path(root.dir, "counts.csv")
write.csv(counts, file=counts.file, append=FALSE, quote=FALSE,
          row.names=TRUE, col.names=TRUE)
genes.file <- file.path(root.dir, "genes.csv")
write.table(hvg, file=genes.file, sep=",", quote=FALSE,
            col.names=FALSE, row.names=FALSE)
