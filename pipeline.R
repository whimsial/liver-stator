## High-resolution mapping of cell states in liver disease
source("properties.R")
source("rnaseq.functions.R")

#' Helper function to load a given Rdata object.
load.rdata <- function(filename) {
    load(filename)
    get(ls()[ls() != "filename"])
}

## Prepare study metadata table which we will use to provide URLs to download
## relevant files (raw metadata is located in liver_studies.tsv)
## -----------------------------------------------------------------------------
## set GEO id
studies.dt <- studies.dt[1:2]
studies.dt[, `GEO ID` := gsub("(.*)_RAW.tar", "\\1",
                              basename(`Data URL`))]

## remove et al and the year from study name
studies.dt[, Study := gsub("^(\\w+).*", "\\1", Study)]

## create URL for the file contining list of files within the study
studies.dt[, `Filelist URL` :=
    file.path(gsub("(.*suppl).*", "\\1", `Data URL`), "filelist.txt")]

## create URL for soft files which contain sample information;
## this file URL can be worked our from the Data URL by walking backwards and
## navigating to `soft` subdirectory
studies.dt[, `Softfile URL` :=
    file.path(dirname(dirname(`Data URL`)), "soft",
              sprintf("%s_family.soft.gz", `GEO ID`))]

## create system paths to where study files are to be downloaded
studies.dt[, study.dir := file.path(root.dir, Study)]
studies.dt[, study.file := file.path(study.dir, basename(`Data URL`))]
studies.dt[, filelist.file := file.path(study.dir, basename(`Filelist URL`))]
studies.dt[, soft.file := file.path(study.dir, basename(`Softfile URL`))]

## Downlaod Ramachandran et al (2019) and Guilliams et al (2022)
## -----------------------------------------------------------------------------
studies <- c(1, 2)
meta.dt <- data.table()
for (this.study in studies) {
    this.meta <- studies.dt[eval(this.study)]
    study.name <- this.meta[, Study]
    study.url <- this.meta[, `Data URL`]
    study.file <- this.meta[, study.file]
    filelist.url <- this.meta[, `Filelist URL`]
    filelist.file <- this.meta[, filelist.file]
    download.study(study.name, study.file, study.url,
                   filelist.file, filelist.url)
    ## read metadata
    meta <- fread(filelist.file)
    meta[, study.name := study.name]
    meta.dt <- rbind(meta.dt, meta)
}

## for Guilliams (2022) we get soft metadata to map samples to disease stages
soft.url <- studies.dt[2, `Softfile URL`]
soft.file <- studies.dt[2, soft.file]
sample.map <- parse.soft(soft.url, soft.file)
sample.map <- sample.map[grep("(Whole Liver)?(Human)", sample.title)]
sample.map[, sample.id := gsub(".*(H\\d+).*", "\\1", sample.title)]

## manually set sample ids for healthy, steatotic and fibrotic livers by
## inspecting sample protocols from the paper
healthy <- "H30"
steatosis <- c("H32", "H35", "H37") ## steatosis > 15%
fibrosis <- c("H33", "H36", "H37", "H38") ## any mention of fibrosis
selected.samples <- sample.map[sample.id %in% c(healthy, steatosis, fibrosis),
                               sample]
## extract nested RAW archive
cmd <- sprintf("tar -xvf %s -C %s", studies.dt[2, study.file],
               studies.dt[2, study.dir])
system(cmd)

## Create data table with metadata for the studies
## -----------------------------------------------------------------------------
sample.files <- list.files(studies.dt[2, study.dir], pattern="*.tsv.gz|*.h5",
                           full.names=TRUE)
if (!all(basename(sample.files) %in% meta$Name)) {
    stop("Files listed in study list were not found in study directory.")
}
meta.dt <- studies.dt[meta.dt, on=c(Study="study.name")]

meta.dt <- meta.dt[, .(study=Study, geo.id=`GEO ID`, study.url=URL,
                       data.url=`Data URL`, filelist.url=`Filelist URL`,
                       softfile.url=`Softfile URL`,
                       study.dir, study.file, filelist.file, soft.file,
                       sample.file=Name)]

meta.dt[, sample := gsub("(GSM\\d+).*", "\\1", sample.file)]
meta.dt[, sample.dir := file.path(study.dir, sample)]

## keep only samples of interest
meta.dt <- meta.dt[!grep(".tar", sample)]
meta.dt <- meta.dt[!grep("blood|mouse", sample.file)]
meta.dt <- meta.dt[study=="Ramachandran" |
                   (study=="Guilliams" & sample %in% selected.samples)]
## filter out other irrelevant files
meta.dt <- meta.dt[!grep("json|png|csv", sample.file)]

## assign file types
meta.dt[study=="Ramachandran", file.type := "mtx"]
meta.dt[study=="Guilliams", file.type := "h5"]

## save metadata to a file
meta.file <- file.path(root.dir, "metadata.txt")
fwrite(meta.dt, file=meta.file, sep="\t")

## Etract sample files into respective subdirectories and run per-sample
## QC steps 1 and 2: detection of empty drops and detection of doublets
## -----------------------------------------------------------------------------
msg(bold, "Extracting sample files into subdirectories")
msg.txt <- sprintf("Extracting sample %s", this.sample)
msg(info, msg.txt)

all.10x.data <- list()
cells <- NULL
for (idx in seq_len(nrow(meta.dt))) {
    this.sample <- meta.dt[idx, sample]
    this.sample.dir <- meta.dt[sample==eval(this.sample), unique(sample.dir)]
    this.archives <- meta.dt[sample==eval(this.sample),
                             file.path(unique(study.dir), sample.file)]

    ## extract sample data from archives
    this.extracts <- extract.samples(this.archives, this.sample.dir)
    meta.dt[sample==eval(this.sample), extract := this.extracts]
    ## read 10X data from extracted archives
    data <- read.10x.data(this.extracts, this.sample.dir)
    all.10x.data[[this.sample]] <- data

    ## perform QC step 1: detection of empty drops
    true.cells <- remove.emptydrops(sce=data$sce, sample=this.sample,
                                    sample.barcodes=data$sample.barcodes$V1)
    cells <- c(cells, true.cells)

    ## perform QC step 2: detection of doublets
    ## first we run jupyter notebook manually to work out suitable threshold
    ## and then run the python script for all samples
    this.file.type <- meta.dt[sample==eval(this.sample), unique(file.type)]
    cmd <- sprintf("python3 doublet.py --sample_dir %s --data_type %s \\
                   --doublet_threshold %f", this.sample.dir, this.file.type,
                   0.15)
    system(cmd)
}

## QC step 3: process all samples, convert to Seurat and merge
## -----------------------------------------------------------------------------
seurat.all <- process.samples.and.merge(meta.dt, output.dir=root.dir)

## QC step 4: filter out genes/cells that do not pass QC thresholds
## -----------------------------------------------------------------------------
seurat.all.filtered <- qc.seurat(seurat.all,
                                 all.samples=meta.dt[, unique(sample)],
                                 output.dir=root.dir)

## QC step 5: bring core genes to HVGs and prepare inputs for Stator
## -----------------------------------------------------------------------------
## load core genes for NAFLD
core.eqtls <- load.rdata(file.path(root.dir,
                                   "nafld.diag.core.genes.eqtls.RData"))
core.pqtls <- load.rdata(file.path(root.dir,
                                   "nafld.diag.core.genes.pqtls.RData"))

core.genes <- c(core.eqtls[pvalue_trans<1E-6, gene_symbol],
                core.pqtls[pvalue_trans<1E-6 & eff.numtransqtls>=5,
                           gene_symbol])

## run QC step 5 to select most variable genes and prepare inputs for Stator
label.genes <- c("ADIPOQ", "ADIPOR1", "ADIPOR2")
seurat.all.filtered.hvg <- process.variable.genes(seurat.all.filtered,
                                                  core.genes=core.genes,
                                                  label.genes=label.genes,
                                                  output.dir=root.dir)
hvg <- seurat.all.filtered.hvg@assays$RNA@var.features

## append ADIPOR1 to HVGs
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
adipor1 <- data.table(gene.id=NULL, gene.symbol="ADIPOR1")
adipor1 <- get.gene.ids(adipor1, ensembl)
hvg <- c(hvg, adipor1$ensembl_gene_id)

## sample 20k cells randomly from all donors
seurat.for.stator <- filter.seurat(seurat.all.filtered.hvg, output.dir=root.dir,
                                   n.cells.keep=20000)

## extract sparse count matrix, convert to dense matrix
count <- seurat.for.stator[["RNA"]]@data
count <- as.matrix(count)
count <- t(count)

## write counts and selected genes to the file
counts.file <- file.path(root.dir, "counts.csv")
write.csv(count, file=counts.file, append=FALSE, quote=FALSE,
          row.names=TRUE, col.names=TRUE)
genes.file <- file.path(root.dir, "genes.csv")
write.table(t(hvg), file=genes.file, sep=",", quote=FALSE, col.names=FALSE,
            row.names=FALSE)

## copy counts and genes files to Eddie
cmd <- sprintf("rsync -chavzP %s eddie:/exports/eddie/scratch/aiakvlie/",
               counts.file)
system(cmd)
cmd <- sprintf("rsync -chavzP %s eddie:/exports/eddie/scratch/aiakvlie/",
               genes.file)
system(cmd)

## Create and save cell annotation data for (Meta_data.csv) which maps cells
## to disease states
## -----------------------------------------------------------------------------
## map samples to diagnoses
meta.dt[sample %in% sample.map[sample.id %in% healthy, sample],
        condition := "healthy"]
meta.dt[grep("healthy", sample.file), condition := "healthy"]
meta.dt[sample %in% sample.map[sample.id %in% steatosis, sample],
        condition := "steatotic"]
meta.dt[sample %in% sample.map[sample.id %in% fibrosis, sample],
        condition := "fibrotic"]
meta.dt[grep("cirrhotic", sample.file), condition := "cirrhotic"]

## map cells to samples
cell.map <- data.table(barcode=rownames(count))
cell.map[, sample := gsub("(\\w+)_[A-Z]+.*", "\\1", barcode)]
cell.map <- cell.map[unique(meta.dt, by="sample")[, .(sample, condition)],
                     on="sample"]

## map Ensembl gene ids to gene symbols
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensembl.ids <- data.table(gene.id=colnames(count))
gene.symbols <- get.gene.symbols(ensembl.ids, ensembl)

## classify cells into cell types using celltype expression map
celltype.exp.file <- file.path(root.dir, "output", "celltypes_map.csv")
celltype.exp.map <- fread(celltype.exp.file, select=c(1, 3:9),
                          col.names=c("celltype", sprintf("gene%d", seq(7))))
gene.symbols[celltype.exp.map, on=c(external_gene_name="gene1")]
celltype.exp.map <- transpose(celltype.exp.map, make.names=1)

count.long <- melt(as.data.table(count, keep.rownames=TRUE),
                   measure.vars=names(count), id.vars="barcode",
                   variable.name="gene.id", value.name="value")
count.long <- count.long[value!=0]
count.long[gene.symbols, on=c(gene.id="ensembl_gene_id"),
           gene.symbol:=external_gene_name]

label.celltype <- function(count.long, celltype.exp.map, label) {
    markers <- celltype.exp.map[get(label)!="", get(label)]
    cutoff <- ceiling(length(markers)/2)
    re <- paste(markers, collapse="|")
    count.long[grep(re, gene.symbol),
               `:=`(candidate=eval(label),
                    labels.expressed=paste(gene.symbol, collapse=","),
                    n.expressed=.N), by="barcode"]
    count.long[is.na(celltype) & !is.na(candidate) & n.expressed >= cutoff,
               celltype := candidate]
    return(count.long)
}

## label MPs
count.long[, celltype := NA]
count.long <- label.celltype(count.long, celltype.exp.map, label="MP")
## label pDC
count.long <- label.celltype(count.long, celltype.exp.map, label="pDC")
## label ILC
count.long <- label.celltype(count.long, celltype.exp.map, label="ILC")
## label T cells
count.long <- label.celltype(count.long, celltype.exp.map, label="T cell")
## label B cells
count.long <- label.celltype(count.long, celltype.exp.map, label="B cell")
## label Plasma cells
count.long <- label.celltype(count.long, celltype.exp.map, label="Plasma cell")
## label Mast cells
count.long <- label.celltype(count.long, celltype.exp.map, label="Mast cell")
## label Endothelial cells
count.long <- label.celltype(count.long, celltype.exp.map, label="Endothelia")
## label Mesenchyme cells
count.long <- label.celltype(count.long, celltype.exp.map, label="Mesenchyme")
## label Hepatocytes cells
count.long <- label.celltype(count.long, celltype.exp.map, label="Hepatocyte")
## label Cholangiocyte cells
count.long <- label.celltype(count.long, celltype.exp.map,
                             label="Cholangiocyte")
## label Cycling cells
count.long <- label.celltype(count.long, celltype.exp.map, label="Cycling")
## label the rest as other

count.long[, `:=`(candidate=NULL, n.expressed=NULL, labels.expressed=NULL)]

## annotate cell ids in the counts matrix by the cell types
cell.map[count.long[!is.na(celltype), .(barcode, celltype)], on="barcode",
         celltype:=celltype]
cell.map[is.na(celltype), celltype := "other"]
