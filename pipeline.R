## High-resolution mapping of cell states in liver disease
source("properties.R")
source("helperfunctions.R")
source("rnaseq.functions.R")

## Process Ramachandran et al (2019) (see full list in liver_studies.tsv)
## -----------------------------------------------------------------------------
study.name <- studies.dt[1,
                         gsub("^(\\w+).*", "\\1", Study)]
study.url <- studies.dt[grep(study.name, Study), `Data URL`]
filelist.url <- studies.dt[grep(study.name, Study),
                           gsub("(.*suppl).*", "\\1", `Data URL`)]

study.meta <- download.study(study.name, study.url, filelist.url, root.dir)

## study variable should be a list with at least 2 elements:
## $dir -> name of the directory containing extract from the GEO archive
## $study.file.list -> content of filelist.txt which contains a list of all
##                     files for the study (downloaded separately from GEO)
sample.files <- list.files(study$dir, pattern="*.tsv.gz", full.names=TRUE)
meta <- fread(study$study.file.list)
stopifnot(all(basename(sample.files) %in% meta$Name))

## create dir path for each sample
meta[, sample := gsub("(GSM\\d+).*", "\\1", Name)]
meta[, sample.dir := file.path(study$dir, sample)]
meta <- meta[!grep(".tar", sample)]
meta <- meta[!grep("blood|mouse", Name)]
meta.file <- file.path(study$dir, "metadata.txt")
fwrite(meta, file=meta.file)

## extract sample files and perform QC step 1: detection of empty drops
source("qc/qc1_emptydrop.R")

## run QC_2RMdoublet.ipynb in Jupyter notebook: detection of duplets

## run QC step 3
seurat.study1 <- process.samples.and.merge(meta)

## run QC step 4 to filter out genes/cells that do not pass the thresholds
seurat.study1.filtered <- qc.seurat(seurat.study1, output.dir=study$dir)

## load core genes for NAFLD
core.eqtls <- load.rdata(file.path(root.dir, "nafld.diag.core.genes.eqtls.RData"))
core.pqtls <- load.rdata(file.path(root.dir, "nafld.diag.core.genes.pqtls.RData"))

core.genes <- c(core.eqtls[pvalue_trans<1E-6, gene_symbol],
                core.pqtls[pvalue_trans<1E-6 & eff.numtransqtls>=5, gene_symbol])

## run QC step 5 to select most variable genes and prepare inputs for Stator
label.genes <- c("ADIPOQ", "ADIPOR1", "ADIPOR2")
seurat.study1.filtered.hvg <- process.variable.genes(seurat.study1.filtered,
                                                     core.genes=core.genes,
                                                     label.genes=label.genes,
                                                     output.dir=study$dir)
hvg <- seurat.object@assays$RNA@var.features

## append ADIPOR1 to HVGs
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
label.genes <- get.gene.ids(data.table(gene.symbol=label.genes), ensembl)
ensembl.adipor1 <- label.genes[external_gene_name=="ADIPOR1", ensembl_gene_id]
hvg <- c(hvg, ensembl.adipor1)
seurat.for.stator.study1 <- filter.seurat(seurat.study1.filtered.hvg,
                                          n.cells.keep=4000)
seurat.for.stator.file <- file.path(study$dir, "seurat.for.stator.Rdata.gz")
save(seurat.for.stator.study1, file=seurat.for.stator.file)

## extract sparse count matrix, convert to dense matrix and save
count <- seurat.for.stator.study1[["RNA"]]@data
count <- as.matrix(count)
count <- t(count)

## write counts and selected genes to the file
counts.file <- file.path(study$dir, "counts.csv")
write.csv(count, file=counts.file, append=F, quote=F, row.names=T, col.names=T)
genes.file <- file.path(study$dir, "genes.csv")
write.table(t(hvg), file=genes.file, sep=",", quote=FALSE, col.names=FALSE,
            row.names=FALSE)

## sopy counts and genes files to Eddie
cmd <- sprintf("rsync -chavzP %s eddie:/exports/eddie/scratch/aiakvlie/",
               counts.file)
system(cmd)
cmd <- sprintf("rsync -chavzP %s eddie:/exports/eddie/scratch/aiakvlie/",
               genes.file)
system(cmd)

## Process Guilliams et al (2022)
## -----------------------------------------------------------------------------
this.study <- 2 ## Guilliams et al (2022)
study <- download.study(this.study, studies.dt, root.dir)
## download soft meta data to map samples
soft.url <- gsub("(suppl.*)", "soft/GSE192742_family.soft.gz",
                 studies.dt[this.study, `Data URL`])
soft.file <- file.path(study$dir, "samples.soft")
sample.map <- parse.soft(soft.url, soft.file)
sample.map <- sample.map[grep("(Whole Liver)?(Human)", sample.title)]
sample.map[, sample.id := gsub(".*(H\\d+).*", "\\1", sample.title)]

## sample ids for healthy, steatosis and fibrosis
healthy <- "H30"
steatosis <- c("H32", "H35", "H37") ## steatosis > 15%
fibrosis <- c("H33", "H36", "H37", "H38") ## any mention of fibrosis
selected.samples <- sample.map[sample.id %in% c(healthy, steatosis, fibrosis),
                               sample]

## extract archive
cmd <- sprintf("tar -xvf %s -C %s", study$study.file, study$dir)
system(cmd)

## study variable should be a list with at least 2 elements:
## $dir -> name of the directory containing extract from the GEO archive
## $study.file.list -> content of filelist.txt which contains a list of all
##                     files for the study (downloaded separately from GEO)
sample.files <- list.files(study$dir, pattern="*.txt.gz", full.names=TRUE)
meta <- fread(study$study.file.list)
stopifnot(all(basename(sample.files) %in% meta$Name))

## create dir path for each sample
meta[, sample := gsub("(GSM\\d+).*", "\\1", Name)]
meta[, sample.dir := file.path(study$dir, sample)]
meta <- meta[!grep(".tar", sample)]
meta <- meta[Type=="H5"]
## keep only files for relevant samples
meta <- meta[sample %in% selected.samples]
meta.study2 <- meta
fwrite(meta, file=file.path(study$dir, "metadata.txt"))

## extract sample files and perform QC step 1
source("qc/qc1_emptydrop.R")
