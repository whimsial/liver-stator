## BiocManager::install("DropletUtils")
## load
library(DropletUtils)

prepare.10x.files <- function(this.meta, study) {
    ## create sample subdirectory and extract sample files
    this.sample.dir <- this.meta[, unique(sample.dir)]
    this.archives <- file.path(study$dir, this.meta$Name)
    this.extracts <- file.path(this.sample.dir,
        gsub(".*(barcodes.tsv|genes.tsv|matrix.mtx).gz", "\\1", this.meta$Name))

    dir.create(this.sample.dir, showWarnings=FALSE)
    cmd <- sprintf("gunzip -c %s > %s", this.archives, this.extracts)
    for (this.cmd in cmd) system(this.cmd)

    this.meta[, extracts := this.extracts]
    return(this.meta)
}

if (FALSE) this.sample <- "GSM4041150" ## for profiling we use only 1 sample

Cell <- NULL
for (this.sample in meta[, unique(sample)]) {
    cat("Analyzing sample", this.sample, "\n")
    this.meta <- meta[sample==eval(this.sample)]

    ## create 1 subdirectory per sample and relevant symlinks
    this.meta <- prepare.10x.files(this.meta, study)

    sample.barcodes.file <- this.meta[grep("barcodes", extracts), extracts]
    sample.barcodes <- fread(sample.barcodes.file, header=FALSE)

    ## read 10x experiment from sample dir
    sce <- read10xCounts(this.meta[, unique(sample.dir)])

    ## run emptyDrops
    ## see https://support.bioconductor.org/p/123554/#123562 in case the matrix
    ## has been already filtered. If it has, this step will be omitted!
    set.seed(100)
    my.count <- counts(sce)
    colnames(my.count) <- paste0(sample, ":", sample.barcodes$V1)
    colnames(my.count) <- gsub(pattern='-1', replacement="x",
                               colnames(my.count))

    cat("Starting emptyDrops\n")
    e.out <- tryCatch(emptyDrops(my.count), error = function(e) {
        warning("The matrix must have been pre-filtered already.")
        return(NA)  # Return NA on error
    })

    if (is.na(e.out)) True_Cell <- colnames(my.count)
    else {
        is.cell <- e.out$FDR <= 0.01
        is.cell[is.na(is.cell)] <- FALSE
        True_Cell <- colnames(my.count)[is.cell]
    }

    Cell <- c(Cell, True_Cell)
    cat("Done.\n")
}

save(Cell, file=file.path(study$dir, "emptydrop.Rdata"))
