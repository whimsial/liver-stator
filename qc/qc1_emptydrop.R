## BiocManager::install("DropletUtils")
## load
library(DropletUtils)

prepare.10x.files <- function(this.meta, study) {
    ## create sample subdirectory and symlink relevant files
    this.sample.dir <- this.meta[, unique(sample.dir)]
    sample.barcodes.file <- file.path(this.sample.dir, "barcodes.tsv")
    sample.genes.file <- file.path(this.sample.dir, "genes.tsv")
    sample.matrix.file <- file.path(this.sample.dir, "matrix.mtx")

    dir.create(this.sample.dir, showWarnings=FALSE)
    if (!any(file.exists(this.meta[, file.path(this.sample.dir, Name)]))) {
        file.symlink(this.meta[, file.path(study$dir, Name)],
                     this.meta[, file.path(this.sample.dir, Name)])
        file.rename(this.meta[grep("barcodes", Name), file.path(sample.dir, Name)],
                    sample.barcodes.file)
        file.rename(this.meta[grep("genes", Name), file.path(sample.dir, Name)],
                    sample.genes.file)
        file.rename(this.meta[grep("matrix", Name), file.path(sample.dir, Name)],
                    sample.matrix.file)
    }

    symlinks <- list(sample.barcodes.file=sample.barcodes.file,
                     sample.genes.file=sample.genes.file,
                     sample.matrix.file=sample.matrix.file)
    return(symlinks)
}

if (FALSE) this.sample <- "GSM4041150" ## for profiling we use only 1 sample

Cell <- NULL
for (this.sample in meta[, unique(sample)]) {
    cat("Analyzing sample", this.sample, "\n")
    this.meta <- meta[sample==eval(this.sample)]

    ## create 1 subdirectory per sample and relevant symlinks
    simlinks <- prepare.10x.files(this.meta, study)

    sample.barcodes <- fread(simlinks$sample.barcodes.file, header=FALSE)

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
