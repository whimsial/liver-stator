## BiocManager::install("DropletUtils")
## load
library(DropletUtils)
library(HDF5Array)

if (FALSE) this.sample <- "GSM4041150" ## for profiling we use only 1 sample
if (FALSE) this.sample <- "GSM5764242" ## for profiling we use only 1 sample

Cell <- NULL
for (this.sample in meta[, unique(sample)]) {
    cat("Analyzing sample", this.sample, "\n")
    this.meta <- meta[sample==eval(this.sample)]
    this.sample.dir <- this.meta[, unique(sample.dir)]
    this.archives <- file.path(study$dir, this.meta$Name)

    if (study$name=="Ramachandran") {
        ## create sample subdirectory and extract sample files
        this.extracts <- file.path(this.sample.dir,
            gsub(".*(barcodes.tsv|genes.tsv|matrix.mtx).gz", "\\1",
                 this.meta$Name))

        dir.create(this.sample.dir, showWarnings=FALSE)
        cmd <- sprintf("gunzip -c %s > %s", this.archives, this.extracts)
        for (this.cmd in cmd) system(this.cmd)

        this.meta[, extracts := this.extracts]

        ## read barcodes
        sample.barcodes.file <- this.meta[grep("barcodes", extracts), extracts]
        sample.barcodes <- fread(sample.barcodes.file, header=FALSE)

        ## read 10x experiment from sample dir
        sce <- read10xCounts(this.meta[, unique(sample.dir)])
    }

    if (study$name=="Guilliams") {
        this.extracts <- file.path(this.sample.dir, this.meta$Name)
        dir.create(this.sample.dir, showWarnings=FALSE)
        file.symlink(this.archives, this.extracts)
        this.meta[, extracts := this.extracts]

        ## read barcodes
        barcode.path <- "matrix/barcodes"
        sample.barcodes <- h5read(this.extracts, barcode.path)
        sample.barcodes <- data.table(V1=sample.barcodes)

        ## read 10x experiment from sample dir
        sce <- read10xCounts(this.extracts)
    }

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
        return(NA)  ## Return NA on error
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
