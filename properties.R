library(data.table)

## install packages for QC
## these packages are only installable with gmm/apps/R/4.0.2.gcc.9.4.0
if (FALSE) {
    ## install in specified user library
    my.lib <- .libPaths()[[1]]
    BiocManager::install("S4Arrays", version=3.12, lib=my.lib, ask=FALSE)
    BiocManager::install("S4Vectors", lib=my.lib, ask=FALSE)
    BiocManager::install("DropletUtils", lib=my.lib, ask=FALSE)
    BiocManager::install("Seurat", lib=my.lib, ask=FALSE)
    install.packages(c("ggExtra", "dplyr", "ggplot2", "xlsx", "stringr",
                       "rjson", "randomcoloR"), lib=my.lib, quiet=TRUE)
    BiocManager::install("sctransform", lib=my.lib, ask=FALSE)
    BiocManager::install("loomR", lib=my.lib, ask=FALSE)
    BiocManager::install("sctransform", lib=my.lib, ask=FALSE)
}

root.dir <- "/exports/igmm/eddie/khamseh-lab/aiakovliev"
codebase.dir <- file.path(root.dir, "liver-stator")

## downloads
studies.list.file <- file.path(codebase.dir, "liver_studies.tsv")
studies.dt <- fread(studies.list.file)
