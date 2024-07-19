library(data.table)

## install packages for QC
if (FALSE) {
    ## manually set C++ standard (make sure gcc v7+ is loaded as module)
    Sys.setenv("PKG_CXXFLAGS"="-std=c++17")

    ## install in specified user library
    my.lib <- .libPaths()[[1]]
    BiocManager::install("rhdf5", lib=my.lib, ask=FALSE)
    BiocManager::install("DropletUtils", lib=my.lib, ask=FALSE)
    BiocManager::install("Seurat", lib=my.lib, ask=FALSE)
    BiocManager::install("S4Vectors", lib=my.lib, ask=FALSE)
    install.packages("ggExtra", lib=my.lib, quiet=TRUE)
    install.packages("dplyr", lib=my.lib, quiet=TRUE)
    install.packages("ggplot2", lib=my.lib, quiet=TRUE)
    install.packages("xlsx", lib=my.lib, quiet=TRUE)
    install.packages("stringr", lib=my.lib, quiet=TRUE)
    BiocManager::install("sctransform", lib=my.lib, ask=FALSE)
    install.packages("rjson", lib=my.lib, quiet=TRUE)
    BiocManager::install("loomR", lib=my.lib, ask=FALSE)
    BiocManager::install("sctransform", lib=my.lib, ask=FALSE)
    install.packages("randomcoloR", lib=my.lib, quiet=TRUE)
}

root.dir <- "/exports/igmm/eddie/khamseh-lab/aiakovliev"
codebase.dir <- file.path(root.dir, "liver-stator")

## downloads
studies.list.file <- file.path(codebase.dir, "liver_studies.tsv")
studies.dt <- fread(studies.list.file)
