library(data.table)

## install packages for QC
## these packages are only installable with gmm/apps/R/4.0.2.gcc.9.4.0
if (FALSE) {
    ## install in specified user library
    my.lib <- .libPaths()[[1]]
    install.packages(c("ggExtra", "dplyr", "ggplot2", "xlsx", "stringr",
                       "rjson", "randomcoloR"), lib=my.lib, quiet=TRUE)
    BiocManager::install("rhdf5", lib=my.lib, ask=FALSE)
    BiocManager::install(c("DropletUtils", "Seurat", "S4Vectors"),
                         lib=my.lib, ask=FALSE)
    ## loomR has to be installed from GitHub
    devtools::install_github(repo="hhoeflin/hdf5r")
    devtools::install_github(repo="mojaveazure/loomR", ref="develop")
    BiocManager::install("sctransform", lib=my.lib, ask=FALSE)
}

# root.dir <- "/exports/igmm/eddie/khamseh-lab/aiakovliev"
root.dir <- "/opt/datastore/aiakovliev/liver"
# codebase.dir <- file.path(root.dir, "liver-stator")
codebase.dir <- "~/repos/liver-stator"

## downloads
studies.list.file <- file.path(codebase.dir, "liver_studies.tsv")
studies.dt <- fread(studies.list.file)
