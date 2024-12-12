install_if_not_installed <- function(package_names, lib) {
    for (package_name in package_names) {
        # Attempt to load the package
        if (!require(package_name, character.only=TRUE)) {
            # Package is not installed, install it
            message(paste("Installing package:", package_name))
            install.packages(package_name, lib)
        }

        # Try loading the package again to confirm it works
        load_result <- tryCatch({
            library(package_name, character.only=TRUE)
            TRUE  # Return TRUE if loading succeeds
        }, error=function(e) {
            FALSE  # Return FALSE if there is an error
        })

        # Check if loading the package was successful
        if (!load_result) {
            # If loading fails, reinstall and try again
            message(paste("Reinstallation needed for package:", package_name))
            install.packages(package_name, lib)

            # Final attempt to load the package
            final_result <- tryCatch({
                library(package_name, character.only=TRUE)
                TRUE
            }, error = function(e) {
                FALSE
            })

            if (!final_result) {
                stop(paste("Failed to properly load package:", package_name))
            }
        } else {
            message(paste("Package loaded successfully:", package_name))
        }
    }
}

## Problem: can't load R4.4
## Add to ~/.bashrc:
##. /etc/profile.d/modules.sh
## source /exports/applications/support/set_qlogin_environment.sh
## module load igmm/apps/R/4.0.2.gcc.9.4.0

## Start new session, then
## module load roslin/R/4.4.0

## Path to my R library (Lana)
# my.lib <- "/exports/igmm/eddie/ponting-lab/sbraich2/Rlibrary"

## When R suggests to load packages locally do this to install packages in the store .libPaths( c(my.lib, .libPaths()))
## my.lib <- .libPaths()[[1]]

## install devtools. If installation throws an error, restart R and install
## htmltools with install.packages("htmltools", lib=my.lib)
## For R4.3 if error persists install.packages("httpuv", lib=my.lib)
install_if_not_installed("devtools", my.lib)
install_if_not_installed("BiocManager", my.lib)
BiocManager::install(version="3.20")
BiocManager::install("biomaRt", lib=my.lib)

## load devtools and BiocManager (will need these later)
library(devtools)
library(BiocManager)

## install other dependencies
install_if_not_installed(c("ggExtra", "dplyr", "ggplot2", "xlsx",
                            "stringr", "rjson"), lib=my.lib)

## on Eddie hdf5r has to be installed first loading hdf5 library
if (grepl("^node*", Sys.info()["nodename"])) {
    system("module load igmm/apps/hdf5/1.8.16")
}
devtools::install_github(repo="hhoeflin/hdf5r", lib=my.lib)

## DropletUtils is installable using Bioconductor
install_if_not_installed("tensor", lib=my.lib)
BiocManager::install(c("DropletUtils", "S4Vectors"),
                        lib=my.lib, ask=FALSE)
## seurat@v4.4.0 is available on GitHub
devtools::install_github("satijalab/seurat@v4.4.0", lib=my.lib, ask=FALSE)

## by default SeuratObject_5.0.2 is installed with Seurat but it changes the
## layout of Seurat objects which makes subsetting tricky. Therefore, we
## install SeuratObject_4.1.4
devtools::install_github("satijalab/seurat-object@v4.1.4", lib=my.lib,
                         ask=FALSE)

## loomR may need to be installed from GitHub if it was not installed with
## Seurat
load_result <- tryCatch({
            library(loomR, character.only=TRUE)
            TRUE  # Return TRUE if loading succeeds
        }, error=function(e) {
            FALSE  # Return FALSE if there is an error
        })
if (!load_result) {
    devtools::install_github(repo="mojaveazure/loomR", ref="develop")
}
## sctransform required above dependencies so we install it last if it has not
## been installed togather with Seurat
load_result <- tryCatch({
            library(sctransform, character.only=TRUE)
            TRUE  # Return TRUE if loading succeeds
        }, error=function(e) {
            FALSE  # Return FALSE if there is an error
        })
if (!load_result) BiocManager::install("sctransform", lib=my.lib, ask=FALSE)

## On some linux distributions this error occurs in Seurat
## https://github.com/satijalab/seurat/issues/8100
## reinstall Matrix and irlba from source and then install seurat
## I did not experience this error on Eddie so the block below is not executed
## by default
if (FALSE) {
    install.packages("Matrix", type="source")
    install.packages("irlba", type="source")
}

## UPDATE: Install scDblFinder for doublet detection - substitution for Scrublet

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scDblFinder")
