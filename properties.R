library(data.table)
library(Seurat) ## install version 4.4.0 and seurat-object 4.1.4
library(stringr)
library(DropletUtils)
library(HDF5Array)
library(biomaRt)
library(ggplot2)
library(cowplot)

## install packages for QC
if (install.dependencies) source("dependencies.R")

# root.dir <- "/exports/igmm/eddie/khamseh-lab/aiakovliev"
root.dir <- "/opt/datastore/aiakovliev/liver"
# codebase.dir <- file.path(root.dir, "liver-stator")
codebase.dir <- "~/repos/liver-stator"

## downloads
studies.list.file <- file.path(codebase.dir, "liver_studies.tsv")
studies.dt <- fread(studies.list.file)

## dublet threshold (inferred from distribution plots in QC step 2)
dublet.threshold <- 0.15

new.run <- FALSE
