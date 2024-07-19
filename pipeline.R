## High-resolution mapping of cell states in liver disease
source("properties.R")
source("helperfunctions.R")

## Process Ramachandran et al (2019) (see full list in liver_studies.tsv)
## -----------------------------------------------------------------------------
this.study <- 1 ## Ramachandran et al (2019)
study <- download.study(this.study, studies.dt, root.dir)

## extract archive
cmd <- sprintf("tar -xvf %s -C %s", study$study.file, study$dir)
system(cmd)

## study variable should be a list with at least 2 elements:
## $dir -> name of the directory containing extract from the GEO archive
## $study.file.list -> content of filelist.txt which contains a list of all
##                     files for the study (downloaded separately from GEO)
sample.files <- list.files(study$dir, pattern="*.tsv.gz", full.names=TRUE)
meta <- fread(study$study.file.list)
stopifnot(all(basename(sample.files) %in% meta$Name))

## create dir path for each sample
meta <- fread(study$study.file.list)
meta[, sample := gsub("(GSM\\d+).*", "\\1", Name)]
meta[, sample.dir := file.path(study$dir, sample)]
meta <- meta[!grep(".tar", sample)]
meta <- meta[!grep("blood|mouse", Name)]

## QC step 1
source("qc/qc1_emptydrop.R")

this.study <- 2 ## Guilliams et al (2022)
study <- download.study(this.study, studies.dt, root.dir)
