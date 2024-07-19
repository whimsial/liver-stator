## High-resolution mapping of cell states in liver disease
source("properties.R")
source("helperfunctions.R")

## download studies with scRNA-seq data (see full list in liver_studies.tsv)
this.study <- 1 ## Ramachandran et al (2019)
study1 <- download.study(this.study, studies.dt, root.dir)

this.study <- 2 ## Guilliams et al (2022)
study2 <- download.study(this.study, studies.dt, root.dir)

## QC step 1
