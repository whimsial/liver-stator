root.dir <- "/exports/igmm/eddie/khamseh-lab/aiakovliev"

## download Ramachandran et al (2019)
ramachandran.url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/GSE136103_RAW.tar"
ramach.file <- "GSE136103_RAW.tar"
curl.cmd <- paste0("curl %s --output %s --retry 100 --retry-delay 2 -s")
system(sprintf(curl.cmd, ramachandran.url, file.path(root.dir, ramach.file)))
