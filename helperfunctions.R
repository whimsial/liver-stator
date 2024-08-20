#' Helper function to load a given Rdata object.
load.rdata <- function(filename) {
    load(filename)
    get(ls()[ls() != "filename"])
}

## download specified study
download.study <- function(this.study, studies.dt, root.dir) {
    this.study.name <- studies.dt[eval(this.study),
                                  gsub("^(\\w+).*", "\\1", Study)]
    this.archive.url <- studies.dt[eval(this.study), `Data URL`]
    this.filelist.url <- studies.dt[eval(this.study),
                                    gsub("(.*suppl).*", "\\1", `Data URL`)]
    this.filelist.url <- file.path(this.filelist.url, "filelist.txt")

    study.dir <- file.path(root.dir, this.study.name)
    dir.create(study.dir, showWarnings=FALSE)

    this.study.file <- file.path(study.dir, basename(this.archive.url))
    this.filelist.file <- file.path(study.dir, basename(this.filelist.url))

    if (!file.exists(this.study.file)) {
        curl.cmd <- paste0("curl %s --output %s --retry 100 --retry-delay 2 -s")
        system(sprintf(curl.cmd, this.archive.url, this.study.file))
        system(sprintf(curl.cmd, this.filelist.url, this.filelist.file))
    }

    study.properties <- list(name=this.study.name, dir=study.dir,
                             study.file=this.study.file,
                             study.file.list=this.filelist.file)
    return(study.properties)
}

## get and parse GEO soft meta data into R data.table
parse.soft <- function(soft.url, soft.file) {
    ## download soft file
    cmd <- sprintf("curl %s --output %s --retry 100 --retry-delay 2 -s",
                   soft.url, paste0(soft.file, ".gz"))
    system(cmd)

    ## extract and read into R
    cmd <- sprintf("gunzip -c %s > %s", paste0(soft.file, ".gz"), soft.file)
    system(cmd)
    tmp <- readLines(soft.file)

    ## grep ^SAMPLE and !Sample_title key words
    dt <- data.table(sample=tmp[grep("\\^SAMPLE", tmp)],
                     sample.title=tmp[grep("\\!Sample_title", tmp)])

    ## cleanup the strings
    dt[, sample := gsub("\\^SAMPLE = ", "", sample)]
    dt[, sample.title := gsub("\\!Sample_title = ", "", sample.title)]
    return(dt)
}
