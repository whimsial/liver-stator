## download specified study
download.study <- function(this.study, studies.dt) {
    this.study.name <- studies.dt[eval(this.study), gsub("^(\\w+).*", "\\1", Study)]
    this.archive.url <- studies.dt[eval(this.study), `Data URL`]
    this.filelist.url <- studies.dt[eval(this.study),
                                    gsub("(.*suppl).*", "\\1", `Data URL`)]
    this.filelist.url <- file.path(this.filelist.url, "filelist.txt")

    study.dir <- file.path(root.dir, this.study.name)
    dir.create(study.dir, showWarnings=FALSE)

    this.study.file <- file.path(study.dir, basename(this.archive.url))
    this.filelist.file <- file.path(study.dir, basename(this.filelist.url))

    curl.cmd <- paste0("curl %s --output %s --retry 100 --retry-delay 2 -s")
    system(sprintf(curl.cmd, this.archive.url, this.study.file))
    system(sprintf(curl.cmd, this.filelist.url, this.filelist.file))
}
