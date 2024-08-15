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

## Get and parse GEO soft meta data into R data.table
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

## count mitochondrial genes measured by scRNA-seq and append to Seurat object
count.mt.genes <- function(seurat.obj) {
    require(Seurat)
    ## if MT- cannot be matched in gene names, then try excluding MT genes which
    ## are matched by Ensembl id
    if (sum(PercentageFeatureSet(seurat.obj, pattern="^MT-|^mt-"))!=0) {
        seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj,
                                                           pattern="^MT-|^mt-")
    } else {
        ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        ## fetch mitochondrial genes using a filter
        mt.genes <- setDT(
            getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                  filters='chromosome_name',
                  values='MT',
                  mart=ensembl))
        ## get all gene names from the Seurat object
        all.genes <- rownames(seurat.obj)

        ## find MT genes in Seurat object
        mt.ensembl <- mt.genes[ensembl_gene_id %in% all.genes, ensembl_gene_id]

        ## save counts of MT genes
        seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj,
                                                           features=mt.ensembl)
    }
    return(seurat.obj)
}

#' Create a Seurat object from SingleCellExperiment data
#'
#' This function reads single-cell data from either .mtx or .h5 files located
#' in a specified directory, creates a Seurat object, and enriches it with
#' metadata including doublet predictions, mitochondrial RNA content, and
#' additional QC plots.
#'
#' @param this.sample A character string specifying the default sample name.
#'                    Default is "GSM4041150".
#' @param this.meta A data.table object containing metadata. This should
#'                  have exactly one row and must include 'sample',
#'                  'sample.dir', and 'Name' columns.
#'
#' @return A Seurat object with additional metadata fields and QC metrics.
#'
#' @details The function handles two specific data formats described in:
#'          Ramachandran et al., which involves .mtx files, and
#'          Guilliams et al., which involves .h5 files.
#'          Depending on the files present in the directory, appropriate reading
#'          functions are called. The function assumes the presence of
#'          'barcodes.tsv' for naming conventions in cases of .mtx files and
#'          additionally handles doublet detection and mitochondrial content as
#'          part of quality control.
#'
#' @examples
#' # Assuming 'metadata' is a data.table with appropriate structure
#' seurat_obj <- create.seurat(this.sample = "Sample123", this.meta = metadata)
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom data.table fread, setDT
#' @importFrom ggplot2 ggplot, geom_violin, geom_jitter, geom_point, labs, theme_minimal
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggsave
create.seurat <- function(this.sample="GSM4041150", this.meta) {
    if (this.meta[, .N] != 1) {
        stop("Metadata should be a data.table with 1 row")
    }

    if (this.meta$sample != this.sample) {
        stop("Sample name is not matched to metadata. Check input.")
    }

    this.sample.dir <- this.meta[, unique(sample.dir)]
    if (!dir.exists(this.sample.dir)) {
        stop("Sample directory listed in metadata does not exist. Check input.")
    }

    ## Ramachandran et al have data for SingleCellExperiment consisting of
    ## matrix.mtx, genes.tsv and barcodes.tsv.
    ## It should be read using read10xCounts and then the sparse matrix of
    ## counts should be extracted and passed to Seurat
    ## Guilliams et al report data in h5 format. These files should be read with
    ## Read10X_h5 function. Thus, we handle both cases here conditional on file
    ## types.
    if (any(grepl("mtx", list.files(this.sample.dir)))) {
        pre <- read10xCounts(this.sample.dir, col.names=TRUE)
        pre <- counts(pre)
    } else if (any(grepl("h5", list.files(this.sample.dir)))) {
        this.h5.file <- this.meta[, Name]
        pre <- Read10X_h5(file.path(this.sample.dir, this.h5.file),
                          use.names=TRUE, unique.features=TRUE)
    } else {
        stop("Either .mtx or .h5 file should be present. Check input directory.")
    }

    pre <- CreateSeuratObject(counts=pre, project=this.sample)

    ## create SC index as <sample_id>:<barcode>
    barcodes.file <- file.path(this.sample.dir, "barcodes.tsv")
    if (file.exists(barcodes.file)) {
        barcodes <- fread(barcodes.file, header=FALSE)$V1
        barcodes <- paste0(this.sample, ":", barcodes)
    } else {
        barcodes <- paste0(this.sample, ":", colnames(pre))
    }
    pre$barcodes <- barcodes
    pre$sample <- this.sample
    pre$orig.ident <- this.sample

    ## read Dublet predictions and scores
    prediction.file <- file.path(this.sample.dir,
                                 "scrublet_EDR0.008_PredictedDoublets.csv")
    dublet.file <- file.path(this.sample.dir,
                             "scrublet_EDR0.008_DoubletScores.csv")
    doublet_prediction <- fread(prediction.file)
    doublet_score <- fread(dublet.file)

    ## a given cell is a duplet if doublet_score > threshold.
    ## original threshold was 0.35, but we set a more stringent threshold
    ## of 0.15 which corresponds to the threshold for dublet detection
    ## from step 2.
    dublets <- data.table(prediction=doublet_prediction$V1,
                          score=doublet_score$V1)
    dublets[, new.prediction := ifelse(score>dublet.threshold, 1, 0)]

    ## append info on duplets and mitochondrial RNA to Seurat metadata object
    pre[["doublet_prediction"]] <- dublets[, new.prediction]
    pre[["doublet_score"]] <- dublets[, score]
    pre <- count.mt.genes(pre)

    ## compute log10 of the proportion of identified genes to counts of all
    ## RNA molecules in each cell to identify anomalies in gene expression due
    ## to mitochondrial or ribosomal RNA.
    pre$log10GenesPerUMI <- log10(pre$nFeature_RNA) / log10(pre$nCount_RNA)

    ## plot numbers of genes, RNA counts, MT RNAs and ratio of counts of genes
    ## to RNAs
    dt <- setDT(pre@meta.data)
    dt[, cell.idx := .I]
    p1 <- ggplot(dt, aes(x=cell.idx, y=nFeature_RNA)) +
        geom_violin(trim=FALSE) +
        geom_jitter(width=0.1, size=0.2, alpha=0.2) +
        labs(x="Cells", y="nFeature_RNA") +
        theme_minimal()
    p2 <- ggplot(dt, aes(x=cell.idx, y=nCount_RNA)) +
        geom_violin(trim=FALSE) +
        geom_jitter(width=0.1, size=0.2, alpha=0.2) +
        labs(x="Cells", y="nCount_RNA") +
        theme_minimal()
    p3 <- ggplot(dt, aes(x=cell.idx, y=log10GenesPerUMI)) +
        geom_violin(trim=FALSE) +
        geom_jitter(width=0.1, size=0.2, alpha=0.2) +
        labs(x="Cells", y="log10GenesPerUMI") +
        theme_minimal()
    p4 <- ggplot(dt, aes(x=cell.idx, y=percent.mt)) +
        geom_violin(trim=FALSE) +
        geom_jitter(width=0.1, size=0.2, alpha=0.2) +
        labs(x="Cells", y="percent.mt") +
        theme_minimal()
    p <- plot_grid(p1, p2, p3, p4, nrow=1, align='v')
    feature.plot <- file.path(this.sample.dir, "feature.plot.png")
    ggsave(file=feature.plot, p)

    ## plot scatter plots of nCount_RNA against percent.mt and nCount_RNA
    ## against nFeature_RNA
    p1 <- ggplot(dt, aes(x=nCount_RNA, y=percent.mt)) +
        geom_point(size=0.2, alpha=0.2) +
        labs(x="nCount_RNA", y="percent.mt") +
        theme_minimal()
    p2 <- ggplot(dt, aes(x=nCount_RNA, y=nFeature_RNA)) +
        geom_point(size=0.2, alpha=0.2) +
        labs(x="nCount_RNA", y="nFeature_RNA") +
        theme_minimal()
    p <- plot_grid(p1, p2, nrow=1, align='v')
    scatter.plot <- file.path(this.sample.dir, "scatter.plot.png")
    ggsave(file=scatter.plot, p)
    return(pre)
}
