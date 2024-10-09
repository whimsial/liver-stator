## count mitochondrial genes measured by scRNA-seq and append to Seurat object
count.mt.genes <- function(seurat.obj) {
    require(Seurat)
    ## if MT- cannot be matched in gene names, then try excluding MT genes which
    ## are matched by Ensembl id
    if (sum(PercentageFeatureSet(seurat.obj, pattern="^MT-|^mt-")) != 0) {
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
#' @param counts A 10x counts matrix.
#' @param this.sample String specifying the name of the sample.
#' @param this.sample.dir Full path to the directory containing 10x data for the
#'        specified sample.
#'
#' @return A Seurat object with additional metadata fields and QC metrics.
#'
#' @details The function assumes the presence of 'barcodes.tsv' for naming
#'          conventions in cases of .mtx files and
#'          additionally handles doublet detection and mitochondrial content as
#'          part of quality control.
#'
#' @examples
#' seurat_obj <- create.seurat(counts, this.sample="GSM4041169", this.sample.dir="GSM4041169")
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom data.table fread, setDT
#' @importFrom ggplot2 ggplot, geom_violin, geom_jitter, geom_point, labs, theme_minimal
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggsave
create.seurat <- function(counts, this.sample, this.sample.dir) {

    pre <- CreateSeuratObject(counts=counts, project=this.sample)

    ## create SC index as <sample_id>:<barcode>
    barcodes.file <- file.path(this.sample.dir, "barcodes.tsv")
    if (file.exists(barcodes.file)) {
        barcodes <- fread(barcodes.file, header=FALSE)$V1
        barcodes <- paste0(this.sample, "_", barcodes)
    } else {
        barcodes <- paste0(this.sample, "_", colnames(pre))
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

process.samples.and.merge <- function(meta) {
    seurat.list <- list()
    all.samples <- meta[, unique(sample)]

    ## Loop through each sample and create a Seurat object
    for (this.sample in all.samples) {
        this.meta <- meta[sample == eval(this.sample)][1]
        cat("Processing sample:", this.sample, "\n")

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
        ## Guilliams et al report data in h5 format. These files should be read
        ## with Read10X_h5 function.
        ## Thus, we handle both cases here conditional on file types.
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

        seurat.list[[this.sample]] <- create.seurat(pre, this.sample,
                                                    this.sample.dir)
    }

    ## Merge all Seurat objects into one
    if (length(seurat.list) > 1) {
        Merge <- merge(x=seurat.list[[1]], y=seurat.list[2:length(seurat.list)],
                       add.cell.ids=all.samples)
    } else {
        Merge <- seurat.list[[1]]
    }

    ## Optional: Add clinical metadata for samples (if any)
    ## Assuming clinical metadata needs to be merged
    ## Merge <- AddClinicalData(Merge, clinical_data)

    ## Save the merged Seurat object
    seurat.file <- file.path(unique(dirname(meta$sample.dir)),
                             "merged.seurat.Rdata.gz")
    save(Merge, file=seurat.file)

    ## Return the merged Seurat object
    return(Merge)
}
