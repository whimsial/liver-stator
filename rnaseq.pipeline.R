## get gene symbols from Ensembl ids
get.gene.symbols <- function(ensembl.ids, ensembl) {
    gene.symbols <- setDT(getBM(attributes=c("ensembl_gene_id",
                                             "external_gene_name"),
                                filters="ensembl_gene_id",
                                values=ensembl.ids$gene.id,
                                mart=ensembl))
    cat(nrow(gene.symbols), "of", nrow(ensembl.ids),
        "genes matched to Ensembl by ID.\n")
    return(gene.symbols)
}

## get Ensembl ids from gene symbols
get.gene.ids <- function(gene.symbols, ensembl) {
    ensembl.ids <- setDT(getBM(attributes=c("ensembl_gene_id",
                                            "external_gene_name"),
                               filters="external_gene_name",
                               values=gene.symbols$gene.symbol,
                               mart=ensembl))
    ensembl.ids <- unique(ensembl.ids, by="external_gene_name")
    cat(nrow(ensembl.ids), "of", nrow(gene.symbols),
        "genes matched to Ensembl by gene name.\n")
    return(ensembl.ids)
}

## plot variable genes
plot.variable.genes <- function(seurat.object, genes.dt, output.file) {
    p <- VariableFeaturePlot(seurat.object, pt.size=0.2)
    p2 <- LabelPoints(plot=p, points=genes.dt$gene.id,
                      labels=genes.dt$gene.symbol, repel=TRUE,
                      xnudge=0, ynudge=0, max.overlaps=Inf, force=7,
                      colour="blue")
    p2 <- p2 + scale_y_log10() + labs(x="Average Expression",
                                      y="Standardized Variance, Log10") +
          theme_minimal()
    cat("Saving plot to:", output.file, "\n")
    ggsave(output.file, plot=p2, width=7, height=7, units="in")
}

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
#' seurat_obj <- create.seurat(counts, this.sample="GSM4041169",
#'                             this.sample.dir="GSM4041169")
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

#' Filter Cells in a Seurat Object
#'
#' This function performs multiple quality control (QC) checks to filter cells
#' in a Seurat object based on gene expression, mitochondrial content, and
#' other metadata-based criteria. It also generates a diagnostic plot to
#' visually assess the quality of the data after filtering.
#'
#' @param seurat.object `Seurat` object with single-cell RNA sequence data.
#' @param output.dir Full path to the directory where the diagnostic plot
#'        will be saved.
#' @param min.features Integer specifying minimum number of features (genes)
#'        per cell. Default is 500.
#' @param min.counts Integer specifying minimum number of transcripts (total
#'        RNA counts) per cell. Default is 1000.
#' @param max.mt Numeric specifying maximum percentage of mitochondrial genes
#'        allowed per cell. Cells with a higher percentage will be removed.
#'        Default is 10.
#' @param n.deviation Integer number of median absolute deviations (MAD)
#'        from the median expression values used to identify and filter
#'        outliers. Default is 3.
#' @param n.cells Integer specifying minimum number of cells in which a gene
#'        must be detected to be retained. This helps filter out
#'        rarely expressed genes. Default is 30.
#'
#' @return A filtered `Seurat` object with cells and genes that pass filtering.
#'
#' @examples
#' # Filter with default parameters
#' filtered.seurat <- filter.seurat(my.seurat.object)
#'
#' # Filter with custom parameters
#' filtered.seurat <- filter.seurat(my.seurat.object,
#'                                  min.features=300,
#'                                  output.dir="/path/to/output")
#'
#' @importFrom Seurat subset
#' @import ggplot2
#' @importFrom data.table as.data.table, setDT
#' @importFrom stats median, mad
filter.seurat <- function(seurat.object, output.dir,
                          min.features=500, min.counts=1000, max.mt=10,
                          n.deviation=3, n.cells=30) {
    require(Seurat)
    require(data.table)
    require(ggplot2)

    ## Remove cells with transcripts/cell counts greater than 3 median absolute
    ## deviation (MAD) away from the median
    all.samples <- meta[, unique(sample)]

    cells.to.remove <- NULL
    for (this.sample in all.samples) {
        meta.data <- setDT(seurat.object@meta.data)[sample==eval(this.sample)]

        count.threshold <- median(meta.data$nCount_RNA) +
                           n.deviation * mad(meta.data$nCount_RNA)
        feature.threshold <- median(meta.data$nFeature_RNA) +
                             n.deviation * mad(meta.data$nFeature_RNA)
        meta.data[, qc.fail := 0]
        meta.data[nCount_RNA > count.threshold | nFeature_RNA > feature.threshold,
                  qc.fail := 1]
        to.remove <- meta.data[qc.fail==1, barcodes]
        cat(length(to.remove),
            "cells will be removed from sample", this.sample,
            "due to feature/RNA count threshold.\n")
        cells.to.remove <- c(cells.to.remove, to.remove)
    }

    if (!all(cells.to.remove %in% colnames(seurat.object)))
        stop("Cells barcodes were not matched to Seurat colnames.
              Check barcode names.")

    cells.to.keep <- colnames(seurat.object)[!colnames(seurat.object)
                                             %in% cells.to.remove]

    df <- seurat.object@meta.data
    df <- df[barcodes %in% cells.to.keep]

    cat(nrow(df), "of", ncol(seurat.object),
        "cells remain after filtering by MAD\n")

    ## Remove cells with N features < min.features
    df <- df[nFeature_RNA>min.features]
    cat(nrow(df), "of", ncol(seurat.object),
        "cells remain after filtering by nFeature_RNA\n")

    ## Remove cells with N transcripts < min.counts (RNA counts)
    df <- df[nCount_RNA>min.counts]
    cat(nrow(df), "of", ncol(seurat.object),
        "cells remain after filtering by nCount_RNA\n")

    ## Remove cells with % MT genes > max.mt%
    df <- df[percent.mt<max.mt]
    cat(nrow(df), "of", ncol(seurat.object),
        "cells remain after filtering by percent.mt\n")

    ## Remove cells with doublet prediction > 0
    df <- df[doublet_prediction==0]
    cat(nrow(df), "of", ncol(seurat.object),
        "cells remain after removing doublets\n")

    ## perform subset and update metadata manually (this is a workaround for the
    ## issue where subset function from Seurat package would not run properly)
    cells.to.keep <- df$barcode
    seurat.object.filtered <- seurat.object[, cells.to.keep]
    seurat.object.filtered@meta.data <- df[barcodes %in% cells.to.keep]

    ## run some checks to ensure consistency between meta data and Seurat object
    if (any(is.na(seurat.object.filtered@meta.data)))
        stop("NAs detected in meta.data.")
    if (!identical(colnames(seurat.object.filtered),
        seurat.object.filtered@meta.data$barcode))
        stop("Cell names in counts matrix do not match meta.data.")
    if (any(!cells.to.keep %in% seurat.object.filtered@meta.data$barcode))
        stop("Cells labeled for removal detected in the filtered Seurat object.")

    ## plot the counts after filtering
    dt <- setDT(seurat.object.filtered@meta.data)
    p <- ggplot(dt, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
        geom_point() +
        scale_colour_gradient(low="gray90", high="black") +
        stat_smooth(method=lm) +
        scale_x_log10() +
        scale_y_log10() +
        labs(x="RNA counts, Log10", y="N transcripts, Log10")
    plot.file <- file.path(output.dir, "feature.plot.png")
    ggsave(file=plot.file, p)

    ## keep only genes which are expressed in n.cells or more cells
    gene.detection.rate <- PercentageFeatureSet(seurat.object.filtered,
                                                pattern="^", assay="RNA")
    counts <- GetAssayData(object=seurat.object.filtered, slot="counts")
    nonzero.counts <- counts > 0
    keep.genes <- Matrix::rowSums(nonzero.counts) >= n.cells
    cat(sum(keep.genes), "of", nrow(seurat.object.filtered),
        "genes are expressed in more than", n.cells,
        "cells and will be retained.\n")
    seurat.object.filtered <- seurat.object.filtered[keep.genes, ]
    seurat.object.filtered@meta.data <- setDF(dt)
    save(seurat.object.filtered,
         file=file.path(output.dir, "merged.filtered.seurat.Rdata.gz"))
    return(seurat.object.filtered)
}

#' Process and visualise variable genes in a Seurat object
#'
#' This function performs a series of operations on a Seurat object to identify
#' and visualise highly variable genes. It normalises data, identifies variable
#' features, maps gene IDs to Ensembl, and optionally plots these genes.
#'
#' The function handles both cases where genes are initially identified by
#' symbols or Ensembl IDs.
#'
#' @param seurat.object Seurat object containing single-cell RNA-seq data.
#' @param output.dir Full path to the directory where output files, such as
#'        plots, will be saved.
#' @param core.genes Optional character vector of core genes of interest for
#'        further analysis. If provided, the function will specifically analyse
#'        these genes alongside the general analysis of variable genes.
#' @param n.top.variable.genes Integer number of top variable genes to process
#'        and visualise.
#' @param label.genes Optional character vector of gene labels to be
#'        specifically marked in the visualisation.
#'
#' @return Seurat object with additional layers containing variable features
#' and normalised counts.
#'
#' @examples
#' seurat_obj <- process.variable.genes(seurat.object=seurat_obj,
#'                                      output.dir="path/to/output",
#'                                      core.genes=c("Gene1", "Gene2"),
#'                                      label.genes=c("GeneX", "GeneY"))
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures
#' @importFrom data.table data.table
#' @importFrom biomaRt useMart getBM
#' @import ggplot2
process.variable.genes <- function(seurat.object, output.dir, core.genes=NULL,
                                   n.top.variable.genes=20,
                                   label.genes=NULL) {
    ## normalize and identify variable features
    seurat.object <- NormalizeData(seurat.object,
                                   normalization.method="LogNormalize",
                                   scale.factor=10000)
    seurat.object <- FindVariableFeatures(seurat.object, selection.method="vst",
                                          nfeatures=1000)

    ## obtain top variable genes
    top.variable.genes <- head(VariableFeatures(seurat.object),
                               n.top.variable.genes)
    top.variable.genes <- data.table(gene.id=top.variable.genes,
                                     gene.symbol=top.variable.genes)

    ## connect to Ensembl via BioMart
    ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

    ## map gene IDs to Ensembl, depending on the format
    if (any(grepl("ENSG", top.variable.genes))) {
        ensembl.ids <- top.variable.genes[grep("ENSG", gene.id)]
        gene.symbols <- get.gene.symbols(ensembl.ids, ensembl)
        top.variable.genes[gene.symbols, on=c(gene.id="ensembl_gene_id"),
                           gene.symbol:=external_gene_name]
    } else {
        gene.symbols <- top.variable.genes[, gene.symbol]
        ensembl.ids <- get.gene.ids(gene.symbols, ensembl)
        top.variable.genes[ensembl.ids, on=c(gene.symbol="external_gene_name"),
                           gene.id:=ensembl_gene_id]
    }

    cat("Printing top", n.top.variable.genes,
        "variable genes with Ensembl ids.\n")
    print(top.variable.genes)

    ## plot most variable genes on Expression-Variance plot
    output.file <- file.path(output.dir, "most.variable.genes.png")
    plot.variable.genes(seurat.object, genes.dt=top.variable.genes, output.file)

    ## process core genes if specified
    if (!is.null(core.genes)) {
        cat(length(core.genes), "genes of interest were provided.\n")
        selected.genes <- data.table(gene.id=core.genes, gene.symbol=core.genes)
        ensembl.ids <- get.gene.ids(selected.genes, ensembl)
        selected.genes[ensembl.ids, on=c(gene.symbol="external_gene_name"),
                       gene.id:=ensembl_gene_id]

        all.hvg <- head(VariableFeatures(seurat.object), 1000)
        select.ids <- intersect(all.hvg, selected.genes$gene.id)
        selected.genes.matched <- selected.genes[gene.id %in% select.ids]
        cat(nrow(selected.genes.matched), "of", nrow(selected.genes),
            "genes matched in 1000 most highly variable genes.\n")
        print(selected.genes.matched)

        ## append labeled genes if specified
        if (!is.null(label.genes)) {
            label.genes.dt <- data.table(gene.id=label.genes,
                                         gene.symbol=label.genes)
            ensembl.ids <- get.gene.ids(label.genes.dt, ensembl)
            label.genes.dt[ensembl.ids, on=c(gene.symbol="external_gene_name"),
                           gene.id:=ensembl_gene_id]
            matched1 <- label.genes.dt[gene.id %in% rownames(seurat.object)]
            matched2 <- label.genes.dt[gene.symbol %in% rownames(seurat.object)]
            label.genes.dt <- rbind(matched1, matched2)
            selected.genes.matched <- rbind(selected.genes.matched, label.genes.dt)
        }

        ## plot core genes and labelled genes
        selected.genes.matched <- unique(selected.genes.matched)
        output.file <- file.path(output.dir, "selected.variable.genes.png")
        plot.variable.genes(seurat.object, genes.dt=selected.genes.matched,
                            output.file)
    }
    return(seurat.object)
}
