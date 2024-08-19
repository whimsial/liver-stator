#' Filter Cells in a Seurat Object
#'
#' This function performs multiple quality control (QC) checks to filter cells
#' in a Seurat object based on gene expression, mitochondrial content, and
#' other metadata-based criteria. It also generates a diagnostic plot to
#' visually assess the quality of the data after filtering.
#'
#' @param seurat.object `Seurat` object with single-cell RNA sequence data.
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
#' @param output.dir Full path to the directory where the diagnostic plot
#'        will be saved. Default is the current working directory (`getwd()`).
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
filter.seurat <- function(seurat.object,
                          min.features=500, min.counts=1000, max.mt=10,
                          n.deviation=3, n.cells=30, output.dir=getwd()) {
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
