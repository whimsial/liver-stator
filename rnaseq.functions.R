#' Download GEO study data and file list
#'
#' This function downloads a specified study and its associated file list to a
#' designated directory. It creates a new directory for the study,
#' downloads the main study file and a file list from provided URLs,
#' and returns a list of properties about the downloaded study.
#'
#' @param study.name Character string specifying the name of the study.
#'        This name is used to create a subdirectory within `root.dir`.
#' @param study.url GEO URL from where the study should be downloaded.
#' @param filelist.url GEO URL containing the URL from where the file list
#'        should be downloaded.
#' @param root.dir Full path to the root directory under which a new directory
#'        will be created for the study.
#'
#' @return A list with the following components:
#'   \itemize{
#'     \item \code{name}: The name of the study.
#'     \item \code{dir}: The directory path where the study and file list are stored.
#'     \item \code{study.file}: The path to the downloaded study file.
#'     \item \code{study.file.list}: The path to the downloaded file list.
#'   }
#'
#' @examples
#' \dontrun{
#'   study.info <- download.study(
#'      "Ramachandran",
#'      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/GSE136103_RAW.tar",
#'      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl",
#'      getwd())
#'   print(study.info)
#' }
#'
download.study <- function(study.name, study.url, filelist.url, root.dir) {
    study.dir <- file.path(root.dir, study.name)
    dir.create(study.dir, showWarnings=FALSE)

    study.file <- file.path(study.dir, basename(study.url))
    filelist.file <- file.path(study.dir, basename(filelist.url))

    if (!file.exists(study.file)) {
        cat("Downloading from", study.url, "to", study.file, "...\n")
        curl.cmd <- "curl %s --output %s --retry 100 --retry-delay 2 -s"
        system(sprintf(curl.cmd, study.url, study.file))
        system(sprintf(curl.cmd, filelist.url, filelist.file))
    }

    ## extract archive
    if (!dir.exists(study.dir)) {
        cmd <- sprintf("tar -xvf %s -C %s", study.file, study.dir)
        system(cmd)
    }

    study.properties <- list(name=study.name, dir=study.dir,
                             study.file=study.file,
                             study.file.list=filelist.file)
    return(study.properties)
}

#' Download and parse GEO Soft Metadata into R Data Table
#'
#' This function downloads a GEO soft file from a specified URL, decompresses it,
#' and parses it to extract sample identifiers and titles, organizing them into
#' a data.table object.
#'
#' @param soft.url GEO URL from where the GEO soft file should be downloaded.
#' @param soft.file Full path with the filename of the file where the
#'        downloaded and extracted data will be stored. The function
#'        will append ".gz" to this filename for the compressed download
#'        and remove it after extraction.
#'
#' @return A `data.table` object containing two columns:
#'   \itemize{
#'     \item \code{sample}: Extracted sample identifiers.
#'     \item \code{sample.title}: Corresponding sample titles.
#'   }
#'
#' @details The function performs the following steps:
#'   * Downloads the .gz compressed soft file from the provided URL.
#'   * Extracts the .gz file into the specified path.
#'   * Reads the file lines and extracts lines starting with "^SAMPLE" for sample
#'     identifiers and lines starting with "!Sample_title" for sample titles.
#'   * Cleans up the extracted strings to remove preceding keywords and equal
#'     signs.
#'
#' @examples
#' \dontrun{
#'   soft.url <- "http://example.com/sample.soft.gz"
#'   soft.file <- "path/to/sample.soft"
#'   sample.data <- parse.soft(soft.url, soft.file)
#'   print(sample.data)
#' }
#'
#' @import data.table
#' @export
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

    ## cleanup strings
    dt[, sample := gsub("\\^SAMPLE = ", "", sample)]
    dt[, sample.title := gsub("\\!Sample_title = ", "", sample.title)]
    return(dt)
}

#' Retrieve gene symbols from Ensembl IDs
#'
#' This function takes a set of Ensembl gene IDs and uses the biomaRt package to
#' retrieve the corresponding gene symbols from the specified Ensembl database.
#'
#' @param ensembl.ids Data.table or data.frame containing the Ensembl gene IDs
#'        in a column named `gene.id`.
#' @param ensembl BioMart database object representing the Ensembl database;
#'        typically this is created using the `useMart` function from the
#'        biomaRt package.
#'
#' @return A data.table with two columns:
#'   \itemize{
#'     \item \code{ensembl_gene_id}: The Ensembl gene IDs queried.
#'     \item \code{external_gene_name}: The corresponding gene symbols for the queried IDs.
#'   }
#'   The function also prints the number of matched and total queried gene IDs.
#'
#' @details The function requires the biomaRt package and expects that the ensembl
#'          object passed to it has been properly initialized with the useMart
#'          function. It uses the `getBM` function to perform the query.
#'
#' @examples
#' \dontrun{
#'   library(data.table)
#'   library(biomaRt)
#'   ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#'   ensembl.ids <- data.table(gene.id = c("ENSG00000139618", "ENSG00000155657"))
#'   gene.symbols <- get.gene.symbols(ensembl.ids, ensembl)
#'   print(gene.symbols)
#' }
#'
#' @importFrom biomaRt getBM
#' @importFrom data.table setDT
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

#' Retrieve Ensembl IDs from gene symbols
#'
#' This function queries an Ensembl database to find Ensembl gene IDs corresponding
#' to a given list of gene symbols using the biomaRt package.
#'
#' @param gene.symbols Data.table or data.frame containing the gene symbols
#'        in a column named `gene.symbol`.
#' @param ensembl BioMart database object representing the Ensembl database;
#'        this is typically initialized with the `useMart` function from
#'        the biomaRt package.
#'
#' @return A unique data.table that includes:
#'   \itemize{
#'     \item \code{ensembl_gene_id}: The Ensembl gene IDs corresponding to the queried gene symbols.
#'     \item \code{external_gene_name}: The gene symbols provided in the query.
#'   }
#'   Additionally, a message is printed showing the number of gene symbols matched
#'   and the total number queried.
#'
#' @details The function makes use of the `getBM` function from the biomaRt package
#'          to perform the query. It is important that the `ensembl` object has been
#'          properly initialized. The results are filtered to provide unique matches
#'          based on gene symbols.
#'
#' @examples
#' \dontrun{
#'   library(data.table)
#'   library(biomaRt)
#'   ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#'   gene.symbols <- data.table(gene.symbol = c("BRCA2", "TP53"))
#'   ensembl.ids <- get.gene.ids(gene.symbols, ensembl)
#'   print(ensembl.ids)
#' }
#'
#' @importFrom biomaRt getBM
#' @importFrom data.table setDT
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

#' Plot variable genes from Seurat object
#'
#' Visualizes variable genes on a scatter plot using data from a Seurat object.
#' This function highlights specific genes and saves the plot to a file.
#'
#' @param seurat.object A Seurat object containing variable gene data and expression
#'        metrics calculated, typically with functions such as `FindVariableFeatures`.
#' @param genes.dt Data.table or data.frame containing the gene identifiers
#'        and symbols to be highlighted; must include columns `gene.id`
#'        for gene identifiers and `gene.symbol` for gene names.
#' @param output.file Full path to the file where the plot should be saved.
#'        Supports any file format compatible with `ggsave`,
#'        including but not limited to png, pdf, and jpeg formats.
#'
#' @details The function first creates a basic variable feature plot using the
#'          `VariableFeaturePlot` function. It then overlays labels using
#'          `LabelPoints` from the Seurat package, which repels labels to avoid
#'          overlap and adjusts aesthetics. The plot is then logged on the y-axis,
#'          labeled, and styled with a minimal theme. Finally, it is saved to the
#'          specified output file using `ggsave`.
#'
#' @return Invisible. The function is called for its side effects: it creates
#'         and saves a plot. It does not return a value.
#'
#' @examples
#' \dontrun{
#'   library(Seurat)
#'   library(data.table)
#'   # Assume `seurat_obj` is a preprocessed Seurat object.
#'   genes.dt <- data.table(gene.id = c("gene1", "gene2"), gene.symbol = c("Gene1", "Gene2"))
#'   output.file <- "variable_genes_plot.pdf"
#'   plot.variable.genes(seurat.object = seurat_obj, genes.dt = genes.dt, output.file = output.file)
#' }
#'
#' @importFrom Seurat VariableFeaturePlot LabelPoints
#' @importFrom ggplot2 ggsave scale_y_log10 labs theme_minimal
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

#' Count Mitochondrial Genes in scRNA-Seq data and append to Seurat object
#'
#' This function calculates the percentage of mitochondrial genes
#' (based on gene symbols starting with "MT-" or "mt-" or by matching Ensembl
#' IDs if no gene symbols match) and appends this as metadata to the Seurat
#' object.
#'
#' @param seurat.obj Seurat object containing single-cell RNA sequencing data.
#'
#' @return A modified Seurat object that includes a new metadata field `percent.mt`
#'         which contains the percentage of mitochondrial gene expression per cell.
#'
#' @details The function first attempts to identify mitochondrial genes by looking
#'          for gene names that start with "MT-" or "mt-". If no such genes are
#'          found, it will then try to identify mitochondrial genes by their
#'          Ensembl gene IDs, specifically those that are located on the
#'          mitochondrial chromosome (MT). This involves querying the Ensembl
#'          database using the `biomaRt` package.
#'
#' @examples
#' \dontrun{
#'   library(Seurat)
#'   # Assuming `seurat_obj` is a Seurat object already created and preprocessed
#'   updated_seurat_obj <- count.mt.genes(seurat_obj)
#'   print(updated_seurat_obj[["percent.mt"]])
#' }
#'
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom biomaRt useMart getBM
#' @importFrom data.table setDT
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

#' Process samples and merge into a single Seurat object
#'
#' This function reads individual sample data from specified directories and file
#' formats, creates Seurat objects for each sample, and merges them into a single
#' Seurat object. The function supports both "mtx" and "h5" file formats.
#'
#' @param meta Data.table containing metadata for each sample.
#'             The metadata table must include the columns `sample`, `sample.dir`,
#'             and potentially `Name` for h5 files. Each row should correspond to
#'             a sample with its directory path where the data files are stored.
#'
#' @return A merged Seurat object containing data from all samples.
#'
#' @details The function processes each sample based on metadata entries. It checks
#'          for the presence of either "mtx" or "h5" files in the directory specified
#'          by `sample.dir`. Depending on the file type, it reads the data using
#'          appropriate methods (`read10xCounts` for "mtx" and `Read10X_h5` for "h5").
#'          Each dataset is converted into a Seurat object which are then merged into
#'          one Seurat object. The function requires the Seurat, data.table, and
#'          DropletUtils packages.
#'
#'          It is assumed that all required libraries are installed and loaded.
#'          Error handling includes checks for valid directory paths, correct sample
#'          naming in metadata, and the presence of necessary data files.
#'
#' @examples
#' \dontrun{
#'   # Assume meta is a data.table with appropriate columns
#'   merged_seurat <- process.samples.and.merge(meta)
#'   print(merged_seurat)
#' }
#'
#' @importFrom Seurat CreateSeuratObject merge
#' @importFrom DropletUtils read10xCounts
#' @importFrom rhdf5 Read10X_h5
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
