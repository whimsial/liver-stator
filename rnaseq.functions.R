#' Report a message to the terminal
#'
#' @param mode One of \sQuote{info} (normal messages), \sQuote{note} (messages
#'        that require some highlighting), \sQuote{warn} (important information
#'        the user should definitely notice).
#' @param ... Strings to be reported.
#' @param LF Whether a newline character should be added at the end of the
#'        message (\code{TRUE} by default).
#'
#' @import crayon
#' @noRd
msg <- function(mode, ..., LF=TRUE) {
    message(mode(...), appendLF=LF)
}
info <- crayon::reset
note <- crayon::green
warn <- crayon::yellow
bold <- crayon::bold
error <- crayon::red

#' Download GEO study data and file list
#'
#' This function downloads a specified study and its associated file list to a
#' designated directory. It creates a new directory for the study and downloads
#' the main study file and a file list from provided URLs.
#'
#' @param study.name Character string specifying the name of the study.
#'        This name is used to create a subdirectory within `root.dir`.
#' @param study.file Full system path to the archive to which the study will be
#'        saved. This file will be created if it does not exist.
#' @param study.url GEO URL from where the study should be downloaded.
#' @param filelist.file Full system path to the file where study file list will
#'        be saved. This file will be created if it does not exist.
#' @param filelist.url GEO URL containing the URL from where the file list
#'        should be downloaded.
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#'   study.info <- download.study(
#'      "Ramachandran",
#'      "~/Ramachandran/GSE136103_RAW.tar",
#'      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/GSE136103_RAW.tar",
#'      "~/Ramachandran/filelist.txt",
#'      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/filelist.txt")
#' }
download.study <- function(study.name, study.file, study.url,
                           filelist.file, filelist.url) {
    study.dir <- dirname(study.file)
    dir.create(study.dir, showWarnings=FALSE)

    msg.txt <- sprintf("Downloading %s study to %s ...",
                       study.name, study.dir)
    msg(note, msg.txt)

    if (!file.exists(study.file)) {
        msg.txt <- sprintf("Downloading from %s to %s ...",
                           study.url, study.file)
        msg(info, msg.txt)
        curl.cmd <- "curl %s --output %s --retry 100 --retry-delay 2 -s"
        system(sprintf(curl.cmd, study.url, study.file))
        system(sprintf(curl.cmd, filelist.url, filelist.file))
    }

    ## extract archive
    if (!dir.exists(study.dir)) {
        cmd <- sprintf("tar -xvf %s -C %s", study.file, study.dir)
        system(cmd)
    }

    msg(note, "Done.")
}

#' Download and parse GEO Soft Metadata into R data.table
#'
#' This function downloads a GEO soft file from a specified URL, decompresses it,
#' and parses it to extract sample identifiers and titles, organizing them into
#' a data.table object.
#'
#' @param soft.url GEO URL from where the GEO soft file should be downloaded.
#' @param soft.file Full path with the filename of the file where the
#'        downloaded and extracted data will be stored.
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
#'   soft.file <- "path/to/sample.soft.gz"
#'   sample.data <- parse.soft(soft.url, soft.file)
#'   print(sample.data)
#' }
#'
#' @import data.table
parse.soft <- function(soft.url, soft.file) {
    msg.txt <- sprintf("Downloading softfile from %s to %s.",
                       soft.url, soft.file)
    msg(info, msg.txt)

    ## download soft file
    cmd <- sprintf("curl %s --output %s --retry 100 --retry-delay 2 -s",
                   soft.url, soft.file)
    system(cmd)

    ## extract and read into R
    cmd <- sprintf("gunzip -c %s > %s", soft.file,
                   gsub(".gz", "", soft.file))
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

#' Extract sample files from GEO archives
#'
#' This function checks for the existence of specified archive files, creates a
#' directory for sample extraction, and extracts relevant files based on file
#' type. Extracted files are either uncompressed or linked depending on the
#' file extension.
#'
#' @param this.archives Character vector containing paths to archive files.
#' @param this.sample.dir Character string indicating the directory where
#'        samples should be extracted.
#'
#' @return A character vector of paths to extracted files.
#'
#' @details
#' The function first checks if all provided archive files exist. If any file
#' does not exist, the function stops and returns an error message.
#'
#' If the archive files exist, the function proceeds to create the specified
#' output directory if it does not already exist.
#' It then checks the types of files contained in the archives:
#' - For compressed `.gz` files containing `barcodes.tsv`, `genes.tsv`, or
#'   `matrix.mtx`, it decompresses them into the specified directory.
#' - For `.h5` files, it creates a symbolic link in the target directory
#' pointing to the original file.
#'
#' Only files matching the specified patterns are processed. Other file types
#' within the archives are ignored.
extract.samples <- function(this.archives, this.sample.dir) {
    if (!all(file.exists(this.archives))) {
        msg.txt <- sprintf("Provided archives do not exist:\n %s \n",
                           this.archives)
        stop(msg(error, msg.txt))
    }
    dir.create(this.sample.dir, showWarnings=TRUE)

    ## catch if archives contain csv and mtx files and extract them
    re <- ".*(barcodes.tsv|genes.tsv|matrix.mtx).gz"
    if (any(grepl(re, this.archives))) {
        this.extracts <- gsub(re, "\\1", basename(this.archives))
        this.extracts <- file.path(this.sample.dir, this.extracts)
        cmd <- sprintf("gunzip -c %s > %s", this.archives, this.extracts)
        for (this.cmd in cmd) system(this.cmd)
    }

    ## otherwise if its an h5 archive we make a symbolic link
    if (any(grepl("h5", this.archives))) {
        this.extracts <- file.path(this.sample.dir, basename(this.archives))
        if (!file.exists(this.extracts))
            file.symlink(this.archives, this.extracts)
    }

    return(this.extracts)
}

#' Read 10x genomics data
#'
#' This function reads 10x Genomics data from provided file paths or directories
#' containing either MTX files (barcodes.tsv, genes.tsv, matrix.mtx) or an H5
#' file.
#'
#' @param this.extracts Character vector containing paths to extracted files or
#'        an H5 file.
#' @param this.sample.dir Character string indicating the directory where the
#'        10x experiment data is located or will be read from.
#'
#' @return A list containing two components:
#' \itemize{
#'  \item{sample.barcodes}{A data table of barcodes read from the provided files.}
#'  \item{sce}{A SingleCellExperiment object containing the expression data read from the files.}
#' }
#'
#' @details
#' The function first checks if MTX files (barcodes.tsv, genes.tsv, matrix.mtx)
#' are present in the file paths provided by `this.extracts`.
#' If these files are present, it reads the barcodes and uses
#' \code{read10xCounts} to read the 10x experiment data from `this.sample.dir`.
#'
#' If an H5 file is provided, it reads the barcodes and the single-cell e
#' xperiment data directly from the H5 file using \code{h5read} and
#' \code{read10xCounts}, respectively.
#'
#' The function returns an error if neither MTX files nor an H5 file are found.
#'
#' @importFrom data.table fread data.table
#' @importFrom rhdf5 h5read
#' @importFrom DropletUtils read10xCounts
read.10x.data <- function(this.extracts, this.sample.dir, mtx.files=c("barcodes.tsv", "genes.tsv", "matrix.mtx")) {
    ## check if MTX files were provided
    # mtx.files <- c("barcodes.tsv", "genes.tsv", "matrix.mtx")
    if (all(mtx.files %in% basename(this.extracts))) {
        ## read barcodes
        sample.barcodes.file <- this.extracts[grep("barcodes", this.extracts)]
        sample.barcodes <- fread(sample.barcodes.file, header=FALSE)

        ## read 10x experiment from sample dir
        sce <- read10xCounts(this.sample.dir)
        data <- list(sample.barcodes=sample.barcodes, sce=sce)
        return(data)
    }

    ## alternatively h5 files should have been provided
    if (grepl("h5", this.extracts)) {
        ## read barcodes
        barcode.path <- "matrix/barcodes"
        sample.barcodes <- h5read(this.extracts, barcode.path)
        sample.barcodes <- data.table(V1=sample.barcodes)

        ## read 10x experiment from sample dir
        sce <- read10xCounts(this.extracts)
        data <- list(sample.barcodes=sample.barcodes, sce=sce)
        return(data)
    }

    msg.txt <- "Files not recognised.
    Either barcodes.tsv/genes.tsv/matrix.mtx or an h5 file are required."
    stop(msg(error, msg.txt))
}

#' Identify non-empty droplets in Single-Cell RNA sequencing data
#'
#' @param sce SingleCellExperiment object containing count data.
#' @param sample Character string that specifies the sample identifier which
#'        will be prefixed to each barcode to create unique column names in the
#'        count matrix.
#' @param sample.barcodes Character vector containing barcode sequences.
#'
#' @return A character vector containing the names of droplets considered to
#' contain cells based on the `emptyDrops` analysis or all droplet names if the
#' matrix was pre-filtered.
#'
#' @details
#' The function retrieves the count matrix from the `sce` object. It modifies
#' the column names of the count matrix to incorporate the sample identifier.
#' The matrix is then subjected to the `emptyDrops` method which identifies
#' likely cell-containing droplets using ambient RNA levels as a reference.
#' If `emptyDrops` fails (likely due to a pre-filtered matrix), a warning is
#' issued, and all columns are considered as containing cells. See
#' https://support.bioconductor.org/p/123554/#123562 for this case.
#'
#' @importFrom DropletUtils emptyDrops
#' @importFrom stats set.seed
#' @importFrom methods tryCatch
#' @references
#' Lun ATL, McCarthy DJ, Marioni JC (2019). A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. F1000Research, 5:2122.
#' @seealso \code{\link[DropletUtils]{emptyDrops}}, \code{\link[S4Vectors]{DataFrame}}
remove.emptydrops <- function(sce, sample, sample.barcodes) {
    set.seed(100)
    my.count <- counts(sce)
    colnames(my.count) <- paste0(sample, "_", sample.barcodes)

    msg(info, "Starting emptyDrops")
    e.out <- tryCatch(
        emptyDrops(my.count),
        error=function(e) {
            warning("The matrix has been alrady pre-filtered.")
            return(NA)  ## Return NA on error
        })

    if (class(e.out)!="DFrame")  {
        true.cells <- colnames(my.count)
    } else {
        is.cell <- e.out$FDR <= 0.01
        is.cell[is.na(is.cell)] <- FALSE
        true.cells <- colnames(my.count)[is.cell]
    }
    msg(info, "Done.")
    return(true.cells)
}

#' Retrieve gene symbols from Ensembl IDs
#'
#' This function takes a set of Ensembl gene IDs and uses the biomaRt package to
#' retrieve the corresponding gene symbols from the specified Ensembl database.
#'
#' @param ensembl.ids Data.table, data.frame or character vector containing the
#'        Ensembl gene IDs. If data.table or data.drame then Ensembl ids should
#'        be provided in a column named `gene.id`.
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
#' @details The function requires the biomaRt package and expects that the
#'          ensembl object passed to it has been properly initialized with the
#'          useMart function. It uses the `getBM` function to perform the query.
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
    if ("data.table" %in% class(ensembl.ids)) {
        values <- ensembl.ids$gene.id
    } else {
        values <- ensembl.ids
    }

    gene.symbols <- setDT(getBM(attributes=c("ensembl_gene_id",
                                             "external_gene_name"),
                                filters="ensembl_gene_id",
                                values=values,
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

#' Helper function to map genes to Ensembl by gene id or gene symbol
map.ensembl <- function(genes.dt, ensembl) {
     ## map gene IDs to Ensembl, depending on the format
    ensembl.ids <- genes.dt[grep("ENSG", gene.id)]
    if (nrow(ensembl.ids)>0) {
        mapped.ids <- get.gene.symbols(ensembl.ids, ensembl)
    } else {
        mapped.ids <- data.table()
    }
    gene.symbols <- genes.dt[grep("ENSG", gene.id, invert=TRUE)]
    if (nrow(gene.symbols)>0) {
        mapped.symbol <- get.gene.ids(gene.symbols, ensembl)
    } else {
        mapped.symbol <- data.table()
    }

    mapped.genes <- rbind(mapped.ids, mapped.symbol)
    genes.dt[mapped.genes, on=c(gene.id="ensembl_gene_id"),
             gene.symbol:=external_gene_name]
    genes.dt[mapped.genes, on=c(gene.symbol="external_gene_name"),
             gene.id:=ensembl_gene_id]
    return(genes.dt)
}

#' Plot variable genes from Seurat object
#'
#' Visualizes variable genes on a scatter plot using data from a Seurat object.
#' This function highlights specific genes and saves the plot to a file.
#'
#' @param seurat.object A Seurat object containing variable gene data and
#'        expression metrics calculated, typically with functions such as
#'        `FindVariableFeatures`.
#' @param points Character vector with data points to label.
#' @param labels Character vector with labels.
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
plot.variable.genes <- function(seurat.object, points, labels, output.file) {
    p <- VariableFeaturePlot(seurat.object, pt.size=0.2)
    p2 <- LabelPoints(plot=p, points=points,
                      labels=labels, repel=TRUE,
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
    found.mt <- PercentageFeatureSet(seurat.obj, pattern="^MT-|^mt-")
    tryCatch({
        msg.txt <- "Requesting gene identifiers from Ensembl to map gene names"
        msg(info, msg.txt)
        ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        ensembl.dt <- setDT(getBM(attributes=c("ensembl_gene_id",
                                               "external_gene_name",
                                               "external_synonym"),
                                  mart=ensembl))
    }, error = function(e) {
        cat("Error requesting data: ", e$message, "\n")
        msg(info, "loading from file")
        ensembl.dt <- fread("ensemblgenes_latest.csv", select=1:2, 
                            col.names=c("ensembl_gene_id", 
                                        "external_gene_name"))
        ## TODO: can't find this, needs fixing
        ensembl.dt[, external_synonym := external_gene_name]
    })
    
    mt.genes <- ensembl.dt[grep("^MT", external_gene_name)]
    all.genes <- rownames(seurat.obj)
    genes <- mt.genes[external_gene_name %in% all.genes, 
                      external_gene_name]
    if (length(genes) > 0) {
        matched.mt <- PercentageFeatureSet(seurat.obj, features=genes)
        all.mt <- found.mt + matched.mt
    } else {
        all.mt <- found.mt
    }

    seurat.obj[["percent.mt"]] <- all.mt

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
#' @param ensembl.dt Data.table containing gene identifiers and gene symbols
#'        from Ensembl. It is used to map transcripts to official gene names.
#' @param doublet.method.scDblFinder a flag to use scDblFinder for doublet 
#'        detection, srcublet is used if FALSE
#'
#' @return A Seurat object with additional metadata fields and QC metrics.
#'
#' @details The function assumes the presence of 'barcodes.tsv' for naming
#'          conventions in cases of .mtx files and
#'          additionally handles doublet detection and mitochondrial content as
#'          part of quality control.
#'
#' @examples
#' ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#' ensembl.dt <- setDT(getBM(attributes=c("ensembl_gene_id",
#'                                        "external_gene_name",
#'                                        "external_synonym"),
#'                           mart=ensembl))
#' seurat_obj <- create.seurat(counts, this.sample="GSM4041169",
#'                             this.sample.dir="GSM4041169",
#'                             ensembl.dt)
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom data.table fread, setDT
#' @importFrom ggplot2 ggplot, geom_violin, geom_jitter, geom_point, labs, theme_minimal
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggsave
create.seurat <- function(true.cells, this.sample, this.sample.dir, ensembl.dt, 
                          doublet.method.scDblFinder=TRUE) {
    ## Studies may have sample data in mtx or h5 format which are to be
    ## read by read10xCounts and Read10X_h5 respectively.
    ## Thus, we handle both cases here conditional on file types.
    sample.files <- list.files(this.sample.dir)
    if (any(grepl("mtx", sample.files))) {
        pre <- read10xCounts(this.sample.dir, col.names=TRUE)
        pre <- counts(pre)
    } else if (any(grepl("h5", sample.files))) {
        this.h5.file <- sample.files[grep("*.h5", sample.files)]
        pre <- Read10X_h5(file.path(this.sample.dir, this.h5.file),
                          use.names=TRUE, unique.features=TRUE)
    } else {
        msg.txt <- "Either .mtx or .h5 file should be present.
        Check input directory."
        stop(msg(error, msg.txt))
    }
    ## filter 10x counts to keep only true cells
    this.true.cells <- remove.sample.from.barcode(true.cells)
    pre <- pre[, colnames(pre) %in% this.true.cells]
    ## Create Seurat object with filtered counts
    pre <- CreateSeuratObject(counts=pre, project=this.sample)

    msg.txt <- sprintf("Mapping transcripts to official gene symbols")
    msg(info, msg.txt)

    transcripts <- data.table(transcript=rownames(pre), gene.symbol="NA")
    transcripts[ensembl.dt, on=c(transcript="external_synonym"),
                gene.symbol := external_synonym]
    transcripts[ensembl.dt, on=c(transcript="external_gene_name"),
                gene.symbol := external_gene_name]
    transcripts[ensembl.dt, on=c(transcript="ensembl_gene_id"),
                gene.symbol := external_gene_name]
    msg.txt <- sprintf("%s transcripts were not matched to Ensembl",
                       transcripts[gene.symbol=="", .N])
    msg(info, msg.txt)
    msg(info, "Using default study identifiers for these.")
    transcripts[gene.symbol=="" | gene.symbol=="NA", gene.symbol := transcript]

    rownames(pre[["RNA"]]@counts) <- transcripts$gene.symbol
    rownames(pre[["RNA"]]@data) <- transcripts$gene.symbol
    pre <- UpdateSeuratObject(object=pre)

    msg(info, "Processing Seurat meta.data.")
    ## pass SC index as <sample_id>_<barcode>
    pre$barcodes <- true.cells
    pre$sample <- this.sample
    pre$orig.ident <- this.sample

    ## update Seurat meta.data
    pre <- update.seurat.meta.data(pre)

    ## read doublet predictions and scores
    if (doublet.method.scDblFinder){
        prediction.file <- file.path(this.sample.dir,
                                     "scDblFinder_EDR0.008_PredictedDoublets.csv")
        doublet.file <- file.path(this.sample.dir,
                                 "scDblFinder_EDR0.008_DoubletScores.csv")
    } else {
        prediction.file <- file.path(this.sample.dir,
                                     "scrublet_EDR0.008_PredictedDoublets.csv")
        doublet.file <- file.path(this.sample.dir,
                                 "scrublet_EDR0.008_DoubletScores.csv")
        ## a given cell is a duplet if doublet_score > threshold.
        ## original threshold was 0.35, but we set a more stringent threshold
        ## of 0.15 which corresponds to the threshold for doublet detection
        ## from step 2.
    }
    doublet_prediction <- fread(prediction.file)
    doublet_score <- fread(doublet.file)
    doublets <- data.table(prediction=doublet_prediction$V1,
                           score=doublet_score$V1)
    if (doublet.method.scDblFinder){
        pre[["doublet_prediction"]] <- doublets[, prediction]
    } else {
        doublets[, new.prediction := ifelse(score>doublet.threshold, 1, 0)]

        ## append info on duplets to Seurat meta.data
        pre[["doublet_prediction"]] <- doublets[, new.prediction]
    }
    pre[["doublet_score"]] <- doublets[, score]

    ## save created Seurat object
    seurat.file <- file.path(this.sample.dir, "seurat.RDS")
    if (!file.exists(seurat.file)) {
        msg.txt <- sprintf("Saving Seurat object to %s", seurat.file)
        msg(info, msg.txt)
        saveRDS(pre, file=seurat.file)
    }

    msg.txt <- sprintf("Saving diagnostic plots to %s", this.sample.dir)
    msg(info, msg.txt)
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

#' Process samples, convert to Seurat objects and merge
#'
#' This function reads individual sample data from specified directories and
#' file ormats, creates Seurat objects for each sample, and merges them into a
#' single Seurat object. The function supports both "mtx" and "h5" file formats.
#'
#' @param meta Data.table containing metadata for each sample.
#'             The metadata table must include the columns `sample`,
#'             `sample.dir`, and potentially `file.type` for h5 files. Each row
#'             should correspond to a sample with its directory path where the
#'             data files are stored.
#' @param true.cells Character vector with ids of true cells identified by
#'        running empty drops algorythm.
#' @param output.dir Full system path to the directory where merged Seurat file
#'        will be saved.
#'
#' @return A merged Seurat object containing data from all samples.
#'
#' @details The function processes each sample based on metadata entries. It
#'          checks for the presence of either "mtx" or "h5" files in the
#'          directory specified by `sample.dir`. Depending on the file type, it
#'          reads the data using appropriate methods (`read10xCounts` for "mtx"
#'          and `Read10X_h5` for "h5"). Each dataset is converted into a Seurat
#'          object which are then merged into one Seurat object. The function
#'          requires the Seurat, data.table, and DropletUtils packages.
#'
#'          It is assumed that all required libraries are installed and loaded.
#'          Error handling includes checks for valid directory paths, correct
#'          sample naming in metadata, and the presence of necessary data files.
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
process.samples.and.merge <- function(meta, true.cells, output.dir) {
    ## read and return the Merged Seurat object if exists
    seurat.file <- file.path(output.dir, "merged.seurat.RDS")
    if (file.exists(seurat.file)) {
        msg.txt <- "Read the detected merged.seurat.RDS file"
        msg(warn, msg.txt)
        Merge <- readRDS(seurat.file)
        return(Merge)
    }
    
    tryCatch({
        msg.txt <- "Requesting gene identifiers from Ensembl to map gene names"
        msg(info, msg.txt)
        ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        ensembl.dt <- setDT(getBM(attributes=c("ensembl_gene_id",
                                               "external_gene_name",
                                               "external_synonym"),
                                  mart=ensembl))
    }, error = function(e) {
        cat("Error requesting data: ", e$message, "\n")
        msg(info, "loading from file")
        ensembl.dt <- fread("ensemblgenes_latest.csv", select=1:2, 
                            col.names=c("ensembl_gene_id", 
                                        "external_gene_name"))
        ## TODO: can't find this, needs fixing
        ensembl.dt[, external_synonym := external_gene_name]
        return(ensembl.dt)
    })
    # browser()
      
    seurat.list <- list()
    all.samples <- meta[, unique(sample)]

    ## Loop through samples and create Seurat objects
    for (this.sample in all.samples) {
        this.meta <- meta[sample == eval(this.sample)][1]
        msg.txt <- sprintf("Processing sample: %s, %s of %s samples...",
                           this.sample,
                           which(all.samples==this.sample),
                           length(all.samples))
        msg(info, msg.txt)

        this.sample.dir <- this.meta[, sample.dir]
        if (!dir.exists(this.sample.dir)) {
            msg.txt <- "Sample directory listed in metadata does not exist."
            stop(msg(error, msg.txt))
        }
        
        seurat.list[[this.sample]] <- create.seurat(true.cells, this.sample,
                                                    this.sample.dir, ensembl.dt)
        msg(info, "Done.\n")
    }

    ## Merge all Seurat objects into one
    if (length(seurat.list) > 1) {
        Merge <- merge(x=seurat.list[[1]], y=seurat.list[2:length(seurat.list)],
                       add.cell.ids=all.samples)
    } else {
        Merge <- seurat.list[[1]]
    }

    ## update and check the meta.data of the merged Seurat object
    Merge <- update.seurat.meta.data(Merge)

    ## Optional: Add clinical metadata for samples (if any)
    ## Assuming clinical metadata needs to be merged
    ## Merge <- AddClinicalData(Merge, clinical_data)

    ## Save the merged Seurat object
    saveRDS(Merge, file=seurat.file)

    ## Return the merged Seurat object
    return(Merge)
}

#' Update metadata in a Seurat Object
#'
#' This function recalculates metadata for a Seurat object, specifically
#' updating `nCount_RNA`, `nFeature_RNA`, and computing the log10 ratio of
#' features to total counts. It also identifies cells with zero counts and adds
#' mitochondrial content to the metadata.
#'
#' @param seurat.obj Seurat object with single-cell RNA sequence data.
#'
#' @return Seurat object with updated metadata.
#'
#' @importFrom Matrix colSums
update.seurat.meta.data <- function(seurat.obj) {
    ## recompute nCount_RNA and nFeature_RNA
    n.counts <- Matrix::colSums(seurat.obj[["RNA"]]@counts)
    if (any(n.counts==0)) {
        msg.txt <- sprintf("%s cells with 0 counts detected.",
                           sum(n.counts==0))
        msg(warn, msg.txt)
        zero.counts <- names(n.counts[n.counts==0])
        print(ifelse(length(zero.counts)>20, zero.counts[1:20], zero.counts))
    }
    seurat.obj@meta.data$nCount_RNA <- n.counts

    n.features <- Matrix::colSums(seurat.obj[["RNA"]]@counts > 0)
    if (any(n.features==0)) {
        msg.txt <- sprintf("%s cells with 0 features detected.",
                           sum(n.features==0))
        msg(warn, msg.txt)
        zero.counts <- names(n.features[n.features==0])
        print(ifelse(length(zero.counts)>20, zero.counts[1:20], zero.counts))
    }
    seurat.obj@meta.data$nFeature_RNA <- n.features

    ## append mitochondrial RNA to Seurat meta.data object
    ## Temporary!
    seurat.obj <- count.mt.genes(seurat.obj)

    ## compute log10 of the proportion of identified genes to counts of all
    ## RNA molecules in each cell to identify anomalies in gene expression due
    ## to mitochondrial or ribosomal RNA.
    seurat.obj@meta.data$log10GenesPerUMI <-
        log10(seurat.obj@meta.data$nFeature_RNA) /
        log10(seurat.obj@meta.data$nCount_RNA)
    return(seurat.obj)
}

#' Filter cells and genes in a Seurat Object
#'
#' This function performs multiple quality control (QC) checks to filter cells
#' in a Seurat object based on gene expression, mitochondrial content, and
#' other metadata-based criteria. It also filteres out genes expressed in small
#' numbers of cells and generates a diagnostic plot to visually assess the
#' quality of the data after filtering.
#'
#' @param seurat.object Seurat object with single-cell RNA sequence data.
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
#' @return A filtered Seurat object with cells and genes that pass filtering.
#'
#' @examples
#' # Filter with default parameters
#' seurat <- qc.seurat(my.seurat.object, output.dir=getwd())
#'
#' # Filter with custom parameters
#' seurat <- qc.seurat(my.seurat.object,
#'                     output.dir=getwd(),
#'                     min.features=300)
#'
#' @require Seurat
#' @require data.table
#' @require ggplot2
#' @require biomaRt
qc.seurat <- function(seurat.object, output.dir,
                      min.features=500, min.counts=1000, max.mt=10,
                      n.deviation=3, n.cells=30) {
    require(Seurat)
    require(data.table)
    require(ggplot2)
    require(biomaRt)

    msg.txt <- "Removing cells with outlier transcripts/cell counts"
    msg(info, msg.txt)

    ## ensure the meta.data is up-to-date
    seurat.object <- update.seurat.meta.data(seurat.object)

    meta.data <- setDT(seurat.object@meta.data)
    ## check cell barcodes against Seurat object
    ## TODO: rewrite this function to apply to multiple samples
    if (!all(meta.data$barcodes %in% colnames(seurat.object))) {
        stop("Cells were not matched to Seurat object. Check barcodes.")
    }

    ## infer cells with outlier transcripts/cell counts
    meta.data[, count.threshold :=
                median(nCount_RNA) + n.deviation * mad(nCount_RNA), by="sample"]
    meta.data[, feature.threshold :=
                median(nFeature_RNA) + n.deviation * mad(nFeature_RNA),
              by="sample"]
    qc.expr <- quote(nCount_RNA > count.threshold |
                     nFeature_RNA > feature.threshold)
    meta.data[, qc.fail := eval(qc.expr)]
    msg.txt <- sprintf("Remove %s cells in %s samples -- outlier transcripts.",
                        meta.data[qc.fail==TRUE, .N],
                        meta.data[qc.fail==TRUE, uniqueN(sample)])
    msg(info, msg.txt)
    meta.data[, `:=`(count.threshold=NULL, feature.threshold=NULL)]

    ## infer which cells fail user specified QC
    qc.expr <- quote(nFeature_RNA < min.features | nCount_RNA < min.counts |
                     percent.mt > max.mt | doublet_prediction != 0)
    meta.data[eval(qc.expr), qc.fail := TRUE]

    msg.txt <- sprintf("Remove %s cells in %s samples -- fail supplied QC.",
                       meta.data[eval(qc.expr), .N],
                       meta.data[eval(qc.expr), uniqueN(sample)])
    msg(info, msg.txt)

    cells.to.keep <- meta.data[qc.fail==FALSE, barcodes]
    msg.txt <- sprintf("Keep %s cells after all filtering steps",
                       length(cells.to.keep))
    msg(info, msg.txt)

    ## filter genes based on how many cells are they expressed in
    counts <- GetAssayData(object=seurat.object, slot="counts")
    ## infer lowly expressed genes on a subset of cell that we keep
    counts <- counts[, cells.to.keep, drop=FALSE]
    nonzero.counts <- counts > 0
    genes.expressed <- Matrix::rowSums(nonzero.counts) >= n.cells
    genes.to.keep <- names(genes.expressed[genes.expressed])

    msg.txt <- sprintf("Keep %s genes which are expressed in >%s cells.",
                       length(genes.to.keep), n.cells)
    msg(info, msg.txt)

    ## subset Seurat objects and check its correctness
    seurat.object <- seurat.object[genes.to.keep, cells.to.keep]
    counts <- GetAssayData(object=seurat.object, slot="counts")
    row.sums <- Matrix::rowSums(counts)
    col.sums <- Matrix::colSums(counts)
    msg.txt <- "Summary of gene counts after QC"
    msg(info, msg.txt)
    print(summary(row.sums))
    msg.txt <- "Summary of cell counts after QC"
    msg(info, msg.txt)
    print(summary(col.sums))

    df <- as.data.frame(meta.data[qc.fail==FALSE])
    ## consistency checks
    if (any(is.na(df))) {
        stop("NAs detected in meta.data.")
    }
    if (!identical(colnames(seurat.object), df$barcode)) {
        stop("Cell names in counts matrix do not match Seurat metadata.")
    }

    seurat.object@meta.data <- df

    ## update the meta.data after subsetting
    seurat.object <- update.seurat.meta.data(seurat.object)

    ## plot the counts after filtering
    plot.features(df, output.dir)

    saveRDS(seurat.object, file=file.path(output.dir, "filtered.seurat.RDS"))

    return(seurat.object)
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

    ## connect to Ensembl via BioMart
    ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

    ## obtain top variable genes
    msg.txt <- sprintf("Selected %s most variable genes", n.top.variable.genes)
    msg(info, msg.txt)

    top.variable.genes <- head(VariableFeatures(seurat.object),
                               n.top.variable.genes)
    top.variable.genes <- data.table(gene.id=top.variable.genes,
                                     gene.symbol=top.variable.genes)

    ## map gene IDs to Ensembl, depending on the format
    top.variable.genes <- map.ensembl(top.variable.genes, ensembl)
    ensembl.ids <- top.variable.genes[grep("ENSG", gene.id)]

    msg.txt <- sprintf("%s most variable genes mapped to Ensembl",
                       n.top.variable.genes)
    msg(info, msg.txt)
    print(top.variable.genes)

    msg.txt <- "Plot most variable genes on Expression-Variance plot"
    msg(info, msg.txt)
    output.file <- file.path(output.dir, "most.variable.genes.png")
    plot.variable.genes(seurat.object, points=top.variable.genes$gene.symbol,
                        labels=top.variable.genes$gene.symbol, output.file)

    ## process core genes if specified
    if (!is.null(core.genes)) {
        msg.txt <- sprintf("%s genes of interest were specified",
                           length(core.genes))
        msg(info, msg.txt)
        selected.genes <- data.table(gene.id=core.genes, gene.symbol=core.genes)
        selected.genes <- map.ensembl(selected.genes, ensembl)
        selected.genes[gene.symbol %in% rownames(seurat.object),
                       matched.symbol := TRUE]

        ## append labeled genes if specified
        if (!is.null(label.genes)) {
            msg.txt <- sprintf("%s additional genes to label",
                               length(label.genes))
            msg(info, msg.txt)
            label.genes.dt <- data.table(gene.id=label.genes,
                                         gene.symbol=label.genes)
            label.genes.dt <- map.ensembl(label.genes.dt, ensembl)
            label.genes.dt[, matched.symbol := FALSE]
            label.genes.dt[gene.symbol %in% rownames(seurat.object),
                           matched.symbol := TRUE]
            selected.genes <- rbind(selected.genes, label.genes.dt)
        }

        ## plot core genes and labelled genes
        selected.genes <- unique(selected.genes[matched.symbol==TRUE],
                                 by="gene.symbol")

        msg.txt <- "Map HVGs to Ensembl"
        msg(info, msg.txt)
        all.hvg <- head(VariableFeatures(seurat.object), 1000)
        all.hvg <- data.table(gene.id=all.hvg, gene.symbol=all.hvg)
        all.hvg <- map.ensembl(all.hvg, ensembl)

        selected.genes[, match.hvg := FALSE]
        selected.genes[all.hvg, on="gene.id", match.hvg := TRUE]
        selected.genes[all.hvg, on="gene.symbol", match.hvg := TRUE]
        selected.genes[label.genes.dt, on="gene.id", match.hvg := TRUE]

        msg.txt <- sprintf("Label %s core genes on 1000 HVGs",
                           selected.genes[match.hvg==TRUE, .N])
        msg(info, msg.txt)
        output.file <- file.path(output.dir, "selected.variable.genes.png")
        plot.variable.genes(seurat.object,
                            points=selected.genes[match.hvg==TRUE, gene.symbol],
                            labels=selected.genes[match.hvg==TRUE, gene.symbol],
                            output.file)
    }
    return(seurat.object)
}

cluster.cells <- function(seurat.object, cell.map, output.dir) {
    seurat.object <- NormalizeData(seurat.object,
                                   normalization.method="LogNormalize",
                                   scale.factor=10000)
    seurat.object <- FindVariableFeatures(seurat.object, selection.method="vst",
                                          nfeatures=2000)
    all.genes <- rownames(seurat.object)
    seurat.object <- ScaleData(seurat.object, features=all.genes)
    seurat.object <- RunPCA(seurat.object,
                            features=VariableFeatures(object=seurat.object))
    seurat.object <- FindNeighbors(seurat.object, dims=1:10)
    seurat.object <- FindClusters(seurat.object, resolution=0.5)
    # seurat.object <- RunUMAP(seurat.object, dims=1:10)
    seurat.object@meta.data$condition <-
        cell.map[match(colnames(seurat.object), barcode), as.factor(condition)]
    p <- DimPlot(seurat.object, reduction="pca", group.by="condition",
                 label=TRUE, combine=FALSE)
    ggsave(p, file=file.path(output.dir, "clusters.png"))
}

sample.seurat <- function(seurat.object, output.dir, cell.map, n.cells.keep) {
    ## map cell counts to samples
    cell.counts <- as.data.table(table(Idents(seurat.object)))
    setnames(cell.counts, old=c("V1", "N"), new=c("sample", "cell.count"))

    ## join with cell.map to bring in conditions
    cell.map <- cell.map[cell.counts, on="sample"]
    cell.map[, cell.count.by.condition := .N, by=condition]

    msg.txt <- "Cells in Seurat object were mapped to the following conditions."
    msg(info, msg.txt)
    print(cell.map[, .(N.cells=unique(cell.count.by.condition)), by=condition])

    ## sample size for each condition to ensure equal representation in the
    ## final sample
    min.cells <- cell.map[, n.cells.keep / uniqueN(condition)]
    cell.map[, sample.size := ifelse(cell.count.by.condition < min.cells,
                                     cell.count.by.condition, min.cells)]

    ## randomly sample cells ensuring proportional representation of samples
    msg.txt <- sprintf("Sampling %s cells from each of %s conditions.",
                       min.cells, n.cells.keep)
    msg(info, msg.txt)
    msg.txt <- sprintf("Default sample: %s cells, smallest sample: %s cells",
                       min.cells,
                       cell.map[, min(cell.count.by.condition)])
    msg(info, msg.txt)
    msg.txt <- sprintf("In total, %s of %s requested cells will be sampled",
                       cell.map[, unique(sample.size), by=condition][, sum(V1)],
                       n.cells.keep)
    msg(info, msg.txt)

    cells.to.keep <- NULL
    for (this.condition in cell.map[, unique(condition)]) {
        this.cell.map <- cell.map[condition==eval(this.condition)]
        this.sample <- this.cell.map[, sample(barcode, size=unique(sample.size),
                                              replace=FALSE)]
        cells.to.keep <- c(cells.to.keep, this.sample)
    }

    ## subset Seurat object
    dt <- setDT(seurat.object@meta.data)
    seurat.object.filtered <- seurat.object[, cells.to.keep]
    dt <- dt[match(colnames(seurat.object.filtered), barcodes)]

    if (!dt[, identical(barcodes, colnames(seurat.object.filtered))])
        stop(msg(error, "Barcodes in counts matrix do not match meta.data."))

    seurat.object.filtered@meta.data <- setDF(dt)

    seurat.object.filtered <- update.seurat.meta.data(seurat.object.filtered)

    ## run some checks to ensure consistency between meta data and Seurat object
    if (any(is.na(seurat.object.filtered@meta.data)))
        stop(msg(error, "NAs detected in meta.data."))
    if (any(!cells.to.keep %in% seurat.object.filtered@meta.data$barcode)) {
        msg.txt <- "Cells labeled for removal are in filtered Seurat object."
        stop(msg(error))
    }

    seurat.for.stator.file <- file.path(output.dir, "seurat.for.stator.RDS")
    saveRDS(seurat.object.filtered, file=seurat.for.stator.file)

    return(seurat.object.filtered)
}

#' Helper function to create a feature plot for the metadata from Seurat object
plot.features <- function(df, output.dir) {
    p <- ggplot(df, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
        geom_point() +
        scale_colour_gradient(low="gray90", high="black") +
        stat_smooth(method=lm) +
        scale_x_log10() +
        scale_y_log10() +
        labs(x="RNA counts, Log10", y="N transcripts, Log10")
    plot.file <- file.path(output.dir, "feature.plot.png")
    ggsave(plot.file, plot=p)
}

#' Helper function to filter out genes which cannot be mappes to Ensembl
filter.ensembl <- function(seurat.object) {
        ## filer all genes which cannot be matched to Ensembl
    ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    genes <- rownames(seurat.object)
    ensembl.ids <- genes[grep("^ENSG", genes)]
    gene.symbols <- genes[!genes %in% ensembl.ids]

    msg.txt <- sprintf("%s Ensembl ids and %s gene symbols detected.",
                       length(ensembl.ids), length(gene.symbols))
    msg(info, msg.txt)

    if (length(gene.symbols) != 0) {
        msg(info, "Converting gene symbols to Ensembl ids.")
        gene.symbols <- data.table(gene.id=NA, gene.symbol=gene.symbols)
        matched.ids <- get.gene.ids(gene.symbols, ensembl)
    }

    all.ensembl.ids <- c(ensembl.ids, matched.ids$ensembl_gene_id)
    msg.txt <- sprintf("%s genes with Ensembl ids will be retained.",
                       length(all.ensembl.ids))
    msg(info, msg.txt)
    seurat.object <- seurat.object[all.ensembl.ids, ]
}

#' Helper function to indentify doublets with scDblFinder
find.doublets.scDblFinder <- function(sce, this.sample.dir) {
    this.sce <- scDblFinder(sce)
    doublet_scores <- this.sce$scDblFinder.score
    doublet_classification <- this.sce$scDblFinder.class
    doublet.dt <- data.table(prediction=doublet_classification, score=doublet_scores)
    doublet.dt[, is.doublet := ifelse(prediction=="doublet", 1, 0)]
    prediction.file <- file.path(this.sample.dir,
                                 "scDblFinder_EDR0.008_PredictedDoublets.csv")
    doublet.file <- file.path(this.sample.dir,
                             "scDblFinder_EDR0.008_DoubletScores.csv")
    fwrite(doublet.dt[, .(is.doublet)], file=prediction.file, col.names=FALSE)
    fwrite(doublet.dt[, .(score)], file=doublet.file, col.names=FALSE)
}

## Function to remove samples names from barcodes
remove.sample.from.barcode <- function(obj, del="_"){
    re <- paste0(".*", del, "(.*)")
    this.obj <- gsub(re, "\\1", obj)
    return(this.obj)
}