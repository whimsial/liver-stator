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

## for each sample read 10X experiment, create Seurat object, append dublet
## info, and create descriptive plots
create.seurat <- function(this.sample="GSM4041150", meta) {
    this.meta <- meta[sample==eval(this.sample)]
    this.sample.dir <- this.meta[, unique(sample.dir)]

    ## read Dublet predictions and scores
    prediction.file <- file.path(this.sample.dir,
                                 "scrublet_EDR0.008_PredictedDoublets.csv")
    dublet.file <- file.path(this.sample.dir,
                             "scrublet_EDR0.008_DoubletScores.csv")
    doublet_prediction <- fread(prediction.file)
    doublet_score <- fread(dublet.file)

    ## Ramachandran et al have data for SingleCellExperiment consisting of
    ## matrix.mtx, genes.tsv and barcodes.tsv.
    ## It should be read using read10xCounts and then the sparse matrix of
    ## counts should be extracted and passed to Seurat
    if (study$name=="Ramachandran") {
        pre <- read10xCounts(this.sample.dir)
        pre <- counts(pre)
    }

    ## Guilliams et al report data in h5 format. These files should be read with
    ## Read10X_h5 function.
    if (study$name=="Guilliams") {
        this.h5.file <- this.meta[, Name]
        pre <- Read10X_h5(file.path(this.sample.dir, this.h5.file),
                          use.names=TRUE, unique.features=TRUE)
    }
    pre <- Seurat::CreateSeuratObject(counts=pre, project=this.sample)

    ## create SC index as <sample_id>:<cell_id>
    pre$barcodes <- rownames(pre@meta.data)
    pre$sample <- this.sample
    scvelo_index <- paste0(pre$sample, ":", pre$barcodes)
    scvelo_index <- gsub(pattern='-1', replacement="x", scvelo_index)

    ## different samples may have the same barcode,
    ## so add the name of the sample to the barcode
    pre <- RenameCells(object=pre, new.names=scvelo_index)

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

    ## plot scatter plots of nCount_RNA against percent.mt and nCount_RNA against
    ## nFeature_RNA
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
