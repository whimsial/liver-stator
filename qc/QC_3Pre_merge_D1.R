process.samples.and.merge <- function(meta) {
    seurat.list <- list()
    all.samples <- meta[, unique(sample)]

    ## Loop through each sample and create a Seurat object
    for (this.sample in all.samples) {
        this.meta <- meta[sample == eval(this.sample)][1]
        cat("Processing sample:", this.sample, "\n")
        seurat.list[[this.sample]] <- create.seurat(this.sample, this.meta)
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
    seurat.file <- unique(meta$surat.file)
    save(Merge, file=seurat.file)

    ## Return the merged Seurat object
    return(Merge)
}
