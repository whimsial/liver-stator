seurat.list <- list()
all.samples <- meta[, unique(sample)]
for (this.sample in all.samples) {
    this.meta <- meta[sample==eval(this.sample)][1]
    cat("Processing sample:", this.sample, "\n")
    seurat.list[this.sample] <- create.seurat(this.sample, this.meta)
}
Merge <- merge(x=seurat.list[[1]], y=seurat.list[2:length(seurat.list)],
               add.cell.ids=all.samples)

## add clinical metadata for samples (if any)

## save the object
seurat.file <- unique(meta$surat.file)
save(Merge, file=seurat.file)

## cleanup to save some memory
rm(seurat.list)
gc()
