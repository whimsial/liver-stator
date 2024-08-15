# liver-stator

Analysis of scRNA-seq from healthy and pathologic livers using Stator


## QC pipeline

### Step 1: detection of empty drops

### Step 2: detection of duplets

This step is implemented in a Jupyter notebook and relies on `scrub_doublets` function from the `scrublet` package which is designed to detect doublets (also known as multiplets) in single-cell RNA sequencing (scRNA-seq) data. Doublets occur when two or more cells are inadvertently sequenced together as a single cell. Identifying and removing doublets is crucial for accurate downstream analysis because they can introduce significant noise and bias into the results.

The `scrub_doublets` function processes scRNA-seq data and returns **doublet scores** for each cell, consisting of:

- a continuous score that represents the likelihood that each cell (or transcriptomic profile) is a doublet. The score is typically a value between 0 and 1, where a higher score indicates a higher probability of the cell being a doublet.
- this function usually also involves setting a threshold score to classify cells as singlets or doublets. This threshold might be determined automatically by the algorithm based on the distribution of doublet scores across all cells, or it can be manually set by the user.

### Step 3: create a Seurat object from SingleCellExperiment data

For each sample, we read single-cell data from either .mtx or .h5 files located
in a specified directory, creates a Seurat object, and enriches it with
metadata including doublet predictions, mitochondrial RNA content, and
additional QC plots.

We then proceed by merging the Seurat objects for all samples within each study to generate a study-wide Seurat object.
