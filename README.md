# liver-stator

Analysis of scRNA-seq from healthy and pathologic livers using [Stator](https://github.com/AJnsm/Stator).

The main analysis pipeline is implemented in [`pipeline.R`](https://github.com/whimsial/liver-stator/blob/main/pipeline.R). It also includes the QC steps written as a collection of R functions in [`rnaseq.functinos.R`](https://github.com/whimsial/liver-stator/blob/main/rnaseq.functions.R) to download and process the raw data from 
[Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) followed by the QC steps described below. 

## Pipeline

Depends on:

- R packages: `data.table`, `Seurat`, `DropletUtils`, `HDF5Array`, `biomaRt` (for full list see [`properties.R`](https://github.com/whimsial/liver-stator/blob/main/properties.R)).

R packages should be installed user R library (see [`properties.R`](https://github.com/whimsial/liver-stator/blob/main/properties.R) for details).


- Python packages: `scanpy`, `scrublet` (for full list see [`doublet.py`](https://github.com/whimsial/liver-stator/blob/main/doublet.py) and [`requirements.txt`](https://github.com/whimsial/liver-stator/blob/main/requirements.txt)).

The Python code can (and probably should) also be run interactively using Jupyter notebook [`doublets.ipynb`](https://github.com/whimsial/liver-stator/blob/main/doublets.ipynb). It is useful to inspect distributions of doublets by eye to come up with a suitable threshold to be used in [`doublet.py`](https://github.com/whimsial/liver-stator/blob/main/doublet.py) (default is 0.15).

I found it useful to run Jupyter notebook inside Python virtual environment on the HPC and connect to it via SSH tunnel:

- start SSH tunnel to HPC mapping one of the open ports (I use port 9999)
- clone this repository on HPC and navigate to it
- start python virtual environment: `python3 -m venv venv`
- activate python virtual env: `source venv/bin/activate`
- install dependencies: `pip install -r requirements.txt` and check everything has been installed: `pip list`
- run the notebook server: `jupyter notebook --no-browser --port=9999`
- connect to the server from a local browser by pasting in the URL printed by the above command

### Step 1: detection of empty drops

Retrieve the count matrix from the Single Cell Experiment object and pass to  the `emptyDrops` method which identifies likely cell-containing droplets using ambient RNA levels as a reference.

If `emptyDrops` fails (likely due to a pre-filtered matrix), a warning will be issued, and all columns will be considered as containing cells.

### Step 2: detection of doublets

The key operation in this step is a call to `scrub_doublets` function from the `scrublet` package which is designed to detect doublets (also known as multiplets) in single-cell RNA sequencing (scRNA-seq) data. Doublets occur when two or more cells are inadvertently sequenced together as a single cell. Identifying and removing doublets is crucial for accurate downstream analysis because they can introduce significant noise and bias into the results.

The `scrub_doublets` function processes scRNA-seq data and returns **doublet scores** for each cell, consisting of:

- a continuous score that represents the likelihood that each cell (or transcriptomic profile) is a doublet. The score is typically a value between 0 and 1, where a higher score indicates a higher probability of the cell being a doublet.
- this function usually also involves setting a threshold score to classify cells as singlets or doublets. This threshold might be determined automatically by the algorithm based on the distribution of doublet scores across all cells, or it can be manually set by the user.

The identification of the doublets cutoff should be performed interactively using Jupyter notebook.  

Having identified the threshold, [`doublet.py`](https://github.com/whimsial/liver-stator/blob/main/doublet.py) can be run for all samples passing in the threshold as a command line argument. An example of running this script from R can be found in [`pipeline.R`](https://github.com/whimsial/liver-stator/blob/main/pipeline.R). 


### Step 3: create a Seurat object from SingleCellExperiment data

For each sample, we read single-cell data from either .mtx or .h5 files located
in a specified directory, create a Seurat object, and enrich it with
metadata including doublet predictions, mitochondrial RNA content, and produce QC plots.

We then proceed by merging the Seurat objects for all samples within specified studies to generate a single analysis-wide Seurat object.

### Step 4: filter out cells/genes not passing QC

In this step, we perform multiple quality control (QC) checks to filter cells
in a Seurat object based on gene expression, mitochondrial content, and
other metadata-based criteria. It also generates a diagnostic plot to
visually assess the quality of the data after filtering.

### Step 5: process and visualize variable genes in a Seurat object

Perform a series of operations on a Seurat object to identify
and visualize highly variable genes. Normalize the data, identify variable
features, map gene IDs to Ensembl, and plot these genes on mean expression vs variance plot. 

In this step we also add and label core genes from GATE analysis and check if they appear as highly variable genes in single cell data.

Finally, we write the counts matrix and a list of highly variable genes to files in a format expected by Stator.


## Known issues

While running this analysis I encountered several issues with the newest release of Seurat 5 R package. These are to do with the new layers introduced to the standard Seurat objects. Downgrading to version 4.4.0 together with seurat-object 4.1.4 solved these issues for me and thus I recommend to run this pipeline with these versions. 

In time I will open an issue on Seurat's GitHub to ask for help with Seurat 5.
