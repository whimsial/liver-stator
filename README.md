# liver-stator

Analysis of scRNA-seq from healthy and pathologic livers using [Stator](https://github.com/AJnsm/Stator).

The main analysis pipeline is implemented in [`pipeline.R`](https://github.com/whimsial/liver-stator/blob/main/pipeline.R), for a newer version of the pipeline developed for ME analysis refer to [`pipeline.me.R`](https://github.com/whimsial/liver-stator/blob/main/pipeline.me.R). It also includes the QC steps written as a collection of R functions in [`rnaseq.functinos.R`](https://github.com/whimsial/liver-stator/blob/main/rnaseq.functions.R) to download and process the raw data from
[Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) followed by the QC steps described below.

Alternatively, if data integration between multiple studies is not needed a
simplified example script [`example.R`](https://github.com/whimsial/liver-stator/blob/main/example.R) can be used.

## Example script

To run the example script on HPC cluster (e.g. EDDIE):

1. login to wildwest1 (ssh node2c15).
2. clone this repository and navigate to it (all subsequent steps should be run from the
root of this repository).
2. in the shell session run [`setup.sh`](https://github.com/whimsial/liver-stator/blob/main/setup.sh) to load R and prepare virtual environment for Python.
3. then start R, and source [`example.R`](https://github.com/whimsial/liver-stator/blob/main/example.R) but first modify it according to your analysis.

**Note**: You will need to install R dependencies so set `install.dependencies<-TRUE` in
the `example.R` script. This should be done once. The script will then source [`dependencies.R`](https://github.com/whimsial/liver-stator/blob/main/dependencies.R) which will try to
install multiple R packages from CRAN, BiocManager, or GitHub. Many of these will have to
be built from source which takes a long time. See [`dependencies.R`](https://github.com/whimsial/liver-stator/blob/main/dependencies.R) for details.

## Pipeline

Depends on:

- R packages: `data.table`, `Seurat`, `DropletUtils`, `HDF5Array`, `biomaRt` (for full list see [`dependencies.R`](https://github.com/whimsial/liver-stator/blob/main/dependencies.R)).

R packages should be installed user R library (see [`dependencies.R`](https://github.com/whimsial/liver-stator/blob/main/dependencies.R) for details).


- Python packages: `scanpy`, `scrublet` (for full list see [`doublet.py`](https://github.com/whimsial/liver-stator/blob/main/doublet.py) and [`requirements.txt`](https://github.com/whimsial/liver-stator/blob/main/requirements.txt)).

The Python code can also be run interactively using Jupyter notebook [`doublets.ipynb`](https://github.com/whimsial/liver-stator/blob/main/doublets.ipynb). It is useful to inspect distributions of doublets by eye to come up with a suitable threshold to be used in [`doublet.py`](https://github.com/whimsial/liver-stator/blob/main/doublet.py) (default is 0.15).

I found it useful to run Jupyter notebook inside Python virtual environment on the HPC and connect to it via SSH tunnel:

- start SSH tunnel to HPC mapping one of the open ports (I use port 9999, add `-R 9999:localhost:9999` to your ssh command)
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

There is an option here chosing between 2 packages `scrublet` and `scDblFinder`. The first one is older and implemented in Python, the second one is more modern, implemented in R and tends to discover more doublets. The `scDblFinder` is implemented directly in the pipeline.
For `scrublet` the key operation in this step is a call to `scrub_doublets` function from the `scrublet` package which is designed to detect doublets (also known as multiplets) in single-cell RNA sequencing (scRNA-seq) data. Doublets occur when two or more cells are inadvertently sequenced together as a single cell. Identifying and removing doublets is crucial for accurate downstream analysis because they can introduce significant noise and bias into the results.

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

Note, to update percentage of mitochondrial DNA removed you just need to update parameter `max.mt=10`, which specifes 10\% by default. This is callable from `qc.seurat` function.

### Step 5: process and visualize variable genes in a Seurat object

Perform a series of operations on a Seurat object to identify
and visualize highly variable genes. Normalize the data, identify variable
features, map gene IDs to Ensembl, and plot these genes on mean expression vs variance plot.

In this step we also add and label core genes from GATE analysis and check if they appear as highly variable genes in single cell data.

Finally, we write the counts matrix and a list of highly variable genes to files in a format expected by Stator.

### Step 6: run Stator on EDDIE (wildwest node2c15)

 - First you need to install Nextflow23 via anaconda. Follow `1-eddie-conda-setup.sh` via this [`Tutorial`](https://gist.github.com/laic/7b23e0fd21685f0527c91378fb45c395) but create environment nextflow23 instead

    `conda create --name nextflow23 bioconda::nextflow=23.04.4 conda-forge::singularity`

    `conda activate nextflow23`

    This needs to be setup only once and then can be used for all Stator runs. 

 - To run Stator you need mRNA expression counts and list of the highly variable genes. The path to the files containing these should be provided in [`stator.params.json`](https://github.com/whimsial/liver-stator/blob/main/stator.params.json), also please make sure to update number of cells and highly variable genes if changed. The rest of the parameters could be saved as default. 
 - Use [`stator.conda.sh`](https://github.com/whimsial/liver-stator/blob/main/stator.conda.sh) to run Stator.


## Known issues

While running this analysis I encountered several issues with the newest release of Seurat 5 R package. These are to do with the new layers introduced to the standard Seurat objects. Downgrading to version 4.4.0 together with seurat-object 4.1.4 solved these issues for me and thus I recommend to run this pipeline with these versions.

In time I will open an issue on Seurat's GitHub to ask for help with Seurat 5.

When running Stator sometimes it may fail on the last step. Exact reasons why are not known, but it could be due to memory overload on the wildwest node. Just try to resbmit the job adding `-resume` flad to the last line in `stator.conda.sh`. 

Also sometimes I get Nextflow error 

```Command exit status:                          140

Command output:
  Modules imported                                                                        Calculating linkage matrix...
  Linkage matrix calculated
  Calculating using 64 cores...

Command error:
  INFO:    Converting SIF file to temporary sandbox...```

Resuming speficying number of cores `--requestedCPU 31` seem to remedy it.
