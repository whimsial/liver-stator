#!/bin/bash

## Bash script to run Stator using nextflow

## load modules
module add roslin/openjdk/13.0.1
# module add igmm/bac/nextflow/24.04.2
module load roslin/nextflow/22.10.7
module load singularity
# module load igmm/apps/singularity/3.9.9
# module load igmm/apps/apptainer/1.3.4

## pull development branch of Stator
# nextflow pull AJnsm/Stator -r develop
# nextflow pull avakhamseh/Stator -r main

## create environment variables
export NXF_SINGULARITY_CACHEDIR=/exports/igmm/eddie/ponting-lab/sbraich2/.singularity
export SINGULARITY_CACHEDIR=/exports/igmm/eddie/ponting-lab/sbraich2/mycontainers
export NXF_HOME=/exports/igmm/eddie/ponting-lab/sbraich2/.nextflow
export NXF_TEMP=/exports/igmm/eddie/ponting-lab/sbraich2/.tmp

## run the pipeline
NXF_VER=22.10.7 nextflow run avakhamseh/Stator -r main -profile eddie_singularity -params-file 1.params.json -resume -work-dir /exports/igmm/eddie/ponting-lab/sbraich2/
