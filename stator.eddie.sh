#!/bin/bash

## Bash script to run Stator using nextflow

## load modules
module add roslin/openjdk/13.0.1
# module add igmm/bac/nextflow/24.04.2
module load roslin/nextflow/22.10.7
module load singularity

## pull development branch of Stator
nextflow pull AJnsm/Stator -r develop

## create environment variables
export NXF_SINGULARITY_CACHEDIR=/exports/igmm/eddie/khamseh-lab/aiakovliev/.singularity
export SINGULARITY_CACHEDIR=/exports/igmm/eddie/khamseh-lab/aiakovliev/mycontainers

## run the pipeline
NXF_VER=22.10.7 nextflow run whimsial/Stator -r develop -profile eddie_singularity -params-file stator.params.json -resume
