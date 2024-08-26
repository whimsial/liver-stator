#!/bin/bash

## Bash script to submit a Stator job from Wild west end nodes on Eddie assuming
## counts.csv and genes.csv are in /exports/eddie/scratch/aiakvlie/

## load NextFlow
module load igmm/bac/nextflow/24.04.2

## pull the latest version of the pipeline
nextflow pull AJnsm/Stator -r main

## run the pipeline
NXF_VER=24.04.2 nextflow run AJnsm/Stator -r main -profile eddie_singularity -params-file stator.params.json
