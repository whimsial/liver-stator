#!/bin/bash

## Bash script to submit a Stator job from Wild west end nodes on Eddie assuming
## counts.csv and genes.csv are in /exports/eddie/scratch/aiakvlie/

## pull the latest version of the pipeline
nextflow pull AJnsm/Stator -r main

## run the pipeline
NXF_VER=23.04.4 nextflow run AJnsm/Stator -r main -profile eddie_singularity -params-file params.json
