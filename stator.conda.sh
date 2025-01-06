#!/bin/bash

module load anaconda
conda activate nextflow23

module add roslin/openjdk/22.0.0
module load singularity

export NXF_SINGULARITY_CACHEDIR=/exports/igmm/eddie/ponting-lab/sbraich2/.singularity
export SINGULARITY_CACHEDIR=/exports/igmm/eddie/ponting-lab/sbraich2/mycontainers
NXF_VER=23.04.4 nextflow run AJnsm/Stator -r main -profile eddie_singularity -params-file 1.params.json
