#!/bin/bash

## To reattach your session if you get bumped of you'll need to get to the same login node that you started from
## to get to login1 you neeed to do: ssh login01-ext.ecdf.ed.ac.uk

## load anaconda
module load anaconda

## Get setup with anaconda: 
## SLP students should have access to this group space:
## /exports/chss/eddie/ppls/groups/lel_hcrc_cstr_students
## So, you can make a directory there called UUN_Firstname_Lastname and work in that
## You can probably just put things in your home directory or in the scratch space (/exports/eddie/scratch/<you-uun>)
## PLEASE NOTE THAT FILES IN THE SCRATCH SPACE GET DELETED AFTER 1 MONTH!

## See this page on setting anaconda directories
## https://www.wiki.ed.ac.uk/display/ResearchServices/Anaconda

CONDADIR=/exports/igmm/eddie/ponting-lab/sbraich2/anaconda
mkdir -p $CONDADIR
mkdir -p $CONDADIR/envs
mkdir -p $CONDADIR/pkgs

conda create --name nextflow23 bioconda::nextflow=23.04.4 conda-forge::singularity
