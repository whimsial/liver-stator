#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

## set graphics parameters
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size = 14)
plt.rcParams['pdf.fonttype'] = 42

## get input parameters
parser = argparse.ArgumentParser(description="Doublet detection.")

parser.add_argument(
    '--sample_dir', type=str,
    help='A directory path for the sample'
)

parser.add_argument(
    '--data_type', type=str, required=True,
    help='Specify type of sample files (mtx or h5)'
)

parser.add_argument(
    '--doublet_threshold', type=float, required=True,
    help='Specify doublet threshold after examining distributions'
)

args = parser.parse_args()

print("You have passed the following directory:")
s = args.sample_dir
print(s)
if os.path.isdir(s):
    print(f"Directory: {s}")
else:
    print(f"Invalid directory: {s}")

## load count matrix
if args.data_type=="h5":
    counts_matrix = sc.read_10x_h5(s)
else:
    counts_matrix = sc.read_10x_mtx(s, var_names='gene_symbols', cache=True)

## doublet rates in 1000 cells normally ~0.8% per 1k cells
edr = 0.008
expected_doublet_rate = float(edr) * counts_matrix.shape[0] / 1000
print('EDR: {}'.format(expected_doublet_rate))

scrub = scr.Scrublet(counts_matrix.X,
                     expected_doublet_rate = float(expected_doublet_rate))

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2,
    min_cells=3,
    min_gene_variability_pctl=60,
    n_prin_comps=30,
    log_transform=True,
    mean_center=True,
    normalize_variance=True,
    synthetic_doublet_umi_subsampling=1)

scrub.call_doublets(threshold=args.doublet_threshold)

outf1 = s + "/scrublet_EDR" + str(edr) + "_DoubletScores.csv"
outf2 = s + "/scrublet_EDR" + str(edr) + "_PredictedDoublets.csv"
np.savetxt(outf1, doublet_scores, delimiter=',')
np.savetxt(outf2, predicted_doublets, delimiter=',')
print("Saved " + outf1)
print("Saved " + outf2)
