#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
import getopt


# In[2]:


#  plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42


# In[3]:


input_dir = '/Volumes/khamseh-lab/ava/ME_scRNA/raw_data/'
sample_names = pd.read_csv('/Volumes/khamseh-lab/ava/ME_scRNA/raw_data/sample_names.csv',header = None)[0].tolist()


# In[4]:


#counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
#genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1))
#print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
#print('Number of genes in gene list: {}'.format(len(genes)))


# In[6]:


sample_names


# In[5]:


for s in sample_names:
    print(s[0:22])


# In[7]:


get_ipython().run_line_magic('matplotlib', 'inline')

for s in sample_names:
    
    # load count matrix
    counts_matrix=sc.read_10x_h5(input_dir+s)
    
    # doublet rates in 1000 cells normally ~0.8% per 1k cells
    edr = 0.008
    expected_doublet_rate = float(edr)*counts_matrix.shape[0]/1000
    print('EDR: {}'.format(expected_doublet_rate))
    
    scrub = scr.Scrublet(counts_matrix.X, expected_doublet_rate = float(expected_doublet_rate))
    
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                      min_cells=3, 
                                                      min_gene_variability_pctl=60, 
                                                      n_prin_comps=30,
                                                      log_transform=True,
                                                      mean_center=True,
                                                      normalize_variance=True,
                                                      synthetic_doublet_umi_subsampling = 1)
    
    # after examining the bimodal distribution by eye, insert the threshold here
    scrub.call_doublets(threshold=0.18)
    scrub.plot_histogram()

    # s[0:22] only takes the first 22 characters, which contains the unique sample name
    outf1 = "/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/" + s[0:22] + "_scrublet_EDR" + str(edr) + "_DoubletScores.csv"
    outf2 = "/Volumes/khamseh-lab/ava/ME_scRNA/qc_analysis/" + s[0:22] + "_scrublet_EDR" + str(edr) + "_PredictedDoublets.csv"
    np.savetxt(outf1, doublet_scores, delimiter=',')
    np.savetxt(outf2, predicted_doublets, delimiter=',')
    print("sample",s,"done")


# In[ ]:





# In[ ]:





# ## Check 2 of the samples that look different

# In[10]:


# load count matrix
counts_matrix=sc.read_10x_h5("/Volumes/khamseh-lab/ava/ME_scRNA/raw_data/GSM6603131_COR-5218-D1_filtered_feature_bc_matrix.h5")
    
# doublet rates in 1000 cells normally ~0.8% per 1k cells
edr = 0.008
expected_doublet_rate = float(edr)*counts_matrix.shape[0]/1000
print('EDR: {}'.format(expected_doublet_rate))
    
scrub = scr.Scrublet(counts_matrix.X, expected_doublet_rate = float(expected_doublet_rate))
    
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=60, 
                                                          n_prin_comps=30,
                                                          log_transform=True,
                                                          mean_center=True,
                                                          normalize_variance=True,
                                                          synthetic_doublet_umi_subsampling = 1)
    
# after examining the bimodal distribution by eye, insert the threshold here
scrub.call_doublets(threshold=0.18)
scrub.plot_histogram()


# In[ ]:




