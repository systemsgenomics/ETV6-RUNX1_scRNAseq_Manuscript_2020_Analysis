# Nalm6 scanpy analysis
### Samples Nalm6 LUC & ER
### No batch correction



## Imports
import numpy as np
import pandas as pd
import scanpy as sc
import os.path
import re
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
plt.switch_backend('agg')



## Settings
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=160, dpi_save=150, vector_friendly=True)  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()
sc.settings.autoshow = False
wd = '/research/groups/allseq/data/scRNAseq/results/juha_wrk/Nalm6_scanpy/'
results_file = wd + 'Nalm6.h5ad'
results_folder = wd + 'data/'
results_folder_raw = results_folder + 'raw/'
n_jobs = 16

if not os.path.exists(wd):
	os.makedirs(wd)

if not os.path.exists(results_folder):
	os.makedirs(results_folder)

if not os.path.exists(results_folder_raw):
	os.makedirs(results_folder_raw)

os.chdir(wd)



## Read data
LUC = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/Nalm6_LUC/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ER = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/Nalm6_ETV6RUNX1/outs/filtered_feature_bc_matrix.h5', genome='hg19')

adata = LUC.concatenate(ER, batch_categories=['LUC', 'ER'], index_unique='_')

del LUC, ER



### Remove useless variables
for i in adata.var.columns.tolist():
	del adata.var[i]



### Log transform and store raw data for further use
adata.raw = sc.pp.log1p(adata, copy=True)



## Preprocessing

### Soft filtering for cells and genes
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=100) # n_obs × n_vars = 4706 × 10370 



### Count some metrics and plot
mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
ribo_genes = [name for name in adata.var_names if name.startswith('RP')]
adata.obs['percent_ribo'] = np.sum(adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], stripplot=False, multi_panel=True, save='_QC.pdf')
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_nCounts_pctMito.pdf')
sc.pl.scatter(adata, x='n_counts', y='percent_ribo', save='_nCounts_pctRibo.pdf')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_nCounts_nGenes.pdf')
sc.pl.scatter(adata, x='n_genes', y='percent_mito', save='_nGenes_pctMito')
sc.pl.scatter(adata, x='n_genes', y='percent_ribo', save='_nGenes_pctRibo')
sc.pl.scatter(adata, x='percent_mito', y='percent_ribo', save='_pctMito_pctRibo')



### Filter low quality cells
adata = adata[adata.obs['percent_mito'] < 0.10] # n_obs × n_vars = 4194 × 10370
adata = adata[adata.obs['n_counts'] < 40000] # n_obs × n_vars = 4086 × 10370
adata = adata[adata.obs['n_genes'] < 5000] # n_obs × n_vars = 4079 × 10370



### Filter genes again
sc.pp.filter_genes(adata, min_cells=500) # n_obs × n_vars = 4079 × 6953



### Normalize (CPM) and get HVG
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.02, max_mean=5, min_disp=1) # sum(adata.var['highly_variable']) >>> 1043
sc.pl.highly_variable_genes(adata, save='.pdf')



### Save
adata.write(results_file)



### Filter by HVG
adata = adata[:, adata.var['highly_variable']]



### Regress out effect of n_counts and pct_mito
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'], n_jobs=n_jobs)
sc.pp.scale(adata, max_value=10)



### Save
#adata.write(results_file)



### Compute PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, save='.png')
#adata.write(results_file)



### Compute neighborhood graph
sc.pp.neighbors(adata)



### Compute UMAP
sc.tl.umap(adata)



### Save
#adata.write(results_file)



### Color markers and QC features (UMAP)
sc.pl.umap(adata, color=['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], save='_QC.png')



### Batch plots
sc.pl.umap(adata, color='batch', save='_batch.png')



### Clustering
sc.tl.louvain(adata)
#adata.write(results_file)



### Cluster plots
sc.pl.umap(adata, color='louvain', save='_clusters.png')
sc.pl.umap(adata, color='louvain', legend_loc='on data', save='_clusters_labels.png')



## Finding marker genes
sc.tl.rank_genes_groups(adata, 'batch', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, save='.png')
sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, use_raw=True, swap_axes=True, save='_top5_markergenes_raw.png')
adata.write(results_file)



### Cell cycle scoring
with open('/research/groups/sysgen/PROJECTS/LEUKEMIA/juha_wrk/regev_lab_cell_cycle_genes_S_phase.txt') as file:
	s_genes = file.read().splitlines()

with open('/research/groups/sysgen/PROJECTS/LEUKEMIA/juha_wrk/regev_lab_cell_cycle_genes_G2M_phase.txt') as file:
	g2m_genes = file.read().splitlines()

s_genes = list(set(s_genes) & set(adata.raw.var_names))
g2m_genes = list(set(g2m_genes) & set(adata.raw.var_names))

adata.raw.obs = adata.obs
adata.raw.obs_names = adata.obs_names

sc.tl.score_genes_cell_cycle(adata.raw, s_genes=s_genes, g2m_genes=g2m_genes) # Use temp
sc.pl.umap(adata, color='phase', save='_cellcycle_phase.png')

adata.write(results_file)



### Write data out
d = pd.DataFrame(adata.raw.X.toarray())
d.to_csv(results_folder_raw + 'X.csv.gz', sep=',', header=False, index=False, compression='gzip')

d = pd.DataFrame(adata.raw.obs.index)
d.to_csv(results_folder_raw + 'obs.txt', sep=',', header=False, index=False)

d = pd.DataFrame(adata.raw.var.index)
d.to_csv(results_folder_raw + 'var.txt', sep=',', header=False, index=False)

del d

adata.write_csvs(dirname=results_folder, skip_data=True)



