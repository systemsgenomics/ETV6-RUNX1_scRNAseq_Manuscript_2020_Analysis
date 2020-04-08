# Nalm6 scanpy analysis
### Samples Nalm6 DMSO +/- RUNX1 overexpression (24h)
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
sc.settings.set_figure_params(dpi=160, dpi_save=160, vector_friendly=True)  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()
sc.settings.autoshow = False
wd = '/research/groups/allseq/data/scRNAseq/results/juha_wrk/Nalm6_RUNX1_overexp_scanpy/'
results_file = wd + 'Nalm6_RUNX1.h5ad'
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
DMSO = sc.read_10x_h5(filename='/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S2/outs/filtered_feature_bc_matrix.h5')
DMSO_gem = pd.read_csv('/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S2/outs/analysis/gem_classification.csv')
DMSO_gem.index = DMSO_gem.barcode
DMSO.obs['hg19_counts'] = DMSO_gem['hg19']
DMSO.obs['mm10_counts'] = DMSO_gem['mm10']
DMSO.obs['call'] = DMSO_gem['call']
DMSO.var_names_make_unique()

lvRUNX1 = sc.read_10x_h5(filename='/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S7/outs/filtered_feature_bc_matrix.h5')
lvRUNX1_gem = pd.read_csv('/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S7/outs/analysis/gem_classification.csv')
lvRUNX1_gem.index = lvRUNX1_gem.barcode
lvRUNX1.obs['hg19_counts'] = lvRUNX1_gem['hg19']
lvRUNX1.obs['mm10_counts'] = lvRUNX1_gem['mm10']
lvRUNX1.obs['call'] = lvRUNX1_gem['call']
lvRUNX1.var_names_make_unique()

adata = DMSO.concatenate(lvRUNX1, batch_categories=['DMSO', 'lvRUNX1'], index_unique='_')



### Remove mouse cells and genes from both data separately and save results (for GEO submission)
DMSO = DMSO[(DMSO.obs['call'] == 'hg19') & (DMSO.obs['mm10_counts'] < 500),:]
DMSO = DMSO[:, DMSO.var_names.str.startswith('hg19')]
DMSO.var_names = [re.sub('hg19_', '', var) for var in DMSO.var_names]

lvRUNX1 = lvRUNX1[(lvRUNX1.obs['call'] == 'hg19') & (lvRUNX1.obs['mm10_counts'] < 500),:]
lvRUNX1 = lvRUNX1[:, lvRUNX1.var_names.str.startswith('hg19')]
lvRUNX1.var_names = [re.sub('hg19_', '', var) for var in lvRUNX1.var_names]

from scipy.io import mmwrite

mmwrite(target='/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S2/outs/Nalm6_DMSO_matrix.mtx', a=DMSO.X)
DMSO.obs.to_csv('/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S2/outs/Nalm6_DMSO_barcode_metadata.tsv.gz', sep='\t', compression='gzip')
DMSO.var.to_csv('/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S2/outs/Nalm6_DMSO_features.tsv.gz', sep='\t', compression='gzip')

mmwrite(target='/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S7/outs/Nalm6_lvRUNX1_matrix.mtx', a=lvRUNX1.X)
lvRUNX1.obs.to_csv('/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S7/outs/Nalm6_lvRUNX1_barcode_metadata.tsv.gz', sep='\t', compression='gzip')
lvRUNX1.var.to_csv('/research/groups/sysgen/PROJECTS/students/scATAC-seq/cellranger/S7/outs/Nalm6_lvRUNX1_features.tsv.gz', sep='\t', compression='gzip')


del DMSO, lvRUNX1



### Remove useless variables
for i in adata.var.columns.tolist():
	del adata.var[i]



### Remove mouse cells and genes
adata = adata[(adata.obs['call'] == 'hg19') & (adata.obs['mm10_counts'] < 500),:]
adata = adata[:, adata.var_names.str.startswith('hg19')]
adata.var_names = [re.sub('hg19_', '', var) for var in adata.var_names]


### Log transform and store raw data for further use
adata.raw = sc.pp.log1p(adata, copy=True)



## Preprocessing

### Soft filtering for cells and genes
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=100) # n_obs × n_vars = 11564 × 13161 



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
adata = adata[adata.obs['percent_mito'] < 0.10] # n_obs × n_vars = 4606 × 13161
adata = adata[adata.obs['n_counts'] < 40000] # n_obs × n_vars = 4596 × 13161
adata = adata[adata.obs['n_genes'] < 5000] # n_obs × n_vars = 4431 × 13161

##NR.filt <- subset(NR, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 20 & mm10counts < 500)


### Filter genes again
sc.pp.filter_genes(adata, min_cells=500) # n_obs × n_vars = 4431 × 7368



### Normalize (CPM) and get HVG
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.02, max_mean=5, min_disp=1) # sum(adata.var['highly_variable']) >>> 1013
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
sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, use_raw=True, swap_axes=True, save='_top20_markergenes_raw.png')
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


markers = ['LAT2', 'CYFIP2', 'IGLL1', 'VPREB1']
sc.pl.violin(adata, keys=markers, groupby='batch', stripplot=False, use_raw=True, save='_markers.pdf')


sc.pl.violin(adata[adata.obs['phase'] == 'G1'], keys=markers, groupby='batch', stripplot=False, use_raw=True, bw=0.5, save='_markers_G1.pdf')

### Write data out
# d = pd.DataFrame(adata.raw.X.toarray())
# d.to_csv(results_folder_raw + 'X.csv.gz', sep=',', header=False, index=False, compression='gzip')

# d = pd.DataFrame(adata.raw.obs.index)
# d.to_csv(results_folder_raw + 'obs.txt', sep=',', header=False, index=False)

# d = pd.DataFrame(adata.raw.var.index)
# d.to_csv(results_folder_raw + 'var.txt', sep=',', header=False, index=False)

# del d

# adata.write_csvs(dirname=results_folder, skip_data=True)



