# ALL scanpy analysis
### Samples ALL1, ALL3, ALL8, ALL9, ALL10, ALL12, ALL10d15, ALL12d15
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
wd = '/research/groups/allseq/data/scRNAseq/results/juha_wrk/ALL_scanpy/new_2019_ver5/'
results_file = wd + 'ALL.h5ad'
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
ALL1 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL1/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL3 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL3/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL8 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL8/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL9 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL9/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL10 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL10/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL12 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL12/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL10d15 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL10d15/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL12d15 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL12d15/outs/filtered_feature_bc_matrix.h5', genome='hg19')

adata = ALL1.concatenate(ALL3, ALL8, ALL9, ALL10, ALL12, ALL10d15, ALL12d15, batch_categories=['ALL1', 'ALL3', 'ALL8', 'ALL9', 'ALL10', 'ALL12', 'ALL10d15', 'ALL12d15'], index_unique='_')

del ALL1, ALL3, ALL8, ALL9, ALL10, ALL12, ALL10d15, ALL12d15



### Remove useless variables
for i in adata.var.columns.tolist():
	del adata.var[i]



## Preprocessing

### Soft filtering for cells and genes
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=100) # n_obs × n_vars = 50747 × 14138



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



### Filter low quality cells
adata = adata[adata.obs['percent_mito'] < 0.10] # n_obs × n_vars = 44826 × 14138
adata = adata[adata.obs['n_counts'] < 50000] # n_obs × n_vars = 44752 × 14138
adata = adata[adata.obs['n_genes'] < 6000] # n_obs × n_vars = 44746 × 14138



### Filter genes again
sc.pp.filter_genes(adata, min_cells=200) # n_obs × n_vars = 44746 × 12643



### Log transform and store raw data for further use
adata.raw = sc.pp.log1p(adata, copy=True)



### Normalize (CPM) and get HVG
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5) # sum(adata.var['highly_variable']) >>> 1425
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
sc.pl.umap(adata, color=['IL7R', 'NKG7', 'GNLY', 'CD3D', 'CD8A'], save='_T_cell_markers.png')
sc.pl.umap(adata, color=['HBA1', 'HBA2', 'HBB', 'HBD'], save='_erythroid_precursor_cell_markers.png')
sc.pl.umap(adata, color=['LYZ', 'CD14', 'FCGR3A', 'MS4A7'], save='_monocyte_cell_markers.png')
sc.pl.umap(adata, color=['MS4A1', 'CD74', 'CD37', 'HLA-DQA1'], save='_mature_B_cell_markers.png')
sc.pl.umap(adata, color=['CD34', 'KIT', 'CD33', 'CD38'], save='_HSC_progenitor_cell_markers.png')
sc.pl.umap(adata, color=['STAG3', 'CYGB', 'SOCS2', 'DNTT', 'TP53INP1', 'FAM111B', 'TYMS', 'SPN', 'PTPRC'], save='_pro-B_cell_markers.png')
sc.pl.umap(adata, color=['IGLL1', 'VPREB1', 'PCDH9', 'SOX4', 'BCL7A', 'NSMCE1', 'ARPP21', 'STMN1', 'TCL1A'], save='_pre-B_cell_markers.png')
sc.pl.umap(adata, color=['MME', 'SOX4', 'IGLL1', 'GNG11', 'CTGF'], save='_ER_markers.png')
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
sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, save='.png')
sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, use_raw=True, swap_axes=True, save='_top5_markergenes_raw.png')
adata.write(results_file)



### Cell cycle scoring
with open('/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/regev_lab_cell_cycle_genes_S_phase.txt') as file:
	s_genes = file.read().splitlines()

with open('/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/regev_lab_cell_cycle_genes_G2M_phase.txt') as file:
	g2m_genes = file.read().splitlines()

s_genes = list(set(s_genes) & set(adata.raw.var_names))
g2m_genes = list(set(g2m_genes) & set(adata.raw.var_names))

adata.raw.obs = adata.obs
adata.raw.obs_names = adata.obs_names

sc.tl.score_genes_cell_cycle(adata.raw, s_genes=s_genes, g2m_genes=g2m_genes) # Use temp
sc.pl.umap(adata, color='phase', save='_cellcycle_phase.png')

adata.write(results_file)



## Calculate MAD (Median Absolute Deviation) for various features and plot
def calculate_MAD(adata, groupby, variable):
	df = pd.DataFrame()
	df[variable] = adata.obs[variable]
	df[groupby] = adata.obs[groupby]
	median = df.groupby(groupby)[variable].median()
	df['groupby_median'] = df[groupby].map(median.to_dict())
	df['groupby_absolute_deviation'] = abs(df[variable] - df['groupby_median'])
	mad = df.groupby(groupby)['groupby_absolute_deviation'].median()
	df['groupby_MAD'] = df[groupby].map(mad.to_dict())
	df['groupby_MAD_diff_x'] = df['groupby_absolute_deviation'] / df['groupby_MAD']
	df['groupby_MAD_diff'] = df['groupby_MAD'] - df['groupby_absolute_deviation']
	return df

n_genes_mad = calculate_MAD(adata=adata, groupby='louvain', variable='n_genes')
pct_mito_mad = calculate_MAD(adata=adata, groupby='louvain', variable='percent_mito')
n_counts_mad = calculate_MAD(adata=adata, groupby='louvain', variable='n_counts')

adata.obs['groupby_n_genes_MAD_diff_x'] = n_genes_mad['groupby_MAD_diff_x']
adata.obs['groupby_n_counts_MAD_diff_x'] = n_counts_mad['groupby_MAD_diff_x']
adata.obs['groupby_percent_mito_MAD_diff_x'] = pct_mito_mad['groupby_MAD_diff_x']

sc.pl.umap(adata, color='groupby_n_genes_MAD_diff_x', save='_MAD_n_genes_diff_x.png')
sc.pl.umap(adata, color='groupby_n_counts_MAD_diff_x', save='_MAD_n_counts_diff_x.png')
sc.pl.umap(adata, color='groupby_percent_mito_MAD_diff_x', save='_MAD_percent_mito_diff_x.png')

sc.pl.violin(adata, 'groupby_n_genes_MAD_diff_x', groupby='louvain', stripplot=False, save='_MAD_n_genes_diff_x.png')
sc.pl.violin(adata, 'groupby_n_counts_MAD_diff_x', groupby='louvain', stripplot=False, save='_MAD_n_counts_diff_x.png')
sc.pl.violin(adata, 'groupby_percent_mito_MAD_diff_x', groupby='louvain', stripplot=False, save='_MAD_percent_mito_diff_x.png')

adata.obs['groupby_n_genes_MAD_diff_x_gt_3'] = (adata.obs['groupby_n_genes_MAD_diff_x'] > 3).astype('category')
adata.obs['groupby_percent_mito_MAD_diff_x_gt_3'] = (adata.obs['groupby_percent_mito_MAD_diff_x'] > 3).astype('category')
adata.obs['groupby_n_counts_MAD_diff_x_gt_3'] = (adata.obs['groupby_n_counts_MAD_diff_x'] > 3).astype('category')
sc.pl.umap(adata, color=['groupby_n_genes_MAD_diff_x_gt_3', 'groupby_n_counts_MAD_diff_x_gt_3', 'groupby_percent_mito_MAD_diff_x_gt_3'], save='_MAD_diff_x_gt_3.png')

adata.obs['groupby_n_genes_MAD_diff_x_gt_5'] = (adata.obs['groupby_n_genes_MAD_diff_x'] > 5).astype('category')
adata.obs['groupby_percent_mito_MAD_diff_x_gt_5'] = (adata.obs['groupby_percent_mito_MAD_diff_x'] > 5).astype('category')
adata.obs['groupby_n_counts_MAD_diff_x_gt_5'] = (adata.obs['groupby_n_counts_MAD_diff_x'] > 5).astype('category')
sc.pl.umap(adata, color=['groupby_n_genes_MAD_diff_x_gt_5', 'groupby_n_counts_MAD_diff_x_gt_5', 'groupby_percent_mito_MAD_diff_x_gt_5'], save='_MAD_diff_x_gt_5.png')



### Filter distinct outliers further with mito% and n_counts MAD difference > 5
ids = adata.obs.index[~(adata.obs['groupby_n_counts_MAD_diff_x_gt_5'] | adata.obs['groupby_percent_mito_MAD_diff_x_gt_5'])]



### Set new paths
wd = '/research/groups/allseq/data/scRNAseq/results/juha_wrk/ALL_scanpy/new_2019_ver5/subset/'
results_file = wd + 'ALL_subset.h5ad'
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



## Read data (again)
ALL1 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL1/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL3 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL3/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL8 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL8/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL9 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL9/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL10 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL10/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL12 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL12/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL10d15 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL10d15/outs/filtered_feature_bc_matrix.h5', genome='hg19')
ALL12d15 = sc.read_10x_h5(filename='/research/groups/allseq/data/scRNAseq/results/juha_wrk/cellranger3/ALL12d15/outs/filtered_feature_bc_matrix.h5', genome='hg19')

adata2 = ALL1.concatenate(ALL3, ALL8, ALL9, ALL10, ALL12, ALL10d15, ALL12d15, batch_categories=['ALL1', 'ALL3', 'ALL8', 'ALL9', 'ALL10', 'ALL12', 'ALL10d15', 'ALL12d15'], index_unique='_')

del ALL1, ALL3, ALL8, ALL9, ALL10, ALL12, ALL10d15, ALL12d15



### Remove useless variables
for i in adata2.var.columns.tolist():
	del adata2.var[i]



### Keep cells passing MAD filter
adata2 = adata2[ids]



## Preprocessing

### Soft filtering for cells and genes
sc.pp.filter_cells(adata2, min_genes=1)

### Log transform and store raw data for further use
adata2.raw = sc.pp.log1p(adata2, copy=True)

d = pd.DataFrame(adata2.X.toarray())
d.to_csv(results_folder_raw + 'X_allgenes.csv.gz', sep=',', header=False, index=False, compression='gzip')

d = pd.DataFrame(adata2.var.index)
d.to_csv(results_folder_raw + 'var_allgenes.txt', sep=',', header=False, index=False)

del d

sc.pp.filter_genes(adata2, min_cells=100) # n_obs × n_vars = 41180 × 13805



### Count some metrics and plot
mito_genes = [name for name in adata2.var_names if name.startswith('MT-')]
adata2.obs['percent_mito'] = np.sum(adata2[:, mito_genes].X, axis=1).A1 / np.sum(adata2.X, axis=1).A1
ribo_genes = [name for name in adata2.var_names if name.startswith('RP')]
adata2.obs['percent_ribo'] = np.sum(adata2[:, ribo_genes].X, axis=1).A1 / np.sum(adata2.X, axis=1).A1
adata2.obs['n_counts'] = adata2.X.sum(axis=1).A1

sc.pl.violin(adata2, ['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], stripplot=False, multi_panel=True, save='_QC.pdf')
sc.pl.scatter(adata2, x='n_counts', y='percent_mito', save='_nCounts_pctMito.pdf')
sc.pl.scatter(adata2, x='n_counts', y='percent_ribo', save='_nCounts_pctRibo.pdf')
sc.pl.scatter(adata2, x='n_counts', y='n_genes', save='_nCounts_nGenes.pdf')
sc.pl.scatter(adata2, x='n_genes', y='percent_mito', save='_nGenes_pctMito')
sc.pl.scatter(adata2, x='n_genes', y='percent_ribo', save='_nGenes_pctRibo')



### Filter genes again
sc.pp.filter_genes(adata2, min_cells=200) # n_obs × n_vars = 41180 × 12483



### Normalize (CPM) and get HVG
sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5) # sum(adata2.var['highly_variable']) >>> 1938
sc.pl.highly_variable_genes(adata2, save='.pdf')



### Save
adata2.write(results_file)



### Filter by HVG
adata2 = adata2[:, adata2.var['highly_variable']]



### Regress out effect of n_counts and pct_mito
sc.pp.regress_out(adata2, ['n_counts', 'percent_mito'], n_jobs=n_jobs)
sc.pp.scale(adata2, max_value=10)



### Save
#adata2.write(results_file)



### Compute PCA
sc.tl.pca(adata2, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata2, log=True, save='.png')
#adata2.write(results_file)



### Compute neighborhood graph
sc.pp.neighbors(adata2)



### Compute UMAP
sc.tl.umap(adata2)



### Save
#adata2.write(results_file)



### Color markers and QC features (UMAP)
sc.pl.umap(adata2, color=['IL7R', 'NKG7', 'GNLY', 'CD3D', 'CD8A'], save='_T_cell_markers.png')
sc.pl.umap(adata2, color=['HBA1', 'HBA2', 'HBB', 'HBD'], save='_erythroid_precursor_cell_markers.png')
sc.pl.umap(adata2, color=['LYZ', 'CD14', 'FCGR3A', 'MS4A7'], save='_monocyte_cell_markers.png')
sc.pl.umap(adata2, color=['MS4A1', 'CD74', 'CD37', 'HLA-DQA1'], save='_mature_B_cell_markers.png')
sc.pl.umap(adata2, color=['CD34', 'KIT', 'CD33', 'CD38'], save='_HSC_progenitor_cell_markers.png')
sc.pl.umap(adata2, color=['STAG3', 'CYGB', 'SOCS2', 'DNTT', 'TP53INP1', 'FAM111B', 'TYMS', 'SPN', 'PTPRC'], save='_pro-B_cell_markers.png')
sc.pl.umap(adata2, color=['IGLL1', 'VPREB1', 'PCDH9', 'SOX4', 'BCL7A', 'NSMCE1', 'ARPP21', 'STMN1', 'TCL1A'], save='_pre-B_cell_markers.png')
sc.pl.umap(adata2, color=['MME', 'SOX4', 'IGLL1', 'GNG11', 'CTGF'], save='_ER_markers.png')
sc.pl.umap(adata2, color=['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], save='_QC.png')



### Batch plots
sc.pl.umap(adata2, color='batch', save='_batch.png')



### Clustering
sc.tl.louvain(adata2)
#adata2.write(results_file)



### Try multiple resolutions with clustering
for i in range(11, 21):
	r = int(i) / 10
	k = 'louvain_' + str(r)
	sc.tl.louvain(adata2, resolution=r, key_added=k)

louvains = [l for l in adata2.obs.columns if l.startswith('louvain')]
sc.pl.umap(adata2, color=louvains, legend_loc='on data', save='_clustering_resolutions.png')



### Choose res 1.4 as clustering (res 1.0 was mixing up come ALL10d15 cells with T cells)
adata2.obs['louvain_1.0'] = adata2.obs['louvain']
adata2.obs['louvain'] = adata2.obs['louvain_1.4']



### Cluster plots
sc.pl.umap(adata2, color='louvain', save='_clusters.png')
sc.pl.umap(adata2, color='louvain', legend_loc='on data', save='_clusters_labels.png')



## Finding marker genes
sc.tl.rank_genes_groups(adata2, 'louvain', method='wilcoxon')
sc.pl.rank_genes_groups(adata2, n_genes=20, save='.png')
sc.pl.rank_genes_groups_stacked_violin(adata2, n_genes=5, use_raw=True, swap_axes=True, save='_top5_markergenes_raw.png')
adata2.write(results_file)



### Cell cycle scoring
with open('/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/regev_lab_cell_cycle_genes_S_phase.txt') as file:
	s_genes = file.read().splitlines()

with open('/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/regev_lab_cell_cycle_genes_G2M_phase.txt') as file:
	g2m_genes = file.read().splitlines()

s_genes = list(set(s_genes) & set(adata2.raw.var_names))
g2m_genes = list(set(g2m_genes) & set(adata2.raw.var_names))

adata2.raw.obs = adata2.obs
adata2.raw.obs_names = adata2.obs_names

sc.tl.score_genes_cell_cycle(adata2.raw, s_genes=s_genes, g2m_genes=g2m_genes) # Use temp
sc.pl.umap(adata2, color='phase', save='_cellcycle_phase.png')

adata2.write(results_file)



### Plot QC violin plots per cluster
sc.pl.violin(adata2, 'n_genes', groupby='louvain', stripplot=False, save='_n_genes_bygroup.png')
sc.pl.violin(adata2, 'n_counts', groupby='louvain', stripplot=False, save='_n_counts_bygroup.png')
sc.pl.violin(adata2, 'percent_mito', groupby='louvain', stripplot=False, save='_percent_mito_bygroup.png')
sc.pl.violin(adata2, 'percent_ribo', groupby='louvain', stripplot=False, save='_percent_ribo_bygroup.png')
sc.pl.scatter(adata2, x='n_counts', y='percent_mito', save='_nCounts_pctMito.pdf')



### Cell types
mapping = {
0:'leukemic',
1:'leukemic',
2:'leukemic',
3:'CD4+ T', # IL7R+, CD3D+, CD8A-
4:'immature B', # MS4A1+, HLA-DQA1+, CD37+
5:'leukemic',
6:'leukemic',
7:'leukemic',
8:'leukemic',
9:'leukemic', #erythoid precursor 1
10:'CD8+ T',
11:'CD8+ T', # CD8A+, CD3D+
12:'erythoid precursor', # NKG7+, CD3D+
13:'erythoid precursor',
14:'CD4+ T', # IL7R+, CD3D+, CD8A-
15:'leukemic',
16:'leukemic',
17:'erythoid precursor',
18:'leukemic', # LYZ+, 
19:'NK', # NKG7+, GNLY+, CD3D-
20:'Monocyte', # IL7R+, CD3D+, CD8A-
21:'CD4+ T',
22:'leukemic',
23:'erythoid precursor',
24:'leukemic',
25:'CD34+ precursor',
26:'Plasma B cell'
}
adata2.obs['celltype'] = adata2.obs['louvain'].astype('int32').map(mapping).astype('category')

sc.pl.umap(adata2, color='celltype', legend_loc='on data', legend_fontsize=8, save='_celltypes.png')

adata2.write(results_file)



### Extract leukemic cells (G1) and plot n_genes vs mito%
#adata2_subset = adata2[np.array(adata2.obs['louvain'].isin(['0','1','2','4','6','7','10','18'])) & np.array(adata.obs['phase'] == 'G1')]

#sc.pl.scatter(adata2_subset, x='n_genes', y='percent_mito', save='_leukemic_G1_nGenes_pctMito')



### barplot of cc prop per cluster
df = adata2.obs.groupby('louvain')['phase'].value_counts(normalize=True).rename(columns={0:'score'}).reset_index().rename(columns={0:'score'})
g = sns.catplot(x='louvain', y='score', hue='phase', data=df, kind='bar')
g.despine(left=True)
g.set_ylabels("%")
g.set_xticklabels(rotation=90)
g.savefig('figures/barplot_louvain_cc_prop.png')



### Plot comparison of splitted cluster (res 1.0 vs. 1.4 clusters, 20 vs 21 & 24)
markergenes = ['IL7R', 'MME', 'GNG11', 'SOX4', 'CTGF']
sc.pl.matrixplot(adata2, var_names=markergenes, groupby='louvain_1.0', save='_res1.0_markergenes.png')
sc.pl.matrixplot(adata2, var_names=markergenes, groupby='louvain_1.4', save='_res1.4_markergenes.png')



### Score genesets
resdir = '/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/'
files = os.listdir(resdir)
files = [i for i in files if 'txt' in i if 'chen' in i]
for i in files:
	file = open(resdir + i, 'r')
	gs = file.read().split('\n')
	file.close()
	n = i[:-4]
	gs = adata2.raw.var_names[adata2.raw.var_names.isin(gs)].tolist()
	if len(gs) == 0:
		continue
	sc.tl.score_genes(adata2, gene_list=gs, ctrl_size=len(gs), score_name=n, use_raw=True)
	savename = '_' + n + '.png'
	sc.pl.umap(adata2, color=n, save=savename)

adata2.write(results_file)



### Write data out
d = pd.DataFrame(adata2.raw.X.toarray())
d.to_csv(results_folder_raw + 'X.csv.gz', sep=',', header=False, index=False, compression='gzip')
del d

d = pd.DataFrame(adata2.raw.obs.index)
d.to_csv(results_folder_raw + 'obs.txt', sep=',', header=False, index=False)
del d

d = pd.DataFrame(adata2.raw.var.index)
d.to_csv(results_folder_raw + 'var.txt', sep=',', header=False, index=False)
del d

adata2.write_csvs(dirname=results_folder, skip_data=True)



### Plot genes indicated in Muschen article
sc.pl.umap(adata2, color=['CD72', 'BLNK', 'SPN', 'PTPN6', 'IL7R', 'KIT', 'CRLF2', 'IGLC2'], save='_muschen_story_genes.png')



### Genes related to ALL relapse
# with open('/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/senft_jeremias_ALL_relapse_enriched_genes.txt') as file:
# 	relapse = file.read().splitlines()

# relapse = adata2.raw.var_names[adata2.raw.var_names.isin(relapse)].tolist()
# sc.tl.score_genes(adata2, gene_list=relapse, ctrl_size=len(relapse), score_name='senft_jeremias_relapse_geneset', use_raw=True)

# sc.pl.umap(adata2, color=relapse, save='_senft_jeremias_relapse_genes.png')
# sc.pl.umap(adata2, color='senft_jeremias_relapse_geneset', save='_senft_jeremias_relapse_geneset_score.png')



### Subset leukemic cells and plot QC 
leukemic = adata2[adata2.obs['celltype'] == 'leukemic'].copy()
sc.pl.violin(leukemic, 'n_genes', stripplot=False, groupby='batch', save='_leukemic_n_genes.pdf')
leukemic_G1 = leukemic[leukemic.obs['phase'] == 'G1']
sc.pl.violin(leukemic_G1, 'n_genes', stripplot=False, groupby='batch', save='_leukemic_G1_n_genes.pdf')
leukemic_cc = leukemic[leukemic.obs['phase'] != 'G1']
sc.pl.violin(leukemic_cc, 'n_genes', stripplot=False, groupby='batch', save='_leukemic_cycling_n_genes.pdf')
leukemic_G2M = leukemic[leukemic.obs['phase'] == 'G2M']
sc.pl.violin(leukemic_G2M, 'n_genes', stripplot=False, groupby='batch', save='_leukemic_G2M_n_genes.pdf')
leukemic_S = leukemic[leukemic.obs['phase'] == 'S']
sc.pl.violin(leukemic_S, 'n_genes', stripplot=False, groupby='batch', save='_leukemic_S_n_genes.pdf')










