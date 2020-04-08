# HCA BM scanpy analysis

## Cellranger v3
## cellranger filtered genes



## Imports
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import os.path
import re
import matplotlib.pyplot as plt
import matplotlib
plt.switch_backend('agg')



## Settings
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=160, dpi_save=150, vector_friendly=True)  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()
sc.settings.autoshow = False
wd = '/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/'
results_file = wd + 'HCA.h5ad'
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
fnames = ['/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/immune_census/cellranger.v3/MantonBM' + str(i) + '_HiSeq_' + str(j) + '/outs/filtered_feature_bc_matrix.h5' for i in range(1,9) for j in range(1,9)]
adata = sc.read_10x_h5(filename=fnames[0], genome='hg19')
name = re.search('MantonBM._HiSeq_.', fnames[0]).group(0)
adata.obs_names = [name + '-' + obs_name for obs_name in adata.obs_names]
for i in range(1,len(fnames)):
	if os.path.isfile(fnames[i]):
		tmp = sc.read_10x_h5(filename=fnames[i], genome='hg19')
		name = re.search('MantonBM._HiSeq_.', fnames[i]).group(0)
		tmp.obs_names = [name + '-' + obs_name for obs_name in tmp.obs_names]
		adata = adata.concatenate(tmp, batch_key=name)



### Fix cell names
adata.obs_names = [re.search(".*[TCGA]-1", i).group(0) for i in adata.obs_names.tolist()]



### Remove useless metadata
for i in adata.var.columns.tolist():
	del adata.var[i]

for i in adata.obs.columns.tolist():
	del adata.obs[i]



### Set batch info
adata.obs['batch'] = [re.sub("_.*", "", x) for x in adata.obs_names.tolist()]



### Save raw counts to .raw slot
adata.raw = sc.pp.log1p(adata, copy=True)



### Save
adata.write(results_file)



## Preprocessing

### Soft filtering for cells and genes
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=100) # n_obs × n_vars = 183632 × 16439



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
adata = adata[adata.obs['percent_mito'] < 0.10] # n_obs × n_vars = 170393 × 16439
adata = adata[adata.obs['n_counts'] < 50000] # n_obs × n_vars = 170106 × 16439
adata = adata[adata.obs['n_genes'] < 6000] # n_obs × n_vars = 170097 × 16439



### Filter genes again
sc.pp.filter_genes(adata, min_cells=400) # n_obs × n_vars = 170097 × 13528



### Normalize (CPM) and get HVG
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5) # sum(adata.var['highly_variable']) >>> 2046
sc.pl.highly_variable_genes(adata, save='.pdf')



### Save
adata.write(results_file)



### Filter by HVG
adata = adata[:, adata.var['highly_variable']]



### Regress out effect of n_counts and pct_mito
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'], n_jobs=n_jobs)
sc.pp.scale(adata, max_value=10)



### Save
adata.write(results_file)



### MNN batch correction
sample1 = adata[adata.obs['batch'] == 'MantonBM1']
sample2 = adata[adata.obs['batch'] == 'MantonBM2']
sample3 = adata[adata.obs['batch'] == 'MantonBM3']
sample4 = adata[adata.obs['batch'] == 'MantonBM4']
sample5 = adata[adata.obs['batch'] == 'MantonBM5']
sample6 = adata[adata.obs['batch'] == 'MantonBM6']
sample7 = adata[adata.obs['batch'] == 'MantonBM7']
sample8 = adata[adata.obs['batch'] == 'MantonBM8']

adata_raw = adata.raw

correction = sce.pp.mnn_correct(sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8, n_jobs=n_jobs, batch_categories=['MantonBM1', 'MantonBM2', 'MantonBM3', 'MantonBM4', 'MantonBM5', 'MantonBM6', 'MantonBM7', 'MantonBM8'])
adata = correction[0]

n_cells = adata.var['n_cells-MantonBM1']
means = adata.var['means-MantonBM1']
dispersions = adata.var['dispersions-MantonBM1']
dispersions_norm = adata.var['dispersions_norm-MantonBM1']

for i in adata.var.columns.tolist():
	del adata.var[i]

adata.var['n_cells'] = n_cells
adata.var['means'] = means
adata.var['dispersions'] = dispersions
adata.var['dispersions_norm'] = dispersions_norm

temp = ad.AnnData(X=adata_raw.X, var=adata_raw.var)
adata.raw = temp

sc.pp.scale(adata, max_value=10)

del correction, adata_raw, temp



### Save
adata.write(results_file)



### Compute PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, save='.png')
adata.write(results_file)



### Compute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=30)



### Compute UMAP
sc.tl.umap(adata)



### Save
adata.write(results_file)



### Color markers and QC features (UMAP)
sc.pl.umap(adata, color=['IL7R', 'NKG7', 'GNLY', 'CD3D', 'CD8A'], save='_T_cell_markers.png')
sc.pl.umap(adata, color=['HBA1', 'HBA2', 'HBB', 'HBD', 'GYPA'], save='_erythroid_precursor_cell_markers.png')
sc.pl.umap(adata, color=['LYZ', 'CD14', 'FCGR3A', 'MS4A7'], save='_monocyte_cell_markers.png')
sc.pl.umap(adata, color=['MS4A1', 'CD74', 'CD37', 'HLA-DQA1'], save='_mature_B_cell_markers.png')
sc.pl.umap(adata, color=['CD34', 'KIT', 'CD33', 'CD38'], save='_HSC_progenitor_cell_markers.png')
sc.pl.umap(adata, color=['STAG3', 'CYGB', 'SOCS2', 'DNTT', 'TP53INP1', 'FAM111B', 'TYMS', 'SPN', 'PTPRC'], save='_pro-B_cell_markers.png')
sc.pl.umap(adata, color=['IGLL1', 'VPREB1', 'PCDH9', 'SOX4', 'BCL7A', 'NSMCE1', 'ARPP21', 'STMN1', 'TCL1A'], save='_pre-B_cell_markers.png')
sc.pl.umap(adata, color=['MME', 'SOX4', 'IGLL1', 'GNG11', 'CTGF'], save='_ER_markers.png')
sc.pl.umap(adata, color=['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], save='_QC.png')



### Markers based on van Galen et al. 2019 Cell
sc.pl.umap(adata, color=['GZMB', 'IGJ', 'MZB1', 'TCF4', 'IRF8', 'PTPRS'], save='_galen_pDC_markers.png')
sc.pl.umap(adata, color=['FCER1A', 'CLEC10A'], save='_galen_cDC_markers.png')
sc.pl.umap(adata, color=['MPO', 'ELANE', 'CTSG', 'AZU1'], save='_galen_GMP_markers.png')
sc.pl.umap(adata, color=['FCER1G', 'FCN1', 'CD14', 'C5AR1', 'MNDA', 'LYZ', 'CEBPD', 'LYST'], save='_galen_promonocyte_markers.png')



### Batch plots
sc.pl.umap(adata, color='batch', save='_batch.png')



### Clustering
sc.tl.louvain(adata)
adata.write(results_file)



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



### Cell type assignment
adata.obs['celltype'] = adata.obs['louvain']
celltypes = ['CD14+ monocyte 2', # LYZ+, CD14+
'CD4+ T 1', # IL7R+, CD3D+, CD8A-
'CD14+ monocyte 1', # LYZ+, CD14+ 
'NK', # NKG7+, GNLY+, CD3D-
'NK T 1', # NKG7+, GNLY+, CD3D+
'NK T 2', # NKG7+, GNLY-, CD3D+
'Late erythoid precursor 1', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
'Immature B', # MS4A1+, HLA-DQA1+
'CD8+ T', # IL7R+, CD3D+, CD8A+
'Promonocyte', # LYZ+, CD14-, FCGR3A-
'Pre B', # IGLL1+, VPREB1+, PCDH9+, SOX4+, TCL1A+
'HSC', # CD34+, SPINK2+
'Conventional dendritic cell', # LYZ+, CD14-, FCGR3A-
'CD4+ T 2', # IL7R+, CD3D+, CD8A-
'CD16+ monocyte', # LYZ+, FCGR3A+ (CD16)
'Pro B', # DNTT+, SOCS2+, CYGB+
'Late erythoid precursor 2', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
'Early erythoid precursor', # HBA1+, HBA2+, HBB+, HBD+, GYPA-
'Plasma B cell', # MZB1+, FKBP11+, SEC11C+, DERL3+, SSR4+
'Plasmacytoid dendritic cell',
'Pre B (cycling)', # IGLL1+, VPREB1+, PCDH9+, SOX4+, TCL1A+
'CD4+ T 3', # IL7R+, CD3D+, CD8A-
'Late erythoid precursor 3', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
'Megakaryocyte', # GNG11+, PPBP+, SDPR+, PF4+, HIST1H2AC+
'Stromal cell like', # C1QA+, C1QB+, SEPP1+, APOE+
'Stromal cell'] # CXCL12+ (stromal cell-derived factor 1 (SDF1))
adata.rename_categories('celltype', celltypes)

adata.obs['celltype_short'] = adata.obs['louvain']
celltypes2 = ['CD14+ Mono 2', # LYZ+, CD14+
'CD4+ T 1', # IL7R+, CD3D+, CD8A-
'CD14+ Mono 1', # LYZ+, CD14+ 
'NK', # NKG7+, GNLY+, CD3D-
'NK T 1', # NKG7+, GNLY+, CD3D+
'NK T 2', # NKG7+, GNLY-, CD3D+
'late Ery 1', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
'Imm B', # MS4A1+, HLA-DQA1+
'CD8+ T', # IL7R+, CD3D+, CD8A+
'Promono', # LYZ+, CD14-, FCGR3A-
'Pre B', # IGLL1+, VPREB1+, PCDH9+, SOX4+, TCL1A+
'HSC', # CD34+, SPINK2+
'cDC', # LYZ+, CD14-, FCGR3A-
'CD4+ T 2', # IL7R+, CD3D+, CD8A-
'CD16+ Mono', # LYZ+, FCGR3A+ (CD16)
'Pro B', # DNTT+, SOCS2+, CYGB+
'late Ery 2', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
'early Ery', # HBA1+, HBA2+, HBB+, HBD+, GYPA-
'Plasma B cell', # MZB1+, FKBP11+, SEC11C+, DERL3+, SSR4+
'pCD',
'Pre B (cc)', # IGLL1+, VPREB1+, PCDH9+, SOX4+, TCL1A+
'CD4+ T 3', # IL7R+, CD3D+, CD8A-
'late Ery 3', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
'Mk', # GNG11+, PPBP+, SDPR+, PF4+, HIST1H2AC+
'Stromal-like', # C1QA+, C1QB+, SEPP1+, APOE+
'Stromal'] # CXCL12+ (stromal cell-derived factor 1 (SDF1))
adata.rename_categories('celltype_short', celltypes2)

sc.pl.umap(adata, color='celltype_short', save='_celltypes_short.png')
sc.pl.umap(adata, color='celltype_short', legend_loc='on data', legend_fontsize=6, save='_celltypes_short_labels.png')



### Plot to UMAP top genes for cluster 24
sc.pl.umap(adata, color=adata.uns['rank_genes_groups']['names']['24'].tolist()[:10], save='_cluster24_10_topGenes.png')
sc.pl.stacked_violin(adata, groupby='celltype_short', var_names=adata.uns['rank_genes_groups']['names']['24'].tolist()[:10], swap_axes=True, save='_cluster24_10_topGenes.png')





### Score genesets
resdir = '/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/'
files = os.listdir(resdir)
files = [i for i in files if 'txt' in i if 'chen' in i]
for i in files:
	file = open(resdir + i, 'r')
	gs = file.read().split('\n')
	file.close()
	n = i[:-4]
	gs = adata.raw.var_names[adata.raw.var_names.isin(gs)].tolist()
	if len(gs) == 0:
		continue
	sc.tl.score_genes(adata, gene_list=gs, ctrl_size=len(gs), score_name=n, use_raw=True)
	savename = '_' + n + '.png'
	sc.pl.umap(adata, color=n, save=savename)

adata.write(results_file)



### Save
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






## Extract HSC and B lineage clusters and redo analysis workflow

### Get relevant cell IDs
celltypes = ['HSC', 'Pro B', 'Pre B (cc)', 'Pre B', 'Imm B']
ids = adata.obs.index[adata.obs['celltype_short'].isin(celltypes)]
ids = [re.sub('-MantonBM.$', '', id) for id in ids]


### Set new paths
wd = '/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/subset/'
results_file = wd + 'HCA_subset.h5ad'
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


### Read data (again)
fnames = ['/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/immune_census/cellranger.v3/MantonBM' + str(i) + '_HiSeq_' + str(j) + '/outs/filtered_feature_bc_matrix.h5' for i in range(1,9) for j in range(1,9)]
adata2 = sc.read_10x_h5(filename=fnames[0], genome='hg19')
name = re.search('MantonBM._HiSeq_.', fnames[0]).group(0)
adata2.obs_names = [name + '-' + obs_name for obs_name in adata2.obs_names]
for i in range(1,len(fnames)):
	if os.path.isfile(fnames[i]):
		tmp = sc.read_10x_h5(filename=fnames[i], genome='hg19')
		name = re.search('MantonBM._HiSeq_.', fnames[i]).group(0)
		tmp.obs_names = [name + '-' + obs_name for obs_name in tmp.obs_names]
		adata2 = adata2.concatenate(tmp, batch_key=name)



### Fix cell names
adata2.obs_names = [re.search(".*[TCGA]-1", i).group(0) for i in adata2.obs_names.tolist()]



### Remove useless metadata
for i in adata2.var.columns.tolist():
	del adata2.var[i]

for i in adata2.obs.columns.tolist():
	del adata2.obs[i]



### Subset HSC and B lineage cells
adata2 = adata2[ids]



### Set batch info
adata2.obs['batch'] = [re.sub("_.*", "", x) for x in adata2.obs_names]



### Save
adata2.write(results_file)



## Preprocessing

### Soft filtering for cells and genes
sc.pp.filter_cells(adata2, min_genes=1)
sc.pp.filter_genes(adata2, min_cells=100) # n_obs × n_vars = 24151 × 12649



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



### Filter low quality cells
#adata2 = adata2[adata2.obs['percent_mito'] < 0.10] # n_obs × n_vars = 170393 × 16439
#adata2 = adata2[adata2.obs['n_counts'] < 50000] # n_obs × n_vars = 170106 × 16439
#adata2 = adata2[adata2.obs['n_genes'] < 6000] # n_obs × n_vars = 170097 × 16439



### Filter genes again
#sc.pp.filter_genes(adata2, min_cells=100) # n_obs × n_vars = 170097 × 13528



### Log transform and store raw data for further use
adata2.raw = sc.pp.log1p(adata2, copy=True)



### Normalize (CPM) and get HVG
sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=1) # sum(adata2.var['highly_variable']) >>> 1226
sc.pl.highly_variable_genes(adata2, save='.pdf')



### Save
adata2.write(results_file)



### Filter by HVG
adata2 = adata2[:, adata2.var['highly_variable']]



### Regress out effect of n_counts and pct_mito
sc.pp.regress_out(adata2, ['n_counts', 'percent_mito'], n_jobs=n_jobs)
sc.pp.scale(adata2, max_value=10)



### Save
adata2.write(results_file)



### MNN batch correction
sample1 = adata2[adata2.obs['batch'] == 'MantonBM1']
sample2 = adata2[adata2.obs['batch'] == 'MantonBM2']
sample3 = adata2[adata2.obs['batch'] == 'MantonBM3']
sample4 = adata2[adata2.obs['batch'] == 'MantonBM4']
sample5 = adata2[adata2.obs['batch'] == 'MantonBM5']
sample6 = adata2[adata2.obs['batch'] == 'MantonBM6']
sample7 = adata2[adata2.obs['batch'] == 'MantonBM7']
sample8 = adata2[adata2.obs['batch'] == 'MantonBM8']

adata2_raw = adata2.raw

correction = sce.pp.mnn_correct(sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8, n_jobs=n_jobs, batch_categories=['MantonBM1', 'MantonBM2', 'MantonBM3', 'MantonBM4', 'MantonBM5', 'MantonBM6', 'MantonBM7', 'MantonBM8'])
adata2 = correction[0]

n_cells = adata2.var['n_cells-MantonBM1']
means = adata2.var['means-MantonBM1']
dispersions = adata2.var['dispersions-MantonBM1']
dispersions_norm = adata2.var['dispersions_norm-MantonBM1']

for i in adata2.var.columns.tolist():
	del adata2.var[i]

adata2.var['n_cells'] = n_cells
adata2.var['means'] = means
adata2.var['dispersions'] = dispersions
adata2.var['dispersions_norm'] = dispersions_norm

temp = ad.AnnData(X=adata2_raw.X, var=adata2_raw.var) # Can't assign AnnData.Raw object to adata2.raw directly. Stupid...
adata2.raw = temp

sc.pp.scale(adata2, max_value=10)

del correction, adata2_raw, temp



### Save
adata2.write(results_file)



### Compute PCA
sc.tl.pca(adata2, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata2, log=True, save='.png')
adata2.write(results_file)



### Compute neighborhood graph
sc.pp.neighbors(adata2, n_neighbors=30)



### Compute UMAP
sc.tl.umap(adata2)



### Save
adata2.write(results_file)



### Color markers and QC features (UMAP)
sc.pl.umap(adata2, color=['IL7R', 'NKG7', 'GNLY', 'CD3D', 'CD8A'], save='_T_cell_markers.png')
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
adata2.write(results_file)



### Cluster plots
sc.pl.umap(adata2, color='louvain', save='_clusters.png')
sc.pl.umap(adata2, color='louvain', legend_loc='on data', save='_clusters_labels.png')



### Finding marker genes
sc.tl.rank_genes_groups(adata2, 'louvain', method='wilcoxon')
sc.pl.rank_genes_groups(adata2, n_genes=20, save='.png')
sc.pl.rank_genes_groups_stacked_violin(adata2, n_genes=5, use_raw=True, swap_axes=True, save='_top5_markergenes_raw.png')
adata2.write(results_file)



### Cluster 7 is CD3D+ and IL32+ indicating T cells
### Cluster 10 is CD3D+, IL32+ and GZMA+ indicating NK T cells
### Cluster 11 is LYZ+, S100A8+, S100A8+
### Cluster 13 is HBA1+, HBA2+, HBB+ and HBD+ indicating erythrocyte precursors
### Remove these and redo analysis

clusters = ['7', '10', '11', '13']
ids = adata2.obs.index[~adata2.obs['louvain'].isin(clusters)]
ids = [re.sub('-MantonBM.$', '', id) for id in ids]



### Set new paths
wd = '/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/subset2/'
results_file = wd + 'HCA_subset2.h5ad'
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



### Read data (again)
fnames = ['/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/immune_census/cellranger.v3/MantonBM' + str(i) + '_HiSeq_' + str(j) + '/outs/filtered_feature_bc_matrix.h5' for i in range(1,9) for j in range(1,9)]
adata2 = sc.read_10x_h5(filename=fnames[0], genome='hg19')
name = re.search('MantonBM._HiSeq_.', fnames[0]).group(0)
adata2.obs_names = [name + '-' + obs_name for obs_name in adata2.obs_names]
for i in range(1,len(fnames)):
	if os.path.isfile(fnames[i]):
		tmp = sc.read_10x_h5(filename=fnames[i], genome='hg19')
		name = re.search('MantonBM._HiSeq_.', fnames[i]).group(0)
		tmp.obs_names = [name + '-' + obs_name for obs_name in tmp.obs_names]
		adata2 = adata2.concatenate(tmp, batch_key=name)



### Fix cell names
adata2.obs_names = [re.search(".*[TCGA]-1", i).group(0) for i in adata2.obs_names.tolist()]



### Remove useless metadata
for i in adata2.var.columns.tolist():
	del adata2.var[i]

for i in adata2.obs.columns.tolist():
	del adata2.obs[i]



### Subset HSC and B lineage cells
adata2 = adata2[ids]



### Set batch info
adata2.obs['batch'] = [re.sub("_.*", "", x) for x in adata2.obs_names]



### Save
adata2.write(results_file)



## Preprocessing

### Soft filtering for cells and genes
sc.pp.filter_cells(adata2, min_genes=1)
sc.pp.filter_genes(adata2, min_cells=100) # n_obs × n_vars = 22324 × 12428



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



### Log transform and store raw data for further use
adata2.raw = sc.pp.log1p(adata2, copy=True)



### Normalize (CPM) and get HVG
sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2, min_mean=0.1, max_mean=3, min_disp=0.5) # sum(adata2.var['highly_variable']) >>> 1974
sc.pl.highly_variable_genes(adata2, save='.pdf')



### Save
adata2.write(results_file)



### Filter by HVG
adata2 = adata2[:, adata2.var['highly_variable']]



### Regress out effect of n_counts and pct_mito
sc.pp.regress_out(adata2, ['n_counts', 'percent_mito'], n_jobs=n_jobs)
sc.pp.scale(adata2, max_value=10)



### Save
adata2.write(results_file)



### MNN batch correction
sample1 = adata2[adata2.obs['batch'] == 'MantonBM1']
sample2 = adata2[adata2.obs['batch'] == 'MantonBM2']
sample3 = adata2[adata2.obs['batch'] == 'MantonBM3']
sample4 = adata2[adata2.obs['batch'] == 'MantonBM4']
sample5 = adata2[adata2.obs['batch'] == 'MantonBM5']
sample6 = adata2[adata2.obs['batch'] == 'MantonBM6']
sample7 = adata2[adata2.obs['batch'] == 'MantonBM7']
sample8 = adata2[adata2.obs['batch'] == 'MantonBM8']

adata2_raw = adata2.raw

correction = sce.pp.mnn_correct(sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8, n_jobs=n_jobs, batch_categories=['MantonBM1', 'MantonBM2', 'MantonBM3', 'MantonBM4', 'MantonBM5', 'MantonBM6', 'MantonBM7', 'MantonBM8'])
adata2 = correction[0]

n_cells = adata2.var['n_cells-MantonBM1']
means = adata2.var['means-MantonBM1']
dispersions = adata2.var['dispersions-MantonBM1']
dispersions_norm = adata2.var['dispersions_norm-MantonBM1']

for i in adata2.var.columns.tolist():
	del adata2.var[i]

adata2.var['n_cells'] = n_cells
adata2.var['means'] = means
adata2.var['dispersions'] = dispersions
adata2.var['dispersions_norm'] = dispersions_norm

temp = ad.AnnData(X=adata2_raw.X, var=adata2_raw.var) # Can't assign AnnData.Raw object to adata2.raw directly. Stupid...
adata2.raw = temp

sc.pp.scale(adata2, max_value=10)

del correction, adata2_raw, temp



### Save
adata2.write(results_file)



### Compute PCA
sc.tl.pca(adata2, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata2, log=True, save='.png')
adata2.write(results_file)



### Compute neighborhood graph
sc.pp.neighbors(adata2, n_neighbors=15, n_pcs=20)



### Compute UMAP
sc.tl.umap(adata2)



### Save
adata2.write(results_file)



### Color markers and QC features (UMAP)
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
adata2.write(results_file)



### Cluster plots
sc.pl.umap(adata2, color='louvain', save='_clusters.png')
sc.pl.umap(adata2, color='louvain', legend_loc='on data', save='_clusters_labels.png')



## Finding marker genes
sc.tl.rank_genes_groups(adata2, 'louvain', method='wilcoxon')
sc.pl.rank_genes_groups(adata2, n_genes=20, save='.png')
sc.pl.rank_genes_groups_stacked_violin(adata2, n_genes=5, use_raw=True, swap_axes=True, save='_top5_markergenes_raw.png')
#adata2.write(results_file)



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



### Plot QC violin plots per cluster
sc.pl.violin(adata2, 'n_genes', groupby='louvain', stripplot=False, save='_n_genes_bygroup.png')
sc.pl.violin(adata2, 'n_counts', groupby='louvain', stripplot=False, save='_n_counts_bygroup.png')
sc.pl.violin(adata2, 'percent_mito', groupby='louvain', stripplot=False, save='_percent_mito_bygroup.png')
sc.pl.violin(adata2, 'percent_ribo', groupby='louvain', stripplot=False, save='_percent_ribo_bygroup.png')
sc.pl.scatter(adata2, x='n_counts', y='percent_mito', save='_nCounts_pctMito.pdf')



### Calculate MAD (Median Absolute Deviation) for various features and plot
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

n_genes_mad = calculate_MAD(adata=adata2, groupby='louvain', variable='n_genes')
pct_mito_mad = calculate_MAD(adata=adata2, groupby='louvain', variable='percent_mito')
n_counts_mad = calculate_MAD(adata=adata2, groupby='louvain', variable='n_counts')

adata2.obs['groupby_n_genes_MAD_diff_x'] = n_genes_mad['groupby_MAD_diff_x']
adata2.obs['groupby_n_counts_MAD_diff_x'] = n_counts_mad['groupby_MAD_diff_x']
adata2.obs['groupby_percent_mito_MAD_diff_x'] = pct_mito_mad['groupby_MAD_diff_x']

sc.pl.umap(adata2, color='groupby_n_genes_MAD_diff_x', save='_MAD_n_genes_diff_x.png')
sc.pl.umap(adata2, color='groupby_n_counts_MAD_diff_x', save='_MAD_n_counts_diff_x.png')
sc.pl.umap(adata2, color='groupby_percent_mito_MAD_diff_x', save='_MAD_percent_mito_diff_x.png')

sc.pl.violin(adata2, 'groupby_n_genes_MAD_diff_x', groupby='louvain', stripplot=False, save='_MAD_n_genes_diff_x.png')
sc.pl.violin(adata2, 'groupby_n_counts_MAD_diff_x', groupby='louvain', stripplot=False, save='_MAD_n_counts_diff_x.png')
sc.pl.violin(adata2, 'groupby_percent_mito_MAD_diff_x', groupby='louvain', stripplot=False, save='_MAD_percent_mito_diff_x.png')

ax1 = sc.pl.violin(adata2, 'groupby_n_genes_MAD_diff_x', groupby='louvain', stripplot=False, show=False)
x, y = ax1.get_xlim()
ax1.hlines(5, x, y, colors='red')
plt.savefig('figures/violinplot_MAD_n_genes_diff_x_hline_5.pdf', dpi=300, bbox_inches='tight')

ax2 = sc.pl.violin(adata2, 'groupby_n_counts_MAD_diff_x', groupby='louvain', stripplot=False, show=False)
x, y = ax2.get_xlim()
ax2.hlines(5, x, y, colors='red')
plt.savefig('figures/violinplot_MAD_n_counts_diff_x_hline_5.pdf', dpi=300, bbox_inches='tight')

ax3 = sc.pl.violin(adata2, 'groupby_percent_mito_MAD_diff_x', groupby='louvain', stripplot=False, show=False)
x, y = ax3.get_xlim()
ax3.hlines(5, x, y, colors='red')
plt.savefig('figures/violinplot_MAD_percent_mito_diff_x_hline_5.pdf', dpi=300, bbox_inches='tight')

adata2.obs['groupby_n_genes_MAD_diff_x_gt_3'] = (adata2.obs['groupby_n_genes_MAD_diff_x'] > 3).astype('category')
adata2.obs['groupby_percent_mito_MAD_diff_x_gt_3'] = (adata2.obs['groupby_percent_mito_MAD_diff_x'] > 3).astype('category')
adata2.obs['groupby_n_counts_MAD_diff_x_gt_3'] = (adata2.obs['groupby_n_counts_MAD_diff_x'] > 3).astype('category')

sc.pl.umap(adata2, color='groupby_n_genes_MAD_diff_x_gt_3', save='_MAD_n_genes_diff_x_gt_3.png')
sc.pl.umap(adata2, color='groupby_n_counts_MAD_diff_x_gt_3', save='_MAD_n_counts_diff_x_gt_3.png')
sc.pl.umap(adata2, color='groupby_percent_mito_MAD_diff_x_gt_3', save='_MAD_percent_mito_diff_x_gt_3.png')
sc.pl.umap(adata2, color=['groupby_n_genes_MAD_diff_x_gt_3', 'groupby_n_counts_MAD_diff_x_gt_3', 'groupby_percent_mito_MAD_diff_x_gt_3'], save='_MAD_diff_x_gt_3.png')

adata2.obs['groupby_n_genes_MAD_diff_x_gt_5'] = (adata2.obs['groupby_n_genes_MAD_diff_x'] > 5).astype('category')
adata2.obs['groupby_percent_mito_MAD_diff_x_gt_5'] = (adata2.obs['groupby_percent_mito_MAD_diff_x'] > 5).astype('category')
adata2.obs['groupby_n_counts_MAD_diff_x_gt_5'] = (adata2.obs['groupby_n_counts_MAD_diff_x'] > 5).astype('category')

sc.pl.umap(adata2, color='groupby_n_genes_MAD_diff_x_gt_5', save='_MAD_n_genes_diff_x_gt_5.png')
sc.pl.umap(adata2, color='groupby_n_counts_MAD_diff_x_gt_5', save='_MAD_n_counts_diff_x_gt_5.png')
sc.pl.umap(adata2, color='groupby_percent_mito_MAD_diff_x_gt_5', save='_MAD_percent_mito_diff_x_gt_5.png')
sc.pl.umap(adata2, color=['groupby_n_genes_MAD_diff_x_gt_5', 'groupby_n_counts_MAD_diff_x_gt_5', 'groupby_percent_mito_MAD_diff_x_gt_5'], save='_MAD_diff_x_gt_5.png')



### Filter distinct outliers further with mito% and n_counts MAD difference > 5
ids = adata2.obs.index[~(adata2.obs['groupby_n_counts_MAD_diff_x_gt_5'] | adata2.obs['groupby_percent_mito_MAD_diff_x_gt_5'])]
ids = [re.sub('-MantonBM.$', '', id) for id in ids]



### Set new paths
wd = '/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/subset3/'
results_file = wd + 'HCA_subset3.h5ad'
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



### Read data (again)
fnames = ['/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/immune_census/cellranger.v3/MantonBM' + str(i) + '_HiSeq_' + str(j) + '/outs/filtered_feature_bc_matrix.h5' for i in range(1,9) for j in range(1,9)]
adata2 = sc.read_10x_h5(filename=fnames[0], genome='hg19')
name = re.search('MantonBM._HiSeq_.', fnames[0]).group(0)
adata2.obs_names = [name + '-' + obs_name for obs_name in adata2.obs_names]
for i in range(1,len(fnames)):
	if os.path.isfile(fnames[i]):
		tmp = sc.read_10x_h5(filename=fnames[i], genome='hg19')
		name = re.search('MantonBM._HiSeq_.', fnames[i]).group(0)
		tmp.obs_names = [name + '-' + obs_name for obs_name in tmp.obs_names]
		adata2 = adata2.concatenate(tmp, batch_key=name)



### Fix cell names
adata2.obs_names = [re.search(".*[TCGA]-1", i).group(0) for i in adata2.obs_names.tolist()]



### Remove useless metadata
for i in adata2.var.columns.tolist():
	del adata2.var[i]

for i in adata2.obs.columns.tolist():
	del adata2.obs[i]



### Subset HSC and B lineage cells
adata2 = adata2[ids]



### Set batch info
adata2.obs['batch'] = [re.sub("_.*", "", x) for x in adata2.obs_names]



### Log transform and store raw data for further use
adata2.raw = sc.pp.log1p(adata2, copy=True)


### Save
adata2.write(results_file)



## Preprocessing

### Soft filtering for cells and genes
sc.pp.filter_cells(adata2, min_genes=1)

d = pd.DataFrame(adata2.X.toarray())
d.to_csv(results_folder_raw + 'X_allgenes.csv.gz', sep=',', header=False, index=False, compression='gzip')

d = pd.DataFrame(adata2.var.index)
d.to_csv(results_folder_raw + 'var_allgenes.txt', sep=',', header=False, index=False)

d = pd.DataFrame(adata2.obs.index)
d.to_csv(results_folder_raw + 'obs_allgenes.txt', sep=',', header=False, index=False)

del d

sc.pp.filter_genes(adata2, min_cells=100) # n_obs × n_vars = 20753 × 12198



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



### Normalize (CPM) and get HVG
sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2, min_mean=0.1, max_mean=3, min_disp=0.75) # sum(adata2.var['highly_variable']) >>> 1218
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



### MNN batch correction
sample1 = adata2[adata2.obs['batch'] == 'MantonBM1']
sample2 = adata2[adata2.obs['batch'] == 'MantonBM2']
sample3 = adata2[adata2.obs['batch'] == 'MantonBM3']
sample4 = adata2[adata2.obs['batch'] == 'MantonBM4']
sample5 = adata2[adata2.obs['batch'] == 'MantonBM5']
sample6 = adata2[adata2.obs['batch'] == 'MantonBM6']
sample7 = adata2[adata2.obs['batch'] == 'MantonBM7']
sample8 = adata2[adata2.obs['batch'] == 'MantonBM8']

adata2_raw = adata2.raw

correction = sce.pp.mnn_correct(sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8, n_jobs=n_jobs, batch_categories=['MantonBM1', 'MantonBM2', 'MantonBM3', 'MantonBM4', 'MantonBM5', 'MantonBM6', 'MantonBM7', 'MantonBM8'])
adata2 = correction[0]

n_cells = adata2.var['n_cells-MantonBM1']
means = adata2.var['means-MantonBM1']
dispersions = adata2.var['dispersions-MantonBM1']
dispersions_norm = adata2.var['dispersions_norm-MantonBM1']

for i in adata2.var.columns.tolist():
	del adata2.var[i]

adata2.var['n_cells'] = n_cells
adata2.var['means'] = means
adata2.var['dispersions'] = dispersions
adata2.var['dispersions_norm'] = dispersions_norm

temp = ad.AnnData(X=adata2_raw.X, var=adata2_raw.var) # Can't assign AnnData.Raw object to adata2.raw directly. Stupid...
adata2.raw = temp

sc.pp.scale(adata2, max_value=10)

del correction, adata2_raw, temp



### Save
#adata2.write(results_file)



### Compute PCA
sc.tl.pca(adata2, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata2, log=True, save='.png')
#adata2.write(results_file)



### Compute neighborhood graph
sc.pp.neighbors(adata2, n_neighbors=15, n_pcs=20)



### Compute UMAP
sc.tl.umap(adata2)



### Save
#adata2.write(results_file)



### Color markers and QC features (UMAP)
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



### Cluster plots
sc.pl.umap(adata2, color='louvain', save='_clusters.png')
sc.pl.umap(adata2, color='louvain', legend_loc='on data', save='_clusters_labels.png')



## Finding marker genes
sc.tl.rank_genes_groups(adata2, groupby='louvain', n_genes=2000, method='wilcoxon')
sc.pl.rank_genes_groups(adata2, n_genes=20, save='.png')
sc.pl.rank_genes_groups_stacked_violin(adata2, n_genes=5, use_raw=True, swap_axes=True, save='_top5_markergenes_raw.png')
#adata2.write(results_file)



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



### Plot QC violin plots per cluster
sc.pl.violin(adata2, 'n_genes', groupby='louvain', stripplot=False, save='_n_genes_bygroup.png')
sc.pl.violin(adata2, 'n_counts', groupby='louvain', stripplot=False, save='_n_counts_bygroup.png')
sc.pl.violin(adata2, 'percent_mito', groupby='louvain', stripplot=False, save='_percent_mito_bygroup.png')
sc.pl.violin(adata2, 'percent_ribo', groupby='louvain', stripplot=False, save='_percent_ribo_bygroup.png')
sc.pl.scatter(adata2, x='n_counts', y='percent_mito', save='_nCounts_pctMito.pdf')



# ### Plot doublet scores calculated with Scrublet
# scrublet = pd.read_csv(wd + 'scrublet/scrublet_results.csv', index_col=0, sep=',')
# scrublet.index = scrublet.index + '-' + adata2.obs['batch'].tolist()
# adata2.obs['scrublet_doublet_score'] = scrublet['doublet_scores']
# adata2.obs['scrublet_predicted_doublet'] = scrublet['predicted_doublets'].astype('category')
# sc.pl.umap(adata2, color=['scrublet_doublet_score', 'scrublet_predicted_doublet'], save='_scrublet_results.png')



# ### Plot doublets identified with DoubletDetection
# dd = pd.read_csv(wd + 'doubletDetection/doubletDetection_results.csv', index_col=0, sep=',')
# dd.index = dd.index + '-' + adata2.obs['batch'].tolist()
# adata2.obs['dd_predicted_doublet'] = dd['predicted_doublet'].astype('bool').astype('category')
# sc.pl.umap(adata2, color='dd_predicted_doublet', save='_DoubletDetection_results.png')


### Plot subset3 clustering to full HCA UMAP
adata = sc.read('/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HCA.h5ad')
adata.obs['louvain_subset3'] = adata2.obs['louvain']
sc.pl.umap(adata, color='louvain_subset3', save='_clusters_in_full_HCA.png')
sc.pl.umap(adata, color='louvain_subset3', legend_loc='on data', save='_clusters_in_full_HCA_labels.png')
del adata



### Cell types
mapping = {
0:'6_preB_G1',
1:'7_immatureB',
2:'7_immatureB',
3:'1_HSC',
4:'7_immatureB',
5:'4_proB_G1',
6:'5_preB_G2MS',
7:'1_HSC',
8:'7_immatureB',
9:'2_earlyLymphoid',
10:'3_proB_G2MS'
}
adata2.obs['celltype'] = adata2.obs['louvain'].astype('int32').map(mapping).astype('category')

sc.pl.umap(adata2, color='celltype', legend_loc='on data', legend_fontsize=8, save='_celltypes.png')



### Save
adata2.write(results_file)



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



### Calculate diffusion map
sc.tl.diffmap(adata2)
sc.pl.diffmap(adata2, color='celltype', components=['1,2','1,3'], save='_celltypes.png')



### Simplify celltypes for label transfer and add specific labels for B lineage
celltypes_labelTransfer = {
0: 'Monocyte', # LYZ+, CD14+
1: 'T cell', # IL7R+, CD3D+, CD8A-
2: 'Monocyte', # LYZ+, CD14+ 
3: 'NK', # NKG7+, GNLY+, CD3D-
4: 'T cell', # NKG7+, GNLY+, CD3D+
5: 'T cell', # NKG7+, GNLY-, CD3D+
6: 'Erythroid precursor cell', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
7: 'B cell', # MS4A1+, HLA-DQA1+
8: 'T cell', # IL7R+, CD3D+, CD8A+
9: 'Monocyte', # LYZ+, CD14-, FCGR3A-
10: 'B cell', # IGLL1+, VPREB1+, PCDH9+, SOX4+, TCL1A+
11: 'HSC', # CD34+, SPINK2+
12: 'Dendritic cell', # LYZ+, CD14-, FCGR3A-
13: 'T cell', # IL7R+, CD3D+, CD8A-
14: 'Monocyte', # LYZ+, FCGR3A+ (CD16)
15: 'B cell', # DNTT+, SOCS2+, CYGB+
16: 'Erythroid precursor cell', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
17: 'Erythroid precursor cell', # HBA1+, HBA2+, HBB+, HBD+, GYPA-
18: 'Plasma B cell', # MZB1+, FKBP11+, SEC11C+, DERL3+, SSR4+
19: 'Dendritic cell',
20: 'B cell', # IGLL1+, VPREB1+, PCDH9+, SOX4+, TCL1A+
21: 'T cell', # IL7R+, CD3D+, CD8A-
22: 'Erythroid precursor cell', # HBA1+, HBA2+, HBB+, HBD+, GYPA+
23: 'Megakaryocyte', # GNG11+, PPBP+, SDPR+, PF4+, HIST1H2AC+
24: 'Stromal-like cell', # C1QA+, C1QB+, SEPP1+, APOE+
25: 'Stromal-like cell' # CXCL12+ (stromal cell-derived factor 1 (SDF1))
} 
adata.obs['celltype_labelTransfer'] = adata.obs['louvain'].astype('int32').map(celltypes_labelTransfer).astype('category')
subset_celltypes = pd.DataFrame({'celltype_labelTransfer': adata2.obs['celltype']})
adata.obs.update(subset_celltypes)

adata.write_csvs(dirname='/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/', skip_data=True)
adata.write('/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HCA.h5ad')



### Plot Hay et al annotations
# hay = pd.read_csv('/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/Hay_and_own_annotation-IDs.txt', sep='\t', index_col=3)
# adata2.obs['Hay_annotation'] = hay['Hay_annotation']
# sc.pl.umap(adata2, color='Hay_annotation', save='_Hay_annotation.pdf')
# sc.pl.umap(adata2, color='Hay_annotation', legend_loc='on data', legend_fontsize=4, save='_Hay_annotation2.pdf')

# adata.obs['Hay_annotation'] = hay['Hay_annotation']
# sc.pl.umap(adata, color='Hay_annotation', save='_Hay_annotation.pdf')
# sc.pl.umap(adata, color='Hay_annotation', legend_loc='on data', legend_fontsize=3, save='_Hay_annotation2.pdf')



### Plot mito% and figure out if it varies too much between celltypes to have an effect which should be removed
sc.pl.umap(adata2, color=['percent_mito', 'louvain', 'celltype'], save='_percent_mito.png')
sc.pl.violin(adata2, 'percent_mito', groupby='celltype', stripplot=False, save='_percent_mito_by_celltype.png')
sc.pl.violin(adata2, 'percent_mito', groupby='louvain', stripplot=False, save='_percent_mito_by_cluster.png')

#adata2.layers['sct_noBatch'] = pd.read_csv('/research/groups/allseq/data/scRNAseq/results/juha_wrk/normalization_DE_tests/vst_noBatch_correctedUMI.csv').T.values
#adata2.layers['sct_Batch'] = pd.read_csv('/research/groups/allseq/data/scRNAseq/results/juha_wrk/normalization_DE_tests/vst_Batch_correctedUMI.csv').T



# Stats by celltype
sc.pl.violin(adata2, keys=['n_counts', 'n_genes', 'percent_mito'], groupby='celltype', stripplot=False, rotation=90, save='_QC_by_celltype.png')



# n_counts & n_genes per phase per celltype
sc.pl.violin(adata2[adata2.obs['celltype'] == '1_HSC'], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_HSC_QC_by_phase.png')
sc.pl.violin(adata2[adata2.obs['celltype'] == '2_earlyLymphoid'], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_earlyLymphoid_QC_by_phase.png')
sc.pl.violin(adata2[adata2.obs['celltype'] == '3_proB_G2MS'], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_proB_G2MS_QC_by_phase.png')
sc.pl.violin(adata2[adata2.obs['celltype'] == '4_proB_G1'], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_proB_G1_QC_by_phase.png')
sc.pl.violin(adata2[adata2.obs['celltype'] == '5_preB_G2MS'], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_preB_G2MS_QC_by_phase.png')
sc.pl.violin(adata2[adata2.obs['celltype'] == '6_preB_G1'], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_preB_G1_QC_by_phase.png')
sc.pl.violin(adata2[adata2.obs['celltype'] == '7_immatureB'], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_immatureB_QC_by_phase.png')
sc.pl.violin(adata2[adata2.obs['celltype'].isin(['3_proB_G2MS', '4_proB_G1'])], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_proB_QC_by_phase.png')
sc.pl.violin(adata2[adata2.obs['celltype'].isin(['5_preB_G2MS', '6_preB_G1'])], keys=['n_counts', 'n_genes'], groupby='phase', stripplot=False, save='_preB_QC_by_phase.png')



# Downsample
# adata_ds = ad.AnnData(X=adata2.raw.X, var=adata2.raw.var, obs=adata2.obs)
# sc.pp.downsample_counts(adata_ds, counts_per_cell=6000) # 6000 is ~80% from median of proB G2M/S n_counts
# d = pd.DataFrame(adata_ds.X.toarray())
# d.to_csv('/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/subset3/data/downsampled/X_n6000.csv.gz', sep=',', header=False, index=False, compression='gzip')
# del d



