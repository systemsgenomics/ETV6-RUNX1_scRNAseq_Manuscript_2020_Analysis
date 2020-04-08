# Run pySCENIC for HCA B lineage cells from scanpy analysis.

## Run with initially downsampling data to include the same number of different cell types
## Add jitter to zeros to see if repressive regulons are detected better.


import os
import glob
import pickle
import pandas as pd
import numpy as np
import anndata as ad
import seaborn as sns

from scipy.stats import pearsonr

from dask.diagnostics import ProgressBar

from distributed import Client, LocalCluster

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize



## Set paths
data_folder = '/research/groups/allseq/data/scRNAseq/results/juha_wrk/SCENIC/HCA_jitter/'
db_folder = os.path.join('/research/work/jmehtone/cisTarget_databases', 'hg19-*.feather')
motif_annotations = os.path.join('/research/work/jmehtone', 'motifs-v9-nr.hgnc-m0.001-o0.0.tbl')
TFs_file = os.path.join('/research/work/jmehtone', 'HGNC_TFs.txt')

if not os.path.exists(data_folder):
	os.makedirs(data_folder)

final_regulons_fname = os.path.join(data_folder, 'final_regulons.p')
final_regulons_aucell_fname = os.path.join(data_folder, 'aucell_final_regulons.csv.gz')

## Load data
adata = ad.read_h5ad('/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/subset3/HCA_subset3.h5ad')
data = pd.DataFrame(adata.raw.X.toarray(), index=adata.obs.index, columns=adata.raw.var.index)
metadata = adata.obs
del adata

## Add jitter
np.random.seed(1)
jitter = np.random.uniform(low=-0.01, high=0.01, size=data.shape)
jitter[data != 0] = 0
data += jitter

## Run parameters
iterations = 10
n_cores = 4
train_pct = 0.7
cor_p_thr = 0.001
memory_limit = 50e9
n_cells = 600
grouping_variable = 'celltype'


if __name__ == '__main__':

	if sum(metadata.columns == grouping_variable) < 1:
		exit('Grouping variable not found in metadata.')

	## Load randing databases
	db_fnames = glob.glob(db_folder)
	def name(fname):
		return os.path.basename(fname).split(".")[0]
	dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
	dbs


	## Initialize cluster
	local_cluster = LocalCluster(n_workers=n_cores, threads_per_worker=1, processes=False, memory_limit=memory_limit)
	custom_client = Client(local_cluster)


	## Load TFs
	tf_names = load_tf_names(TFs_file)


	## Collect here regulons passing correlation filter
	cortest_passed_regulons = []


	for i in range(0, iterations):

		## Split to train and test
		data[grouping_variable] = metadata[grouping_variable]
		data_sampled = data.groupby(grouping_variable).apply(lambda x: x.sample(n=min(n_cells, len(x)), random_state=i))
		data_sampled.index = data_sampled.index.get_level_values(1)
		data_train = data_sampled.groupby(grouping_variable).apply(lambda x: x.sample(n=round(len(x) * train_pct), random_state=i))
		data_train.index = data_train.index.get_level_values(1)
		data_test = data_sampled[~data_sampled.index.isin(data_train.index)]
		del data_train[grouping_variable], data_test[grouping_variable]



		## Paths for results
		data_folder_iter = os.path.join(data_folder, 'iter_' + str(i))
		network_fname = os.path.join(data_folder_iter, 'network.csv.gz')
		modules_fname = os.path.join(data_folder_iter, 'modules.p')
		motifs_fname = os.path.join(data_folder_iter, 'motifs.csv')
		regulons_fname = os.path.join(data_folder_iter, 'regulons.p')
		aucell_train_fname = os.path.join(data_folder_iter, 'aucell_train_scores.csv.gz')
		aucell_test_fname = os.path.join(data_folder_iter, 'aucell_test_scores.csv.gz')
		if not os.path.exists(data_folder_iter):
			os.makedirs(data_folder_iter)

		os.chdir(data_folder_iter)



		## Run GRNBoost2 (faster equivalent of GENIE3) from arboreto to infer co-expression modules
		if not os.path.isfile(network_fname):
			adjacencies = grnboost2(data_train, tf_names=tf_names, verbose=True, client_or_address=custom_client, seed=i)
			adjacencies.to_csv(network_fname, sep=',', header=True, index=False, compression='gzip')
		else:
			adjacencies = pd.read_csv(network_fname)



		## Derive potential regulons from co-expression modules
		if not os.path.isfile(modules_fname):
			modules = list(modules_from_adjacencies(adjacencies, data_train, keep_only_activating=False))
			pickle.dump(modules, open(modules_fname, 'wb'))
		else:
			modules = pickle.load(open(modules_fname, 'rb'))
		
		del adjacencies



		## Prune modules for targets with cis regulatory footprints (aka RcisTarget)

		### Calculate a list of enriched motifs and the corresponding target genes for all modules.
		if not os.path.isfile(motifs_fname):
			df = prune2df(dbs, modules, motif_annotations, num_workers=n_cores)
			df.to_csv(motifs_fname)
		else:
			df = pd.read_csv(motifs_fname)
		
		del modules



		### Create regulons from this table of enriched motifs.
		if not os.path.isfile(regulons_fname):
			regulons = df2regulons(df)
			pickle.dump(regulons, open(regulons_fname, 'wb'))
		else:
			regulons = pickle.load(open(regulons_fname, 'rb'))

		del df



		## Cellular regulon enrichment matrices
		if not os.path.isfile(aucell_train_fname):
			auc_train = aucell(data_train, regulons, num_workers=n_cores)
			auc_train.to_csv(aucell_train_fname, sep=',', header=True, index=True, compression='gzip')
		else:
			auc_train = pd.read_csv(aucell_train_fname, index_col=0)

		if not os.path.isfile(aucell_test_fname):
			auc_test = aucell(data_test, regulons, num_workers=n_cores)
			auc_test.to_csv(aucell_test_fname, sep=',', header=True, index=True, compression='gzip')
		else:
			auc_test = pd.read_csv(aucell_test_fname, index_col=0)



		## Filter regulons with low correlation between train and test
		auc_train[grouping_variable] = metadata[grouping_variable]
		auc_test[grouping_variable] = metadata[grouping_variable]

		auc_train_mean = auc_train.groupby(grouping_variable).mean()
		auc_test_mean = auc_test.groupby(grouping_variable).mean()

		cortest = pd.Series(np.NaN, index=auc_train_mean.columns)
		for reg in cortest.index:
			cortest[reg] = pearsonr(auc_train_mean[reg], auc_test_mean[reg])[1]

		print('Filtering {} out of {} regulons for not passing correlation p-value threshold {}!'.format(sum(cortest >= cor_p_thr), len(cortest), cor_p_thr))
		regulon_names_filtered = cortest[cortest < cor_p_thr].index

		passed_regulons = [x for x in regulons for y in regulon_names_filtered if x.name == y]

		cortest_passed_regulons.extend(passed_regulons)



	## Save all regulons to file
	pickle.dump(cortest_passed_regulons, open(final_regulons_fname, 'wb'))



	## Calculate AUCell for whole dataset with all found regulons
	adata = ad.read_h5ad('/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/subset3/HCA_subset3.h5ad')
	adata2 = ad.read_h5ad('/research/groups/allseq/data/scRNAseq/results/juha_wrk/ALL_scanpy/new_2019_ver5/subset/ALL_subset.h5ad')
	adata2 = adata2[adata2.obs['celltype'] == 'leukemic']
	X1 = pd.DataFrame(adata.raw.X.toarray(), index=adata.obs.index, columns=adata.raw.var.index)
	X2 = pd.DataFrame(adata2.raw.X.toarray(), index=adata2.obs.index, columns=adata2.raw.var.index)
	data = pd.concat([X1, X2], join='inner')

	auc_mtx = aucell(data, cortest_passed_regulons, num_workers=n_cores)
	auc_mtx.to_csv(final_regulons_aucell_fname, sep=',', header=True, index=True, compression='gzip')



	## Calculate mean score per regulon from multiple iterations
	auc_mtx = auc_mtx.T
	auc_mtx['regulon'] = [re.sub('[.].*', '', i) for i in auc_mtx.index]
	auc_mtx_mean = auc_mtx.groupby('regulon').mean().T
	auc_mtx_mean.to_csv(final_regulons_aucell_means_fname, sep=',', header=True, index=True, compression='gzip')



	## Binarize mean scores
	bin_mtx, _ = binarize(auc_mtx_mean, num_workers=n_cores)
	bin_mtx.to_csv(final_regulons_aucell_means_bin_fname, sep=',', header=True, index=True, compression='gzip')










