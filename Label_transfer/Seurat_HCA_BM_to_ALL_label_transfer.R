# Label transfer HCA BM --> ALL

options(stringsAsFactors = F)
library(Seurat)


# ALL data
X <- data.table::fread(file = "/research/groups/allseq/data/scRNAseq/results/juha_wrk/ALL_scanpy/new_2019_ver5/subset/data/raw/X.csv.gz", data.table = F)
colnames(X) <- readLines("/research/groups/allseq/data/scRNAseq/results/juha_wrk/ALL_scanpy/new_2019_ver5/subset/data/raw/var.txt")
metadata.X <- read.csv("/research/groups/allseq/data/scRNAseq/results/juha_wrk/ALL_scanpy/new_2019_ver5/subset/data/obs.csv", stringsAsFactors = F)
rownames(X) <- rownames(metadata.X) <- metadata.X$index <- gsub("-1", "", metadata.X$index)
metadata.X$louvain <- paste0("ALL.", metadata.X$louvain)
metadata.X <- metadata.X[, c("batch", "louvain", "celltype", "phase")]

# HCA data
Y <- data.table::fread(file = "/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/raw/X.csv.gz", data.table = F)
colnames(Y) <- readLines("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/raw/var.txt")
metadata.Y <- read.csv("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/obs.csv", stringsAsFactors = F)
rownames(Y) <- rownames(metadata.Y) <- metadata.Y$index
metadata.Y$louvain <- paste0("HCA.", metadata.Y$louvain)
metadata.Y <- metadata.Y[, c("batch", "louvain", "celltype_labelTransfer", "phase")]
colnames(metadata.Y) <- c("batch", "louvain", "celltype", "phase")
metadata.Y <- metadata.Y[! metadata.Y$celltype %in% c("B cell", "HSC"), ] # Remove some cells which were filtered out from final HSC + B lineage cell set
Y <- Y[rownames(metadata.Y), ]

X <- as(t(X), "dgCMatrix")
Y <- as(t(Y), "dgCMatrix")

HCA <- CreateSeuratObject(counts = Y, meta.data = metadata.Y)
HCA <- NormalizeData(object = HCA, verbose = T)
HCA <- FindVariableFeatures(object = HCA, selection.method = "vst", nfeatures = 2000, verbose = T)

ALL <- CreateSeuratObject(counts = X, meta.data = metadata.X)
ALL.list <- SplitObject(object = ALL, split.by = "batch")

for(i in 1:length(ALL.list)) {
  ALL.list[[i]] <- NormalizeData(object = ALL.list[[i]], verbose = T)
  ALL.list[[i]] <- FindVariableFeatures(object = ALL.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = T)
}

for(i in 1:length(ALL.list)) {
  batchname <- unique(ALL.list[[i]]$batch)
  fname <- paste0("ALL_Seurat_obj_", batchname, ".rds")
  anchors <- FindTransferAnchors(reference = HCA, query = ALL.list[[i]], dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = HCA$celltype, dims = 1:30)
  ALL.list[[i]] <- AddMetaData(object = ALL.list[[i]], metadata = predictions)
  saveRDS(ALL.list[[i]], file = fname)
}
