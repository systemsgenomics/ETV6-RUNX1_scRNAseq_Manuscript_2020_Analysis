# Filtering and visualization of SCENIC HCA B lineage + ALL result 

options(stringsAsFactors = F)

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggpubr)
#library(ggrastr)

setwd("/research/groups/allseq/data/scRNAseq/results/juha_wrk/SCENIC/SCENIC_analysis/FINAL/HCA_ALL_FINAL/")



# Helper functions

## Mean across variables based label vector
regulonMeanScore <- function(data, regulons) {
  df <- data.frame(t(data), regulons, check.names = F)
  df <- df %>% group_by(regulons) %>% summarise_all(mean)
  df <- as.data.frame(df)
  rownames(df) <- df[, 1]
  df[, -1]
}


## Linear model fitting through subsampling. Return
linearFitPermutation <- function(data, labels, n_cells = 600, iterations = 100) {
  
  df <- data.frame(data, label = labels, check.names = F)
  
  res <- do.call(cbind, lapply(1:iterations, function(i) {
    set.seed(i)
    df.sampled <- df %>% group_by(label) %>% sample_n(n_cells)
    df.sampled <- as.data.frame(df.sampled)
    labels.sampled <- df.sampled$label
    df.sampled <- df.sampled[, colnames(df.sampled) != "label"]
    Rsqrt <- apply(df.sampled, 2, function(x) {
      dat.df <- data.frame(val = x, label = labels.sampled)
      summary(lm(val~label, data = dat.df))$r.squared
    })
    return(Rsqrt)
  }))
  
  return(res)
}


## Filter regulons by score compared to TF expression level
filterRegulons <- function(aucell, gexp, celltypes, regulon_activity_threshold = 0.5, gene_expressed_threshold = 0.05, aucell.mean = F) {
  gexp2 <- gexp[, gsub("[(].*", "", colnames(aucell))]
  gexp2 <- (gexp2 > 0) * 1
  gexp2 <- regulonMeanScore(gexp2, celltypes)
  gexp.expressed <- gexp2 > gene_expressed_threshold
  if (aucell.mean) {
    aucell2 <- aucell
  } else {
    aucell2 <- regulonMeanScore(t(aucell), celltypes)
  }
  aucell2 <- apply(aucell2, 2, function(x) {
    x <- x - min(x)
    x <- x / max(x)
    x
  })
  aucell.active <- aucell2 > regulon_activity_threshold
  fail <- colSums(aucell.active & (! gexp.expressed)) > 0
  return(aucell[, ! fail])
}





####################################################
####################################################

# Prepare data for analysis

## TF annotation from Garcia-Alonso et al. 2019 Genome Research https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/
TF_classification <- readxl::read_xlsx("/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/GarciaAlonso_supplemental_table_S1.xlsx", skip = 1, na = "-")
TF_classification <- as.data.frame(TF_classification)

## TF annotation from Heinäniemi et al. 2013 Nature Methods https://www.nature.com/articles/nmeth.2445
TF_classification2 <- readxl::read_xls("/research/groups/allseq/data/scRNAseq/results/juha_wrk/resources/shmulevichTableS5.xls", skip = 2, na = "na")
TF_classification2 <- as.data.frame(TF_classification2)
mapping <- read.delim("/research/groups/biowhat_share/anno_data/uniq_homosapiens_entrez2refseq2symbol.txt", header = F)
mapping <- unique(mapping[, -2])
colnames(mapping) <- c("EntrezID", "Gene")
TF_classification2$TF <- mapping$Gene[match(TF_classification2$EntrezID, mapping$EntrezID)]
TF_classification2$Category[TF_classification2$Category == "NA"] <- NA

## SCENIC results from HCA only and HCA + ALL runs
HCA_ALL_jitter <- data.frame(data.table::fread(file = "/research/groups/allseq/data/scRNAseq/results/juha_wrk/SCENIC/HCA_ALL_jitter/aucell_final_regulons_means.csv.gz"), row.names = 1, check.names = F)
HCA_ALL_jitter.orig <- HCA_ALL_jitter

## Metadata
load("/research/groups/allseq/data/scRNAseq/results/ER_project/new_2019/data/HCA_ALL_combined.RData")
metadata$celltype[metadata$phase != "G1" & metadata$celltype == "leukemic"] <- "leukemic_cycling"
metadata$celltype2 <- metadata$celltype
metadata$celltype2[metadata$celltype == "leukemic"] <- paste0("leukemic_", metadata$batch[metadata$celltype == "leukemic"])
metadata$celltype3 <- metadata$celltype
metadata$celltype3[metadata$louvain == "HCA.4"] <- "6_preB-II_G1"
metadata$celltype4 <- metadata$celltype2
metadata$celltype4[metadata$louvain == "HCA.4"] <- "6_preB-II_G1"
metadata <- metadata[rownames(HCA_ALL_jitter), ]
metadata.orig <- metadata

## Remove I samples from data
HCA_ALL_jitter <- HCA_ALL_jitter[! metadata$batch %in% c("ALL10-d15", "ALL12-d15"), ]
X <- X[! metadata$batch %in% c("ALL10-d15", "ALL12-d15"), ]
metadata <- metadata[! metadata$batch %in% c("ALL10-d15", "ALL12-d15"), ]

## Calculate mean scores for cell types
HCA_ALL_jitter.means <- regulonMeanScore(t(HCA_ALL_jitter), metadata$celltype)
HCA_ALL_jitter.means2 <- regulonMeanScore(t(HCA_ALL_jitter), metadata$celltype4)

## Remove immature B cells not used in scDD analysis (for linear fit)
HCA_ALL_jitter <- HCA_ALL_jitter[! metadata$louvain %in% c("HCA.1", "HCA.2", "HCA.8"), ]
X2 <- X[! metadata$louvain %in% c("HCA.1", "HCA.2", "HCA.8"), ]
metadata2 <- metadata[! metadata$louvain %in% c("HCA.1", "HCA.2", "HCA.8"), ]

####################################################
####################################################

# Fit a linear model to HCA + ALL result and plot

HCA_ALL.lm <- linearFitPermutation(HCA_ALL_jitter, metadata2$celltype, iterations = 100)

HCA_ALL.TF_annotation <- TF_classification[match(gsub("[(].*", "", colnames(HCA_ALL_jitter)), TF_classification$TF), ]
HCA_ALL.TF_annotation$regulon <- colnames(HCA_ALL_jitter)
HCA_ALL.TF_annotation$regulon_type <- gsub(".*[(]|[)]", "", HCA_ALL.TF_annotation$regulon)
HCA_ALL.TF_annotation$regulon_type <- plyr::mapvalues(HCA_ALL.TF_annotation$regulon_type, c("+", "-"), c("activators", "repressors"))
HCA_ALL.TF_annotation$correct <- "False"
HCA_ALL.TF_annotation$correct[(HCA_ALL.TF_annotation$regulon_type == "activators" & HCA_ALL.TF_annotation$mode_of_regulation %in% c("activators", "dual"))
                              | (HCA_ALL.TF_annotation$regulon_type == "repressors" & HCA_ALL.TF_annotation$mode_of_regulation %in% c("repressors", "dual"))] <- "True"
HCA_ALL.TF_annotation$correct[is.na(HCA_ALL.TF_annotation$mode_of_regulation)] <- NA

HCA_ALL_jitter.filt <- filterRegulons(HCA_ALL_jitter.means, X, metadata$celltype, 0.7, 0.04, aucell.mean = T)
HCA_ALL.TF_annotation$activity_filter <- plyr::mapvalues(colnames(HCA_ALL_jitter.means) %in% colnames(HCA_ALL_jitter.filt), c(T, F), c("pass", "fail"))

HCA_ALL.TF_annotation$TF_category <- TF_classification2$Category[match(gsub("[(].*", "", colnames(HCA_ALL_jitter)), TF_classification2$TF)]
HCA_ALL.TF_annotation$TF_category[is.na(HCA_ALL.TF_annotation$TF_category)] <- "NA"

HCA_ALL.TF_annotation$Rsquared_mean <- rowMeans(HCA_ALL.lm)

HCA_ALL_jitter.means.scaled <- t(scale(HCA_ALL_jitter.means2))

filter <- HCA_ALL.TF_annotation$Rsquared_mean > 0.5 &
  ! is.na(rowSums(HCA_ALL.lm)) &
  HCA_ALL.TF_annotation$activity_filter == "pass" &
  ! apply(HCA_ALL_jitter.means, 2, which.max) %in% c(1:3, 5)

mat <- HCA_ALL_jitter.means.scaled[filter, ]
ann <- HCA_ALL.TF_annotation[filter, ]
lmr <- HCA_ALL.lm[filter, ]

mat[mat > 3] <- 3
mat[mat < -3] <- -3

col_fun = circlize::colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

TF_category_colors <- RColorBrewer::brewer.pal(n = length(unique(ann$TF_category)), name = "Set3")
names(TF_category_colors) <- unique(ann$TF_category)

hmr <- Heatmap(mat,
               col = col_fun,
               #row_split = ann$activity_filter,
               cluster_rows = T,
               clustering_distance_rows = "pearson",
               show_row_dend = T,
               cluster_columns = F,
               show_row_names = T,
               show_column_names = T,
               row_names_gp = gpar(fontsize = 6),
               right_annotation = rowAnnotation(#"regulon_type_SCENIC" = ann$regulon_type,
                                                #"TF_type_annotated" = ann$mode_of_regulation,
                                                #"correct_annotation" = ann$correct,
                                                "TF_category" = ann$TF_category,
                                                col = list("TF_type_annotated" = c("activators" = "#1B9E77", "repressors" = "#D95F02", "dual" = "#7570B3"),
                                                           "regulon_type_SCENIC" = c("activators" = "#1B9E77", "repressors" = "#D95F02"),
                                                           "correct_annotation" = c("True" = "#66A61E", "False" = "#E7298A"),
                                                           "TF_category" = TF_category_colors)),
               show_heatmap_legend = T,
               heatmap_legend_param = list(title = "z-score"),
               use_raster = F
)

name_prefix <- "HCA_ALL_jitter_sampledLinearFit_Rsqrt_gt_05_activityFilter"

pdf(paste0(name_prefix, ".pdf"), width = 6, height = 12)
draw(hmr)
graphics.off()

write.table(HCA_ALL.TF_annotation, paste0(name_prefix, "2.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)






# Violin plots of regulons' scores
pdf(paste0(name_prefix, "_violinPlots.pdf"), width = 7, height = 4)
for(i in rownames(mat)) {
  df.tmp <- data.frame(score = HCA_ALL_jitter.orig[, i], celltype = factor(metadata.orig$celltype4))
  p <- ggplot(df.tmp, aes(x = celltype, y = score)) +
    geom_violin() +
    labs(title = i, x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(p)
}
graphics.off()








