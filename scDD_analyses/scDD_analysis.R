# scDD analysis with vst
#
# - HCA data: HSC + B lineage cells
# - ALL data: leukemic cells from ALL1, ALL3, ALL8, ALL9, ALL10, ALL10d15, ALL12 and ALL12d15
# - Nalm6 cell line: LUC and E/R induced
# - Use vst normalized data
# - Custom scDD installation devtools::install_github("juhaa/scDD")
#   - adds jitter to counts of genes which cause numerical problems with mclust
#



# Load libraries, setup options
source("scDD_wrapper_scripts.R")
future::plan(strategy = 'multicore', workers = 2)
options(future.globals.maxSize = 10 * 1024 ^ 3, stringsAsFactors = F)
setwd("/research/groups/allseq/data/scRNAseq/results/ER_project/new_2019/scDD/vst")



# Set up HCA + ALL data
load("/research/groups/allseq/data/scRNAseq/results/ER_project/new_2019/data/HCA_ALL_combined.RData")
X <- t(X)
metadata$log_umi <- log10(colSums(X))
metadata$n_counts <- colSums(X)
metadata$celltype_leukemic <- metadata$celltype
metadata$celltype_leukemic[! metadata$batch %in% paste0("MantonBM", 1:8)] <- metadata$batch[! metadata$batch %in% paste0("MantonBM", 1:8)]
metadata$celltype_leukemic_noI <- metadata$celltype
metadata$celltype_leukemic_noI[metadata$batch %in% c("ALL10d15", "ALL12d15")] <- "I"
metadata$celltype_leukemic_noI_noCC <- metadata$celltype_leukemic_noI
metadata$celltype_leukemic_noI_noCC[metadata$phase != "G1"] <- "CC"
metadata$celltype_leukemic_noI_CC <- metadata$celltype_leukemic_noI
metadata$celltype_leukemic_noI_CC[metadata$phase == "G1"] <- "G1"
metadata$celltype2 <- metadata$celltype
metadata$celltype2[metadata$celltype == "leukemic"] <- paste0("leukemic_", metadata$batch[metadata$celltype == "leukemic"])



# Set up Nalm6 data
nalm <- data.table::fread("/research/groups/allseq/data/scRNAseq/results/juha_wrk/Nalm6_scanpy/data/raw/X.csv.gz", data.table = F)
colnames(nalm) <- readLines("/research/groups/allseq/data/scRNAseq/results/juha_wrk/Nalm6_scanpy/data/raw/var.txt")
nalm.anno <- read.csv("/research/groups/allseq/data/scRNAseq/results/juha_wrk/Nalm6_scanpy/data/obs.csv", stringsAsFactors = F, row.names = 1)
rownames(nalm) <- rownames(nalm.anno) <- gsub("-1", "", rownames(nalm.anno))
nalm <- t(nalm)
nalm.anno$log_umi <- log10(colSums(nalm))



# vst
set.seed(42)
HCA_ALL <- sctransform::vst(umi = X, cell_attr = metadata, batch_var = "batch", return_corrected_umi = T)
saveRDS(HCA_ALL, "vst_HCABlineage_ALL_log_counts_batch.rds")
write.csv(t(as.matrix(HCA_ALL$umi_corrected)), gzfile("vst_HCABlineage_ALL_log_counts_batch_UMIcorrected.csv.gz"))
#HCA_ALL <- readRDS("vst_HCABlineage_ALL_log_counts_batch.rds")

set.seed(42)
nalm.vst <- sctransform::vst(umi = nalm, cell_attr = nalm.anno, batch_var = "batch", return_corrected_umi = T)
saveRDS(nalm.vst, "vst_batch_Nalm6_log_counts.rds")
write.csv(t(as.matrix(nalm.vst$umi_corrected)), gzfile("vst_batch_Nalm6_log_counts_UMIcorrected.csv.gz"))
#nalm.vst <- readRDS("vst_batch_Nalm6_log_counts.rds")


# scDD runs
metadata$n_counts_vst <- Matrix::colSums(HCA_ALL$umi_corrected)

scDD_run(HCA_ALL$umi_corrected, metadata$celltype, grep("1_|2_", metadata$celltype, value = T), n_cores = 5, prefix = "HCA_Blineage_ALL_sct_batch_2_earlyLymphoid_vs_1_HSC", ref = "2_earlyLymphoid", reuse = F, min.size = 30)
scDD_run(HCA_ALL$umi_corrected, metadata$celltype, grep("2_|3_", metadata$celltype, value = T), n_cores = 5, prefix = "HCA_Blineage_ALL_sct_batch_3_proB_G2MS_vs_2_earlyLymphoid", ref = "3_proB_G2MS", reuse = F, min.size = 30)
scDD_run(HCA_ALL$umi_corrected, metadata$celltype, grep("3_|5_", metadata$celltype, value = T), n_cores = 20, prefix = "HCA_Blineage_ALL_sct_batch_5_preB_G2MS_vs_3_proB_G2MS", ref = "5_preB_G2MS", reuse = F, min.size = 30)
scDD_run(HCA_ALL$umi_corrected[, metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500], metadata$celltype[metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500], grep("4_|6_", metadata$celltype, value = T), n_cores = 10, prefix = "HCA_Blineage_ALL_sct_batch_6_preB_G1_vs_4_proB_G1_minncounts3000_maxncounts3500", ref = "6_preB_G1", reuse = F, min.size = 30)
scDD_run(HCA_ALL$umi_corrected, metadata$louvain, c("HCA.0", "HCA.4"), n_cores = 10, prefix = "HCA_Blineage_ALL_sct_batch_cluster0_vs_cluster4", ref = "HCA.4", reuse = F, min.size = 30)
scDD_run(HCA_ALL$umi_corrected[, metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500], metadata$celltype_leukemic_noI_noCC[metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500], c("leukemic", "4_proB_G1"), n_cores = 20, prefix = "HCA_Blineage_ALL_sct_batch_4_proB_G1_vs_leukemic_G1_noI_minncounts3000_maxncounts3500", ref = "leukemic", reuse = F, min.size = 30)
scDD_run(HCA_ALL$umi_corrected, metadata$celltype_leukemic_noI_CC, grep("leukemic|3_proB_G2MS", metadata$celltype_leukemic_noI_CC, value = T), n_cores = 20, prefix = "HCA_Blineage_ALL_sct_batch_3_proB_G2MS_vs_leukemic_noG1_noI", ref = "leukemic", reuse = F, min.size = 30)

scDD_run(HCA_ALL$umi_corrected[, metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500 & metadata$batch != "ALL12"],
         metadata$celltype_leukemic_noI_noCC[metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500 & metadata$batch != "ALL12"],
         c("leukemic", "4_proB_G1"),
         n_cores = 10,
         prefix = "HCA_Blineage_ALL_sct_batch_4_proB_G1_vs_leukemic_G1_noI_no_ALL12_minncounts3000_maxncounts3500",
         ref = "leukemic",
         reuse = F,
         min.size = 30)

scDD_run(HCA_ALL$umi_corrected[, metadata$batch != "ALL12"],
         metadata$celltype_leukemic_noI_CC[metadata$batch != "ALL12"],
         c("leukemic", "3_proB_G2MS"),
         n_cores = 10,
         prefix = "HCA_Blineage_ALL_sct_batch_3_proB_G2MS_vs_leukemic_noG1_noI_noALL12",
         ref = "leukemic",
         reuse = F,
         min.size = 30)

scDD_run(HCA_ALL$umi_corrected[, metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500 & metadata$phase == "G1"],
         metadata$celltype2[metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500 & metadata$phase == "G1"],
         c("leukemic_ALL12", "4_proB_G1"),
         n_cores = 10,
         prefix = "HCA_Blineage_ALL_sct_batch_4_proB_G1_vs_leukemic_ALL12_G1_minncounts3000_maxncounts3500",
         ref = "leukemic_ALL12",
         reuse = F,
         min.size = 30)

scDD_run(HCA_ALL$umi_corrected[, metadata$phase != "G1"],
         metadata$celltype2[metadata$phase != "G1"],
         c("leukemic_ALL12", "3_proB_G2MS"),
         n_cores = 10,
         prefix = "HCA_Blineage_ALL_sct_batch_3_proB_G2MS_vs_leukemic_ALL12_noG1",
         ref = "leukemic_ALL12",
         reuse = F,
         min.size = 30)




scDD_run(nalm.vst$umi_corrected, nalm.anno$batch, nalm.anno$batch, n_cores = 20, prefix = "Nalm6_vst_batch_LUC_vs_ER", ref = "ER", reuse = F, min.size = 30)



# Density plots of selected genes
data <- as.matrix(HCA_ALL$umi_corrected)
proB <- data[, (metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500 & metadata$celltype_leukemic_noI_noCC == "4_proB_G1") | metadata$celltype == "3_proB_G2MS"]
leukemic <- data[, metadata$celltype == "leukemic"]

df <- data.frame(log(t(cbind(proB, leukemic)) + 1),
                 celltype = factor(c(rep("proB", ncol(proB)), rep("leukemic", ncol(leukemic))), levels = c("leukemic", "proB")))
colnames(df) <- gsub("[.]", "-", colnames(df))

mdata <- metadata[rownames(df), ]

genes <- sort(c("TERF2", "INSR", "TGFB1", "HLA-DOB", "HLA-E", "LY6E", "ISG20", "IGLL1", "LAT2", "CYFIP2", "VPREB1"))

library(ggplot2)

pdf("proB_all_vs_leukemic_all_densityPlots.pdf", width = 5, height = 3)

for(i in genes) {
  df.tmp <- data.frame(gexp = df[, i], celltype = df$celltype)
  p <- ggplot(df.tmp, aes(x = gexp, fill = celltype, colour = celltype)) +
    geom_density(alpha = 0.4, adjust = 2) +
    scale_fill_brewer(palette = "Dark2") +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = i) +
    theme_classic()
  plot(p)
  
}

graphics.off()


df2 <- df[mdata$phase == "G1", ]

pdf("proB_G1_vs_leukemic_G1_densityPlots.pdf", width = 5, height = 3)

for(i in genes) {
  df.tmp <- data.frame(gexp = df2[, i], celltype = df2$celltype)
  p <- ggplot(df.tmp, aes(x = gexp, fill = celltype, colour = celltype)) +
    geom_density(alpha = 0.4, adjust = 2) +
    scale_fill_brewer(palette = "Dark2") +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = i) +
    theme_classic()
  plot(p)
  
}

graphics.off()




# Density plots of selected genes per sample
data <- as.matrix(HCA_ALL$umi_corrected)
proB <- data[, (metadata$n_counts_vst > 3000 & metadata$n_counts_vst < 3500 & metadata$celltype_leukemic_noI_noCC == "4_proB_G1") | metadata$celltype == "3_proB_G2MS"]

genes <- sort(c("TERF2", "INSR", "TGFB1", "HLA-DOB", "HLA-E", "LY6E", "ISG20", "IGLL1", "LAT2", "CYFIP2", "VPREB1"))

for(j in unique(metadata$batch[metadata$celltype == "leukemic"])) {
  
  leukemic <- data[, metadata$celltype == "leukemic" & metadata$batch == j]
  
  df <- data.frame(log(t(cbind(proB, leukemic)) + 1),
                   celltype = factor(c(rep("proB", ncol(proB)), rep("leukemic", ncol(leukemic))), levels = c("leukemic", "proB")))
  colnames(df) <- gsub("[.]", "-", colnames(df))
  
  mdata <- metadata[rownames(df), ]
  
  pdf(paste0("proB_all_vs_leukemic_", j, "_all_densityPlots.pdf"), width = 5, height = 3)
  
  for(i in genes) {
    df.tmp <- data.frame(gexp = df[, i], celltype = df$celltype)
    p <- ggplot(df.tmp, aes(x = gexp, fill = celltype, colour = celltype)) +
      geom_density(alpha = 0.4, adjust = 2) +
      scale_fill_brewer(palette = "Dark2") +
      scale_colour_brewer(palette = "Dark2") +
      labs(title = i) +
      theme_classic()
    plot(p)
    
  }
  
  graphics.off()
  
}





for(j in unique(metadata$batch[metadata$celltype == "leukemic"])) {
  
  leukemic <- data[, metadata$celltype == "leukemic" & metadata$batch == j]
  
  df <- data.frame(log(t(cbind(proB, leukemic)) + 1),
                   celltype = factor(c(rep("proB", ncol(proB)), rep("leukemic", ncol(leukemic))), levels = c("leukemic", "proB")))
  colnames(df) <- gsub("[.]", "-", colnames(df))
  
  mdata <- metadata[rownames(df), ]
  
  df <- df[mdata$phase == "G1", ]
  
  pdf(paste0("proB_G1_vs_leukemic_", j, "_G1_densityPlots.pdf"), width = 5, height = 3)
  
  for(i in genes) {
    df.tmp <- data.frame(gexp = df[, i], celltype = df$celltype)
    p <- ggplot(df.tmp, aes(x = gexp, fill = celltype, colour = celltype)) +
      geom_density(alpha = 0.4, adjust = 2) +
      scale_fill_brewer(palette = "Dark2") +
      scale_colour_brewer(palette = "Dark2") +
      labs(title = i) +
      theme_classic()
    plot(p)
    
  }
  
  graphics.off()
  
}






for(j in unique(metadata$batch[metadata$celltype == "leukemic"])) {
  
  leukemic <- data[, metadata$celltype == "leukemic" & metadata$batch == j]
  
  df <- data.frame(log(t(cbind(proB, leukemic)) + 1),
                   celltype = factor(c(rep("proB", ncol(proB)), rep("leukemic", ncol(leukemic))), levels = c("leukemic", "proB")))
  colnames(df) <- gsub("[.]", "-", colnames(df))
  
  mdata <- metadata[rownames(df), ]
  
  df <- df[mdata$phase != "G1", ]
  
  pdf(paste0("proB_G2MS_vs_leukemic_", j, "_G2MS_densityPlots.pdf"), width = 5, height = 3)
  
  for(i in genes) {
    df.tmp <- data.frame(gexp = df[, i], celltype = df$celltype)
    p <- ggplot(df.tmp, aes(x = gexp, fill = celltype, colour = celltype)) +
      geom_density(alpha = 0.4, adjust = 2) +
      scale_fill_brewer(palette = "Dark2") +
      scale_colour_brewer(palette = "Dark2") +
      labs(title = i) +
      theme_classic()
    plot(p)
    
  }
  
  graphics.off()
  
}






# Density plots for selected genes from diag vs d15 results
genes_both <- sort(c("RUNX1", "TCF3", "LEF1", "SOX4", "SMAD1", "JUNB"))
genes_ALL12 <- sort(c("ERG", "IKZF1", "ARID5A", "JUN", "ETS2", "POU4F1", "TCF4", "ELK3"))
genes_ALL10 <- sort(c("IRF8", "CEBPD", "POU2F2", "KLF2", "KLF6", "AFF3", "SPIB", "MS4A1"))
genes_all <- c(genes_both, genes_ALL10, genes_ALL12)

df <- data.frame(log(t(data[, metadata$batch %in% c("ALL10", "ALL12", "ALL10d15", "ALL12d15")] + 1)),
                 sample = metadata$batch[metadata$batch %in% c("ALL10", "ALL12", "ALL10d15", "ALL12d15")])


pdf("densityPlots_ALL10_ALL10d15.pdf", width = 5, height = 3)

for(i in genes_all) {
  df.tmp <- data.frame(gexp = df[df$sample %in% c("ALL10", "ALL10d15"), i], sample = df$sample[df$sample %in% c("ALL10", "ALL10d15")])
  p <- ggplot(df.tmp, aes(x = gexp, fill = sample, colour = sample)) +
    geom_density(alpha = 0.4, adjust = 2) +
    scale_fill_brewer(palette = "Dark2") +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = i) +
    theme_classic()
  plot(p)
  
}

graphics.off()


pdf("densityPlots_ALL12_ALL12d15.pdf", width = 5, height = 3)

for(i in genes_all) {
  df.tmp <- data.frame(gexp = df[df$sample %in% c("ALL12", "ALL12d15"), i], sample = df$sample[df$sample %in% c("ALL12", "ALL12d15")])
  p <- ggplot(df.tmp, aes(x = gexp, fill = sample, colour = sample)) +
    geom_density(alpha = 0.4, adjust = 2) +
    scale_fill_brewer(palette = "Dark2") +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = i) +
    theme_classic()
  plot(p)
  
}

graphics.off()



library(dplyr)
df2 <- df[, c(genes_all, "sample")] %>% group_by(sample) %>% summarise_all(mean)
mat <- as.data.frame(df2)
rownames(mat) <- mat[, 1]
mat <- mat[, -1]
mat <- t(scale(mat))
mat <- mat[, c(1, 3, 2, 4)]

hmr <- Heatmap(mat,
               cluster_rows = F,
               cluster_columns = F,
               show_row_names = T,
               show_column_names = T,
               row_names_gp = gpar(fontsize = 4),
               show_heatmap_legend = T,
               use_raster = F
)

pdf("heatmap_diag_d15.pdf", width = 5, height = 12)
draw(hmr)
graphics.off()



# Ridge plots from selected genes (Fig5D)
# Note: Add CD20 represented as protein level (FACS)
TFs <- c("RUNX1", "POU2F2", "ERG",
         "TCF3", "KLF6", "ELK3",
         "SOX4", "AFF3", "MS4A1",
         "SMAD1", "SPIB")

labels <- c(TFs, "CD20")

CD20 <- read.delim("/research/groups/allseq/data/scRNAseq/results/ER_project/data/CD20_FACS.txt")
load("/research/groups/allseq/data/scRNAseq/results/ER_project/new_2019/data/HCA_ALL_combined.RData")
HCA_ALL <- readRDS("vst_HCABlineage_ALL_log_counts_batch.rds")
data <- as.matrix(HCA_ALL$umi_corrected)[TFs]
df <- data.frame(log(t(data[, metadata$batch %in% c("ALL10", "ALL12", "ALL10d15", "ALL12d15") & metadata$phase == "G1"] + 1)),
                 batch = metadata$batch[metadata$batch %in% c("ALL10", "ALL12", "ALL10d15", "ALL12d15") & metadata$phase == "G1"])

df$day <- "d0"
df$day[df$batch %in% c("ALL10d15", "ALL12d15")] <- "d15"
df$donor <- "ALL10"
df$donor[df$batch %in% c("ALL12", "ALL12d15")] <- "ALL12"
df$day <- factor(df$day)
df$donor <- factor(df$donor)

library(ggridges)
library(ggpubr)

pl <- lapply(TFs, function(x) {
  ggplot(df, aes_string(x = x, y = "donor", fill = "day")) + 
    geom_density_ridges(alpha = .5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_ridges()
})

CD20.df <- rbind.data.frame(data.frame(CD20 = CD20$ALL10_d0_CD20[! is.na(CD20$ALL10_d0_CD20)], donor = "ALL10", day = "d0"),
                            data.frame(CD20 = CD20$ALL10_d15_CD20[! is.na(CD20$ALL10_d15_CD20)], donor = "ALL10", day = "d15"),
                            data.frame(CD20 = CD20$ALL12_d0_CD20[! is.na(CD20$ALL12_d0_CD20)], donor = "ALL12", day = "d0"),
                            data.frame(CD20 = CD20$ALL12_d15_CD20[! is.na(CD20$ALL12_d15_CD20)], donor = "ALL12", day = "d15"))

pl[[12]] <- ggplot(CD20.df, aes(x = CD20, y = donor, fill = day)) + 
  geom_density_ridges(alpha = .5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_ridges()


pb <- ggarrange(pl[[1]], pl[[2]], pl[[3]],
                pl[[4]], pl[[5]], pl[[6]],
                pl[[7]], pl[[8]], pl[[9]],
                pl[[10]], pl[[11]], pl[[12]], nrow = 4, ncol = 3, common.legend = T)


pdf("Fig5_TFs_ridgePlot.pdf")
plot(pb)
graphics.off()





# Fix bandwidth of gexp to 0.2

pl <- lapply(TFs, function(x) {
  ggplot(df, aes_string(x = x, y = "donor", fill = "day", color = "day")) + 
    stat_density_ridges(alpha = .5, bandwidth = .2) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_ridges()
})


pl[[12]] <- ggplot(CD20.df, aes(x = CD20, y = donor, fill = day, color = day)) + 
  geom_density_ridges(alpha = .5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_ridges()


pb <- ggarrange(pl[[1]], pl[[2]], pl[[3]],
                pl[[4]], pl[[5]], pl[[6]],
                pl[[7]], pl[[8]], pl[[9]],
                pl[[10]], pl[[11]], pl[[12]], nrow = 4, ncol = 3, common.legend = T)


pdf("Fig5_TFs_ridgePlot_bw_adjusted.pdf")
plot(pb)
graphics.off()


