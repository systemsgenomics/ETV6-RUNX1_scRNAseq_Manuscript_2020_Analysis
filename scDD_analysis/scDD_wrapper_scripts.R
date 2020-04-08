# Custom scripts for scDD analysis

options(stringsAsFactors = F)
library(scDD) # version 1.6.1 (devtools::install_github("juhaa/scDD"))
library(SingleCellExperiment) # version 1.4.1

#' findMaxClusterMeanPerCondition
#' 
#' Find highest average cluster mean or standard deviation
#' 
#' @param data Data matrix (genes in rows, cells in columns)
#' @param Zhat Clustering annotation as matrix (genes in rows, cells in columns)
#' @param what What to calculate. One of "mean" or "sd"
#' @param min.size Minimum size of cluster to take it into account when finding maximum cluster mean/sd
#' 
#' @return A vector of maximum mean/sd per gene
#' 
findMaxClusterMeanPerCondition <- function(data, Zhat, what, min.size = 3) {
  data <- data[, colnames(Zhat)]
  maxmean <- sapply(rownames(data), function(g) {
    D <- data[g, ]
    cl <- Zhat[g, ]
    cl.u <- table(cl)
    cl.u <- cl.u[cl.u > min.size & names(cl.u) != "0"]
    cl.u <- as.numeric(names(cl.u))
    if (length(cl.u) == 0) { return(0) }
    if (what == "mean") {
      return(max(sapply(cl.u, function(i) { mean(D[cl == i]) } )))
    }
    if (what == "sd") {
      return(max(sapply(cl.u, function(i) { sqrt(var(D[cl == i])) } )))
    }
  })
  return(maxmean)
}

#' getResults
#' 
#' Calculate some additional statistics for scDD results
#' 
#' @param data Data matrix (genes in rows, cells in columns)
#' @param anno A data.frame where column \code{condition} contains grouping info used in the scDD run
#' @param scdd \code{metadata} slot from the scDD results object
#' @param ref Reference condition
#' @param min.size Minimum size of cluster to take it into account when finding maximum cluster mean/sd
#' 
#' @return "Enchanced" scDD results
#'  
getResults <- function(data, anno, scdd, ref, min.size = 3) {
  r <- scdd$Genes
  r$c1.celltype <- as.character(anno$condition[which(colnames(data) %in% colnames(scdd$Zhat.c1)[1])])
  r$c2.celltype <- as.character(anno$condition[which(colnames(data) %in% colnames(scdd$Zhat.c2)[1])])
  r$c1.maxmean <- findMaxClusterMeanPerCondition(data = data, Zhat = scdd$Zhat.c1, what = "mean", min.size = min.size)
  r$c2.maxmean <- findMaxClusterMeanPerCondition(data = data, Zhat = scdd$Zhat.c2, what = "mean", min.size = min.size)
  r$c1.maxsd <- findMaxClusterMeanPerCondition(data = data, Zhat = scdd$Zhat.c1, what = "sd", min.size = min.size)
  r$c2.maxsd <- findMaxClusterMeanPerCondition(data = data, Zhat = scdd$Zhat.c2, what = "sd", min.size = min.size)
  r$c1.maxcv <- r$c1.maxsd / r$c1.maxmean
  r$c2.maxcv <- r$c2.maxsd / r$c2.maxmean
  r$c1.sd <- apply(data[, anno$condition == r$c1.celltype[1]], 1, var)
  r$c2.sd <- apply(data[, anno$condition == r$c2.celltype[1]], 1, var)
  if (ref == r$c1.celltype[1]) {
    r$sign <- r$c1.maxmean > r$c2.maxmean
    r$relative_difference <- r$c1.maxmean - r$c2.maxmean
  } else {
    r$sign <- r$c2.maxmean > r$c1.maxmean
    r$relative_difference <- r$c2.maxmean - r$c1.maxmean
  }
  r$sign <- plyr::mapvalues(x = r$sign, from = c(T, F), to = c("up", "down"))
  r$c1.n_zeros <- rowSums(scdd$Zhat.c1 == 0)
  r$c2.n_zeros <- rowSums(scdd$Zhat.c2 == 0)
  r$c1.zero_proportion <- r$c1.n_zeros / ncol(scdd$Zhat.c1)
  r$c2.zero_proportion <- r$c2.n_zeros / ncol(scdd$Zhat.c2)
  if (ref == r$c1.celltype[1]) {
    r$zero_proportion_difference <- r$c1.zero_proportion - r$c2.zero_proportion
  } else {
    r$zero_proportion_difference <- r$c2.zero_proportion - r$c1.zero_proportion
  }
  r$zero_has_lower_pvalue.adj <- r$zero.pvalue.adj < r$nonzero.pvalue.adj
  r$new_DDcategory <- r$DDcategory
  r$new_DDcategory[r$zero_has_lower_pvalue.adj & r$DDcategory != "NS"] <- "DZ"
  r$zero_proportion_difference_sign <- ifelse(r$zero_proportion_difference >= 0, "up", "down")
  r$sign_by_zero <- plyr::mapvalues(x = r$zero_proportion_difference_sign, from = c("up", "down"), to = c("down", "up"))
  
  return(r)
}

#' scDD_run
#' 
#' Wrapper for running scDD
#' 
#' @param data Data matrix (genes in rows, cells in columns) as CPM (or as raw counts if \code{scran_norm = T})
#' @param annotation Annotation vector for \{data}
#' @param groups Two group names from annotation to compare with scDD
#' @param n_cores Number of cores to use in parallel calculation
#' @param prefix Prefix filename for saving results
#' @param ref Reference group
#' @param min.size Minimum size of cluster to take it into account when finding maximum cluster mean/sd
#' @param min.nonzero.prop Minimum proportion of nonzero cells allowed per group for calculating nonzero statistics with scDD
#' @param scran_norm Boolean on whether to normalize using scran (done with \code{Preprocess}-function via scDD). Note: requires raw counts!
#' @param clean.genes Boolean for cleaning *possibly* problematic genes.
#' @param jitter Numeric vector of length two to define limits to jitter.
#' @param reuse Boolean on whether to use existing scDD results and calculate only the additional statistics
#' 
#' @return Saves scDD results as \code{paste0(prefix, "_scDD.rds")} and enhanced results as \code{paste0(prefix, "_scDD_results.tsv")}
#'   
scDD_run <- function(data, annotation, groups, n_cores, prefix, ref, min.size = 3, min.nonzero.prop = 0.05, scran_norm = F, clean.genes = F, jitter = c(-0.1, 0.1), reuse = F) {
  
  if (! ref %in% groups) stop("Reference '", ref, "' not found in groups.")
  if (! all(groups %in% annotation)) stop("All groups not found from annotation.")
  if (anyNA(data)) stop("Input data should not contain NA values.")
  if (! all(data >= 0)) stop("Input data should be nonnegative.")
  
  scdd1 <- paste0(prefix, "_scDD.rds")
  scdd2 <- paste0(prefix, "_scDD_results.tsv")
  
  data <- as.matrix(data[, annotation %in% groups])
  anno <- data.frame(condition = annotation[annotation %in% groups])
  
  # Check and remove all-zero cells & genes
  cs <- colSums(data)
  rs <- rowSums(data)
  if (any(cs == 0)) {
    data <- data[, cs > 0, drop = F]
    anno <- anno[cs > 0, , drop = F]
    message("Notice: Removed ", sum(cs == 0), " all-zero cells.")
  }
  if (any(rs == 0)) {
    data <- data[rs > 0, , drop = F]
    message("Notice: Removed ", sum(rs == 0), " all-zero genes.")
  }
  
  if (clean.genes) {
    # Remove genes where either condition is only zeros
    #rm <- (rowSums(data[, anno$condition == unique(groups)[1]]) == 0 & apply(data[, anno$condition == unique(groups)[2]], 1, function(x) length(unique(x))) < 3) |
    #  (rowSums(data[, anno$condition == unique(groups)[2]]) == 0 & apply(data[, anno$condition == unique(groups)[1]], 1, function(x) length(unique(x))) < 3)
    rm <- apply(data[, anno$condition == unique(groups)[1]], 1, function(x) length(unique(x))) < 3 & apply(data[, anno$condition == unique(groups)[2]], 1, function(x) length(unique(x))) < 3
    data <- data[! rm, , drop = F]
    message("Notice: Removed ", sum(rm), " genes that have low counts and may cause numerical problems during clustering.")
  }
  
  # If scran_norm = T, assume raw counts (not CPM) and use preprocess from scDD to normalize using scran.
  if (scran_norm) {
    sce <- SingleCellExperiment(assays = list(counts = data), colData = anno)
    sce <- preprocess(SCdat = sce, zero.thresh = 1, scran_norm = T)
  } else {
    sce <- SingleCellExperiment(assays = list(normcounts = data), colData = anno)
  }
  
  if (reuse & file.exists(scdd1)) {
    cat("Found existing scDD file. Loading it...\n")
    scdd <- readRDS(scdd1)
  } else {
    min.nonzero <- round(min(table(anno$condition)) * min.nonzero.prop)
    if (n_cores > 1) {
      scdd <- scDD(SCdat = sce, testZeroes = T, param = BiocParallel::MulticoreParam(workers = n_cores), min.size = min.size, min.nonzero = min.nonzero, jitter = jitter)
    } else {
      scdd <- scDD(SCdat = sce, testZeroes = T, param = BiocParallel::SerialParam(), min.size = min.size, min.nonzero = min.nonzero, jitter = jitter)
    }
    
    scdd <- scdd@metadata
    saveRDS(scdd, scdd1)
  }
  r <- getResults(data = data, anno = anno, scdd = scdd, ref = ref, min.size = min.size)
  write.table(r, scdd2, quote = F, col.names = T, row.names = F, sep = "\t")
  invisible(NULL)
}




