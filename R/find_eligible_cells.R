#' Find cells that meet certain cluster-wise criteria
#' 
#' By default, to be eligible, a cell must come from a sample whose presence 
#' in a cluster appears to be non-noise (i.e. not the "soup" part of a mixture),
#' and after pruning out samples that are all-one-cluster (likely artifacts) 
#' as well as clusters that are all-one-sample (likely artifacts), the cluster
#' from which the cells come must have at least `mincells` cells left in it.
#' The default for `mincells` is 10, so this isn't as harsh as it might seem, 
#' but it's not something you want to be relying upon for SmartSeq-like runs.
#' One useful application of find_eligible_cells is to feed plot_eligible_cells,
#' which can help determine whether to raise or lower `maxcells` and `mincells`
#' in downsample_txis and label_cells. 
#' 
#' @param txis      SingleCellExperiment where !is.null(colLabels(txis))
#' @param mincells  min cells per cluster per sample (see Details) (10)
#' 
#' @import mclust
#' 
#' @export
find_eligible_cells <- function(txis, mincells=10) { 

  combos <- .lex_to_num(as.matrix(table(colLabels(txis), txis$sample)))
  if (nrow(combos) > 50) {
    message("You have a LOT of clusters. This could be a problem downstream.")
  }

  # if all or nearly all of a particular sample's cells are in one cluster,
  # or a cluster has less than mincells total, drop the sample, cluster, or 
  # sample*cluster combination (see .plausible for details) 
  plausible <- .plausible(combos, mincells=mincells)
  samplethis <- matrix(0, nrow=nrow(combos), ncol=ncol(combos))
  dimnames(samplethis) <- dimnames(combos)

  # second pass: zero out combinations in the "soup"
  if (ncol(combos) < 10) message("You have few samples; mixture fits may fail.")
  classified <- fit_mixtures(plausible)
  samplethis[rownames(classified), colnames(classified)] <- classified 

  # all set (e.g. to plot) 
  attr(classified, "eligible") <- combos * samplethis 
  return(classified) 

}
  

# order clusters numerically rather than lexically
.lex_to_num <- function(mat) { 

  prefix <- .lcs(rownames(mat))
  numbers <- as.integer(sub(prefix, "", rownames(mat)))
  mat[order(numbers), ]  

}


# longest common substring (since Biostrings::lcprefix is experimental) 
.lcs <- function(..., sep="") {

  paste(Reduce(intersect, strsplit(unlist(...), sep)), collapse = sep)

}


# helper fn -- 
.keepcombos <- function(combos, mincells=10) { 

  csc <- colSums(combos) 
  pct <- round(sweep(combos, 2, csc, `/`), 2)
  keepcols <- names(which(colSums(pct > 0) > 1))

  kept <- combos[, keepcols] 
  keeprows <- names(which(rowSums(kept) >= mincells))

  combos[keeprows, keepcols]

} 


# helper fn 
.plausible <- function(combos, mincells=10) { 

  iters <- 1
  kept_0 <- .keepcombos(combos, mincells=mincells)
  kept_1 <- .keepcombos(kept_0, mincells=mincells)

  # iterate to a stable set 
  while (!identical(dim(kept_0), dim(kept_1))) { 
    iters <- iters + 1
    kept_0 <- kept_1
    kept_1 <- .keepcombos(kept_0, mincells=mincells)
  } 

  message("Converged after ", iters, " iteration", ifelse(iters > 1, "s.", "."))
  message("Dropped ", ncol(combos) - ncol(kept_1), " samples.")
  message("Dropped ", nrow(combos) - nrow(kept_1), " clusters.")
  return(kept_1)

} 
