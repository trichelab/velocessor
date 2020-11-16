#' blockwise downsampling: try to preserve balanced clusters across samples
#' 
#' Originally this was a helper function for label_cells but it is useful 
#' in its own right, so it is exposed as an exported function now. Note that 
#' `maxcells` is a MAXIMUM, i.e. if there are 30 cells in a cluster for a sample
#' and maxcells == 50, then (obviously?) only 30 cells will be returned for that
#' particular combination of cluster and sample.
#' 
#' @param txis      SingleCellExperiment where !is.null(colLabels(txis))
#' @param maxcells  max cells per cluster per sample (see Details) (50)
#' @param mincells  min cells per cluster per sample (see Details) (10)
#' @param ...       additional arguments to accomodate bootstrapping (not yet)
#'
#' @return          colnames(txis) satisfying the sampling scheme (see Details)
#' 
#' @details
#' Especially when using the default Louvain clustering approach, there will 
#' be samples without any cells in a cluster, and vice versa. To avoid having
#' a bunch of artifacts, when sample==TRUE, we fit a mixture model to the number
#' of cells in each cluster, and exclude samples with few or no cells in that 
#' cluster from block sampling. Don't use this on SmartSeq-type data.
#' 
#' Note that attr(downsample_txis(txis), "scheme") is a list with elements 
#' 'mincells', 'maxcells', and 'eligible'. 'mincells' and 'maxcells' are 
#' integers, while 'eligible' is an integer matrix with counts of cells 
#' post-filtering (i.e., subject to `mincells` and per-cluster mixture fits). 
#'
#' The mixture fits assume that a two-component mixture model on either 
#' log(1+cells) or directly on cell number per cluster will remove "noise" 
#' elements. This may be false; the user will have to investigate if so. 
#' 
#' @import mclust
#' 
#' @export
downsample_txis <- function(txis, maxcells=50, mincells=10, ...) {

  if (length(colLabels(txis)) < ncol(txis)) stop("colLabels() is empty!")
  if (is.null(txis$sample)) stop("txis$sample is NULL!")
  
  # tabulate sample-by-cluster combinations
  combos <- .lex_to_num(as.matrix(table(colLabels(txis), txis$sample)))
  if (nrow(combos) > 50) message("You have a LOT of clusters; this may fail.")

  # if all or nearly all of a particular sample's cells are in one cluster,
  # or a cluster has less than mincells total, drop the sample, cluster, or 
  # sample*cluster combination (see .plausible for details) 
  plausible <- .plausible(combos, mincells=mincells)
  samplethis <- matrix(0, nrow=nrow(combos), ncol=ncol(combos))
  dimnames(samplethis) <- dimnames(combos)

  # second pass: zero out combinations in the "soup"
  if (ncol(combos) < 10) message("You have few samples; mixture fits may fail.")
  classified <- .fit_mixtures(plausible)
  samplethis[rownames(classified), colnames(classified)] <- classified 
  eligible <- combos * samplethis # simple dot product 

  # pool names:
  sat <- list()
  tosample <- apply(classified, 2, function(x) names(which(x==1)))
  for (samp in names(tosample)) {
    for (clust in tosample[[samp]]) {
      sat[[paste(samp, clust, sep="_")]] <- .whichcells(txis, samp, clust)
    }
  }

  # return a named vector of sampled colnames
  res <- .nv(do.call(c, lapply(sat, .getcells)))
  scheme <- list(mincells=mincells, 
                 maxcells=maxcells, 
                 eligible=eligible)
  attr(res, "scheme") <- scheme
  return(res) 

}


# for block sampling
.getcells <- function(cells, maxcells=50) {

  maxcells <- min(length(cells), maxcells)
  sample(cells, maxcells)

}


# for block sampling 
.whichcells <- function(txis, samp, clust) {

  picked <- which(txis$sample == samp & colLabels(txis) == clust)
  colnames(txis)[picked]

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


# helper fn -- fit mixture models on log1p(cells) with fallback onto I(cells)
.fit_mixtures <- function(cells) { 

  logcells <- log1p(cells)
  cellfits <- apply(cells, 1, densityMclust, G=1:2, verbose=FALSE)
  logcellfits <- apply(logcells, 1, densityMclust, G=1:2, verbose=FALSE)

  modelnames <- .nv(c("cells", "logcells"))
  comps <- data.frame(cells=unname(sapply(cellfits, `[[`, "G")),
                      logcells=unname(sapply(logcellfits, `[[`, "G")))
  bics <- data.frame(cells=unname(sapply(cellfits, `[[`, "bic")),
                     logcellcomps=unname(sapply(logcellfits, `[[`, "bic")))
  rownames(comps) <- rownames(bics) <- names(cellfits)

  fits <- list() 
  for (cl in rownames(cells)) {
    fits[[cl]] <- list(cells=cellfits[[cl]], logcells=logcellfits[[cl]])
  }

  # default just returns 1 
  model <- rep("none", nrow(cells))
  names(model) <- rownames(cells) 
  res <- matrix(0, nrow=nrow(cells), ncol=ncol(cells)) 
  dimnames(res) <- dimnames(cells) 
  for (cl in rownames(comps)) {
    if (all(comps[cl, ] > 1)) {
      model[cl] <- modelnames[which.max(bics[cl,])] 
      res[cl, ] <- (fits[[cl]][[model[cl]]][["classification"]] - 1)
    } else if (any(comps[cl, ] > 1)) {
      model[cl] <- modelnames[which.max(comps[cl,])] 
      res[cl, ] <- (fits[[cl]][[model[cl]]][["classification"]] - 1)
    } else { 
      model[cl] <- "none"
      res[cl, ] <- rep(1, ncol(res))
    }
  }
  return(res)

}


# mclust helper
.plotDens <- function(fits) { 

  par(mfrow=c(1,2))
  plotDensityMclust1(fits[[1]], data=fits[[1]]$data, xlab="cells")
  plotDensityMclust1(fits[[2]], data=fits[[2]]$data, xlab="log1p(cells)")

} 


# order clusters numerically rather than lexically
.lex_to_num <- function(mat) { 

  prefix <- .lcs(rownames(mat))
  numbers <- as.integer(sub(prefix, "", rownames(mat)))
  mat[order(numbers), ]  

}


# longest common prefix (since Biostrings::lcprefix is experimental) 
.lcs <- function(..., sep="") {

  paste(Reduce(intersect, strsplit(unlist(...), sep)), collapse = sep)

}


# helper fn
.nv <- function(x) { 

  names(x) <- x
  x

}
