#' blockwise downsampling: try to preserve balanced clusters across samples
#' 
#' Originally this was a helper function for label_cells but it is useful 
#' in its own right, so it is exposed as an exported function now. Note that 
#' `maxcells` is a MAXIMUM, i.e. if there are 30 cells in a cluster for a sample
#' and maxcells == 50, then (obviously?) only 30 cells will be returned for that
#' particular combination of cluster and sample.
#' 
#' @param txis      SingleCellExperiment where !is.null(colLabels(txis))
#' @param maxcells  max cells per cluster per sample (see Details) (20)
#' @param mincells  min cells per cluster per sample (see Details) (10)
#' @param ret       whether to return colnames ("colnames") or (default) "sce"
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
#' Note that attr(downsample_txis(txis, ret="colnames"), "scheme") is a list 
#' with elements 'mincells', 'maxcells', and 'eligible'. 'mincells' & 'maxcells'
#' are integers, while 'eligible' is an integer matrix with counts of cells 
#' post-filtering (i.e., subject to `mincells` and per-cluster mixture fits). 
#'
#' The mixture fits assume that a two-component mixture model on either 
#' log(1+cells) or directly on cell number per cluster will remove "noise" 
#' elements. This may be false; the user will have to investigate if so. 
#' 
#' @seealso find_eligible_cells
#' @seealso label_cells
#' 
#' @import mclust
#' 
#' @export
downsample_txis <- function(txis, maxcells=20, mincells=10, ret=c("sce","colnames"), ...) {

  if (length(colLabels(txis)) < ncol(txis)) stop("colLabels() is empty!")
  if (is.null(txis$sample)) stop("txis$sample is NULL!")
  ret <- match.arg(ret)

  classified <- find_eligible_cells(txis, mincells=10) 
  tosample <- apply(classified, 2, function(x) names(which(x==1)))
  sat <- list()
  for (i in names(tosample)) {
    for (j in tosample[[i]]) {
      combo <- paste(i, j, sep="_")
      sat[[combo]] <- .whichcells(txis, i, j)
    }
  }
  res <- .nv(do.call(c, lapply(sat, .getcells, maxcells=maxcells)))
  scheme <- list(mincells=mincells, maxcells=maxcells, 
                 eligible=attr(classified, "eligible"))

  attr(res, "scheme") <- scheme

  if (ret == "sce") { 
    metadata(txis)$downsampling <- res
    return(txis[, res])
  } else { 
    return(res) 
  } 

}


# for block sampling
.getcells <- function(cells, maxcells=20) {

  maxcells <- min(length(cells), maxcells)
  sample(cells, maxcells)

}


# for block sampling 
.whichcells <- function(txis, samp, clust) {

  picked <- which(txis$sample == samp & colLabels(txis) == clust)
  colnames(txis)[picked]

}


# helper fn
.nv <- function(x) { 

  names(x) <- x
  x

}
