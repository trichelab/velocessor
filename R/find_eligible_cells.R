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
  classified <- .fit_mixtures(plausible)
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


# helper fn -- fit mixture models on log1p(cells) with fallback onto I(cells)
.fit_mixtures <- function(cells) { 

  logcells <- log1p(cells)
  cellfits <- apply(cells, 1, densityMclust, G=1:2, verbose=FALSE)
  logcellfits <- apply(logcells, 1, densityMclust, G=1:2, verbose=FALSE)

  modelnames <- c("cells", "logcells")
  names(modelnames) <- modelnames
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
