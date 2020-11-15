#' convenience function for merging datasets after QC (e.g. doublet finding) 
#' 
#' This function assumes that harmonized UMAP and maybe scVelo objects exist
#' (the latter in metadata(x) and/or metadata(y)), but also assumes that the 
#' user may not want to brute-force merge them, so they are saved into the 
#' metadata slot of the new merged object for future reference. scDblFinder 
#' columns in rowRanges are not, however; they cause problems with cbind(). 
#' 
#' @param x   the first SCE
#' @param y   the second SCE 
#' 
#' @return    a merged SCE
#' 
#' @import SingleCellExperiment
#'
#' @export
merge_sces <- function(x, y) { 

  stopifnot(is(x, "SingleCellExperiment"))
  stopifnot(is(y, "SingleCellExperiment"))

  ncd_x <- names(colData(x))
  ncd_y <- names(colData(y))
  ncd_z <- union(ncd_x, ncd_y)

  for (i in setdiff(ncd_z, ncd_x)) colData(x)[, i] <- rep(NA, ncol(x))
  for (j in setdiff(ncd_z, ncd_y)) colData(y)[, i] <- rep(NA, ncol(y))

  rn_z <- intersect(rownames(x), rownames(y))
  if (length(rn_z) == 0) stop("x & y have no shared genes!") 
  if (length(rn_z) < 5000) message("x & y have fewer than 5000 shared genes!") 

  x <- x[rn_z, ] 
  y <- y[rn_z, ]

  nmc_x <- names(mcols(x))
  nmc_y <- names(mcols(y))
  nmc_z <- intersect(nmc_y, nmc_y)
  nmc_z <- nmc_z[!grepl("scDblFinder", nmc_z)]
  
  mcols(x) <- mcols(x)[, nmc_z]
  mcols(y) <- mcols(y)[, nmc_z]

  umap_x <- reducedDim(x, "UMAP")
  umap_y <- reducedDim(y, "UMAP")
  velo_x <- metadata(x)$scVelo
  velo_y <- metadata(y)$scVelo

  z <- scater::runPCA(cbind(x, y))
  z <- cluster_velo_txis(z, ret="sce") 
  z <- harmonize_velo_txis(z, ret="sce")
  z <- scater::runUMAP(z, dimred="HARMONY", ncomp=3)
  metadata(z) <- list(UMAP_x=umap_x, UMAP_y=umap_y, 
                      velo_x=velo_x, velo_y=velo_y)
  return(z) 

}
