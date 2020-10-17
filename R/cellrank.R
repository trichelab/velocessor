#' Cell trajectories with \pkg{CellRank}
#'
#' Use RNA velocity to find root or final states via the \pkg{CellRank} package.
#' Most of this is a straight-up ripoff of velociraptor so that I can learn how
#' to operate Python via basilisk for CellRank and SignatureFinder-GPU to run. 
#'
#' @param x   The output of velociraptor::scvelo(). 
#'
#' @details
#' This function uses the \pkg{CellRank} Python package (which can be found at 
#' \url{https://pypi.org/project/cellrank/}) for RNA velocity-pseudotime 
#' informed root and leaf finding when constructing trajectories and the like. 
#' This is done via the \pkg{basilisk} package - see the documentation for that package for trouble-shooting.
#'
#' @return
#' A \linkS4class{SingleCellExperiment} is returned containing the output.
#'
#' @examples
#' # Using mock data to demonstrate the process:
#' library(scuttle)
#' sce1 <- mockSCE()
#' sce2 <- mockSCE()
#'
#' spliced <- counts(sce1)
#' unspliced <- counts(sce2)
#'
#' mats <- list(X=spliced, spliced=spliced, unspliced=unspliced)
#' velo <- velociraptor::scvelo(mats)
#' cellrank <- cellrank(velo)
#'
#' @author Tim Triche, Jr. 
#' @name cellrank
NULL


#' @importFrom S4Vectors DataFrame
#' @importFrom scuttle normalizeCounts
#' @importFrom Matrix t
.cellrank <- function(velo, cellrank.params=list()) { 

    mode <- match.arg(mode)
    output <- basiliskRun(env=cellrank.env, fun=.run_cellrank, velo=velo)
    return(output) 
}


#' @importFrom reticulate import
#' @importFrom DelayedArray is_sparse t
#' @importFrom zellkonverter AnnData2SCE
.run_cellrank <- function(velo, cellrank.params=list()) {

  ad <- import("anndata")
  cr <- import("cellrank")
  adata <- ad$AnnData(velo, 
                      layers=list(spliced=spliced, unspliced=unspliced))
  adata$obs_names <- rownames(velo)
  adata$var_names <- colnames(spliced)
  # dimred <- .make_np_friendly(dimred)

  do.call(cr$pp$moments, c(list(data=adata), cellrank.params$moments))
  AnnData2SCE(adata)

}


#' @export
#' @rdname cellrank
setGeneric("cellrank", function(x, ...) standardGeneric("cellrank"))


#' @export
#' @rdname cellrank
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
setMethod("cellrank", "SingleCellExperiment", function(x, ...) {
    .cellrank(list(X=x))
})
