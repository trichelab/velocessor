#' Bayesian NMF and ARD with \pkg{SignatureAnalyzer}
#'
#' Use the fast, scalable, GPU-aware Bayesian NMF implementation in the Broad's
#' SignatureAnalyzer to factorize matrices and tensors for various purposes.
#'
#' @param x   A matrix. 
#'
#' @details
#' This function uses the \pkg{SignatureAnalyzer} Python package (found at 
#' \url{https://pypi.org/project/signatureanalyzer/}) for ARD and Bayesian NMF.
#' This is done via the \pkg{basilisk} package.
#'
#' @return
#' A list containing two matrices (W and H) is returned.
#'
#' @author Tim Triche, Jr. 
#' @name bnmf
NULL


#' @importFrom Matrix t
.bnmf <- function(x) { 

    mode <- match.arg(mode)
    output <- basiliskRun(env=bnmf.env, fun=.run_sa_bnmf, x=x)
    return(output) 
}


#' @importFrom reticulate import
#' @importFrom DelayedArray is_sparse t
#' @importFrom zellkonverter AnnData2SCE
.run_sa_bnmf <- function(x, sa.params=list()) {

  cr <- import("signatureanalyzer")
  cr$run_matrix(matrix=x)

}


#' @export
#' @rdname bnmf
setGeneric("bnmf", function(x, ...) standardGeneric("bnmf"))


#' @export
#' @rdname bnmf
setMethod("bnmf", "ANY", .bnmf)


#' @export
#' @rdname bnmf
#' @importFrom BiocGenerics sizeFactors
setMethod("bnmf", "SummarizedExperiment", function(x, ...) {
  .bnmf(list(matrix=assays(x)[[1]]))
})


#' @export
#' @rdname bnmf
setMethod("bnmf", "SingleCellExperiment", function(x, ...) {
  callNextMethod(x, ...)
})
