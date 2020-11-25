#' Get or set rowData(object)[, biotype] columns in a SummarizedExperiment
#' 
#' @name  biotypes 
#' @rdname  biotypes 
#' @aliases gene_biotype
#' @aliases transcript_biotype
#' 
#' @param   object  a SummarizedExperiment
#' @param   value   value for rowData(object)[, biotype] (replaces existing)
#' @param   ...     other arguments, as appropriate (which they rarely are)
#'
#' @return          result of getting/setting rowData(object)[, biotype]
#'
#' @details GET_RD and SET_RD are aped from SingleCellExperiment.
#'
#' @import  SummarizedExperiment
NULL

# aped from SingleCellExperiment
GET_RD <- function(name, ...) {
  (name) # To ensure evaluation
  function(object, ...) rowData(object)[, name]
}

# stolen from SingleCellExperiment
SET_RD <- function(name, ...) {
  (name) # To ensure evaluation
  function(object, ..., value) {
    rowData(object)[, name] <- value
    object
  }
}

#' @rdname biotypes
#' @export
setGeneric("gene_biotype", 
           function(object, ...) standardGeneric("gene_biotype"))

#' @rdname biotypes
#' @export
setMethod("gene_biotype", "SummarizedExperiment",
          GET_RD("gene_biotype"))

#' @rdname biotypes
#' @export
setGeneric("gene_biotype<-", 
           function(object, ..., value) standardGeneric("gene_biotype<-"))

#' @rdname biotypes
#' @export
setReplaceMethod("gene_biotype", c("SummarizedExperiment", "ANY"), 
                 SET_RD("gene_biotype"))

#' @rdname biotypes
#' @export
setGeneric("transcript_biotype", 
           function(object, ...) standardGeneric("transcript_biotype"))

#' @rdname biotypes
#' @export
setMethod("transcript_biotype", "SummarizedExperiment", 
          GET_RD("transcript_biotype"))

#' @rdname biotypes
#' @export
setGeneric("transcript_biotype<-", 
           function(object, ..., value) standardGeneric("transcript_biotype<-"))

#' @rdname biotypes
#' @export
setReplaceMethod("transcript_biotype", c("SummarizedExperiment", "ANY"), 
                 SET_RD("transcript_biotype"))
