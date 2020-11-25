#' Get or set assays particular to velocitized SingleCellExperiments
#' 
#' Specifically, `spliced` and `unspliced` are made into generics. 
#'
#' @name    velocity_assays
#' @rdname  velocity_assays
#' @aliases unspliced
#' @aliases spliced
#' 
#' @param   object  a SingleCellExperiment
#' @param   value   the value for the assay in question 
#' @param   ...     other arguments, as appropriate (which they rarely are)
#'
#' @return          the result of getting or setting the assay
#'
#' @details GET and SET are stolen directly from SingleCellExperiment.
#' 
#' @seealso SingleCellExperiment
#'
#' @import SingleCellExperiment
NULL

# stolen from SingleCellExperiment
GET <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ...) {
    assay(object, i=exprs_values, ...)
  }
}

# stolen from SingleCellExperiment
SET <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ..., value) {
    assay(object, i=exprs_values, ...) <- value
    object
  }
}

#' @rdname velocity_assays
#' @export
setGeneric("spliced", function(object, ...) standardGeneric("spliced"))

#' @rdname velocity_assays
#' @export
setMethod("spliced", "SingleCellExperiment", GET("spliced"))

#' @rdname velocity_assays
#' @export
setGeneric("spliced<-", 
           function(object, ..., value) standardGeneric("spliced<-"))

#' @rdname velocity_assays
#' @export
setReplaceMethod("spliced", c("SingleCellExperiment","ANY"), SET("spliced"))

#' @rdname velocity_assays
#' @export
setGeneric("unspliced", function(object, ...) standardGeneric("unspliced"))

#' @rdname velocity_assays
#' @export
setMethod("unspliced", "SingleCellExperiment", GET("unspliced"))

#' @rdname velocity_assays
#' @export
setGeneric("unspliced<-", 
           function(object, ..., value) standardGeneric("unspliced<-"))

#' @rdname velocity_assays
#' @export
setReplaceMethod("unspliced", c("SingleCellExperiment","ANY"), SET("unspliced"))
