#' Get or set metadata items for velocitized SingleCellExperiments
#' 
#' Specifically, `scVelo` and `embedded` have their own generics. 
#'
#' @name    velocity_metadata
#' @rdname  velocity_metadata
#' @aliases embedded
#' @aliases scVelo
#' 
#' @param   object  a SingleCellExperiment
#' @param   value   value to replace existing metadata(object)[[name]]
#' @param   ...     other arguments, as appropriate (which they rarely are)
#'
#' @return          result of getting or setting metadata(object)[[name]]
#'
#' @import SingleCellExperiment
NULL

# aped from SingleCellExperiment
GET_MD <- function(name, ...) {
  (name) # To ensure evaluation
  function(object, ...) metadata(object)[[name]]
}

# stolen from SingleCellExperiment
SET_MD <- function(name, ...) {
  (name) # To ensure evaluation
  function(object, ..., value) {
    metadata(object)[[name]] <- value
    object
  }
}

#' @rdname  velocity_metadata
#' @export
setGeneric("scVelo", function(object, ...) standardGeneric("scVelo"))

#' @rdname  velocity_metadata
#' @export
setMethod("scVelo", "SingleCellExperiment", GET_MD("scVelo"))

#' @rdname  velocity_metadata
#' @export
setGeneric("scVelo<-", 
           function(object, ..., value) standardGeneric("scVelo<-"))

#' @rdname  velocity_metadata
#' @export
setReplaceMethod("scVelo", c("SingleCellExperiment","ANY"),
                 SET_MD("scVelo"))

#' @rdname  velocity_metadata
#' @export
setGeneric("embedded", function(object, ...) standardGeneric("embedded"))

#' @rdname  velocity_metadata
#' @export
setMethod("embedded", "SingleCellExperiment", GET_MD("embedded"))

#' @rdname  velocity_metadata
#' @export
setGeneric("embedded<-", 
           function(object, ..., value) standardGeneric("embedded<-"))

#' @rdname  velocity_metadata
#' @export
setReplaceMethod("embedded", c("SingleCellExperiment","ANY"), 
                 SET_MD("embedded"))
