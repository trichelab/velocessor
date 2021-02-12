#' find nonzero entries in a sparse (dgcMatrix) matrix
#' 
#' @param   x   the matrix
#' 
#' @return      for `nonzero`, a list with elements "nz", "nrow", and "ncol"
#' 
#' @details
#' `nz` is borrowed from the DAPAR package and is the base function.
#' `nz_rownames` returns the rownames for nonzero rows of `x`.
#' `nz_rows` just wraps unique(x@i) for a dgCMatrix `x`.
#' `nz_colnames` returns the colnames for nonzero columns of `x`.
#' `nz_cols` finds nonzero columns of a dgCMatrix `x`.
#' `nonzero` retains additional information about `x` alongside `nz`.
#' These are handy for (e.g.) computing TPMs and the like in 3'/5' data.
#' 
#' @import Matrix
#' @rdname nonzero 
#' 
#' @export
nonzero <- function(x, agg=NA) {

  stopifnot(inherits(x, "dgCMatrix"))
  res <- list(nz=nz(x), nrow=nrow(x), ncol=ncol(x))
  return(res)

}


#' @rdname nonzero
#' @export
nz <- function(x) { 

  stopifnot(inherits(x, "dgCMatrix"))
  if (all(x@p == 0)) {
    nz <- matrix(0, nrow=0, ncol=2, 
                 dimnames=list(character(0), c("row","col")))
  } else {
    nz <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
  }
  colnames(nz) <- c("row", "col")
  nz[x@x != 0, , drop = FALSE]

}


#' @rdname nonzero
#' @export
nz_rows <- function(x) sort(unique(x@i))


#' @rdname nonzero
#' @export
nz_rownames <- function(x) rownames(x)[nz_rows(x)]


#' @rdname nonzero
#' @export
nz_byrow <- function(x) {
  
  nzx <- nz(x)
  res <- rowsum(rep(1, nrow(nzx)), nzx[, "row"]) 
  rownames(res) <- rownames(x)[as.integer(rownames(res))]
  return(res[,1])

}


#' @rdname nonzero
#' @export 
nz_cols <- function(x) sort(unique(nz(x)[, "col"]))


#' @rdname nonzero
#' @export
nz_colnames <- function(x) colnames(x)[nz_cols(x)]


#' @rdname nonzero
#' @export
nz_bycol <- function(x) { 

  nzx <- nz(x)
  res <- rowsum(rep(1, nrow(nzx)), nzx[, "col"])
  rownames(res) <- colnames(x)[as.integer(rownames(res))]
  return(res[,1])

}
