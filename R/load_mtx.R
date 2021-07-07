#' Load a triple of gene names, cell names, gene x cell counts as from bustools
#' 
#' The returned matrix is a dgCMatrix as used by (e.g.) SingleCellExperiment.
#' This function is useful when processing data from (e.g.) kallisto|bustools, 
#' for example, CITE-seq data processed with KITE, or when comparing uniquely 
#' mapped vs. multimapped counts, as one might do for repetitive elements.
#' 
#' @param   path      where the counts files (features, barcodes, mtx) (".")
#' @param   verbose   be verbose? (TRUE)
#' @param   frags     optional patterns to scan for features, barcodes, matrix
#' @param   splt      split by feature type? (FALSE)
#' @param   spltcol   if splitting by feature type, which column to use? (3)
#' 
#' @return            a dgCMatrix 
#'
#' @details
#' This function works fine on gzipped files (tested on a 180K x 2761K matrix). 
#'
#' @import  Matrix
#' 
#' @export
load_mtx <- function(path=".", verbose=TRUE, frags=NULL, splt=FALSE, spltcol=3){

  if (is.null(frags)) {
    frags <- c(features="(features|genes).(txt|tsv)", 
               barcodes="(cells|barcodes).(txt|tsv)",
               mat=".mtx")
  }
  stopifnot(all(c("features", "barcodes", "mat") %in% names(frags)))

  files <- lapply(frags, function(fr) grep(fr, list.files(path), value=TRUE))
  files <- lapply(files, function(fl) file.path(path, fl))
  dat <- with(files, .readSparseMat(mat))

  rows <- with(files, read.table(features)[,1])
  if (length(rows) == ncol(dat)) dat <- as(t(dat), "CsparseMatrix")
  if(verbose) message("Loaded ", nrow(dat), " rows x ", ncol(dat), " columns.")
  stopifnot(length(rows) == nrow(dat))
  rownames(dat) <- rows
  if(verbose) message("Labeled ", nrow(dat), " rows.")

  cols <- with(files, read.table(barcodes)[,1])
  stopifnot(length(cols) == ncol(dat))
  colnames(dat) <- cols
  if(verbose) message("Labeled ", ncol(dat), " columns.")

  if (splt) {
    rowsplit <- with(files, read.table(features)[, spltcol])
    dat <- split(dat, rowsplit)
  } 
  return(dat)

}


# helper fn, refactored out of the above for ease of repurposing
.readSparseMat <- function(mat) as(t(readMM(mat)), "CsparseMatrix")
