#' Load a triple of gene names, cell names, gene x cell counts as from bustools
#' 
#' The returned matrix is a dgCMatrix as used by (e.g.) SingleCellExperiment.
#' This function is useful when processing data from (e.g.) kallisto|bustools, 
#' for example, CITE-seq data processed with KITE, or when comparing uniquely 
#' mapped vs. multimapped counts, as one might do for repetitive elements.
#' 
#' @param   path      where the counts files are (genes, barcodes, mtx)
#' @param   verbose   be verbose? (TRUE) 
#' 
#' @return            a dgCMatrix 
#' 
#' @import  Matrix
#' 
#' @export
load_mtx <- function(path, verbose=TRUE) {

  frags <- c(mat=".mtx", genes=".genes.txt", barcodes=".barcodes.txt")
  files <- lapply(frags, function(fr) grep(fr, list.files(path), value=TRUE))
  files <- lapply(files, function(fl) file.path(path, fl))
  dat <- as(with(files, t(readMM(mat))), "CsparseMatrix")
  if(verbose) message("Loaded ", nrow(dat), " rows x ", ncol(dat), " columns.")

  rows <- with(files, read.table(genes)[,1])
  stopifnot(length(rows) == nrow(dat))
  rownames(dat) <- rows
  if(verbose) message("Labeled ", nrow(dat), " rows.")

  cols <- with(files, read.table(barcodes)[,1])
  stopifnot(length(barcodes) == ncol(dat))
  colnames(dat) <- cols
  if(verbose) message("Labeled ", ncol(dat), " columns.")

  return(dat)

}
