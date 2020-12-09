#' Extract the sample/batch name for each cell from a SAMPLENAME_BARCODE scheme
#' 
#' This function assumes that process_velo_txis or similar has renamed the cell 
#' barcodes (colnames) in a SingleCellExperiment to indicate whence they came, 
#' and pops the last element of the split vector off the end.  This is handy 
#' for running Harmony as a first pass to get rid of batch effects, although 
#' the right way to do this while respecting spliced-unspliced relationships 
#' is currently an open research problem. 
#' 
#' @param txis  SingleCellExperiment or vector with entries like FOO_BAR_AATGCT
#' @param sep   the separator character used to paste the column names ("_")
#' 
#' @return      a character vector of sample names, e.g. FOO_BAR
#' 
#' @import      SingleCellExperiment
#' 
#' @export
get_sample_from_barcode <- function(txis, sep="_") {
  if (is(txis, "SingleCellExperiment")) .getSampleNames(colnames(txis), sep=sep)
  else .getSampleNames(txis, sep=sep)
}


# helper fn
.getSampleNames <- function(x, sep="_") { 
  vapply(strsplit(x, sep), .popOff, sep=sep, character(1))
}


# helper fn
.popOff <- function(x, sep="_") {
  paste(x[-length(x)], collapse=sep)
}
