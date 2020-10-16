#' extract the sample/batch name for each cell from a SAMPLENAME_BARCODE scheme
#' 
#' This function simply assumes that process_velo_txis has renamed the cell 
#' barcodes in a SingleCellExperiment to indicate where they came from, and 
#' pops the last element of the split vector off the end.  This is handy for 
#' automatically running Harmony as a first pass to get rid of batch effects, 
#' though the right way to handle that task while respecting spliced-unspliced
#' relationships is currently an open research problem. 
#' 
#' @param txis  a SingleCellExperiment with colnames(txis) like FOO_BAR_AATGCT
#' @param sep   the separator character used to paste the column names ("_")
#' 
#' @return      a character vector of sample names, e.g. FOO_BAR
#' 
#' @import SingleCellExperiment
#' 
#' @export
get_sample_from_barcode <- function(txis, sep="_") {
  vapply(strsplit(colnames(txis), sep), .popOff, character(1))
}


# helper fn
.popOff <- function(x, sep="_") paste(x[-length(x)], collapse=sep)
