#' convenience function for fetching ENSEMBL annotations (genes & txs)
#' 
#' @param version   the ENSEMBL version
#' @param species   the ENSEMBL species ("Homo sapiens")
#' 
#' @return          an ensdb object
#'
#' @examples
#'   fetch_ensembl("101", "Homo sapiens")
#'   get_ensembl_transcripts("101", "Mus musculus")
#'   get_gencode_transcripts("m24")
#'   get_gencode_genes("gencode33")
#' 
#' @import AnnotationHub
#' @import GenomicFeatures
#'
#' @aliases get_ensembl_transcripts
#' @aliases get_ensembl_genes
#'
#' @export 
fetch_ensembl <- function(version, species) {

  ah <- AnnotationHub()
  version <- sub("ensembl", "", version, ignore=TRUE)
  ah[[query(ah, c("ensdb", version, species))$ah_id]]

}
