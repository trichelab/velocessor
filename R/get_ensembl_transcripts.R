#' Get transcripts grouped by gene name. 
#'
#' @param   version   the ENSEMBL release version
#' @param   species   the ENSEMBL species ("Homo sapiens")
#' 
#' @return            a GRangesList of transcripts grouped by gene
#' 
#' @importFrom ensembldb transcripts
#' @importFrom ensembldb transcriptsBy
#' 
#' @export 
get_ensembl_transcripts <- function(version, species="Homo sapiens") { 

  ensdb <- fetch_ensembl(version, species=species) 
  transcriptsBy(ensdb, "gene", columns=c("symbol", "gene_biotype"))

}
