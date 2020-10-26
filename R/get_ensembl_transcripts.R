#' @import ensembldb
#' 
#' @export 
get_ensembl_transcripts <- function(version, species="Homo sapiens") { 

  ensdb <- fetch_ensembl(version, species=species) 
  transcriptsBy(ensdb, "gene", columns=c("symbol", "gene_biotype"))

}
