#' @import ensembldb
#' 
#' @export 
get_ensembl_genes <- function(version, species="Homo sapiens") { 

  ensdb <- fetch_ensembl(version=version, species=species)
  genes(ensdb, columns=c("symbol", "gene_biotype"))

}
