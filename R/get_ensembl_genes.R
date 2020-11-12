#' Simple wrapper to annotate genes/transcripts with symbols and so forth 
#' 
#' @param version   ENSEMBL version (from tximeta or gencode mapping)
#' @param species   ENSEMBL species (from tximeta or gencode release)
#' 
#' @return          a GRanges of genes with symbols and biotypes
#' 
#' @importFrom ensembldb genes
#' 
#' @export 
get_ensembl_genes <- function(version, species="Homo sapiens") { 

  ensdb <- fetch_ensembl(version=version, species=species)
  genes(ensdb, columns=c("symbol", "gene_biotype"))

}
