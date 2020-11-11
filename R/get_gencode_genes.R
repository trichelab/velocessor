#' A wrapper for get_ensembl_genes, to handle Gencode release mapping
#'
#' Warning: this function is horrifically poorly tested 
#' 
#' @param  release   GENCODE release (can be got from tximeta json)
#' 
#' @return           a GRanges of genes with symbols and biotypes
#' 
#' @export 
get_gencode_genes <- function(release) { 

  get_ensembl_genes(version=get_gencode_ensembl(release), 
                    species=get_gencode_species(release))

}
