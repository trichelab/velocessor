#' Just like get_ensembl_transcripts, but for gencode. 
#' 
#' @param   release   the GENCODE release string
#' @return            a GRangesList of transcripts grouped by gene
#' 
#' @export 
get_gencode_transcripts <- function(release) {

  get_ensembl_transcripts(version=get_gencode_ensembl(release), 
                          species=get_gencode_species(release))

}
