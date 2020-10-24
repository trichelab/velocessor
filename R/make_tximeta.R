#' downstream processor for make_spliced_index results (post-salmon-index)
#'
#' @param params    the results from make_spliced_index
#' @param source    where the annotations came from originally (ENSEMBL)
#' @param release   what release of the transcriptome it is (NULL, don't know)
#' @param organism  which organism this index is for (guess from genome)
#' 
#' @return          a new entry in the bfc and a new linkedTxome
#' 
#' @import tximeta
#'
#' @export
make_tximeta <- function(params, source="ENSEMBL", release="", organism=NULL){
 
  if (is.null(organism)) {
    if (grepl("(GRCh|hg)", params$genome)) organism <- "Homo sapiens"
    else if (grepl("(GRCm|mm)", params$genome)) organism <- "Mus musculus"
    else if (grepl("(GRCz|dr)", params$genome)) organism <- "Danio rerio"
    else stop("Can't figure out what organism this came from, specify please.")
  } 

  jsonFile <- sub("gtf", "json", params$expandedGtf)

  with(params, 
       tximeta::makeLinkedTxome(indexDir = sidx,
                                source = source,
                                genome = genome,
                                organism = organism, 
                                release = release,
                                fasta = fa,
                                gtf = gtf,
                                write = TRUE,
                                jsonFile = json))

}
