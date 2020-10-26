#' convenience functions for fetching ENSEMBL/Gencode annotations (genes & txs)
#' 
#' Gencode genes are based on ENSEMBL builds. When AnnotationHub's version of 
#' Gencode lags behind current (as often happens), it is possible to retrieve 
#' ENSEMBL builds instead. For example, Gencode 33 is ENSEMBL 99; GENCODE 35 is
#' ENSEMBL 101; and so forth. Moreover, given an ENSEMBL build and a species,
#' a minimum subset of information (gene id, gene name, gene biotype) is easily
#' retrieved for annotation purposes. So by default we just grab ENSEMBL genes.
#' 
#' @param version   the ENSEMBL (preferred) or Gencode (as gencodeXYZ) version
#' @param species   the species to fetch (default is "Homo sapiens") 
#' 
#' @return          ensdb object (fetch_ensembl) or GRanges[List] of genes/[txs]
#'
#' @examples
#'   fetch_ensembl("101", "Homo sapiens")
#'   get_ensembl_transcripts("101", "Mus musculus")
#'   get_ensembl_genes("gencode33")
#' 
#' @import ensembldb
#' @import AnnotationHub
#' @import GenomicFeatures
#'
#' @aliases get_ensembl_transcripts
#' @aliases get_ensembl_genes
#'
#' @export 
fetch_ensembl <- function(version, species) {

  ah <- AnnotationHub()

  if (grepl("gencode", version)) {
    g <- .gencode_mappings(gencode)
    version <- g$ensembl
    species <- g$species
  } else { 
    version <- sub("ensembl", "", version, ignore=TRUE)
  }

  ah[[query(ah, c("ensdb", version, species))$ah_id]]

}


# helper fn
.gencode_mappings <- function(gencode) {

  gencode <- sub("gencode", "", ignore=TRUE, gencode)
  
  # gencode35 is human ensembl 101; M25 is mouse ensembl 101
  gencode_to_ensembl <- c(as.character(101 - seq_len(35) + 1),  # human
                          as.character(101 - seq_len(25) + 1))  # mouse
  names(gencode_to_ensembl) <- c(as.character(seq(1, 35)),      # human
                                 paste0("M", 1:25))             # mouse
  gencode_to_species <- c(rep("Homo sapiens", 35),              # human
                          rep("Mus musculus", 25))              # mouse 
  names(gencode_to_species) <- c(as.character(seq(1, 35)),      # human
                                 paste0("M", 1:25))             # mouse

  res <- list(gencode=paste0("gencode", gencode),
              ensembl=gencode_to_ensembl[gencode],
              species=gencode_to_species[gencode])
  return(res) 

}
