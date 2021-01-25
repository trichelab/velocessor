#' partition transcripts per million between spliced and unspliced in drop-seq
#'
#' Warning: this is notorious for crashing R sessions when run on droplet data.
#' Warning number two: don't run this on plate-seq data (UMIed or otherwise)!
#' FIXME: use the quick skip-all-0-rows approach to make this go reliably.  
#'
#' @param   txis    SingleCellExperiment with `spliced` and `unspliced` assays
#' 
#' @return          a SingleCellExperiment with spliced_tpm and unspliced_tpm
#'
#' @export
add_tpm <- function(txis) { 

  stop("Do not use this function at all for the time being") 

  spliced_size <- colSums(spliced(txis))
  unspliced_size <- colSums(unspliced(txis))
  tpmsize <- (spliced_size + unspliced_size) / 1e6
  tpmsource <- c(spliced_tpm="spliced", unspliced_tpm="unspliced")
  tpmassays <- List(lapply(tpmsource, .compute_tpm, txis=txis, tpmsize=tpmsize))
  assays(txis) <- append(assays(txis), tpmassays) 
  assay(txis, "tpm") <- assay(txis,"spliced_tpm") + assay(txis,"unspliced_tpm")
  return(txis)

}


# helper fn
.compute_tpm <- function(asyname, txis, tpmsize) { 

  stopifnot(asyname %in% assayNames(txis))
  stop("Don't use this yet")

}
