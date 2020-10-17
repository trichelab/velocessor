#' velocity wrapper, modularizes pipelined loading
#' 
#' @param txis    A SingleCellExperiment with assays 'spliced' and 'unspliced'
#' @param ...     Additional arguments passed to velociraptor::scvelo 
#'
#' @return        A SingleCellExperiment with scVelo output in metadata() 
#' 
#' @import scran
#' @import velociraptor
#' @import BiocGenerics
#' @import SingleCellExperiment
#' 
#' @export
compute_velocity <- function(txis, ...) { 

  message("Removing dead cells and low-variance genes for velocity")
  mt <- names(subset(rowRanges(txis), seqnames %in% c("chrM", "chrMT", "MT")))
  txis$mtPercent <- (colSums(counts(txis[mt,])) / colSums(counts(txis))) * 100
  mtCut <- max(quantile(txis$mtPercent, 0.95), 10)
  live <- colnames(txis)[txis$mtPercent < mtCut]
  dec <- modelGeneVar(txis[, live])
  HVGs <- scran::getTopHVGs(dec, n=1000)

  message("Adding velocity...") 
  metadata(txis)$scVelo <- 
    scvelo(txis, subset.row=HVGs, assay.X="spliced", mode="stochastic", ...) 

  message("Adding velocity pseudotime...")
  txis$velocity_pseudotime <- metadata(txis)$scVelo$velocity_pseudotime
    
  message("Embedding velocity onto UMAP coordinates...")
  embedded <- embedVelocity(reducedDim(txis, "UMAP"), metadata(txis)$scVelo)
  metadata(txis)$embedded <- embedded

  message("Added scVelo (stochastic mode) output to metadata(txis)$scVelo")
  return(txis) 

}
