#' velocity wrapper, modularizes pipelined loading
#' 
#' @param txis    A SingleCellExperiment with assays 'spliced' and 'unspliced'
#' @param mode    scVelo mode (default "stochastic"; "dynamical" also available)
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
compute_velocity <- function(txis, mode="stochastic", ...) { 

  message("Removing dead cells and low-variance genes for velocity")
  mt <- names(subset(rowRanges(txis), seqnames %in% c("chrM", "chrMT", "MT")))
  txis$mtPercent <- (colSums(counts(txis[mt,])) / colSums(counts(txis))) * 100
  mtCut <- max(quantile(txis$mtPercent, 0.95), 10, na.rm=TRUE)
  live <- colnames(txis)[txis$mtPercent < mtCut]

  dimred <- ifelse("HARMONY" %in% reducedDimNames(txis), "HARMONY", "PCA")
  if (!dimred %in% reducedDimNames(txis)) txis <- runPCA(txis)
  
  dec <- modelGeneVar(txis[, live])
  HVGs <- scran::getTopHVGs(dec, n=1000)

  an <- c(X="spliced", spliced="spliced", unspliced="unspliced")
  sf <- lapply(an, function(x) scuttle::librarySizeFactors(txis, assay.type=x))

  message("Adding velocity...") 
  metadata(txis)$scVelo <- 
    velociraptor::scvelo(txis, subset.row=HVGs, dimred=dimred, mode=mode,
                         assay.X="spliced", sf.X=sf$spliced, 
                         sf.spliced=sf$spliced, sf.unspliced=sf$unspliced, ...) 

  message("Adding velocity_pseudotime...")
  txis$velocity_pseudotime <- metadata(txis)$scVelo$velocity_pseudotime
 
  # need to 1) fix this and 2) add cellrank support pronto 
  if (mode == "dynamical") message("Adding latent_time...")
  if (mode == "dynamical") txis$latent_time <- metadata(txis)$scVelo$latent_time
    
  message("Embedding velocity onto UMAP coordinates...")
  if (!"UMAP" %in% reducedDimNames(txis) | ncol(reducedDims(txis)$UMAP) < 3) { 
    if ("HARMONY" %in% reducedDimNames(txis)) {
      txis <- scater::runUMAP(txis, ncomponents=3, dimred="HARMONY")
    } else { 
      txis <- scater::runUMAP(txis, ncomponents=3)
    } 
  } 
  embedded <- 
    velociraptor::embedVelocity(reducedDims(txis)$UMAP, metadata(txis)$scVelo)
  metadata(txis)$embedded <- embedded

  message("Added scVelo (", mode, " mode) output to metadata(txis)$scVelo")
  return(txis) 

}
