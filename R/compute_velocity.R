#' velocity wrapper, modularizes pipelined loading
#' 
#' @param txis    A SingleCellExperiment with assays 'spliced' and 'unspliced'
#' @param scvmode scVelo mode (default "stochastic"; "dynamical" also available)
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
compute_velocity <- function(txis, scvmode="stochastic", ...) { 

  message("Removing dead cells and low-variance genes for velocity")
  mt <- names(subset(rowRanges(txis), seqnames %in% c("chrM", "chrMT", "MT")))
  txis$mtPercent <- (colSums(counts(txis[mt,])) / colSums(counts(txis))) * 100
  mtCut <- max(quantile(txis$mtPercent, 0.95), 10, na.rm=TRUE)
  live <- colnames(txis)[txis$mtPercent < mtCut]

  dimred <- ifelse("HARMONY" %in% reducedDimNames(txis), "HARMONY", "PCA")
  if (!dimred %in% reducedDimNames(txis)) txis <- runPCA(txis)
  
  dec <- modelGeneVar(txis[, live])
  HVGs <- scran::getTopHVGs(dec, n=1000)

  # mandatory for velociraptor::scvelo:
  an <- c(spl="spliced", unspl="unspliced")
  sf <- lapply(an, function(x) scuttle::librarySizeFactors(txis, assay.type=x))
  sfs <- do.call(cbind, sf)
  kept <- names(which(apply(sfs, 1, function(x) all(x > 0))))

  message("Adding velocity...") 
  metadata(txis)$scVelo <- velociraptor::scvelo(txis[, kept], 
                                                subset.row=HVGs, 
                                                dimred=dimred, 
                                                assay.X="spliced", 
                                                sf.X=sfs[kept, "spl"],
                                                sf.spliced=sfs[kept, "spl"],
                                                sf.unspliced=sfs[kept, "unspl"],
                                                mode=scvmode, 
                                                ...) 

  message("Adding velocity_pseudotime...")
  txis$velocity_pseudotime <- NA
  with(metadata(txis),
       colData(txis)[kept, "velocity_pseudotime"] <- scVelo$velocity_pseudotime)
 
  # need to 1) fix this and 2) add cellrank support pronto 
  if (scvmode == "dynamical") {
    message("Adding latent_time...")
    txis$latent_time <- NA
    with(metadata(txis),
         colData(txis)[kept, "latent_time"] <- scVelo$latent_time)
  }
    
  message("Embedding velocity onto UMAP coordinates...")
  if (!"UMAP" %in% reducedDimNames(txis) | ncol(reducedDim(txis, "UMAP") < 3)) {
    if ("HARMONY" %in% reducedDimNames(txis)) {
      txis <- scater::runUMAP(txis, ncomponents=3, dimred="HARMONY")
    } else { 
      txis <- scater::runUMAP(txis, ncomponents=3)
    } 
  } 
  embedded <- 
    velociraptor::embedVelocity(reducedDim(txis[, kept], "UMAP"), 
                                metadata(txis)$scVelo)
  metadata(txis)$embedded <- embedded

  message("Added scVelo (", scvmode, " mode) output to metadata(txis)$scVelo")
  return(txis) 

}
