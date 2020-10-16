#' Load some Alevin runs quantified with velocity information 
#' 
#' The default assay (which is then split into spliced and unspliced matrices)
#' is rescaled TPM, which is to say, abundance in transcripts per million, 
#' scaled by the number of fragments/transcripts observed in each cell.
#'
#' We have found when comparing single-cell protocols against bulk and against
#' each other that this allows for both correction of sampling bias, and also 
#' the use of count-like representations instead of purely proportional TPMs. 
#' If you don't have a good reason to pass something other than 'abundance' in
#' the optional arguments to 'asy', we would suggest leaving it as is. 
#'
#' The results of this function can be further QC'ed and/or batch corrected, 
#' for example with Harmony or sva::ComBat, then passed to velociraptor::scvelo
#' and cellrank::cellrank for data-driven trajectory/cell-fate inference. We
#' are investigating whether a factorization approach may allow for sensible 
#' batch effect correction while respecting the relationships between spliced
#' and unspliced transcript abundance, but have not yet finalized it. 
#' 
#' See also plotVelocity for a relatively straightforward exploratory plot 
#' based on UMAP coordinates and a downsampled subset of cells from the 
#' SingleCellExperiment this function and velociraptor::scvelo yields. 
#' 
#' @param runs    a named vector of paths (if unnamed, you'll have problems)
#' @param txstub  the shared portion of the transcriptome annotation files 
#' @param qm      the name and path of the Alevin quants matrix file in 'runs' 
#' @param anno    optional GRanges linking e.g. ENSG transcripts to HGNC symbols
#' @param QC      do rudimentary quality control on each Alevin dataset? (TRUE)
#' @param HARMONY run Harmony on merged dataset? (FALSE, as harmony is non-BioC)
#' @param CLUSTER apply simple Louvain clustering to the result? (TRUE) 
#' @param DEDUPE  de-dupe genes and cells (doublets)? (FALSE unless SCVELO=TRUE)
#' @param SCVELO  run scVelo on the (perhaps harmonized) result? (FALSE, slow)
#' @param ...     other arguments to pass to velociraptor::scvelo, if any 
#' @param BPPARAM if running in parallel, provide a MulticoreParam or similar
#' 
#' @return        a SingleCellExperiment with 
#' 
#' @import scds
#' @import eisaR
#' @import scran
#' @import tximeta
#' @import jsonlite
#' @import fishpond
#' @import basilisk
#' @import scDblFinder
#' @import velociraptor
#' @import BiocParallel 
#' 
#' @export
process_velo_txis <- function(runs, txstub, anno=NULL, qm="alevin/quants_mat.gz", QC=TRUE, HARMONY=FALSE, CLUSTER=TRUE, DEDUPE=FALSE, SCVELO=FALSE, BPPARAM=SerialParam()){

  stopifnot(!is.null(names(runs)))
  txome <- file.path(txpath, paste(txstub, "json", sep="."))
  tximeta::loadLinkedTxome(jsonFile=txome)
  
  feats <- file.path(txpath, paste(txstub, "features", "tsv", sep="."))
  cg <- read.delim(feats, header=TRUE, as.is=TRUE)
  colnames(cg)[colnames(cg) == "intron"] <- "unspliced"

  gtf <- jsonlite::fromJSON(file=txome)[[1]][["gtf"]]
  stopifnot(file.exists(gtf)) # so tximeta doesn't puke

  cdat <- data.frame(names=names(runs),
                     files=paste(runs, qm, sep="/"),
                     stringsAsFactors=FALSE)
  rownames(cdat) <- cdat$names

  # to avoid (some) issues
  if (any(duplicated(cdat$names))) message("Warning: autorename may fail later")

  # can parallelize here; do QC post-merge
  by_sample <- split(cdat, cdat$names)
  txis <- do.call(cbind, bplapply(by_sample, import_velo_txis,
                                  cg=cg, QC=QC, BPPARAM=BPPARAM))
  txis$sample <- get_sample_from_barcode(txis) # better than auto

  # add annotation? 
  if (!is.null(anno) & all(rownames(txis) %in% names(anno))) {
    rowRanges(txis) <- anno[rownames(txis)] 
  } 

  # cluster?
  if (CLUSTER | DEDUPE | SCVELO) {
    txis$cluster <- cluster_velo_txis(txis)
  }
  
  # dedupe? 
  if (DEDUPE | SCVELO) {
    txis <- dedupe_velo_txis(txis)
    txis <- scDblFinder(txis, clusters=txis$cluster, samples=txis$sample)
    txis <- txis[, txis$scDblFinder.class == "singlet"]
  }

  # batch correct?
  if (HARMONY) txis <- harmonize_velo_txis(txis)

  # find velocity?
  if (SCVELO) {
    message("Removing dead cells and low-variance genes for velocity")
    mt <- names(subset(rowRanges(txis), seqnames %in% c("chrM", "chrMT", "MT")))
    txis$mtPercent <- (colSums(counts(txis[mt,])) / colSums(counts(txis))) * 100
    mtCut <- max(quantile(txis$mtPercent, 0.95), 10)
    live <- colnames(txis)[txis$mtPercent < mtCut]
    dec <- modelGeneVar(txis[, live])
    HVGs <- getTopHVGs(dec, n=1000)
    metadata(txis)$scVelo <- 
      scvelo(txis, subset.row=HVGs, assay.X="spliced", mode="stochastic", ...) 

    message("Adding velocity pseudotime...")
    txis$velocity_pseudotime <- metadata(txis)$scVelo$velocity_pseudotime
    
    message("Embedding velocity onto UMAP coordinates...")
    embedded <- embedVelocity(reducedDim(txis, "UMAP"), metadata(txis)$scVelo)
    metadata(txis)$embedded <- embedded

    message("Added scVelo (stochastic mode) output to metadata(txis)$scVelo")
  }

  # done
  return(txis)

}
