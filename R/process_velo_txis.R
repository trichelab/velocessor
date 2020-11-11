#' Load some Alevin runs quantified with velocity information 
#' 
#' This is evolving rapidly to handle microwell-seq and plate-seq data.
#'
#' The results of this function can be further QC'ed and/or batch corrected, 
#' for example with Harmony or sva::ComBat, then passed to velociraptor::scvelo
#' and/or velocessor::cellrank for data-driven cell fate inference. We are 
#' investigating various factorization approaches that may allow for sensible 
#' batch effect correction while respecting the relationships between spliced
#' and unspliced transcript abundance, but have not yet settled on a default.
#' 
#' See also plot_velo for a relatively straightforward exploratory plot 
#' based on UMAP coordinates and the downsampled subset of cells from the 
#' SingleCellExperiment that this function and velociraptor::scvelo yields. 
#' 
#' @param runs    a named vector of paths (if unnamed, you'll have problems)
#' @param txstub  the shared portion of the transcriptome annotation files 
#' @param qm      the name and path of the Alevin quants matrix file in 'runs' 
#' @param QC      do rudimentary quality control on each Alevin dataset? (TRUE)
#' @param HARMONY run Harmony on merged dataset? (FALSE, as harmony is non-BioC)
#' @param CLUSTER apply simple Louvain clustering to the result? (TRUE) 
#' @param DEDUPE  de-dupe genes and cells (doublets)? (FALSE unless SCVELO=TRUE)
#' @param SCVELO  run scVelo on the (perhaps harmonized) result? (FALSE, slow)
#' @param meta    use tximeta for import_velo_txis? (FALSE; avoid bfc issues)
#' @param ...     other arguments to pass to velociraptor::scvelo, if any 
#' @param BPPARAM if running in parallel, provide a MulticoreParam or similar
#' 
#' @return        a SingleCellExperiment with 
#' 
#' @import eisaR
#' @import scran
#' @import tximeta
#' @import jsonlite
#' @import fishpond
#' @import basilisk
#' @import velociraptor
#' @import BiocParallel 
#' 
#' @export
process_velo_txis <- function(runs, txstub, qm="alevin/quants_mat.gz", QC=TRUE, HARMONY=FALSE, CLUSTER=TRUE, DEDUPE=FALSE, SCVELO=FALSE, meta=FALSE, ..., BPPARAM=SerialParam()) {

  # sanity checking prior to processing  
  stopifnot(!is.null(names(runs)))
 
  txome <- paste(txstub, "json", sep=".")
  stopifnot(file.exists(txome))
  tximeta::loadLinkedTxome(jsonFile=txome)
  txjson <- jsonlite::fromJSON(txome)
  if (tolower(txjson[["source"]]) == "gencode") {
    anno <- get_gencode_genes(release=txjson[["release"]])
  } else { 
    anno <- get_ensembl_genes(version=txjson[["release"]],
                              species=txjson[["organism"]])
  }
    
  gtf <- txjson[1, "gtf"]
  stopifnot(file.exists(gtf))
  
  # sanity checking prior to processing  
  feats <- paste(txstub, "features", "tsv", sep=".")
  stopifnot(file.exists(feats))

  # this can feed either 
  cdat <- data.frame(names=names(runs),
                     files=paste(runs, qm, sep="/"),
                     txome=txome)
  rownames(cdat) <- cdat$names

  # to avoid (some) issues
  if (any(duplicated(cdat$names))) message("Warning: autorename may fail later")

  # can parallelize here; do QC post-merge
  txis <- do.call(cbind, 
                  bplapply(split(cdat, cdat$names),
                           import_velo_txis, 
                           quant="alevin",
                           meta=meta, 
                           QC=QC, 
                           BPPARAM=BPPARAM))

  txis$sample <- get_sample_from_barcode(txis) # better than auto
  txis$imported_by <- factor(paste("velocessor", packageVersion("velocessor")))

  # useful downstream 
  txis <- runPCA(txis)

  # add annotation? 
  rowData(txis)$ENSG <- sapply(strsplit(rownames(txis), "\\."), `[`, 1)
  if (!is.null(anno) & all(rowData(txis)$ENSG %in% names(anno))) {
    rowRanges(txis) <- anno[rowData(txis)$ENSG]
  } 
  
  # optional embellishments
  if (CLUSTER | DEDUPE | SCVELO) txis <- .recluster(txis) 
  if (HARMONY) txis <- .recluster(harmonize_velo_txis(txis))
  if (DEDUPE | SCVELO) txis <- drop_doublets(dedupe_velo_txis(txis))
  if (DEDUPE & HARMONY) txis <- harmonize_velo_txis(txis)
  if (SCVELO) txis <- compute_velocity(txis, ...) 

  # done
  return(txis)

}


# helper fn
.recluster <- function(txis) cluster_velo_txis(txis, ret="sce")  
