#' barebones cell labeling function
#' 
#' Uses singleR to label celltypes. More sophisticated users may want to work 
#' directly with SingleR input and output for (e.g.) specific gene sets.
#' 
#' @param txis        a SingleCellExperiment
#' @param species     what species is this from? ("Homo sapiens" -- automate?) 
#' @param ret         what kind of object to return
#' @param downsample  downsample? (if ncol(txis) > 20000, TRUE, else FALSE) 
#' @param maxcells    maximum number of cells per cluster per sample (50)
#' @param ...         other arguments passed to SingleR
#'
#' @return depending on the value of 'ret', either an SCE or a set of labels
#' 
#' @import SingleR 
#' @import celldex
#'
#' @export
label_cells <- function(txis, species=c("Homo sapiens", "Mus musculus"), ret=c("sce", "labels"), downsample=NULL, maxcells=50, ...) {

  ret <- match.arg(ret)
  species <- match.arg(species) 
  training <- switch(species, 
                     "Homo sapiens"=celldex::HumanPrimaryCellAtlasData(), 
                     "Mus musculus"=celldex::ImmGenData())
  rows <- intersect(rownames(txis), rownames(training))
  trained <- trainSingleR(training[rows,], training$label.main)

  # downsample?
  cols <- colnames(txis)
  if (is.null(downsample)) downsample <- (ncol(txis) > 20000)
  if (downsample == TRUE) {
    message("Downsampling prior to cell labeling. Some labels will be NA.")
    cols <- sample_by_cluster_and_source(txis, maxcells=maxcells)
  }

  # now label -- note that this is a ham-fisted default 
  pred <- SingleR(txis[cols,], ref, labels=ref$label, ...)$pruned.labels

  # accommodate downsampling
  sampled <- colnames(txis) %in% cols
  names(sampled) <- colnames(txis)
  colData(txis)[, "celltype.sampled"] <- sampled

  label <- rep(NA, ncol(txis))
  names(label) <- colnames(txis)
  label[cols] <- pred[cols]

  colData(txis)[, "celltype.label"] <- celltype
  return(switch(ret,
                sce=txis, 
                labels=txis$celltype.label))

}
