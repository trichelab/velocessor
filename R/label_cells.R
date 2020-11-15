#' barebones cell labeling function
#' 
#' Uses singleR to label celltypes. More sophisticated users may want to work 
#' directly with SingleR input and output for (e.g.) specific gene sets.
#' 
#' @param txis        a SingleCellExperiment
#' @param species     what species is this from? (autodetect Hs/Mm) 
#' @param ret         what kind of object to return
#' @param downsample  downsample? (if ncol(txis) > 20000, TRUE, else FALSE) 
#' @param maxcells    maximum number of cells per cluster per sample (50)
#' @param ...         other arguments passed to SingleR
#'
#' @return depending on the value of 'ret', either an SCE or a set of labels
#' 
#' @details
#' Autodetection of species will fail if the genome for the SingleCellExperiment
#' does not contain "GRCm", "mm", "GRCh", or "hg". 
#' 
#' @import SingleR 
#'
#' @export
label_cells <- function(txis, species=NULL, ret=c("sce", "labels"), downsample=NULL, maxcells=50, ...) {

  stopifnot(is(txis, "SingleCellExperiment"))
  if (is.null(species)) {
    if (grepl("^(GRCh|hg)", unique(genome(txis)))) {
      species <- "Homo sapiens"
    } else if (grepl("^(GRCm|mm)", unique(genome(txis)))) {
      species <- "Mus musculus"
#    } else if (grepl("^(GRCz|dr)", unique(genome(txis)))) {
#      species <- "Danio rerio"
    } else { 
      stop("Cannot automatically determine species to label cells. Exiting.") 
    } 
  } 

  cols <- colnames(txis)
  if (is.null(downsample)) downsample <- (ncol(txis) > 20000)
  if (downsample == TRUE) {
    message("Downsampling prior to cell labeling. Some labels will be NA.")
    cols <- sample_by_cluster_and_source(txis, maxcells=maxcells)
  }

  ret <- match.arg(ret)
  training <- switch(species, 
                     "Homo sapiens"=celldex::HumanPrimaryCellAtlasData(), 
                     "Mus musculus"=celldex::ImmGenData())
  rows <- intersect(rownames(txis), rownames(training))
  ref <- trainSingleR(training, training$label.main)

  # now label -- note that this is a ham-fisted default 
  pred <- SingleR(txis[, cols], ref, labels=ref$label, ...)$pruned.labels

  # accommodate downsampling
  colData(txis)[, "celltype.sampled"] <- colnames(txis) %in% cols

  label <- rep(NA, ncol(txis))
  names(label) <- colnames(txis)
  label[cols] <- pred[cols]

  colData(txis)[, "celltype.label"] <- label
  return(switch(ret,
                sce=txis, 
                labels=txis$celltype.label))

}
