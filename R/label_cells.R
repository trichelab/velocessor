#' barebones cell labeling function
#' 
#' Uses singleR to label celltypes. More sophisticated users may want to work 
#' directly with SingleR input and output for (e.g.) specific gene sets.
#' 
#' @param txis        a SingleCellExperiment
#' @param species     what species is this from? (autodetect Hs/Mm) 
#' @param ret         what kind of object to return
#' @param downsample  downsample? (if ncol(txis) > 20000, TRUE, else FALSE) 
#' @param maxcells    maximum number of cells _per_sample_ per cluster (20)
#' @param mincells    minimum number of cells per cluster (10)
#' @param label       which label to use ("label.main")
#' @param ref         a reference SummarizedExperiment (default: HPCD/ImmGen)
#' @param ...         other arguments passed to SingleR
#'
#' @return depending on the value of 'ret', either an SCE or a set of labels
#' 
#' @details
#' Autodetection of species will fail if the genome for the SingleCellExperiment
#' does not contain "GRCm", "mm", "GRCh", or "hg". Once reasonable reference 
#' datasets are available for GRCz/dr genomes, we will support those too.
#' 
#' Twenty cells per sample per cluster may be an awful lot of cells if you have
#' a large number of samples (i.e. GEMcode preps). Setting `maxcells` as low as 
#' 10, or perhaps even lower (see upcoming cluster-by-sample plots), if the 
#' cells comprising a cluster come from across many samples. The quick way to 
#' find this out is with a plot, 
#' 
#' @import SingleR 
#'
#' @export
label_cells <- function(txis, species=NULL, ret=c("sce", "labels"), downsample=NULL, maxcells=10, mincells=10, label="label.main", ref=NULL, ...) {

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

  ret <- match.arg(ret)
  cols <- colnames(txis)
  if (is.null(downsample)) downsample <- (ncol(txis) > 20000)
  if (downsample == TRUE) {
    message("Downsampling prior to cell labeling. Some labels will be NA.")
    cols <- downsample_txis(txis=txis, maxcells=maxcells, mincells=mincells)
    sclusts <- sum(apply(attr(cols, "scheme")$eligible > 0, 1, any))
    message(length(cols), " cells sampled from ", sclusts, " clusters.")
  }

  if (is.null(ref)) { 
    ref <- switch(species, 
                  "Homo sapiens"=celldex::HumanPrimaryCellAtlasData(), 
                  "Mus musculus"=celldex::ImmGenData())
  }
  stopifnot(label %in% names(colData(ref)))
  labelings <- colData(ref)[, label]
  rows <- intersect(rownames(txis), rownames(ref))
  pred <- SingleR(txis[rows, cols], ref[rows,], labels=labelings, ...)

  # accommodate downsampling
  colData(txis)[, "celltype.sampled"] <- colnames(txis) %in% cols

  celllabel <- rep("unlabeled", ncol(txis))
  names(celllabel) <- colnames(txis)
  celllabel[cols] <- pred[cols, "pruned.labels"]
  colData(txis)[, "celltype.label"] <- celllabel
  metadata(txis)$celltypePredictions <- pred

  return(switch(ret, sce=txis, labels=celllabel))

}
