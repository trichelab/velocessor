#' barebones cell labeling function
#' 
#' Uses singleR to label celltypes. More sophisticated users may want to work 
#' directly with SingleR input and output for (e.g.) specific gene sets.
#' 
#' @param txis        a SingleCellExperiment
#' @param species     what species is this from? (autodetect Hs/Mm) 
#' @param ret         what kind of object to return
#' @param downsample  downsample? (if ncol(txis) > 20000, TRUE, else FALSE) 
#' @param maxcells    maximum number of cells per cluster per sample (20)
#' @param mincells    minimum number of cells per cluster (10)
#' @param byclust     label cells by cluster? (FALSE, label individually)
#' @param label       which label to use ("label.main")
#' @param ...         other arguments passed to SingleR
#'
#' @return depending on the value of 'ret', either an SCE or a set of labels
#' 
#' @details
#' Autodetection of species will fail if the genome for the SingleCellExperiment
#' does not contain "GRCm", "mm", "GRCh", or "hg". Once reasonable reference 
#' datasets are available for GRCz/dr genomes, we will support those too. 
#' 
#' @import SingleR 
#'
#' @export
label_cells <- function(txis, species=NULL, ret=c("sce", "labels"), downsample=NULL, maxcells=20, mincells=10, label=c("label.main", "label.fine"), byClust=FALSE, ...) {

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
  label <- match.arg(label) 
  if (is.null(downsample)) downsample <- (ncol(txis) > 20000)
  if (downsample == TRUE) {
    message("Downsampling prior to cell labeling. Some labels will be NA.")
    cols <- downsample_txis(txis=txis, maxcells=maxcells, mincells=mincells)
    sclusts <- sum(apply(attr(cols, "scheme")$eligible > 0, 1, any))
    message(length(cols), " cells sampled from ", sclusts, " clusters.")
  }

  clusters <- NULL 
  ref <- switch(species, 
                "Homo sapiens"=celldex::HumanPrimaryCellAtlasData(), 
                "Mus musculus"=celldex::ImmGenData())
  stopifnot(label %in% names(colData(ref)))
  labelings <- colData(ref)[, label]
  if (byClust) clusters <- factor(colLabels(txis)[cols]) 
  rows <- intersect(rownames(txis), rownames(ref))
  pred <- SingleR(txis[rows, cols], ref[rows,], labels=labelings, 
                  clusters=clusters, ...)$pruned.labels

  # accommodate downsampling
  colData(txis)[, "celltype.sampled"] <- colnames(txis) %in% cols

  label <- rep(NA, ncol(txis))
  names(label) <- colnames(txis)
  label[cols] <- pred[cols]
  colData(txis)[, "celltype.label"] <- label

  return(switch(ret, sce=txis, labels=label))

}
