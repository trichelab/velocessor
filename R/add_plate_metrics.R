#' some metrics than have come in handy for plate-seq experiments
#' 
#' Prerequisite: annotate the damned gene symbols (with e.g. get_gencode_genes)
#' 
#' @param txis      a SingleCellExperiment, usually from import_plate_txis
#' @param ...       other parameters, currently just passed to izar_transform
#' 
#' @return          a SingleCellExperiment with additional metrics in metadata
#'
#' @seealso         make_plate_plot
#'
#' @import          scater
#' @import          uwot  
#' 
#' @export
add_plate_metrics <- function(txis, ...) { 
  
  stopifnot("symbol" %in% names(mcols(txis)))
  stopifnot("gene_biotype" %in% names(mcols(txis)))

  txis <- dedupe_velo_txis(txis)
  rownames(txis) <- mcols(txis)$symbol
  txis <- izar_transform(txis, ...) 
  mcols(txis)$expressed <- (rowSums2(assay(txis, "izar") >= 1) > 0)
  by_type <- as.matrix(with(mcols(txis), table(gene_biotype, expressed)))
  biotypes <- data.frame(expressed=by_type[,2], unexpressed=by_type[,1])
  add <- median(biotypes$expressed)
  biotypes$est <- with(biotypes, round(expressed/(expressed+unexpressed+add),2))
  metadata(txis)$expression_by_biotype <- biotypes
  txis$NumGenesExpressed <- colSums(assay(txis, "izar") >= 1)
  txis <- scater::runPCA(txis, exprs_values="izar") # plate 
  txis <- cluster_velo_txis(txis, use="PCA") # plate
  txis <- add_umap_embedding(txis, rdn="PCA") # plate 
  return(txis)

}
