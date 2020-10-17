#' stupid-simple Louvain clustering. Doesn't currently care about unspliced.
#' 
#' If a reduced dimension representation named "HARMONY" exists, the graph will
#' be constructed on that, otherwise "PCA" will be used (or added if not found)
#' unless use.dimred is explicitly specified by the user. 
#' 
#' @param txis        a SingleCellExperiment
#' @param k           number of nearest neighbors (20) 
#' @param type        what type of SNNGraph to construct ("rank")
#' @param use.dimred  which reducedDimension embedding to use (autodetermine)
#' @param ...         optional arguments to pass to igraph::cluster_louvain()
#' 
#' @return            Louvain cluster assignments for each cell. 
#' 
#' @import scran
#' @importFrom igraph cluster_louvain
#' 
#' @export 
cluster_velo_txis <- function(txis, k=20, type="rank", use.dimred=NULL, ...) { 

  if (is.null(use.dimred)) {
    use.dimred <- "PCA"
    if (!"PCA" %in% reducedDimNames(txis)) txis <- runPCA(txis)
    if ("HARMONY" %in% reducedDimNames(txis)) use.dimred <- "HARMONY"
  }

  stopifnot(use.dimred %in% reducedDimNames(txis))
  message("Note: at the moment, we do not use unspliced information here.") 
  graf <- scran::buildSNNGraph(txis, k=k, use.dimred=use.dimred, type=type)
  louvain <- igraph::cluster_louvain(graf, ...)
  clusters <- paste0("cluster", louvain$membership) 
  names(clusters) <- colnames(txis)
  return(clusters) 

}
