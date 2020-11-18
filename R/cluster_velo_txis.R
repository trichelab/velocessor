#' Simple Louvain, density, or k-means clustering. Doesn't use unspliced (yet).
#' 
#' If a reduced dimension representation named "HARMONY" exists, the graph will
#' be constructed on that, otherwise "PCA" will be used (or added if not found)
#' unless use.dimred is explicitly specified by the user or "how" == "density".
#' For mbkmeans, k is the number of clusters and "HARMONY" is used if found.
#' For densityClust, the default representation to cluster upon is "UMAP". 
#' 
#' It is not yet clear whether an approach like that taken in `Rphenoannoy` 
#' would allow for faster and therefore bootstrap-able clustering, and this is 
#' something that will be of interest going forwards. In some cases, using 
#' `rank` or `jaccard` distances with Louvain clustering produces overly fine-
#' grained clusters, to the extent that something like SPADE may be needed to
#' reconcile the results. The 'phenograph' algorithm is itself just Louvain
#' clustering on 'jaccard' distance ("type"), something that 'rank' improves.
#' 
#' The default number of PCAs or Harmony components "k" is set equal to that 
#' which Seurat uses by default, not for any particularly good reason. Note
#' that `k` means something different (cluster number K) for minibatch K-means,
#' but we leave the default the same (because YOLO, and you can always merge).
#' (Which is to say, also not for any particularly good reason.)
#' 
#' It would probably be a good idea to assess stability with `bluster` here.
#' 
#' @param txis        a SingleCellExperiment
#' @param k           number of nearest neighbors, or clusters for mbkmeans (20)
#' @param type        what type of SNNGraph to construct ("rank")
#' @param use         which reducedDimension embedding to use (autodetermine)
#' @param how         cluster via "louvain" (default), "density", or "mbkmeans"?
#' @param ret         return "clust" assignments (default) or "sce"?
#' @param ...         optional arguments to pass to clustering function(s) 
#' 
#' @return            Cluster labels (if ret=="clust") or a now-labeled "sce"
#' 
#' @references \
#' Louvain clustering: \
#' Blondel VD, Guillaume JL, Lambiotte R, Lefebvre E (2008). \
#' Fast unfolding of communities in large networks. \
#' Journal of Statistical Mechanics: Theory and Experiment. 2008 (10): P10008. \
#' doi:10.1088/1742-5468/2008/10/P10008 \
#' \
#' Density clustering: \
#' Rodriguez A, Laio A. (2014). \
#' Clustering by fast search and find of density peaks. \
#' Science, 344(6191), 1492-1496. \
#' doi:10.1126/science.1242072 \
#' \ 
#' Minibatch K-means: \
#' Sculley D. (2010). \
#' Web-Scale K-Means Clustering. \
#' WWW 2010, Raleigh, North Carolina, USA, ACM 978-1-60558-799-8/10/04. \
#' https://doi.org/10.1145/1772690.1772862 \
#' \
#' Hicks SC, Liu R, Ni Y, Purdom E, Risso D. (2020). \
#' mbkmeans: fast clustering for single cell data using mini-batch k-means. \
#' bioRxiv 2020.05.27.119438 \
#' doi: https://doi.org/10.1101/2020.05.27.119438
#'
#' @importFrom  scran         buildSNNGraph
#' @importFrom  igraph        cluster_louvain
#' @importFrom  densityClust  densityClust
#' @importFrom  densityClust  findClusters
#' @import      mbkmeans
#'
#' @export 
cluster_velo_txis <- function(txis, k=20, type="rank", use=NULL, how=c("louvain", "density", "mbkmeans"), ret=c("clust", "sce"), ...) { 

  how <- match.arg(how)
  ret <- match.arg(ret)
  if (is.null(use) & how %in% c("louvain", "mbkmeans")) {
    use <- ifelse("HARMONY" %in% reducedDimNames(txis), "HARMONY", "PCA")
    if (use == "PCA" & !"PCA" %in% reducedDimNames(txis)) txis <- runPCA(txis)
  } else if (how == "density") {
    use <- "UMAP" 
    src <- ifelse("HARMONY" %in% reducedDimNames(txis), "HARMONY", "PCA")
    if (src == "PCA" & !"PCA" %in% reducedDimNames(txis)) txis <- runPCA(txis)
    if (!"UMAP" %in% reducedDimNames(txis)) txis <- runUMAP(txis, dimred=src)
  }
  if (!use %in% reducedDimNames(txis)) stop(use," not in reducedDimNames(txis)")
  message("Note: at present, we ignore unspliced information while clustering.")

  # default resembles Seurat
  clusters <- switch(how,
                     louvain=.clust_louvain(txis, k=k, type=type, use=use, ...),
                     mbkmeans=.clust_mbkmeans(txis, k=k, use=use, ...),
                     density=.clust_density(txis, use=use, ...))

  # switch to using colLabels pervasively downstream
  colLabels(txis) <- factor(clusters[colnames(txis)])
  return(switch(ret, sce=txis, clust=colLabels(txis)))

}


# helper fn 
.clust_louvain <- function(txis, k=20, type="rank", use="HARMONY", ...){
  
  graf <- scran::buildSNNGraph(txis, k=k, use.dimred=use, type=type)
  louvain <- igraph::cluster_louvain(graf, ...) # nb. jaccard may be better
  clusters <- paste0("louvain", louvain$membership) 
  names(clusters) <- colnames(txis)
  return(clusters)

}


# helper fn
.clust_density <- function(txis, gaussian=TRUE, use="UMAP", ...) { 

  dp <- densityClust(dist(reducedDim(txis, use)), gaussian=gaussian)
  if ((!exists("rho") | is.null(rho) | !exists("delta") | is.null(delta))) {
    message("Use 0.95 of the delta and 0.95 of the rho as the cutoff")
    message("for assigning density peaks and clusters")
    rho <- quantile(dp$rho, probs = 0.95)
    delta <- quantile(dp$delta, probs = 0.95)
  }
  dc <- findClusters(dp, rho=rho, delta=delta)
  clusters <- paste0("density", dc$clusters)
  names(clusters) <- colnames(txis) 
  return(clusters)

}


# helper fn 
.clust_mbkmeans <- function(txis, k=20, use="HARMONY", ...){
 
  mbkm <- mbkmeans::mbkmeans(txis, clusters=k, reduceMethod=use, ...)
  clusters <- paste0("mbkmeans", mbkm$Clusters)
  names(clusters) <- colnames(txis)
  return(clusters) 

}
