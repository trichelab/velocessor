#' run UMAP with `uwot` and store the embedding model in metadata(x)$uwot
#' 
#' @param x     a SingleCellExperiment or a slice of one
#' @param nc    how many components to keep for `uwot`'s UMAP embedding (3)
#' @param ret   return the modified "sce" (default) or "uwot" model? ("sce")
#' @param rdn   reducedDimName to use ("HARMONY" if present, "PCA" otherwise)
#' @param ...   additional parameters for uwot::umap (e.g. `y` for SDR)
#' 
#' @return      depending on `ret`, an embellished sce or a uwot model
#'
#' @seealso     uwot::umap
#' 
#' @import      uwot
#' 
#' @export
add_umap_embedding <- function(x, nc=3L, ret=c("sce","uwot"), rdn=NULL, ...) {

  ret <- match.arg(ret)
  rdns <- reducedDimNames(x)
  if (is.null(rdn)) rdn <- ifelse("HARMONY" %in% rdns, "HARMONY", "PCA")
  if (rdn == "PCA" & !"PCA" %in% rdns) x <- runPCA(x)
  rdns <- reducedDimNames(x)
  stopifnot(rdn %in% rdns)

  metadata(x)$uwot <- uwot::umap(X=reducedDim(x, rdn), n_components=nc, ..., 
                                 metric="cosine", ret_model=TRUE, verbose=TRUE)
  reducedDim(x, "UMAP") <- metadata(x)$uwot$embedding
  return(switch(ret, sce=x, uwot=metadata(x)$uwot))

}
