#' given a `uwot` model in metadata(x)$uwot, project y onto those coordinates
#' 
#' @param x     a SingleCellExperiment where !is.null(metadata(x)$uwot)
#' @param y     a SingleCellExperiment with samples to map onto x's model
#' @param rdn   reducedDimName to use ("HARMONY" if present, "PCA" otherwise)
#' @param ...   additional parameters for uwot::umap_transform
#' 
#' @return      embedding coordinates for y onto x
#'
#' @seealso     uwot::umap_transform
#' 
#' @import      uwot
#' 
#' @export
project_onto_umap <- function(x, y, rdn=NULL, ...) {

  rdns <- reducedDimNames(y)
  if (!is.null(rdn)) stopifnot(rdn %in% rdns)
  if (is.null(rdn)) rdn <- ifelse("HARMONY" %in% rdns, "HARMONY", "PCA")
  if (rdn == "PCA" & !"PCA" %in% rdns) y <- runPCA(y)
  if(!"uwot" %in% names(metadata(x))) x <- add_umap_embedding(x)
  uwot::umap_transform(X=reducedDim(y, rdn), model=metadata(x)$uwot, ...)

}
