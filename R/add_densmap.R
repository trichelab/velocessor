#' run densMAP and store the result in reducedDim(x, "DENSMAP")
#' 
#' @param x       a SingleCellExperiment or a slice of one
#' @param nc      how many components to keep for the densmap embedding (6)
#' @param rdn     reducedDimName to use ("HARMONY" if present, "PCA" otherwise)
#' @param ...     other parameters to pass to `densvis::densmap`
#' 
#' @return        a SingleCellExperiment with reducedDim(x, "DENSMAP")
#'
#' @seealso       densmap
#'
#' @import        densvis
#' 
#' @export
add_densmap <- function(x, nc=3L, rdn=NULL, ...) { 

  rdns <- reducedDimNames(x)
  if (is.null(rdn)) rdn <- ifelse("HARMONY" %in% rdns, "HARMONY", "PCA")
  if (rdn == "PCA" & !"PCA" %in% rdns) x <- runPCA(x)
  rdns <- reducedDimNames(x)
  stopifnot(rdn %in% rdns)

  reducedDim(x, "DENSMAP") <- densvis::densmap(x=reducedDim(x, rdn), 
                                               n_components=nc, ...) 
  return(x)

}
