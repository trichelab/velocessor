#' some equivalents to Seurat standards, with or without velocity vectors 
#' 
#' Mostly for joint density plots via Nebulosa with velocity overlays.
#' 
#' @param txis      a SingleCellExperiment
#' @param features  names of features (rows) to plot
#' @param joint     for multiple features, plot joint density? (TRUE)
#' @param ...       other arguments for Nebulosa::plot_density
#'
#' @import SingleCellExperiment
#' @importFrom Nebulosa plot_density
#'
#' @export
plot_features <- function(txis, features, joint=TRUE, ...) {

  Nebulosa::plot_density(txis, features, joint=joint, ...)

}
