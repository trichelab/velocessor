#' another helper that grew into a real live function (albeit a pathetic one)
#' 
#' Plot densities of competing mixture fits for diagnostics.
#' 
#' @seealso fit_mixtures
#' 
#' @param fits  the fits from fit_mixtures
#' 
#' @import mclust
#' 
#' @export
plot_densities <- function(fits) { 

  par(mfrow=c(1,2))
  plotDensityMclust1(fits$raw, data=fits$raw$data, xlab="cells")
  plotDensityMclust1(fits$log, data=fits$log$data, xlab="log1p(cells)")

} 
