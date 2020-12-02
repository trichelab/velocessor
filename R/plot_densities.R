#' another helper that grew into a real live function (albeit a pathetic one)
#' 
#' Plot densities of competing mixture fits for diagnostics.
#' 
#' @seealso fit_mixtures
#' 
#' @param fits  the fits from fit_mixtures
#' @param what  what it is that was fitted
#' 
#' @import mclust
#' 
#' @export
plot_densities <- function(fits, what="cells") { 

  if (length(fits) > 1) par(mfrow=c(1, length(fits)))

  # the default
  plotDensityMclust1(fits$raw, data=fits$raw$data, xlab=what)

  # generalize?
  if ("log" %in% names(fits)) {
    plotDensityMclust1(fits$log, data=fits$log$data, 
                       xlab=paste0("log1p(", what, ")"))
  }

} 
