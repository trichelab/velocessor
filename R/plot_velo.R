#' Generate interactive 3D UMAP plot with cones for velocity. Work in progress.
#' 
#' The opacity of the cones is manipulated so that dot colors can be seen. 
#' 
#' @param txis        a SingleCellExperiment with metadata(txis)$scVelo
#' @param colr        column to use for colors (else just velocity_pseudotime)
#' @param ...         optional arguments to pass to plotly::plot_ly()
#' 
#' @return            a plotly plot
#' 
#' @import plotly
#' @import viridis
#' @import Polychrome
#' 
#' @export 
plot_velo <- function(txis, colr="velocity_pseudotime", ...) {
 
  if (!"scVelo" %in% names(metadata(txis))) stop("Cannot find scVelo output!")
  if (!"velocity_pseudotime" %in% names(colData(txis))) {
    colData(txis)[, "velocity_pseudotime"] <- 
      metadata(txis)$scVelo$velocity_pseudotime
  }
  if (!colr %in% names(colData(txis))) stop("Cannot find txis$", colr, "!")
  
  if ("embedded" %in% names(metadata(txis))) { 
    embedded <- metadata(txis)$embedded
  } else {
    embedded <- embedVelocity(reducedDim(txis, "UMAP"), metadata(txis)$scVelo)
  } 

  grid.df <- gridVectors(reducedDim(txis, "UMAP"), embedded)
  rownames(grid.df) <- colnames(txis)[as.numeric(rownames(grid.df))]
  colnames(grid.df) <- sub("\\.1", "X",
                         sub("\\.2", "Y",
                             sub("\\.3", "Z",
                                 colnames(grid.df))))
  
  # for cones
  grid.df$U <- grid.df$endX - grid.df$startX
  grid.df$V <- grid.df$endY - grid.df$startY
  grid.df$W <- grid.df$endZ - grid.df$startZ

  # for colors
  grid.df$COLOR <- colData(txis)[rownames(grid.df), colr]
  grid.df$SAMPLE <- as.factor(colData(txis)[rownames(grid.df), "sample"])

  seed <- c("#ff0000", "#00ff00", "#0000ff") # R, G, B 
  if (colr == "velocity_pseudotime") {
    grid.df$COLOR <- paste0(round(grid.df$COLOR * 20) * 5, "%")
    pal <- viridis(20)
  } else {
    pal <- createPalette(nlevels(as.factor(grid.df$COLOR)), seed, prefix="colr")
  }
  
  p <- plot_ly(grid.df,
               type = "cone",
               showscale = FALSE,
               colors = pal, 
               opacity = 0.5,
               x = ~startX,
               y = ~startY,
               z = ~startZ,
               u = ~U,
               v = ~V,
               w = ~W,
               text = ~SAMPLE,
               color = ~COLOR,
               colorscale = "Greys",
               alpha_stroke = I(0.5),
               alpha = I(0.5),
               ...)
  p <- add_markers(p)
  p <- layout(p,
              scene = list(xaxis = list(title = 'UMAP1'),
              yaxis = list(title = 'UMAP2'),
              zaxis = list(title = 'UMAP3')))
  return(p)

}
