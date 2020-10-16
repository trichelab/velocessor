#' Generate interactive 3D UMAP plot with cones for velocity. Work in progress.
#' 
#' The opacity of the cones is manipulated so that dot colors can be seen. 
#' 
#' @param txis        a SingleCellExperiment with metadata(txis)$scVelo
#' @param ...         optional arguments to pass to plotly::plot_ly()
#' 
#' @return            a plotly plot
#' 
#' @import plotly
#' @import Polychrome
#' 
#' @export 
plot_velo <- function(txis, ...) {
 
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
  grid.df$CLUSTER <- as.factor(colData(txis)[rownames(grid.df), "cluster"])
  grid.df$SAMPLE <- as.factor(colData(txis)[rownames(grid.df), "sample"])

  seed <- c("#ff0000", "#00ff00", "#0000ff")
  pal <- createPalette(nlevels(grid.df$CLUSTER), seed, prefix="cluster")

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
               color = ~CLUSTER,
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
