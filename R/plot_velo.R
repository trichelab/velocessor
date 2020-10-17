#' Generate interactive 3D UMAP plot with cones for velocity. Work in progress.
#' 
#' The opacity of the cones is manipulated so that dot colors are visible.
#' The provided SingleCellExperiment `txis` MUST have an element `scVelo` in 
#' its metadata, and it would be a good idea for 
#' 
#' @param txis        a SingleCellExperiment with metadata(txis)$scVelo
#' @param embed       which reducedDims() element to use for plotting? (UMAP)
#' @param replace     replace any existing metadata element `embedded`? (FALSE)
#' @param colr        column to use for colors (else just velocity_pseudotime)
#' @param ...         optional arguments to pass to plotly::plot_ly()
#' 
#' @return            a plotly plot
#' 
#' @import plotly
#' @import viridis
#' @import Polychrome
#' @import velociraptor
#' 
#' @export 
plot_velo <- function(txis, embed="UMAP", replace=FALSE, colr="velocity_pseudotime", ...) {

  # sanity checking  
  if (!"scVelo" %in% names(metadata(txis))) stop("Cannot find scVelo output!")
  if (!"velocity_pseudotime" %in% names(metadata(txis)$scVelo)) {
    stop("You need to have actual scVelo output in your SCE metadata.")
  }
  if (!"velocity_pseudotime" %in% names(colData(txis))) {
    colData(txis)[, "velocity_pseudotime"] <- 
      metadata(txis)$scVelo$velocity_pseudotime
  }
  if (!colr %in% names(colData(txis))) {
    stop("Cannot find ", colr, "in names(colData(txis)). Exiting.")
  }
  if ("embedded" %in% names(metadata(txis)) & !replace) { 
    embedded <- metadata(txis)$embedded
  } else {
    embedded <- velociraptor::embedVelocity(reducedDim(txis, embed), 
                                            metadata(txis)$scVelo)
  }

  # now actually get to work 
  dat <- velociraptor::gridVectors(reducedDim(txis, embed), embedded)
  rownames(dat) <- colnames(txis)[as.numeric(rownames(dat))]
  swap <- c(".1"="X", ".2"="Y", ".3"="Z", "start"="")
  for (s in names(swap)) names(dat) <- sub(s, swap[s], fixed=TRUE, names(dat))
  
  # for cones
  starts <- c("U"="X", "V"="Y", "W"="Z")
  ends <- paste0("end", starts)
  for (n in names(ends)) dat[, n] <- dat[, ends[n]] - dat[, starts[n]]

  # for colors and point labels 
  dat$COLOR <- colData(txis)[rownames(dat), colr]
  if (colr == "velocity_pseudotime") {
    dat$COLOR <- round(dat$COLOR * 100) # percent
    pal <- inferno(100) # colorful continuous scale
  } else {
    dat$COLOR <- as.factor(dat$COLOR)
    seed <- c("#ff0000", "#00ff00", "#0000ff") # R, G, B 
    pal <- createPalette(nlevels(dat$COLOR), seed, prefix="color")
  }
  dat$SAMPLE <- as.factor(colData(txis)[rownames(dat), "sample"])

  # for axis labels and formatting
  axes <- lapply(list(xaxis=1, yaxis=2, zaxis=3), 
                 function(i) list(title=paste0(embed, i), 
                                  showticklabels = FALSE,
                                  showline = FALSE,
                                  showgrid = FALSE))

  # build the plot
  p <- plot_ly(dat,
               x = ~X,
               y = ~Y,
               z = ~Z,
               u = ~U,
               v = ~V,
               w = ~W,
               colors = pal, 
               type = "cone",
               opacity = 0.5,
               text = ~SAMPLE,
               color = ~COLOR,
               hoverinfo = "text",
               showscale = FALSE,
               colorscale = "Greys",
               alpha_stroke = I(0.5),
               alpha = I(0.5),
               ...)
  p <- config(p, displayModeBar = FALSE)
  layout(add_markers(p), scene = axes)

}
