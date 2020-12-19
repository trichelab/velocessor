#' Refactoring plot_velo into multiple functions moving forwards 
#' 
#' @param txis        a SingleCellExperiment
#' @param embed       which reducedDims() element to use for plotting? (UMAP)
#' @param colr        column(s) to use for colors (default: colLabels(txis))
#' @param p           an existing plotly object to add the trace to (NULL)
#' @param ...         optional arguments to pass to plotly::plot_ly()
#' 
#' @return            a plotly trace
#' 
#' @import plotly
#' @import viridis
#' @import Polychrome
#' @import viridisLite
#' 
add_static_trace <- function(txis, embed="UMAP", colr=NULL, p=NULL, ...) {

  stop("Placeholder function at this time")

  if (is.null(colr)) {
    message("Using colLabels for default color assignment") 
    txis$col_label <- colLabels(txis)
    colr <- "col_label"
    # don't forget to delete it 
  }

  # now actually get to work 
  message("Refactor tihs")
  dat <- velociraptor::gridVectors(reducedDims(txis)[[embed]], embedded)
  rownames(dat) <- colnames(txis)[as.integer(rownames(dat))]
  swap <- c(".1"="X", ".2"="Y", ".3"="Z", "start"="")
  for (s in names(swap)) names(dat) <- sub(s, swap[s], fixed=TRUE, names(dat))
  
  # for cones
  starts <- c("U"="X", "V"="Y", "W"="Z")
  ends <- paste0("end", starts)
  names(ends) <- names(starts)
  for (n in names(ends)) dat[, n] <- dat[, ends[n]] - dat[, starts[n]]

  # for colors and point labels 
  dat$SAMPLE <- colData(txis)[rownames(dat), "sample"]
  dat$COLORING <- colData(txis)[rownames(dat), colr]
  if (colr != "velocity_pseudotime") dat$COLORING <- factor(dat$COLORING)

  # for colors and point labels 
  if (colr == "velocity_pseudotime") {
    dat$COLORING <- round(dat$COLORING * 100) # percent
    dat$LABEL <- paste0(dat$SAMPLE, ": ", dat$COLORING, "%") 
    pal <- inferno(100) # colorful continuous scale
  } else {
    dat[, colr] <- colData(txis)[rownames(dat), colr] 
    dat[, "LABEL"] <- paste0(dat$SAMPLE, ": ", dat[, colr])
    seed <- c("#ff0000", "#00ff00", "#0000ff") # R, G, B 
    pal <- createPalette(nlevels(dat$COLORING), seed, prefix="color")
    names(pal) <- levels(dat$COLORING)
  }
  dat$SAMPLE <- colData(txis)[rownames(dat), "sample"]

  # for axis labels and formatting
  axes <- lapply(list(xaxis=1, yaxis=2, zaxis=3), .clean_axis, embed=embed)

  # FIXME1: make it so the cones (velocity) and points (cells) can be turned
  #         on/off and/or relabeled/recolored arbitrarily during interaction.
  #
  message("FIXME: add velocity and cells as two separate traces.") 

  # FIXME2: provide a way for users to squash a 3D UMAP down into a 2D UMAP,
  #         ideally resulting in something like the output of plot_features().
  #
  message("FIXME: add 3D-to-2D squashing.") 

  # FIXME3: provide a way for users to switch subsetting/coloring on the fly.
  #
  message("FIXME: add dynamic subsetting for cells (and velocity?).")

  # build the plot
  # add the velocity trace, unless requested to ditch it
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
               text = ~LABEL,
               color = ~COLORING,
               hoverinfo = "text",
               showscale = FALSE,
               colorscale = "Greys",
               alpha_stroke = I(0.5),
               alpha = I(0.75),
               ...)
  p <- config(p, displayModeBar = FALSE)

  # if only one scheme: 
  if (length(colr) < 2) { 
    layout(add_markers(p), scene = axes)
  } else { 
    # add switch button 
    message("FIXME: add swapping buttons!")
    layout(add_markers(p), scene = axes)
  }

}
