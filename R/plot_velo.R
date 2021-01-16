#' Generate interactive 3D UMAP plot with cones for velocity. Work in progress.
#' 
#' The opacity of the cones is manipulated so that dot colors are visible.
#' The provided SingleCellExperiment `txis` MUST have an element `scVelo` in 
#' its metadata, and that element (a SingleCellExperiment) needs to have a 
#' colData column "velocity_pseudotime" unless an argument 'colr' is provided.
#' (Hence the references internally to "COLORING" rather than "COLOR") 
#' 
#' @param txis        a SingleCellExperiment with metadata(txis)$scVelo
#' @param embed       which reducedDims() element to use for plotting? (UMAP)
#' @param replace     replace any existing metadata element `embedded`? (FALSE)
#' @param colr        column(s) to use for colors (default: velocity_pseudotime)
#' @param flip        flip the velocity arrows (ghetto CellRank)? (FALSE)
#' @param ...         optional arguments to pass to plotly::plot_ly()
#' 
#' @return            a plotly plot
#' 
#' @import plotly
#' @import viridis
#' @import Polychrome
#' @import viridisLite
#' @import velociraptor
#' 
#' @export 
plot_velo <- function(txis, embed="UMAP", replace=FALSE, colr="velocity_pseudotime", flip=FALSE, ...) {

  # sanity checking 
  if (!.velo_ok(txis)) stop("Velocity information missing. Cannot proceed.") 
  if (!"velocity_pseudotime" %in% names(colData(txis))) txis <- .add_vpt(txis)
  if (!.colrs_ok(colr, txis)) stop("Missing colData columns. Cannot proceed.") 
  embedded <- .get_embedding(txis, embed=embed, replace=replace) 

  # now actually get to work 
  stopifnot(nrow(reducedDims(txis)[[embed]]) == nrow(embedded))
  stopifnot(ncol(reducedDims(txis)[[embed]]) >= ncol(embedded))
  cols <- seq_len(ncol(embedded))
  dat <- velociraptor::gridVectors(reducedDims(txis)[[embed]][,cols], embedded)
  rownames(dat) <- colnames(txis)[as.integer(rownames(dat))]
  swap <- c(".1"="X", ".2"="Y", ".3"="Z", "start"="")
  for (s in names(swap)) names(dat) <- sub(s, swap[s], fixed=TRUE, names(dat))
  
  # for cones (refactor this out) 
  starts <- c("U"="X", "V"="Y", "W"="Z")
  ends <- paste0("end", starts)
  names(ends) <- names(starts)
  for (n in names(ends)) dat[, n] <- dat[, ends[n]] - dat[, starts[n]]

  # ghetto cellrank
  if (flip) for (n in names(ends)) dat[, n] <- dat[, starts[n]] - dat[, ends[n]]

  # for colors and point labels 
  dat$SAMPLE <- colData(txis)[rownames(dat), "sample"]
  dat$COLORING <- colData(txis)[rownames(dat), colr]
  if (colr != "velocity_pseudotime") dat$COLORING <- factor(dat$COLORING)

  # for colors and point labels 
  if (colr == "velocity_pseudotime") {
    dat$COLORING <- round(dat$COLORING * 100) # percent
    if (flip) dat$COLORING <- 100 - dat$COLORING # reverse it
    dat$LABEL <- paste0(dat$SAMPLE, ": ", dat$COLORING, "%") 
    pal <- inferno(100) # colorful continuous scale
  } else {
    dat[, colr] <- colData(txis)[rownames(dat), colr] 
    dat[, "LABEL"] <- paste0(dat$SAMPLE, ": ", dat[, colr])
    seed <- c("#ff0000", "#0000ff", "#00ff00") # R, B, G 
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


# helper fn 
.velo_ok <- function(txis) { 

  ret <- TRUE 
  if (!"scVelo" %in% names(metadata(txis))) {
    message("Cannot find scVelo output in metadata(txis).")
    ret <- FALSE
  } 
  if (!"velocity_pseudotime" %in% names(colData(metadata(txis)$scVelo))) {
    message("metadata(txis)$scVelo does not appear to be velociraptor output.") 
    ret <- FALSE
  }
  return(ret)

} 


# helper fn 
.colrs_ok <- function(colr, txis) { 

  ret <- TRUE 
  if (!all(colr %in% names(colData(txis)))) {
    for (missed in setdiff(colr, names(colData(txis)))) {
      message("Cannot find ", missed, "in names(colData(txis)).")
    } 
    ret <- FALSE
  } 
  return(ret)

} 


# helper fn 
.add_vpt <- function(txis) { 

  scVelo <- metadata(txis)$scVelo
  colData(txis)[,"velocity_pseudotime"] <- scVelo$velocity_pseudotime
  return(txis) 

} 


# helper fn 
.get_embedding <- function(txis, embed="UMAP", replace=FALSE) { 

  if ("embedded" %in% names(metadata(txis)) & !replace) { 
    return(metadata(txis)$embedded) 
  } else {
    rd <- reducedDims(txis)[[embed]] 
    return(velociraptor::embedVelocity(rd, metadata(txis)$scVelo))
  }

}


# helper fn 
.make_color_menu <- function(txis, colr, ...) { 


} 


# helper fn 
.make_velo_menu <- function(...) { 

  velo_size <- list(
    type = "buttons",
    direction = "right",
    xanchor = "center",
    yanchor = "top",
    pad = list('r'= 0, 't'= 10, 'b' = 10),
    x = 0.5,
    y = 1.27,
  buttons = list(

    list(method = "restyle",
         args = list("type", "heatmap"),
         label = "Heatmap"),

    list(method = "restyle",
         args = list("type", "contour"),
         label = "Contour"),

    list(method = "restyle",
         args = list("type", "surface"),
         label = "Surface")
  ))

} 


# helper fn 
.make_annotations <- function(updatemenus, ...) { 

}


# helper fn
.clean_axis <- function(i, embed="axis") {
   
  list(title=paste0(embed, i), 
       showticklabels = FALSE,
       zeroline = FALSE,
       showline = FALSE,
       showgrid = FALSE)

}


