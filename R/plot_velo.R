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
  if (!.velo_ok(txis)) stop("Velocity information missing. Cannot proceed.") 
  if (!"velocity_pseudotime" %in% names(colData(txis))) txis <- .add_vpt(txis)
  if (!.colrs_ok(colr, txis)) stop("Missing colData columns. Cannot proceed.") 
  embedded <- .get_embedding(txis, embed=embed, replace=replace) 

  # now actually get to work 
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
  dat$COLORING <- colData(txis)[rownames(dat), colr]

  # for colors and point labels 
  if (colr == "velocity_pseudotime") {
    dat$COLORING <- round(dat$COLORING * 100) # percent
    dat$LABEL <- paste0(dat$SAMPLE, ": ", dat$COLORING, "%") 
    pal <- inferno(100) # colorful continuous scale
  } else {
    dat$COLORING <- dat$COLORING
    seed <- c("#ff0000", "#00ff00", "#0000ff") # R, G, B 
    pal <- createPalette(nlevels(dat$COLORING), seed, prefix="color")
    dat$LABEL <- paste0(dat$SAMPLE, ": ", dat[, colr])
  }
  dat$SAMPLE <- colData(txis)[rownames(dat), "sample"]

  # for axis labels and formatting
  axes <- lapply(list(xaxis=1, yaxis=2, zaxis=3), .clean_axis, embed=embed)

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
               text = ~LABEL,
               color = ~COLORING,
               hoverinfo = "text",
               showscale = FALSE,
               showlegend = FALSE,
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


