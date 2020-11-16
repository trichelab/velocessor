#' Plot a heatmap of cells x samples that fit certain sampling criteria. 
#' 
#' See the documentation for find_eligible_cells regarding these criteria. 
#' This can be a useful diagnostic when setting `maxcells` for label_cells.
#' 
#' @param txis      SingleCellExperiment where !is.null(colLabels(txis))
#' @param mincells  min cells per cluster per sample (see Details) (10)
#' @param logcells  use log1p(cells) instead of raw cell numbers? (FALSE) 
#' @param ...       additional arguments to pass to ComplexHeatmap::Heatmap
#' 
#' @seealso label_cells
#' @seealso downsample_txis
#' 
#' @importFrom ComplexHeatmap Heatmap
#' 
#' @export
plot_eligible_cells <- function(txis, mincells=10, logcells=FALSE, ...) { 

  res <- find_eligible_cells(txis=txis, mincells=mincells)
  eligible <- attr(res, "eligible")
  name <- "cells"

  if (logcells == TRUE) { 
    eligible <- log10(eligible + 1) 
    name <- "log10(cells)"
  }

  if (nrow(eligible) > ncol(eligible)) {
    eligible <- t(eligible) 
    row_title <- "source"
    column_title <- "cluster"
  } else { 
    row_title <- "cluster"
    column_title <- "source"
  } 

  Heatmap(as(eligible, "matrix"), name=name,
          row_title=row_title, column_title=column_title,
          show_row_dend=FALSE, show_column_dend=FALSE) 

}
