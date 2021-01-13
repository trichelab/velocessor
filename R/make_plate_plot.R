#' plot outcomes of 96-, 384-, or 1536-well plate sequencing runs
#' 
#' @param txis    a SingleCellExperiment with either $well or $row and $column
#' @param column  which column to use for the data to plot? ($NumGenesExpressed)
#' @param ...     additional arguments for platetools::raw_map
#'
#' @return        a ggplot object 
#' 
#' @examples
#' 
#'  # helper fn
#'  revglue <- function(x) gsub(" ", "", paste(rev(x), collapse=""))
#'  wells <- toupper(apply(expand.grid(1:24, letters[1:16]), 1, revglue))
#' 
#'  # mock up an SCE
#'  maxgenes <- 100
#'  library(scDblFinder) 
#'  sce <- mockDoubletSCE(ncells=c(184, 150, 50), ngenes=maxgenes)[, 1:384]
#'  sce$NumGenesExpressed <- colSums(counts(sce) > 0)
#'  sce$well <- wells
#' 
#'  # plot it with $well
#'  make_plate_plot(sce)
#' 
#'  # plot it using colnames
#'  colnames(sce) <- sce$well
#'  sce$well <- NULL # drop
#'  # ensure case is ignored for now
#'  changecase <- sample(seq_len(ncol(sce)), round(ncol(sce)/2))
#'  colnames(sce)[changecase] <- toupper(colnames(sce)[changecase])
#'  colnames(sce)[-changecase] <- tolower(colnames(sce)[-changecase])
#'  make_plate_plot(sce)
#' 
#' @seealso platetools::raw_map 
#' @seealso platetools::b_map 
#' 
#' @import platetools
#' 
#' @export 
make_plate_plot <- function(txis, column="NumGenesExpressed", ...) { 

  if ("well" %in% names(colData(txis))) {
    well <- .tidy_well_names(txis$well, pad=TRUE)
  } else if (all(c("row","column") %in% names(colData(txis)))) {
    well <- .tidy_well_names(paste0(txis$row, txis$column), pad=TRUE) 
  } else if (all(.colnames_are_wells(txis))) {
    well <- colnames(txis)
  } else { 
    stop("Don't know how to find the well or row or column for runs.")
  }

  if ("plate" %in% names(colData(txis)) & length(unique(txis$plate)) > 1) { 
    stop("Plate-grid plots are not supported yet. Use platetools directly.")
  }

  # choose the number of wells in the obvious way so that I don't have to think
  plt <- ifelse(length(well) < 97, 96, ifelse(length(well) < 384, 384, 1536))
  p <- raw_map(data=colData(txis)[, column], well=well, plate=plt, ...)
  p + legend_title(gsub("([a-z])([A-Z])", "\\1\n\\2", column))

}


# helper fn
.tidy_well_names <- function(well, pad=TRUE)  {
  
  letter <- toupper(substr(well, 1, 1))
  number <- as.integer(substr(well, 2, max(nchar(well))))
  fmt <- paste0("%s%", ifelse(pad, "0.", ""), max(nchar(number)), "i")
  gsub(" ", "", sprintf(fmt, letter, number)) 

}


# helper fn
.colnames_are_wells <- function(txis) {

  grepl("^[a-z]+[0-9]+", ignore.case=TRUE, colnames(txis))

}
