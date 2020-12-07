#' check that an SCE will be suitable for muscat, scvelo, plotting, etc.
#' 
#' @param sce     the SingleCellExperiment
#' @param syms    check that genes/transcripts have unique symbols? (TRUE) 
#' @param velo    check for plot_velo compatibility? (TRUE) 
#' @param muscat  check for muscat compatibility (TRUE) 
#' @param ...     other parameters, currently ignored 
#' 
#' @return        logical: whether all the checks passed
#' 
#' @export
check_sce <- function(sce, syms=TRUE, velo=TRUE, muscat=TRUE, ...) {

  why <- c()

  if (!is(sce, "SingleCellExperiment")) {
    message("Can't check something that isn't a SingleCellExperiment!") 
    return(FALSE) 
  }

  if (syms) { 
    if (!"symbol" %in% names(mcols(sce))) why <- c(why, "syms: no $symbol")
    else if (any(duplicated(rowData(sce)$symbol))) why <- c(why, "syms: dupes")
  }  
  if (velo) { 
    if (!"scVelo" %in% names(metadata(sce))) why <- c(why, "velo: no scVelo")
  }
  if (muscat) { 
    if (is.null(colLabels(sce))) why <- c(why, "muscat: no colLabels")
    if (!"sample" %in% names(colData(sce))) why <- c(why, "muscat: no $sample")
    if (!"group" %in% names(colData(sce))) why <- c(why, "muscat: no $group")
  }

  if (length(why) > 0) {
    message("check_sce failed due to the following errors:")
    for (i in seq_along(why)) message(why[i])
    return(FALSE) 
  }

  # no problems:  
  return(TRUE) 

}
