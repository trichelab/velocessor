#' grab ADT/HTO matrices from an h5 mRNA+ADT file
#' 
#' @param   h5_fname  the h5 filename (required)
#' @param   base      the base name group in the h5 file ("matrix") 
#' @param   stub      a column name stub, if using one (NULL)
#' 
#' @return            a matrix (i.e., as.matrix(subset of a TENxMatrix))
#' 
#' @import  HDF5Array
#' @import  rhdf5
#' 
#' @export
adts_from_hdf5 <- function(h5_fname, base="matrix", stub=NULL) { 
  
  feats <- paste0(base, "/features") 
  message("Loading data matrix from ", h5_fname, "...")
  dat <- TENxMatrix(h5_fname, group=base)
  message("Loaded ", ncol(dat), " cells.")
  message("Loading row names...")
  rnames <- as.character(h5read(h5_fname, feats)["id"][[1]]) 
  message("Found ", nrow(dat), " features.")
  stopifnot(nrow(dat) == length(rnames))
  message("Adding row names...")
  rownames(dat) <- rnames
  message("Subsetting ADTs...")
  ADTs <- adt_seqs_from_hdf5(h5_fname, base=base, fa_fname=NULL) # don't write
  message("Found ", nrow(ADTs), " ADTs.")
  res <- as.matrix(dat[ADTs$id, ])
  if (!is.null(stub)) res <- .tidy_colnames(res, stub=stub)
  return(res) 

}


# helper fn
.tidy_colnames <- function(mat, stub, sep="_") { 
  colnames(mat) <- paste(stub, gsub("\\-[123]$", "", colnames(mat)), sep=sep)
  return(mat) 
}
