#' drop cells with an unspliced size factor of 0 (which will crash scvelo). 
#' 
#' Also attempts to use scran::cleanSizeFactors to recover dubious cells, 
#' if NumGenesExpressed is a column in colData(txis).  This may or may not be 
#'
#' @param txis          a SingleCellExperiment
#' @param clean         fix negative size factors with cleanSizeFactors? (TRUE)
#' 
#' @return              a SingleCellExperiment with any useless cells omitted
#'
#' @importFrom scuttle  librarySizeFactors
#' @importFrom scran    cleanSizeFactors
#' 
#' @export
drop_uninformative_cells <- function(txis, clean=TRUE) { 

  an <- c(spl="spliced", unspl="unspliced")
  sf <- do.call(cbind,
                lapply(an, function(x) 
                       scuttle::librarySizeFactors(txis, assay.type=x)))
  kept <- names(which(apply(sf, 1, function(x) all(x > 0))))
  suspicious <- ncol(txis) - length(kept)
  if (suspicious > 0) message(suspicious, " suspicious cells.")
  
  if (clean & "NumGenesExpressed" %in% names(colData(txis))) {
    sf$num <- txis$NumGenesExpressed
    message("Cleaning up negative size factors...") 
    sf <- do.call(cbind,
                  lapply(an, function(x) 
                         scran::cleanSizeFactors(sf[,an], sf[,"num"])))
    sf$num <- NULL
    kept <- names(which(apply(sf, 1, function(x) all(x > 0))))
  }
  dropt <- ncol(txis) - length(kept) # hopefully less than suspicious

  txis$sf.unspl <- sf[colnames(txis), "unspl"]
  txis$sf.spl <- sf[colnames(txis), "spl"]
  
  if (dropt > 0) message("Dropped ", dropt, " uninformative cells.") 
  
  return(txis[, kept])

}
