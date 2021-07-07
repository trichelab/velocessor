#' fix or drop cells with an size factor of 0 (which will crash scvelo). 
#' 
#' Also attempts to use scuttle::cleanSizeFactors to recover dubious cells, 
#' if NumGenesExpressed is a column in colData(txis).  This may or may not be 
#' the world's greatest idea
#'
#' @param txis          a SingleCellExperiment
#' @param clean         fix negative size factors with cleanSizeFactors? (FALSE)
#' 
#' @return              a SingleCellExperiment with any useless cells omitted
#'
#' @importFrom scuttle  librarySizeFactors
#' @importFrom scuttle  cleanSizeFactors
#' 
#' @export
fix_cells <- function(txis, clean=TRUE) { 

  an <- c(spl="spliced", unspl="unspliced")
  sf <- data.frame(do.call(cbind,
                           lapply(an, function(x) 
                                  librarySizeFactors(txis, assay.type=x))))
  rownames(sf) <- colnames(txis)
  kept <- names(which(apply(sf, 1, function(x) all(x > 0))))
  dropt <- suspicious <- ncol(txis) - length(kept)
  if (suspicious > 0) message(suspicious, " suspicious cells.")
  
  if (clean & "NumGenesExpressed" %in% names(colData(txis))) {
    an2 <- names(sf)
    names(an2) <- an2 
    sf$num <- colData(txis)$NumGenesExpressed
    message("Cleaning up negative size factors...")
    sf <- data.frame(do.call(cbind,
                             lapply(an2, function(x) 
                                    cleanSizeFactors(sf[, x], sf[,"num"]))))
    rownames(sf) <- colnames(txis)
    kept <- names(which(apply(sf, 1, function(x) all(x > 0))))
    dropt <- ncol(txis) - length(kept) # hopefully less than suspicious
  }

  txis$sf.unspl <- sf[colnames(txis), "unspl"]
  txis$sf.spl <- sf[colnames(txis), "spl"]
  
  if (dropt > 0) message("Dropped ", dropt, " cells.") 
  
  return(txis[, kept])

}
