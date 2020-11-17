#' drop cells with an unspliced size factor of 0 (will crash scvelo). 
#' 
#' @param txis          a SingleCellExperiment
#' 
#' @return              a SingleCellExperiment with any useless cells omitted
#'
#' @importFrom scuttle  librarySizeFactors
#' 
#' @export
drop_uninformative_cells <- function(txis) { 

  an <- c(spl="spliced", unspl="unspliced")
  sf <- do.call(cbind,
                lapply(an, function(x) 
                       scuttle::librarySizeFactors(txis, assay.type=x)))
  kept <- names(which(apply(sf, 1, function(x) all(x > 0))))

  dropt <- ncol(txis) - length(kept)
  if (dropt > 0) message("Dropped ", dropt, " uninformative cells.") 
  
  txis$sf.unspl <- sf[colnames(txis), "unspl"]
  txis$sf.spl <- sf[colnames(txis), "spl"]
  return(txis[, kept])

}
