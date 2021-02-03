#' strip silly ENSEMBL gene/transcript version IDs 
#' 
#' tired of typing this out 
#' 
#' @param   txis  an SCE
#' @param   stub  optional leader stub ("^ENS")
#' @param   sep   optional separator ("\\.")
#' @param   idx   which fragment to return (1) 
#' 
#' @return        same SCE but without goofy rownames 
#' 
#' @import        SingleCellExperiment
#' 
#' @export 
fix_rownames <- function(txis, stub="^ENS", sep="\\.", idx=1) { 

  fixable <- grep(stub, rownames(txis))
  if (length(fixable) < 1) { 
    message("No fixable rownames. Returning unaltered.")
    return(txis)
  }
  rownames(txis)[fixable] <- 
    sapply(strsplit(rownames(txis)[fixable], sep), `[`, idx)
  message("Fixed ", length(fixable), " rownames.") 
  return(txis) 

}
