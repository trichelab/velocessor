#' Mark duplicated mcols(txis)$symbol in a SingleCellExperiment object `txis`
#' 
#' For large experiments, it may not be possible to run rowsum on a Matrix in
#' limited memory. And, on occasion, it may not be necessary -- if only one of 
#' the ENSEMBL gene IDs mapping to an HGNC or mouse symbol is nonzero, we don't
#' need to sum them at all; dropping the near-zero rows is sufficient. This 
#' function allows looping through the ENSEMBL IDs for each "duplicated" symbol
#' and determining whether rowsum() needs to be called at all. 
#' 
#' @param txis        a SingleCellExperiment
#' @param mincount    minimum count to register as nonzero (1)
#' @param mincells    how many cells must have at least mincount to retain? (1)
#' 
#' @return            data.frame with symbol, ensembl, index, cutoff, cells
#'
#' @details 
#' 
#' If there are no duplicated symbols, an empty data.frame results.
#' 
#' @import Matrix
#' 
#' @export 
get_dupe_txis <- function(txis, mincount=1, mincells=1) {

  empty <- rep(NA, 0) 
  if (!all(c("symbol","gene_id") %in% names(mcols(txis)))) {
    stop("mcols(txis) needs both 'symbol' and 'gene_id'. Cannot proceed.")
  } else if (!any(duplicated(mcols(txis)$symbol))) {
    message("No duplicated symbols found.")
    syms <- ensg <- idxs <- cuts <- cells <- empty
  } else { 
    duped <- sort(mcols(txis)$symbol[which(duplicated(mcols(txis)$symbol))])
    ensg <- subset(rownames(txis), mcols(txis)$symbol %in% duped)
    syms <- mcols(txis)[ensg, "symbol"]
    idxs <- match(ensg, rownames(txis))
    cuts <- rep(mincount, length(idxs))
    cells <- apply(counts(txis)[idxs,], 1, function(x) sum(x > mincount))
  }

  ret <- data.frame(symbol=syms,
                    ensembl=ensg,
                    index=idxs,
                    cuts=cuts,
                    cells=cells)
  return(ret)  

}
