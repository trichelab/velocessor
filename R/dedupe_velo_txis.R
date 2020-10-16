#' deduplicate symbols in a SingleCellExperiment object with multiple assays
#' 
#' rowData(txis) *must* have columns 'symbol' and 'gene_id' for this to work.
#' Perhaps not coincidentally, this is what ensembldb returns for a given build.
#' 
#' @param txis        a SingleCellExperiment
#' 
#' @return            txis, but with duplicate representations of symbols summed
#'
#' @import Matrix
#' 
#' @export 
dedupe_velo_txis <- function(txis) {

  if (!all(c("symbol","gene_id") %in% names(rowData(txis)))) {
    message("rowData(txis) does not contain 'symbol' and 'gene_id'.")
    message("Returning untouched txis with any duplicates intact.")
    return(txis) 
  } else if (!any(duplicated(rowData(txis)$symbol))) {
    message("No duplicated symbols found. Returning txis unsummed.")
    rownames(txis) <- rowData(txis)$symbol
    return(txis) 
  }

  syms <- rowData(txis)$symbol
  uniquesyms <- unique(syms)
  rownames(txis) <- make.unique(syms)

  for (a in assayNames(txis)) {
    message("Deduplicating symbols in `", a, "`")
    assay(txis, a, withDimnames=FALSE) <- 
      Matrix(rowsum(assays(txis)[[a]], syms), sparse=TRUE)[syms, ]
  }
  return(txis)[uniquesyms, ] 

} 
