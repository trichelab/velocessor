#' deduplicate symbols in a SingleCellExperiment object with multiple assays
#' 
#' rowData(txis) *must* have columns 'symbol' and 'gene_id' for this to work.
#' Perhaps not coincidentally, this is what ensembldb returns for a given build.
#' To avoid potentially pointless calculations, this function first drops all
#' completely-zero ENS[MUS]G rows, then calls velocessor::find_dupe_txis, and
#' only works on the surviving dupes (if any). 
#' 
#' @param txis        a SingleCellExperiment
#' @param mincount    minimum count to register as nonzero (1)
#' @param mincells    how many cells must have at least mincount to retain? (1)
#' 
#' @return            txis, but with duplicate representations of symbols summed
#'
#' @import Matrix
#' 
#' @export 
dedupe_velo_txis <- function(txis, mincount=1, mincells=1) {

  if (!all(c("symbol","gene_id") %in% names(rowData(txis)))) {
    message("rowData(txis) does not contain 'symbol' and 'gene_id'.")
    message("Returning untouched txis with any duplicates intact.")
    return(txis) 
  } else if (!any(duplicated(rowData(txis)$symbol))) {
    message("No duplicated symbols found. Returning txis unsummed.")
    rownames(txis) <- rowData(txis)$symbol
    return(txis) 
  }

  duped <- find_dupe_txis(txis, mincount=mincount, mincells=mincells)
  todrop <- rownames(subset(duped), cells <  mincells)
  tokeep <- setdiff(rownames(txis), todrop) 
  txis <- txis[tokeep, ] 

  stillduped <- find_dupe_txis(txis, mincount=mincount) 
  if (nrow(stillduped) > 0) {

    syms <- rowData(txis)$symbol
    uniquesyms <- unique(syms)
    rownames(txis) <- make.unique(syms)

    for (a in assayNames(txis)) {
      message("Deduplicating symbols in `", a, "`")
      # FIXME: there must be a better way to do this on large Matrix objects.
      # For very large datasets, we may have to just drop the lesser gene :-O
      assay(txis, a, withDimnames=FALSE) <- 
        Matrix(rowsum(assays(txis)[[a]], syms), sparse=TRUE)[syms, ]
    }
    return(txis)[uniquesyms, ] 

  } else { 

    return(txis)

  }

} 
