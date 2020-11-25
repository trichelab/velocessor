#' massage a (possibly velocitized) SingleCellExperiment into an infercnv object
#' 
#' @details
#' We do not require or import infercnv, because the requirement for JAGS can 
#' be a showstopper for HPC installations. This function will stop() if the
#' infercnv package namespace cannot be attached, however, for obvious reasons. 
#' I have no idea why the `infercnv` function name was not taken by `infercnv`, 
#' but since it wasn't, I stole it. 
#' 
#' @param x           a SingleCellExperiment
#' @param anno_col    the colData column containing cell annotations ("status")
#' @param ref_groups  annotations of the cells to use as references ("normal")
#' @param ...         any other arguments for infercnv::CreateInfercnvObject
#'
#' @return            an infercnv object
#' 
#' @import            GenomeInfoDb
#'
#' @export
infercnv <- function(x, anno_col="status", ref_groups=c("normal"), ...) { 

  # needed for coordinates 
  stopifnot(length(unique(genome(x))) > 0) 

  # needed for inferring CNV (duh)
  if (!requireNamespace("infercnv")) {
    message("You can install infercnv with BiocManager::install(\"infercnv\").")
    stop("You cannot create infercnv objects without first installing infercnv")
  } 

  # needed for my sanity
  seqlevelsStyle(x) <- "UCSC"
  x <- sort(keepStandardChromosomes(x))

  # needed to not waste a lot of time
  stopifnot(anno_col %in% names(colData(x)))
  stopifnot(any(colData(x)[, anno_col] %in% ref_groups))
  
  # now we should have all we need to proceed:
  gene_order <- as.data.frame(rowRanges(x))[, c("seqnames","start","end")]
  anno_df <- as.data.frame(colData(x)[, anno_col, drop=FALSE])
  raw_counts <- counts(x) # spliced(x)?
  infercnv::CreateInfercnvObject(raw_counts_matrix=raw_counts, 
                                 gene_order_file=gene_order,
                                 annotations_file=anno_df,
                                 ref_group_names=ref_groups, 
                                 ...)

}
