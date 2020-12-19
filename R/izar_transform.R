#' generalizes the transformation from Izar et al 2020 for plate-seq data
#' 
#' In order to compare the results of plate-seq and droplet-seq preparations,
#' Izar and coauthors proposed a transformation of transcripts per million 
#' (TPM) for transcript `i` in cell `j` from a non-UMI plate-seq library as:
#' 
#' E[i,j] = log2( (TPM[i,j] / 10) + 1 ) 
#' 
#' Here we generalize this slightly by allowing a variable scaling and offset.
#' The goal is to make non-UMI plate-seq and UMI droplet-seq data comparable.
#' We caution the user that 'comparable' is a subjective term here.
#' 
#' @param txis        a SingleCellExperiment
#' @param orig        the name of the assay to transform ("counts") 
#' @param ret         return "sce" (default) or "matrix" of transformed values?
#' @param dupe_rate   an estimated duplication rate (default is 10, per Izar)
#' @param pseudocount a pseudocount to add to TPM/dupe_rate (default 1, ibid)
#' 
#' @details 
#'
#' Most plate-seq protocols do not add unique molecular indices (UMIs) to 
#' each template fragment in a library (plate-seq protocols WITH UMIs include 
#' Quartz-Seq2, SMART-Seq2, and STORM-UMI). In contrast, almost all droplet 
#' protocols apply both a cell barcode (CB) and UMI barcode (UB) to each 
#' template fragment, disambiguating whether two fragments that map to the 
#' same reference sequence are from the same template molecule or not. In 
#' order to compare the results of plate-seq and droplet-seq preparations,
#' Izar and colleagues proposed this transformation. One may additionally 
#' apply methods such as `sctransform` to the resulting estimate, harmonize
#' the resulting matrix, or similar shenanigans. Alternatively, one may use
#' a UMI-enabled plate-seq preparation to elide this transformation. 
#' 
#' @return            a SingleCellExperiment with `izar` or matrix (per `ret`)
#'
#' @references 
#' Izar B, Tirosh I, Stover EH, et al. A single-cell landscape of high-grade 
#' serous ovarian cancer. Nat Med 26, 1271–1279 (2020).
#' https://doi.org/10.1038/s41591-020-0926-0
#' 
#' @import            SingleCellExperiment
#' 
#' @export 
izar_transform <- function(txis, ret=c("sce","matrix"), orig="counts", dupe_rate=10, pseudocount=1) {
  
  ret <- match.arg(ret)
  stopifnot(dupe_rate > 0)
  assay(txis, "izar") <- log2( (assay(txis, orig) / dupe_rate) + pseudocount )
  metadata(txis)$izar_dupe_rate <- dupe_rate
  metadata(txis)$izar_pseudocount <- pseudocount
  if (ret == "sce") return(txis)
  else return(assay(txis, "izar"))

}
