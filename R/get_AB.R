#' a shim function
#'
#' binarize and TF-IDF, then call boundaries
#' 
#' @param   x     something that inherits from RangedSummarizedExperiment
#' @param   asy   the name of the assay to summarize down ("counts") 
#' @param   chr   what chromosome to work on ("chr3") 
#' @param   res   resolution for compartmap bins, in base pairs (1e5)
#' @param   ...   additional parameters to feed to getATACABsignal 
#'
#' @return        AB signal 
#'
#' @import compartmap
#' 
get_AB <- function(x, asy="counts", chr="chr3", res=1e5, ...) { 

  gen <- unique(genome(x))[1]
  if (is.null(gen) | is.na(gen)) gen <- "hg19"
  assay(x, paste0("raw", asy)) <- assay(x, asy)
  assay(x, asy) <- t(compartmap::transformTFIDF(assay(x, asy)))
  AB <- getATACABsignal(x, res = res, chr = chr, genome = gen, 
                        num.bootstraps = 100, bootstrap = TRUE, 
                        group = FALSE, boot.parallel = TRUE, boot.cores = 8)
  # what's this do? 
  #
  # k562_scatac_group_chr14_100boot.fix <- 
  #   fixCompartments(k562_scatac_group_chr14_100boot)
  #
  # plotAB(k562_scatac_group_chr14_100boot.fix, 
  #   chr = "chr14", what = "flip.score", with.ci = T, median.conf = T)

  return(AB) 

}
