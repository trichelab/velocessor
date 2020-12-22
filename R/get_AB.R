#' a shim function, mostly to feed plot_ABC()
#'
#' binarize & TF-IDF SCE, then call boundaries
#' requires a genome for the SCE or it will fail 
#' requires coordinates for the SCE or it will fail
#' don't run this on thousands of cells ungrouped
#' 
#' @param   x     something that inherits from RangedSummarizedExperiment
#' @param   asy   assay to feed compartmap ("counts"; see Details)
#' @param   chr   what chromosome to work on ("chr3") 
#' @param   res   resolution for compartmap bins, in base pairs (1e5)
#' @param   boot  how many bootstraps to run for compartmap (10) 
#' @param   minct minimum number of counts to consider "nonzero" (see Details)
#' @param   mincl minimum number of cells to consider "nonzero" (see Details)
#' @param   ...   parameters to feed to getATACABsignal (e.g. group, cores, etc)
#'
#' @return        AB signal
#'
#' @details 
#' 
#' The value "ratio" for `asy` means "unspliced / spliced" and is computed on
#' the fly, for the purpose of identifying "unstable" regions of chromatin. It
#' is not (repeat, NOT) ready for prime time, not least due to TF-IDF vs. LSI.
#' That said, the value of `min` (which defaults to 1) plays an important role
#' in determining which transcripts to consider: any transcript that does not 
#' have at least `minct` counts unspliced in at least `mincl` cells will not
#' be considered in the resulting downstream analyses (which are subsetted to
#' "nonzero" regions of chromatin for conformational modeling purposes). 
#' 
#' @seealso       compartmap::getATACABsignal
#'
#' @import        compartmap
#'
#' @export
get_AB <- function(x, asy="counts", chr="chr3", res=1e5, boot=10, minct=1, mincl=1, ...) { 

  gen <- .get_genome(x)
  seqlevelsStyle(x) <- "UCSC"
  x <- keepStandardChromosomes(x, pruning.mode="coarse") 

  if (asy == "ratio") { 
    message("This feature is experimental.  Use at your own risk.")
    x <- x[which(rowSums2(unspliced(x) >= minct) >= mincl), ] # i.e., nonzero
    assay(x, paste0("raw", asy)) <- assay(x, "unspliced") / assay(x, "spliced")
  } else if (!asy %in% assayNames(x)) { 
    stop("The assay named ", asy, " could not be found in your object.")
  } else { 
    assay(x, paste0("raw", asy)) <- assay(x, asy)
  }
 
  # this does bad things sometimes  
  assay(x, asy) <- t(compartmap::transformTFIDF(as.matrix(assay(x, asy))))
  AB <- getATACABsignal(x, res = res, chr = chr, genome = gen, 
                        num.bootstraps = boot, bootstrap = TRUE, ...) 
  AB <- fixCompartments(AB)

  # what's this do? 
  #
  # plotAB(k562_scatac_group_chr14_100boot.fix, 
  #   chr = "chr14", what = "flip.score", with.ci = T, median.conf = T)
  #
  return(AB) 

}


# helper fn
.get_genome <- function(sce) { 

  gen <- unique(genome(sce))[1]
  if (length(gen) == 0) stop("Cannot proceed without genomic coordinates.")
  else if (tolower(gen) %in% c("grch38", "hg38")) return("hg38") 
  else if (tolower(gen) %in% c("grch37", "hg19")) return("hg38") 
  else if (tolower(gen) %in% c("grcm38", "mm10")) return("mm10") 
  else if (tolower(gen) %in% c("grcm37", "mm9")) return("mm9") 
  else stop("Unsupported genome \"", gen, "\". Exiting.")

}
