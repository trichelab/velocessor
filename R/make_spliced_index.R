#' aped from https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity
#' 
#' Given a fasta and a gtf, make a spliced index out of it. 
#' 
#' @param gtf           a GTF file, usually from ENSEMBL
#' @param fa            a FASTA file, usually a genome (same as .gtf but .fa)
#' @param flankLength   how many bases of intron (90)
#' @param verbose       squawk? (TRUE) 
#' 
#' @return      a list with information about the index and tximeta created
#'
#' @import eisaR
#' @import tximeta
#' @import Biostrings
#' @import GenomicFeatures
#'
#' @export 
make_spliced_index <- function(gtf, fa=NULL, flankLength=90L, verbose=TRUE) {

  # we don't actually care whether the GTF is bgzipped & tabixed
  if (!file.exists(gtf)) stop("Cannot find the GTF file ", gtf)
  if (is.null(fa)) fa <- sub("\\.gz", "", sub("gtf", "fa", gtf))
  if (!file.exists(fa)) stop("Cannot find the FASTA file ", fa)

  # extract the ranges for the cDNA seqs
  grl <- eisaR::getFeatureRanges(gtf=gtf, 
                                 featureType=c("spliced","intron"), 
                                 intronType="separate", 
                                 flankLength=flankLength, 
                                 joinOverlappingIntrons=FALSE,
                                 verbose=TRUE) 

  # extract the cDNA sequences for each tx
  genome <- Biostrings::readDNAStringSet(fa)
  names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
  seqs <- GenomicFeatures::extractTranscriptSeqs(x=genome, transcript=grl)

  # write out an annotated version for salmon to index
  expandedFasta <- sub("\\.gz", "", sub("gtf", "annotation.expanded.fa", gtf))
  Biostrings::writeXStringSet(seqs, filepath=expandedFasta)

  # write out an annotated version of the GTF for tximeta 
  expandedGtf <- sub("\\.gz", "", sub("gtf", "annotation.expanded.gtf", gtf))
  eisaR::exportToGtf(grl, filepath=expandedGtf)

  # write out the feature table 
  feat <- sub("\\.gz", "", sub("gtf", "annotation.expanded.features.tsv", gtf))
  write.table(metadata(grl)$corrgene, row.names=FALSE, col.names=TRUE, 
              file=feat, quote=FALSE, sep="\t")
  
  # write out the transcript to gene table 
  t2g <- sub("\\.gz", "", sub("gtf", "annotation.expanded.tx2gene.tsv", gtf))
  df <- eisaR::getTx2Gene(grl, filepath=t2g)

  # instructions
  stub <- sub("\\.fa$", "", fa)
  chrnames <- paste(stub, "chrnames", "txt", sep=".")
  sidx <- sub("\\.tx2gene\\.tsv", "sidx", t2g)
  message("To index with Salmon:")
  message("")
  message("grep ^\\> ", fa, " | cut -d \\> -f 2 | cut -d \\  -f 1 > ", chrnames)
  message("")
  message("salmon index \\")
  message("  -t <(cat ", expandedFasta, " ", fa, ") \\")
  message("  -i ", sidx, " \\")
  message("  -d ", chrnames, " \\")
  message("  -p 32")

  res <- list(gtf=gtf,
              fa=fa,
              expandedFasta=expandedFasta,
              expandedGtf=expandedGtf,
              genome=stub,
              feat=feat,
              t2g=t2g,
              sidx=sidx)

  return(res) 

}
