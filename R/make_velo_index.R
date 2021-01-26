#' aped from https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity
#' 
#' Given a fasta and a gtf, make a spliced index out of it. Like duh, the idea
#' is to eventually be able to create linkedTxomes that are really pangenomes.
#' 
#' @param gtf           a GTF file, usually from ENSEMBL
#' @param fa            a FASTA file, usually a genome (same as .gtf but .fa)
#' @param flank         how many bases of intron (90)
#' @param intron        include unspliced transcripts? (TRUE; for debugging)
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
make_velo_index <- function(gtf, fa=NULL, flank=90L, intron=TRUE, verbose=TRUE){

  # we don't actually care whether the GTF is bgzipped & tabixed
  if (!file.exists(gtf)) stop("Cannot find the GTF file ", gtf)
  if (is.null(fa)) fa <- sub("\\.gz", "", sub("gtf", "fa", gtf))
  if (!file.exists(fa)) stop("Cannot find the FASTA file ", fa)

  # spliced and unspliced, or just spliced?
  feats <- c("spliced")
  if (intron) feats <- c("spliced","intron")

  # extract the ranges for the cDNA seqs
  grl <- eisaR::getFeatureRanges(gtf=gtf, 
                                 featureType=feats,
                                 intronType="separate", 
                                 flankLength=flank, 
                                 joinOverlappingIntrons=FALSE,
                                 verbose=TRUE) 

  # spliced-only vs. velocity-aware
  expansion <- c("static", "velo")[length(feats)] # somewhat foolproof

  # extract the cDNA sequences for each tx
  genome <- Biostrings::readDNAStringSet(fa)
  names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
  GenomeInfoDb::seqlevelsStyle(names(genome)) <- 
    GenomeInfoDb::seqlevelsStyle(grl) # somewhat awkward
  
  # warn if there are transcripts on contigs not in the FASTA, and prune them
  if (any(!seqlevels(grl) %in% names(genome))) {
    warning("Your GTF has transcripts on contigs not present in your FASTA.")
    grl <- keepSeqlevels(grl, names(genome), pruning.mode="coarse")
    warning("Transcriptome pruned to only reference contigs in FASTA.")
  }
  # in the past, this usually meant the seqlevelsStyle was out of sync

  # proceed to extract the appropriate transcript sequences
  seqs <- GenomicFeatures::extractTranscriptSeqs(x=genome, transcript=grl)

  # write out an annotated version for salmon to index
  faSub <- paste(expansion, "fa", sep=".")
  expandedFasta <- sub("\\.gz", "", sub("gtf", faSub, gtf))
  Biostrings::writeXStringSet(seqs, filepath=expandedFasta)

  # write out an annotated version of the GTF for tximeta 
  gtfSub <- paste(expansion, "gtf", sep=".")
  expandedGtf <- sub("\\.gz", "", sub("gtf", gtfSub, gtf))
  eisaR::exportToGtf(grl, filepath=expandedGtf)

  # write out the feature table 
  featureSub <- paste(expansion, "features", "tsv", sep=".")
  feat <- sub("\\.gz", "", sub("gtf", featureSub, gtf))
  write.table(metadata(grl)$corrgene, row.names=FALSE, col.names=TRUE, 
              file=feat, quote=FALSE, sep="\t")
  
  # write out the transcript to gene table 
  t2gSub <- paste(expansion, "tx2gene", "tsv", sep=".")
  t2g <- sub("\\.gz", "", sub("gtf", "expanded.tx2gene.tsv", gtf))
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
  message("")
  message("Then run make_tximeta(x), where x is results from make_velo_index()")
  message("")
  
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
