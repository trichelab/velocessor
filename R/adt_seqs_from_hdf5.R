#' dump ADT/HTO sequences from an h5 mRNA+ADT file for [re]indexing
#' 
#' The resulting FASTA can then be fed to Salmon, Kallisto, STAR, etc.
#' 
#' @param   h5_fname  the h5 filename (required)
#' @param   fa_fname  the fa filename to write ("ADTs.fa"; NULL elides FASTA)
#' @param   base      the base name group in the h5 file ("matrix") 
#' 
#' @return            invisibly, a DataFrame of the ADT features and sequences
#' 
#' @details 
#' 
#' Depending on the dataset, the base name for the features may differ from 
#' the usual "matrix".  This can be detected using rhdf5::h5ls(h5_fname). 
#' 
#' @import  GenomicRanges
#' @import  Biostrings
#' @import  rhdf5
#' 
#' @export
adt_seqs_from_hdf5 <- function(h5_fname, fa_fname="ADTs.fa", base="matrix") { 

  library(Biostrings) 
  library(GenomicRanges) 
  bname <- paste0(base, "/features")
  ADTs <- subset(as.data.frame(h5read(h5_fname, name=bname)[-1]),
                 feature_type == "Antibody Capture")
  if (any(duplicated(ADTs$name))) stop("Duplicate ADT names!")
  rownames(ADTs) <- ADTs$name
  seqs <- DNAStringSet(ADTs$sequence)
  names(seqs) <- ADTs$name
  ADTs <- DataFrame(ADTs)
  ADTs$sequence <- seqs
  if (!is.null(fa_fname)) writeXStringSet(seqs, fa_fname)
  invisible(ADTs[, c("id","sequence")]) 

}
