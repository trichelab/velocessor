#' import plate-seq runs from Kallisto into a SingleCellExperiment
#'
#' Automated (more or less) importation of STORM-seq and SMART-seq[123] data 
#' for downstream attempts at merging, index visualization, etc.  Assumes that
#' users will quantify against ENSEMBL transcripts instead of ridiculous crap 
#' like GENCODE and similar. Default is ERCC spikes, ENSembl transcripts, and
#' if found, intronic quantifications of the same. Everything else gets labeled
#' as 'viruses_and_repeats' since in our transcriptomes, that's what it is. If
#' altExps is set to FALSE, no splitting will be attempted. 
#' 
#' @param   runs    a vector of .tsv file paths (usually abundance.tsv) w/names
#' @param   what    columns to import (default: tpm, est_counts, eff_length)
#' @param   altExps split spliced, unspliced, spikes, viruses & repeats? (TRUE)
#' @param   BPPARAM optional BiocParallel parameter object for bplapply()
#' @param   ...     additional parameters to pass to .split_altExps
#' 
#' @return          A SingleCellExperiment
#'
#' @import          SingleCellExperiment
#' @import          BiocParallel
#' 
#' @export
import_kallisto_txis <- function(runs, 
                                 what=c(tpm="tpm",
                                        counts="est_counts",
                                        eff_length="eff_length"),
                                 altExps=TRUE, 
                                 BPPARAM=bpparam(),
                                 ...) {

  stopifnot(!is.null(names(runs)))
  message("Processing ", length(runs), " runs.")

  message("Importing first run as a prototype to check transcript names.")
  proto <- .import_kallisto(runs[1])[, c("est_counts", "eff_length", "tpm")]
  
  # yes this would be faster in C/C++.  No, I don't care 
  message("Importing all ", length(runs), " runs in parallel...")
  imported <- bplapply(runs, .import_kallisto, proto=proto, BPPARAM=BPPARAM)
  mats <- lapply(what, .extract_columns, imported=imported, BPPARAM=BPPARAM)

  message("Constructing SingleCellExperiment...")
  sce <- SingleCellExperiment(mats) 
  if (altExps) sce <- .split_altExps(sce, ...)
  return(sce)

}


# helper fn
.import_kallisto_tsv <- function(x, ...) { 
  
  read.table(x, sep="\t", header=TRUE, row.names=1, 
             colClasses=c(target_id="character", 
                          length="integer", 
                          eff_length="numeric",
                          est_counts="numeric", 
                          tpm="numeric"))

}


# helper fn
.import_kallisto <- function(x, proto=NULL) { 

  message("Importing ", x, " ... ", appendLF=FALSE)
  res <- .import_kallisto_tsv(x)
  message("done.")
  if (!is.null(proto)) return(res[rownames(proto), colnames(proto)])
  else return(res)

}


# helper fn
.extract_columns <- function(imported, what="tpm", BPPARAM=bpparam()) {
  
  getcol <- function(x) x[, what, drop=FALSE] 
  message("Extracting ", what, " for each run...")
  res <- do.call(cbind, bplapply(imported, getcol, BPPARAM=BPPARAM))
  rownames(res) <- rownames(imported[[1]])
  colnames(res) <- names(imported)
  return(res)

}


# helper fn 
.split_altExps <- function(sce, 
                           tx="ENS",
                           sep="\\-", 
                           unspliced="\\-I",
                           spikes="ERCC", 
                           t2g=NULL) {

  message("Categorizing transcripts...")
  txs <- grep(unspliced, invert=TRUE, value=TRUE, 
              grep(tx, rownames(sce), value=TRUE))
  txus <- grep(unspliced, value=TRUE, 
               grep(tx, rownames(sce), value=TRUE)) 
  txu2tx <- sapply(strsplit(txus, sep), `[`, 1)
  names(txu2tx) <- txus

  rowData(sce)$type <- "viruses_and_repeats"
  rowData(sce[txs, ])$type <- "spliced"
  rowData(sce[txus, ])$type <- "unspliced"
  rowData(sce[grep(spikes, rownames(sce), value=TRUE),])$type <- "spikes"
  splitAltExps(sce, rowData(sce)$type, "spliced")

}
