#' import plate-seq velocity information (from Salmon) to a SingleCellExperiment
#'
#' Automated (more or less) importation of STORM-seq and SMART-seq[123] data for
#' downstream attempts at merging with droplet data, index sort visualization, 
#' and so forth. 
#' 
#' FIXME: Add Kallisto support and Arkas style txome/repeatome/spikeome support.
#' FIXME: Add in ability to support bootstraps/Gibbs samples from salmon
#'
#' @param   quants  where the quant.sf files are
#' @param   t2g     required file with tx2gene information for TPM calcs
#' @param   type    What type of quantifications are these? ("salmon") 
#' @param   gtf     Where an expanded GTF lives if not collapsed to gene-level
#' @param   qc      Whether QC should be performed to toss low-quality cells
#' @param   fixrn   Fix goofy versioned ENSEMBL gene names? (TRUE) 
#' @param   ...     additional parameters to pass to tximport, if any
#' 
#' @return          A SingleCellExperiment with 'spliced' & 'unspliced' assays.
#'
#' @import scuttle
#' @import tximeta
#' @import tximport
#' @import jsonlite
#' @import rtracklayer
#' @import SingleCellExperiment
#' 
#' @export
import_plate_txis <- function(quants, t2g=NA, type="salmon",
                              gtf=NULL, qc=TRUE, fixrn=TRUE, ...) {

  if (is.na(t2g)) {
    message("No transcript to gene mapping file provided. Cannot compute TPM.")
    message("(Usually a name like mm10.ens101.annotation.expanded.tx2gene.tsv)")
    stop("This can be guessed from the index, but that is usually a poor idea.")
  } else if (!file.exists(t2g)) {
    message(t2g, " does not exist or cannot be read by this process.")
    stop("Cannot calculate per-gene TPM without tx2gene mappings.")
  } else { 
    tx2gene <- read.delim(t2g, sep="\t", head=FALSE)
  }

  if (!all(file.exists(quants))) stop("Some of your quant files don't exist.")
  cmds <- unlist(lapply(quants, FUN = function(x) gsub("quant\\.*sf", "cmd_info.json", x)))
  if (!all(file.exists(cmds))) stop("Some quants don't have cmd_info.json")
  if (!is.null(gtf)) {
    stopifnot(file.exists(gtf))
    gtfs <- gtf
  } else {
    ## this means they were collapsed to gene level already
    gtfs <- sapply(cmds, .getGTF) 
  }
  if (length(unique(gtfs)) > 1) stop("Some quants use different GTF files.") 
  gtf <- unique(gtfs)  
  stopifnot(file.exists(gtf))
  # because the next step can take a LONG time:

  message("Processing ", quants[1], " through ", quants[length(quants)], "...")
  mats <- tximport(quants, type="salmon", txIn=TRUE, tx2gene=tx2gene, ...)
  asys <- mats[ c("counts", "abundance") ]
  names(asys) <- c("counts", "tpm")
  
  message("Reading annotations from ", gtf, "...")
  rr <- .rowRanges(gtf, asys)
  
  message("Constructing sample annotations...")
  cd <- DataFrame(sample=make.unique(quants))
  if (!is.null(names(quants))) cd <- DataFrame(sample=names(quants))
  txi <- SingleCellExperiment(asys, rowRanges=rr, colData=cd)
  
  ## catch and remove empty cells
  txi <- .removeEmptyCells(txi)

  message("Splitting...")
  feats <- sub("\\.gtf", ".features.tsv", gtf)
  cg <- read.delim(feats, header=TRUE, as.is=TRUE)
  if ("intron" %in% colnames(cg)) {
    message("WARNING: using a full-length protocol and only using intron pieces may miss relevant reads/unspliced transcripts.")
    message("Ideally, you want to create the reference with 'unspliced' instead of 'introns' with eisaR.")
    colnames(cg)[colnames(cg) == "intron"] <- "unspliced"
  } else {
    ## check to make sure 'unspliced' is there in the features tsv
    stopifnot("unspliced" %in% colnames(cg))
  }
  txis <- tximeta::splitSE(txi, cg, assayName="counts")
  assay(txis, "counts") <- assays(txis)[[1]] + assays(txis)[[2]]
  txis2 <- tximeta::splitSE(txi, cg, assayName="tpm")
  assayNames(txis2) <- paste(assayNames(txis2), "tpm", sep="_")
  assay(txis2, "tpm") <- assays(txis2)[[1]] + assays(txis2)[[2]]
  for (an in assayNames(txis2)) assay(txis, an) <- assay(txis2, an)
  txis <- as(txis, "SingleCellExperiment") 
  colnames(txis) <- make.unique(colnames(txis)) # else many things will fail
  rm(txis2)

  message("adding NumGenesExpressed...")
  txis$NumGenesExpressed <- colSums(counts(txis) > 0)
  
  if (qc) {
    message("quick QC...")
    txis <- .quickQC(txis)
  }

  if (fixrn) {
    message("Fixing goofy row names...")
    rownames(txis) <- fix_rownames(txis)
  }

  message("adding logNormCounts...")
  txis <- scuttle::logNormCounts(txis)

  message("Done.")
  metadata(txis)$origin <- "plate"
  return(txis) 

}


# helper fn
.getGTF <- function(cmd) jsonlite::fromJSON(cmd)$geneMap


# helper fn
.rowRanges <- function(gtf, asys) { 

  gxs <- subset(rtracklayer::import(gtf), type=="gene")
  names(gxs) <- gxs$gene_id
  granges(gxs)[rownames(asys$counts)]

}

# helper fn
.removeEmptyCells <- function(se) {
  message("Removing failed cells.")
  se[,colSums(assay(se)) > 0]
}

# helper fn
.quickQC <- function(se) {
  qcstats <- scuttle::perCellQCMetrics(se)
  qcfilter <- scuttle::quickPerCellQC(qcstats)
  se[,!qcfilter$discard]
}
