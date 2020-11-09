#' import velocity information (usually from Alevin) into a SingleCellExperiment
#'
#' This function is usually called by process_velo_txis(). The 'cdat' argument
#' takes a single-row data.frame containing 'names' and 'files' as its columns. 
#' 'quant' is the type of quantifications; at present, only alevin is supported.
#' 'sep' allows deduplicating barcodes by sample using cbind without incident 
#' once the results of import_velo_txis are returned to process_velo_txis.
#' 
#' @param   cdat    A row of the colData prepared by process_velo_txis.
#' @param   quant   What type of quantifications are these? ("alevin")
#' @param   meta    Use tximeta for imports? (FALSE; can fail unexpectedly)
#' @param   QC      Perform rudimentary quality control? (TRUE) 
#' @param   sep     What string separates sample name from barcode? ("_")
#' 
#' @return          A SingleCellExperiment with 'spliced' & 'unspliced' assays.
#' 
#' @import scran 
#' @import scater 
#' @import scuttle
#' @import tximeta
#' @import jsonlite
#' @import fishpond
#' @import rtracklayer
#' @import MatrixGenerics
#' @import SingleCellExperiment
#' 
#' @export
import_velo_txis <- function(cdat, quant=c("alevin"), meta=FALSE, QC=TRUE, sep="_") {

  quant <- match.arg(quant) 
  stopifnot(length(cdat[['names']]) == 1)
  stopifnot(all(c("names", "files") %in% names(cdat)))

  root <- strsplit(cdat[['files']], "/")[[1]][1]
  alevin <- file.path(root, "alevin")
  index <- jsonlite::fromJSON(file.path(root, "cmd_info.json"))$index

  txome <- cdat[['txome']]
  gtf <- jsonlite::fromJSON(txome)[1, "gtf"]
  stopifnot(file.exists(gtf)) # so tximeta doesn't puke

  txstub <- sub("\\.sidx$", "", index)
  feats <- paste(txstub, "features", "tsv", sep=".")
  stopifnot(file.exists(feats))

  message("Processing ", cdat[['names']], "...")

  message("Importing...")
  # this is infuriating -- it ditches entire columns of counts.  Fuck that.
  if (meta) { 

    message("Note: if the following import fails, try setting `meta=FALSE`")
    txi <- tximeta(coldata=cdat, type="alevin")
  
  } else { 
  
    message("meta==FALSE; working directly from Alevin counts and GTF file...")
    asys <- list(counts=.readFish(alevin))
    rr <- .rowRangesFromGtf(gtf)[rownames(asys$counts)]
    cd <- DataFrame(attr(asys$counts, "stats"))[colnames(asys$counts),]
    txi <- SingleCellExperiment(asys, rowRanges=rr, colData=cd)

  } 

  message("Renaming...")
  colnames(txi) <- paste(cdat['names'], colnames(txi), sep=sep)
  colData(txi)$numTranscripts <- colSums(assay(txi, "counts")) # for TPM'ing

  # needs to be sped up 
  message("Splitting...")
  cg <- read.delim(feats, header=TRUE, as.is=TRUE)
  colnames(cg)[colnames(cg) == "intron"] <- "unspliced"
  print(system.time(txis <- tximeta::splitSE(txi, cg, assayName="counts")))

  message("Converting...")
  print(system.time(txis <- as(txis, "SingleCellExperiment")))
  print(assayNames(txis))
  assays(txis) <- list(counts=assays(txis)[["spliced"]],
                       spliced=assays(txis)[["spliced"]],
                       unspliced=assays(txis)[["unspliced"]])

  # stupid simple QC
  if (QC) {
    message("QCing...")
    print(system.time(txis <- .doQC(txis)))
  } else { 
    message("adding logNormCounts...")
    print(system.time(txis <- scuttle::logNormCounts(txis)))
  } 

  message("Done.")
  return(txis) 

}


# helper fn
.readFish <- function(alevin, tier=TRUE) { 
 
  cells <- read.table(file.path(alevin, "quants_mat_rows.txt"))[,1]
  ncells <- length(cells)
  genes <- read.table(file.path(alevin, "quants_mat_cols.txt"))[,1]
  ngenes <- length(genes)
  qm <- file.path(alevin, ifelse(tier, "quants_tier_mat.gz", "quants_mat.gz"))
  counts <- fishpond::readEDS(ngenes, ncells, qm, tierImport=tier)
  stats <- read.table(file.path(alevin, "featureDump.txt"), head=TRUE, row=1)
  attr(counts, "stats") <- stats
  rownames(counts) <- genes
  colnames(counts) <- cells
  return(counts) 

}


# helper fn
.rowRangesFromGtf <- function(gtf) { 

  gxs <- subset(rtracklayer::import(gtf), type=="gene")
  names(gxs) <- gxs$gene_id
  granges(gxs) 

}


# helper fn
.doQC <- function(txis, verbose=FALSE) { 

  if (verbose) message("Starting QC...")
  qcstats <- scuttle::perCellQCMetrics(txis)
  qcfilter <- scuttle::quickPerCellQC(qcstats)
  txis <- txis[, !qcfilter$discard]
  if (verbose) message(length(qcfilter$discard), " cells dropped")
  txis <- scuttle::logNormCounts(txis)
  colData(txis)$QCed <- TRUE
  return(txis)

}
