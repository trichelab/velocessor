#' import droplet RNAseq data (usually from Alevin) into a SingleCellExperiment
#'
#' This function is usually called by process_velo_txis(). `quants` is, like in
#' import_plate_txis, a bunch of directories with quantifications in them. `QC`
#' does not exist in import_plate_txis, but it is a flag to perform rudimentary
#' quality control for droplet data (and is enabled by default). `sep` is the 
#' separator character between the assigned sample name and the trimmed barcode.
#' 
#' FIXME: Add Kallisto support and Arkas style txome/repeatome/spikeome support.
#' 
#' @param   quants  Directories containing droplet quantification results 
#' @param   feats   optional file with intron-exon-gene mappings (guess) 
#' @param   type    What type of quantifications are these? ("alevin") 
#' @param   QC      Perform rudimentary quality control? (TRUE) 
#' @param   tpms    compute TPMs? (FALSE; can create a CHOLMOD error)
#' @param   su_tpms compute spliced/unspliced TPMs? (FALSE; as above)
#' @param   sep     What string separates sample name from barcode? ("_")
#' @param   fixrn   Fix goofy versioned ENSEMBL gene names? (TRUE) 
#' @param   ...     additional parameters to pass to tximport, if any
#' 
#' @details 
#' TPM calculation is disabled by default due to a tendency to crash sessions.
#' The function velocessor::add_tpm handles the grunt work for this at present.
#' It would probably be a good idea to do it more efficiently (e.g. in C++). 
#' 
#' @return          A SingleCellExperiment with 'spliced' & 'unspliced' assays.
#' 
#' @import scran 
#' @import scater 
#' @import scuttle
#' @import tximeta
#' @import tximport
#' @import jsonlite
#' @import fishpond
#' @import MatrixGenerics
#' @import SingleCellExperiment
#' 
#' @export
import_droplet_txis <- function(quants, feats=NULL, type=c("alevin"), QC=TRUE, tpms=FALSE, su_tpms=FALSE, sep="_", fixrn=TRUE,...) {

  params <- data.frame() 
  type <- match.arg(type) 
  if (type == "alevin") {

    files <- file.path(quants, "alevin", "quants_mat.gz")
    stats <- file.path(quants, "alevin", "featureDump.txt")

    indices <- sapply(file.path(quants, "cmd_info.json"), 
                      function(x) jsonlite::fromJSON(x)$index)
    if (length(unique(indices)) > 1) stop("Your quants use different indices!")

    gtfs <- sapply(file.path(quants, "cmd_info.json"), function(x) 
                   sub("tx2gene.tsv", "gtf", jsonlite::fromJSON(x)$tgMap))
    if (length(unique(gtfs)) > 1) stop("Your quants use different GTFs!")
    stopifnot(file.exists(unique(gtfs))) # so tximeta doesn't puke
    if (is.null(feats)) feats <- sub("\\.gtf", ".features.tsv", unique(gtfs))
    stopifnot(file.exists(feats))

    params <- data.frame(quants=quants, 
                         files=files, 
                         index=indices,
                         gtf=gtfs,
                         stats=stats, 
                         feats=feats)
    if (!is.null(names(quants))) {
      rownames(params) <- names(quants)
    } else {
      message("names(quants) is NULL; sample names may be confusing")
      rownames(params) <- sub("_alevin", "", basename(quants))
    }

  }  # add kallisto etc. here in future

  stopifnot("feats" %in% names(params)) 
  qnames <- rownames(params)
  message("Processing ", paste(qnames, collapse=" and "))
  names(qnames) <- qnames 

  # could presumably mclapply or bplapply here
  txi <- do.call(cbind, lapply(qnames, .import_alevin, params=params, sep=sep))

  # needs to be sped up 
  message("Splitting...")
  cg <- read.delim(unique(params$feats), header=TRUE, as.is=TRUE)
  colnames(cg)[colnames(cg) == "intron"] <- "unspliced"
  print(system.time(txis <- tximeta::splitSE(txi, cg, assayName="counts")))

  message("Converting...")
  print(system.time(txis <- as(txis, "SingleCellExperiment")))
  assays(txis) <- list(counts=assays(txis)[["spliced"]],
                       spliced=assays(txis)[["spliced"]],
                       unspliced=assays(txis)[["unspliced"]])
 
  if (tpms) {  
    message("Calculating TPMs...")
    txis <- add_tpm(txis)
  }

  # stupid simple QC
  if (QC) {
    message("QCing...")
    print(system.time(txis <- .do_droplet_QC(txis)))
  } else { 
    message("adding logNormCounts...")
    print(system.time(txis <- scuttle::logNormCounts(txis)))
  } 

  if (fixrn) {
    message("Fixing goofy row names...")
    rownames(txis) <- fix_rownames(txis)
  }

  message("Done.")
  metadata(txis)$origin <- "droplet"
  return(txis) 

}


# helper fn
.import_alevin <- function(qname, params, sep="_") {
  
  param <- params[qname, , drop=FALSE] 
  message("Importing ", param[, "files"], "...")
  asys <- with(param, tximport(files, type="alevin"))["counts"] # stay as list
  stats <- with(param, read.table(stats, head=TRUE, row=1))
  rr <- with(param, .rowRanges_from_gtf(gtf)[rownames(asys$counts)])
  cd <- DataFrame(stats)[colnames(asys$counts),]
  txi <- SingleCellExperiment(asys, rowRanges=rr, colData=cd)
  colnames(txi) <- 
    paste(qname, sapply(strsplit(colnames(txi),"\\-"), `[`, 1), sep=sep)
  return(txi)

}


# helper fn
.rowRanges_from_gtf <- function(gtf) { 

  gxs <- subset(rtracklayer::import(gtf), type=="gene")
  names(gxs) <- gxs$gene_id
  granges(gxs) 

}


# helper fn
.do_droplet_QC <- function(txis, verbose=FALSE) { 

  if (verbose) message("Starting QC...")
  qcstats <- scuttle::perCellQCMetrics(txis)
  qcfilter <- scuttle::quickPerCellQC(qcstats)
  txis <- txis[, !qcfilter$discard]
  if (verbose) message(length(qcfilter$discard), " cells dropped")
  txis <- scuttle::logNormCounts(txis)
  colData(txis)$QCed <- TRUE
  return(txis)

}
