#' massage a (possibly velocitized) SingleCellExperiment into an infercnv object
#'
#' @name    infercnv 
#' @rdname  infercnv
#' 
#' @details
#' We do not require or import infercnv, because the requirement for JAGS can 
#' be a showstopper for HPC installations. This function will stop() if the
#' infercnv package namespace cannot be attached, however, for obvious reasons. 
#' I have no idea why the `infercnv` function name was not taken by `infercnv`, 
#' but since it wasn't, I stole it. If `NumGenesExpressed` is among the colData
#' columns for the SingleCellExperiment, and `run`==TRUE, then this function 
#' will use the former to guess a sensible value for the latter. The velocessor
#' package also registers a halfway decent `show` method for infercnv objects.
#' 
#' @param   x           a SingleCellExperiment
#' @param   group_col   the colData column with cell annotations ("cnv_group")
#' @param   ref_prefix  prefix for reference sample clusters ("normal_")
#' @param   obs_prefix  prefix for tumor sample clusters ("tumor_")
#' @param   downsample  downsample? (TRUE, `maxcells` per cluster per sample)
#' @param   maxcells    how many cells per sample per cluster to grab? (100) 
#' @param   run         perform infercnv run with "sensible" defaults? (TRUE)
#' @param   cutoff      expression cutoff for infercnv::run (default: guess)
#' @param   ...         other arguments for infercnv::run, e.g. `num_threads`
#'
#' @return              an infercnv object
#' 
#' @import              GenomeInfoDb
#'
#' @export
infercnv <- function(x, group_col="cnv_group", ref_prefix="normal_", obs_prefix="tumor_", downsample=TRUE, maxcells=100, run=TRUE, cutoff=NULL, ...) { 

  # needed for inferring CNV (duh)
  if (!group_col %in% names(colData(x))) {
    stop("You need to specify groups for infercnv in colData(x)$", group_col)
  } else if (!any(grepl(ref_prefix, colData(x)[, group_col]))) {
    stop("You need some normals (", ref_prefix, "*) in colData(x)$", group_col)
  } else if (!any(grepl(obs_prefix, colData(x)[, group_col]))) { 
    stop("You need some tumors (", obs_prefix, "*) in colData(x)$", group_col)
  }
  if (!requireNamespace("infercnv")) {
    message("Cannot create infercnv objects without having installed infercnv.")
    stop("You can install infercnv by calling BiocManager::install('infercnv')")
  } 

  # needed for my sanity
  seqlevelsStyle(x) <- "UCSC"
  stopifnot(length(unique(genome(x))) > 0) 
  labelpat <- paste0("(", ref_prefix, "|", obs_prefix, ")")
  labeled <- which(grepl(labelpat, colData(x)[, group_col]))
  x <- sort(keepStandardChromosomes(x[, labeled]))
  
  # now we should have all we need to proceed:
  origcells <- ncol(x)
  origgenes <- nrow(x)
  if (downsample) x <- downsample_txis(x, maxcells=maxcells, ret="sce")
  gene_order <- as.data.frame(rowRanges(x))[, c("seqnames","start","end")]
  anno_df <- as.data.frame(colData(x)[, group_col, drop=FALSE])
  ref_groups <- unique(grep(ref_prefix, colData(x)[, group_col], value=TRUE)) 
  
  res <- infercnv::CreateInfercnvObject(raw_counts_matrix=counts(x),
                                        gene_order_file=gene_order,
                                        annotations_file=anno_df,
                                        ref_group_names=ref_groups)

  res@options$origcells <- origcells
  res@options$origgenes <- origgenes
  res@options$downsample <- downsample 
  res@options$maxcells <- maxcells

  if (run) {

    if (is.null(cutoff)) { 
      if ("NumGenesExpressed" %in% names(colData(x))) { 
        if (median(colData(x)[, "NumGenesExpressed"], na.rm=TRUE) > 5000) { 
          message("This looks like plate-seq data, setting `cutoff` to 1...")
          # 10 ** 0
          cutoff <- 1
        } else { 
          message("This looks like droplet data, setting `cutoff` to 0.1...")
          # 10 ** -1
          cutoff <- 0.1
        }
      } else { 
        # Split the difference -- 10**-0.5
        message("No cutoff provided, `NumGenesExpressed` missing, use default:")
        # 10 ** -0.5
        cutoff <- 0.3162 
      }
    }
  
    # report findings
    message("Genes with mean(counts) < ", cutoff, " will be dropped.")
    message("InferCNV output will be written to:")
    message("  ", file.path(getwd(), "infercnv_output"))
    res <- infercnv::run(res, 
                         cutoff=cutoff,
                         plot_steps=FALSE, 
                         prune_outliers=TRUE, 
                         cluster_by_groups=TRUE,
                         out_dir="infercnv_output", 
                         noise_logistic=TRUE,
                         sd_amplifier=3,
                         denoise=TRUE, 
                         HMM=TRUE, 
                         HMM_type="i6", 
                         HMM_report_by="cell",
                         ...)

    # apply median filter
    message("Applying median filter to smooth output...")
    res <- infercnv::apply_median_filtering(res)
    message("(See github.com/broadinstitute/inferCNV/wiki/De-noising-Filters )")

  }
  
  return(res)

}


# helper for .show_infercnv
.countcells <- function(name, cells) paste0(name, ": ", cells, " cells")


# infercnv lacks a `show` method...!
.show_infercnv <- function(object) {
   
  # generic information
  genes <- nrow(object@gene_order)
  chroms <- length(unique(object@gene_order[,1]))
  refs <- sapply(object@reference_grouped_cell_indices, length) 
  refcounts <- mapply(.countcells, names(refs), refs)
  obs <- sapply(object@observation_grouped_cell_indices, length)
  obscounts <- mapply(.countcells, names(obs), obs)
  cellgroups <- c(refs, obs)
  cells <- sum(cellgroups) 
  
  # summarize the above:
  dims <- paste(genes, "genes on", chroms, "chromosomes from", cells, "cells")

  # generic options
  opt <- with(object@options, 
              list(mode=paste0(analysis_mode, ", ",
                               ifelse(HMM == "TRUE", HMM_type, "no"), " HMM",
                               ifelse(denoise == "TRUE", ", denoised", "")),
                   percell=paste(min_max_counts_per_cell[1],"reads/cell"),
                   pergene=paste(cutoff, "reads/gene/cell")))

  # velocessor-specific options
  origdims <- ""
  extras <- c("origcells", "origgenes", "downsample")
  if (all(extras %in% names(object@options))) {
    origdims <- with(object@options, 
                     paste(ifelse(as.logical(downsample) == TRUE, 
                                  "(downsampled", "("), "from", 
                           origgenes, "genes on", origcells, "cells"))
  }

  # print a tidy summary of `object` 
  cat("An", class(object), "object with", dims, "\n")
  if (nchar(origdims) > 0) cat(origdims, "\n")
  cat("Cutoffs:", opt$percell, "and", opt$pergene, "minimum\n")
  cat("Analysis mode:", opt$mode, "\n\n")
  S4Vectors::coolcat("Reference groups (%d): %s\n", names(refs))
  for (ref in names(refcounts)) cat("  ", refcounts[ref], "\n")
  S4Vectors::coolcat("Observation groups (%d): %s\n", names(obs))
  for (obs in names(obscounts)) cat("  ", obscounts[obs], "\n")
  cat("Output path:", object@options$out_dir, "\n\n")

}


#' @export
suppressWarnings(setMethod("show", "infercnv", .show_infercnv))
