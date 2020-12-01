#' this is kind of a helper function, but useful enough to split out 
#' 
#' Given compartment estimates A and B, compute estimates of their differences.
#'
#' @param A     a GRanges
#' @param B     a GRanges 
#' @param conf  what quantile to use for confidence intervals (.75, quartiles)
#' @param steps how many bins to use for the weighted smoother (1; don't smooth)
#' @param what  what assay fed the original A and B calls? (counts; see Details)
#' @param stubs what stubs to use for the difference variables? (c("A","B"))
#' 
#' @return      a GRanges with differences between A and B and damped bounds
#' 
#' @details
#' `what` is left open for future use comparing unspliced/spliced ratios. It 
#' has no effect upon the execution of this function at the moment. Note that 
#' if using the weighted smoother (`steps` > 1), then `steps` MUST be odd.
#' 
#' @import      GenomicRanges
#' @import      compartmap
#'
#' @export
diff_AB <- function(A, B, conf=0.75, steps=1, what="counts", stubs=c("A","B")) {
                 
  if ((steps %% 2) == 0) stop("You must specify an odd number for `steps`.")
  
  stubA <- stubs[1]
  stubB <- stubs[2]
  ABC <- disjoin(c(A, B))
  ABC <- .addColumns(ABC, .renormalize(A, conf=conf, steps=steps), stub=stubA)
  ABC <- .addColumns(ABC, .renormalize(B, conf=conf, steps=steps), stub=stubB)
  ABC <- .addDiffColumns(ABC, stubA=stubA, stubB=stubB)
  ABC$dampdiffstep <- Rle(steps)

  return(ABC) 

}


# helper fns
.squeeze <- function(p0, sqz=1e-6) ((p0 - 0.5) * (1 - sqz)) + 0.5
.unsqueeze <- function(p1, sqz=1e-6) ((p1 - 0.5) / (1 - sqz)) + 0.5
# all(p0 - .unsqueeze(.squeeze(p0)) < sqz) == TRUE, for all (sqz < 1 & sqz >= 0)


# helper fn
.status <- function(x) c("closed", "open")[ceiling((sign(x) + 2)/2)]


# helper fn -- go back and forth via CDF
.renormalize <- function(A, sqz=1e-6, conf=0.75, steps=3) { 

  stopifnot(sqz < 1)
  stopifnot(conf > 0 & conf < 1)
  stopifnot((steps %% 2) > 0)

  hi <- conf
  lo <- 1 - conf
  A$re <- qnorm(.squeeze((A$score - min(A$score))/(max(A$score)-min(A$score))))
  A$mid <- (pnorm(A$re) - 0.5) * 2 
  A$lower <- (pnorm(qnorm(lo, mean=A$re, sd=sqrt(1/A$conf.est))) - .5) * 2
  A$upper <- (pnorm(qnorm(hi, mean=A$re, sd=sqrt(1/A$conf.est))) - .5) * 2
  A$status <- .status(A$re) # open or closed, pre-damping
  A$weight <- A$conf.est * width(A) # usually width will be constant, but...
  A$damped <- .weighted_smooth(A$mid, w=A$conf.est, span=steps)
  A$dampedlo <- (pnorm(qnorm(lo, mean=A$damped, sd=sqrt(1/A$conf.est))) - .5)*2
  A$dampedhi <- (pnorm(qnorm(hi, mean=A$damped, sd=sqrt(1/A$conf.est))) - .5)*2
  A$dampedstatus <- .status(A$damped) # open or closed, post-damping
  A$coord <- (start(A) + end(A)) / 2

  return(A)

}


# helper fn
.weighted_smooth <- function(x, w, span=3, delta=1e-6) {

  stopifnot(length(x) == length(w))
  k <- span - 2 # see ?compartmap::meanSmoother
  compartmap::meanSmoother(x, k=k, iter=2, delta=delta, w=w, na.rm=TRUE)

}


# helper fn
.addColumns <- function(ABC, A, stub="A") {

  dstub <- paste0("damped", stub)
  nameABC <- c(stub,                    # e.g. "A"
               paste0(stub, "lo"),      # e.g. "Alo"
               paste0(stub, "hi"),      # e.g. "Ahi"
               paste0(stub, "status"),  # e.g. "Astatus"
               dstub,                   # e.g. "dampedA"
               paste0(dstub, "lo"),     # e.g. "dampedAlo"
               paste0(dstub, "hi"),     # e.g. "dampedAhi"
               paste0(dstub, "status")) # e.g. "dampedAstatus"

  nameA <- c("mid", 
             "lower", 
             "upper", 
             "status",
             "damped", 
             "dampedlo", 
             "dampedhi", 
             "dampedstatus")
  names(nameA) <- nameABC
  ol <- findOverlaps(ABC, A) 
  for (i in nameABC) ABC <- .addColumn(ABC, A, ol, i, nameA[i])
  return(ABC) 

}


# helper fn
.addColumn <- function(ABC, A, ol, nameABC, nameA) {

  mcols(ABC)[, nameABC] <- NA
  mcols(ABC[queryHits(ol)])[, nameABC] <- mcols(A[subjectHits(ol)])[, nameA]
  return(ABC)

}


# helper fn
.addDiffColumns <- function(ABC, stubA="A", stubB="B") {

  stopifnot(stubA %in% names(mcols(ABC)))
  dstubA <- paste0("damped", stubA)

  stopifnot(stubB %in% names(mcols(ABC)))
  dstubB <- paste0("damped", stubB)

  ABC$diff <- 
    mcols(ABC)[, stubA] - mcols(ABC)[, stubB]
  ABC$diffmin <- 
    mcols(ABC)[, paste0(stubA, "lo")] - mcols(ABC)[, paste0(stubB, "hi")]
  ABC$diffmax <- 
    mcols(ABC)[, paste0(stubA, "hi")] - mcols(ABC)[, paste0(stubB, "lo")]
  
  ABC$dampdiff <- 
    mcols(ABC)[, dstubA] - mcols(ABC)[, dstubB]
  ABC$dampdiffmin <-
    mcols(ABC)[, paste0(dstubA, "lo")] - mcols(ABC)[, paste0(dstubB, "hi")]
  ABC$dampdiffmax <-
    mcols(ABC)[, paste0(dstubA, "hi")] - mcols(ABC)[, paste0(dstubB, "lo")]

  return(ABC) 

}
