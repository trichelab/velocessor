#' plot compartments from samples A & B, perhaps by diffing matrices
#'
#' FIXME: use an exponential smoother instead of the mean smoother
#' 
#' @param   A     a GRanges, usually from compartmap or get_AB
#' @param   B     a GRanges, usually from compartmap or get_AB (but see Details)
#' @param   diff  plot their differences, too? (TRUE, and again, see Details)
#' @param   with  plot with confidence intervals? (TRUE; yet again, see Details)
#' @param   what  what assay fed the original A and B calls? (also, see Details)
#' @param   steps damping step if using mean smoother -- default is 1, no smooth
#' @param   nameA name for A in the plots ("A")
#' @param   nameB name for B in the plots ("B")
#' @param   mat   difference correlation matrices instead? (FALSE, experimental)
#' @param   ...   more parameters to feed to diff_AB
#'
#' @details 
#' There are two nonstandard options here, neither of which should be used just
#' yet (if it breaks, you get to keep the pieces, but I'm happy to accept PRs). 
#' The first is if `what` is set to "ratio" -- the assumption is that we are
#' looking at "chromatin potential" and, accordingly, the function will assume
#' that `A`, `B`, or both were generated from both spliced and unspliced counts.
#' The second is if both `diff` and `with` are specified; then the CIs will be
#' treated as upper and lower bounds for the differences in a third panel. The
#' behavior of this function is undefined if `diff` and `with` are TRUE along 
#' with `what` being set to "ratio".  The behavior of this function is also 
#' undefined if `mat` is set to true, not because it's particularly hard, but
#' because the difference detection is a work in progress. 
#' 
#' @seealso       compartmap::plotAB
#'
#' @import        compartmap
#'
#' @export
plot_ABC <- function(A, B, diff=TRUE, with=TRUE, what="counts", steps=1L, nameA="A", nameB="B", ...){

  if (what == "ratio") stop("This feature is not yet implemented.") # will be
  ABC <- diff_AB(A, B, what=what, steps=steps)

  if (mat) { 
    
    stop("Matrix plotting is not yet implemented.") 

  } else { 

    ABC <- diff_AB(A, B, what=what, ...) 
    rows <- ifelse(diff, 2, 3) 
    columns <- ifelse(what == "ratio", 
                      ifelse(steps < 2, 1, 2), 
                      ifelse(steps < 2, 2, 4)) 

    par(mfrow=c(rows, columns))

    plotAB(ABC, what=nameA, main=paste0(nameA, ", undamped"))
    plotAB(ABC, what=nameB, main=paste0(nameB, ", undamped"))

    if (diff) plotAB(ABC, what="diff", 
                     main=paste0(nameA, " - ", nameB, ", undamped"))

    plotAB(ABC, what=paste0("damped", nameA), main=paste0(nameA, ", damped"))
    plotAB(ABC, what=paste0("damped", nameB), main=paste0(nameB, ", damped"))

    if (diff) plotAB(ABC, what="dampdiff", 
                     main=paste0(nameA, " - ", nameB, ", damped"))

    message("FIXME: add confidence bounds for A, B, and A-B!") 

  }

}
