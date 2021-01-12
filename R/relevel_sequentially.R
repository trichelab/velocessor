#' convenience function for reordering factors without thinking very hard 
#'
#' I wrote this a million years ago for TARGET and have used it ever since.
#' It's not terribly clever but it can help A LOT with plotting to reorder.
#' 
#' @param x     a factor, or something that can become one
#' @param lvls  a vector with level names that include all levels of x
#' 
#' @return      a releveled factor version of x
#'
#' @export
relevel_sequentially <- function(x, lvls) { 

  if (!is(x, "factor")) {
    message("First argument `x` is not a factor; turning it into one.")
    x <- factor(x) 
  }

  if (length(lvls) < nlevels(x) | !all(levels(x) %in% lvls)) { 
    stop("Not all factor levels of `x` are contained within `lvls`. Exiting.")
  }

  lvls <- intersect(lvls, levels(x))
  for(y in rev(lvls)) x <- relevel(x, which(levels(x) == y))
  return(x)

}
