#' Fit mixtures. Turns out this was useful outside of downsampling, so...
#' 
#' Fit a Gaussian mixture model to x and log1p(x) with 1-2 components (default).
#' Whichever has the lowest BIC wins. Nothing fancy! If nothing wins, then put 
#' all of the observations into the same group and leave them all eligible. For
#' prescreening of potential cluster-by-sample combinations, see find_eligible.
#'
#' @seealso     find_eligible
#' 
#' @param x     the thing to fit a 1D mixture to 
#' @param G     vector of component numbers (default is 1:2)
#' @param log   take log1p of x? (TRUE) 
#' @param keep  keep the model fits as an attribute? (FALSE) 
#' 
#' @return      list of fits, length 1 if log == FALSE 
#' 
#' @import      mclust
#' 
#' @export
fit_mixtures <- function(x, G=1:2, log=TRUE, keep=FALSE) { 

  cellfits <- apply(x, 1, densityMclust, G=G, verbose=FALSE)
  
  if (log) logx <- log1p(x)
  if (log) logcellfits <- apply(logx, 1, densityMclust, G=1:2, verbose=FALSE)

  # possible model names
  modelnames <- c("raw")
  if (log) modelnames <- c("raw", "log")
  names(modelnames) <- modelnames

  # comparisons 
  comps <- data.frame(raw=unname(sapply(cellfits, `[[`, "G")))
  if (log) comps$log <- unname(sapply(logcellfits, `[[`, "G"))
  bics <- data.frame(raw=unname(sapply(cellfits, `[[`, "bic")))
  if (log) bics$log <- unname(sapply(logcellfits, `[[`, "G"))
  rownames(comps) <- rownames(bics) <- names(cellfits)

  # the fits
  fits <- list() 
  for (cl in rownames(x)) {
    if (log) fits[[cl]] <- list(raw=cellfits[[cl]], log=logcellfits[[cl]])
    else fits[[cl]] <- list(raw=cellfits[[cl]])
  }

  # default just returns 1 
  model <- rep("none", nrow(x))
  names(model) <- rownames(x) 
  res <- matrix(0, nrow=nrow(x), ncol=ncol(x)) 
  dimnames(res) <- dimnames(x) 
  for (cl in rownames(comps)) {
    if (all(comps[cl, ] > 1)) {
      model[cl] <- modelnames[which.max(bics[cl,])] 
      res[cl, ] <- (fits[[cl]][[model[cl]]][["classification"]] - 1)
    } else if (any(comps[cl, ] > 1)) {
      model[cl] <- modelnames[which.max(comps[cl,])] 
      res[cl, ] <- (fits[[cl]][[model[cl]]][["classification"]] - 1)
    } else { 
      model[cl] <- "none"
      res[cl, ] <- rep(1, ncol(res))
    }
  }

  if (keep) attr(res, "fits") <- fits
  attr(res, "model") <- model
  return(res)

}
