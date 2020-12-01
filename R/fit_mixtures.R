#' Fit mixtures. Turns out this was useful outside of downsampling, so...
#' 
#' Fit a Gaussian mixture model to x and log1p(x) with 1-2 components (default).
#' Whichever has the lowest BIC wins. Nothing fancy! If nothing wins, then put 
#' all of the observations into the same group and leave them all eligible. For
#' prescreening of potential cluster-by-sample combinations, see find_eligible.
#'
#' @seealso find_eligible
#' 
#' @param x the thing to fit a 1D mixture to 
#' @param G a vector of component numbers (default is 1:2)
#' 
#' @return  a list of fits
#' 
#' @import  mclust
#' 
#' @export
fit_mixtures <- function(x, G=1:2) { 

  logx <- log1p(x)
  cellfits <- apply(x, 1, densityMclust, G=G, verbose=FALSE)
  logcellfits <- apply(logx, 1, densityMclust, G=1:2, verbose=FALSE)

  modelnames <- c("raw", "log")
  names(modelnames) <- modelnames
  comps <- data.frame(raw=unname(sapply(cellfits, `[[`, "G")),
                      log=unname(sapply(logcellfits, `[[`, "G")))
  bics <- data.frame(raw=unname(sapply(cellfits, `[[`, "bic")),
                     log=unname(sapply(logcellfits, `[[`, "bic")))
  rownames(comps) <- rownames(bics) <- names(cellfits)

  fits <- list() 
  for (cl in rownames(x)) {
    fits[[cl]] <- list(raw=cellfits[[cl]], log=logcellfits[[cl]])
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
  attr(res, "model") <- model
  return(res)

}
