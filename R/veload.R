#' world's dumbest velocity loader, currently only for Alevin output 
#' 
#' @param dir        where the output folders are (".")
#' @param suffix     what the suffix of an output folder is ("_alevin")
#' 
#' @export
veload <- function(dir=".", suffix="_alevin", ...) {

  if (dir != ".") {
    stopifnot(dir.exists(dir))
    oldwd <- getwd()
    setwd(dir)
  }

  runs <- list.files(pattern=paste0(suffix, "$"))
  names(runs) <- sub(suffix, "", runs)
  txome <- list.files(patt=".*annotation.*json$")
  stopifnot(file.exists(txome))
  txstub <- sub("\\.json$", "", txome)
  res <- process_velo_txis(runs, txstub, ...)

  if (dir != ".") setwd(oldwd)
  return(res)

}
