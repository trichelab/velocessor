#' Gencode releases do not map 1:1 onto ENSEMBL releases, hence this function.
#'
#' We don't even try to be complete because it's a mess before gencode 19. 
#' This ought to use a table in system.file("data", package="velocessor").
#' 
#' @param release   the Gencode release (e.g. M24 or 30) 
#'
#' @return          the ENSEMBL version (e.g. 99 or 96)
#' 
#' @seealso fix_gencode
#' 
#' @export
get_gencode_ensembl <- function(release) { 

         #   human     mouse
  ens <- c("19"=75, 
           "20"=76, 
           "22"=79, 
           "23"=81, 
           "24"=84,  "M9"=84,
           "25"=85,           
                     "M11"=86,
           "26"=88,           
                     "M14"=89,
           "27"=90,  "M15"=90,
                     "M16"=91,
           "28"=92,  "M17"=92,
                     "M18"=93,
           "29"=94,  "M19"=94,
                     "M20"=95,
           "30"=96,  "M21"=96,
           "31"=97,  "M22"=97,
           "32"=98,  "M23"=98,
           "33"=99,  "M24"=99,
           "34"=100, "M25"=100,
           "35"=101            )

  return(ens[fix_gencode(release)])

}
