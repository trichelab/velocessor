#' automates fetching of Valentine Svenssen's single-cell database
#' 
#' scdb() tidies up the database a bit, with Date, URL, Cells, and Velocity
#'
#' @details
#' The database itself lives at http://www.nxn.se/single-cell-studies/data.tsv
#' 
#' @references 
#' Svensson V, Beltrame EdV, Pachter L (2020), A curated database reveals 
#' trends in single-cell transcriptomics, Database, 28 November 2020
#' https://doi.org/10.1093/database/baaa073
#'
#' @param URL   the URL if it changes 
#' @param ...   anything else to tell utils::read.delim
#' 
#' @return      a data.frame, with the Date column as an actual Date
#'
#' @examples 
#' 
#'   db <- scdb()
#'
#'   library(useful)
#'   library(tidyverse)
#'   subset(db, Velocity == TRUE) %>% corner
#'
#'   columns <- c("Title", "Data.location", "URL")
#'   subset(db, Tissue == "Blood" & Cells > 50000)[, columns]
#'   subset(db, grepl("CITE", Technique))[, columns]
#'   subset(db, grepl("Amit", Authors))[, columns]
#' 
#'   a_year_ago <- (Sys.Date() - 365)
#'   subset(db, Date > a_year_ago) %>% corner
#'
#'   five_years_ago <- (Sys.Date() - (5*365.25))
#'   p1 <- ggplot(subset(db, Date > five_years_ago), 
#'                aes(x = Date, fill = Velocity)) + 
#'     geom_histogram() + 
#'     scale_x_date() + 
#'     xlab("year") + 
#'     scale_y_log10() + 
#'     ylab("papers") + 
#'     NULL 
#' 
#'   suppressWarnings(p1 + theme_minimal() + ggtitle("Velocity in papers"))
#' 
#'   p2 <- ggplot(subset(db, 
#'                       Organism %in% c("Human", "Mouse", "Zebrafish") & 
#'                       Date > five_years_ago), 
#'                aes(x=Date, y=Cells, color=Organism, fill=Organism)) + 
#'     geom_point() + 
#'     geom_smooth() + 
#'     scale_y_log10() + 
#'     NULL
#' 
#'   suppressWarnings(p2 + facet_grid(vars(Organism)) + theme_minimal())
#' 
#' @importFrom utils read.delim
#' 
#' @export 
scdb <- function(URL="http://www.nxn.se/single-cell-studies/data.tsv", ...) {

  db <- read.delim(URL, ...)
  db$URL <- with(db, paste0("http://dx.doi.org/", DOI))
  db$Date <- with(db, as.Date(as.character(Date), "%Y%m%d"))
  db$Cells <- with(db, as.integer(gsub(",", "", Reported.cells.total)))
  db$Velocity <- with(db, RNA.Velocity == "Yes")
  return(db)

}
