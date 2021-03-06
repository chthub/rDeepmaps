## code to prepare `DATASET` dataset goes here

#' load_zeisel_2015
#' @description Zeisel 2015: Mouse Brain, 7 cell types, 3005 cells
#' @return NULL
#'

load_zeisel_2015 <- function() {
  library(fst)
  expr <- read_fst('data-raw/Zeisel_expression.fst')
  meta <- read_csv('data-raw/Zeisel_index_label.csv')
  zeisel_2015 <- list(expr=expr, meta=meta)
  usethis::use_data(zeisel_2015, overwrite = TRUE)
}


#' load_yan_2013
#' @description Yan 2013: Human embryo, 7 cell types, 90 cells
#' @return null
#'
load_yan_2013 <- function() {
  expr <- read_csv('data-raw/Yan_2013_expression.csv')
  meta <- read_csv('data-raw/Yan_2013_label.csv')
  yan_2013 <- list(expr=expr, meta=meta)
  usethis::use_data(yan_2013, overwrite = TRUE)
}


