## code to prepare `DATASET` dataset goes here

#' load_zeisel_2015
#' @description Zeisel 2015: Mouse Brain, 7 cell types, 3005 cells
#' @return NULL
#'

load_zeisel_2015 <- function() {
  expr <- fst::read_fst("data-raw/Zeisel_expression.fst")
  meta <- readr::read_csv("data-raw/Zeisel_index_label.csv")
  rownames(expr) <- NULL
  rownames(meta) <- NULL
  expr <- tibble::column_to_rownames(expr, "X1")
  meta <- tibble::column_to_rownames(meta, "Cell")
  zeisel_2015 <- list(expr = expr, meta = meta)
  usethis::use_data(zeisel_2015, overwrite = TRUE)
}


#' load_yan_2013
#' @description Yan 2013: Human embryo, 7 cell types, 90 cells
#' @return null
#'
load_yan_2013 <- function() {
  expr <- readr::read_csv("data-raw/Yan_2013_expression.csv")
  meta <- readr::read_csv("data-raw/Yan_2013_label.csv")
  rownames(expr) <- NULL
  rownames(meta) <- NULL
  expr <- tibble::column_to_rownames(expr, "Gene_ID")
  meta <- tibble::column_to_rownames(meta, "Cell_type")
  yan_2013 <- list(expr = expr, meta = meta)
  usethis::use_data(yan_2013, overwrite = TRUE)
}
