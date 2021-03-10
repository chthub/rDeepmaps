#' Get all Idents names from object
#'
#' @return
#' @export
get_all_idents <- function() {
  return(colnames(e1$obj@meta.data))
}

#' Get all gene names/ids from object
#'
#' @return
#' @export
get_all_genes <- function() {
  return(rownames(e1$obj))
}
