#' Get all Idents names from object
#'
#' @return array all idents names
#' @export
get_all_idents <- function() {
  return(colnames(e1$obj@meta.data))
}

#' Get all gene names/ids from object
#'
#' @return array all gene names
#' @export
get_all_genes <- function() {
  return(rownames(e1$obj))
}


#' Set Seurat Idents by name
#' @param req request payload
#' @param name string idents name
#' @return array levels of new idents
#' @export
#'
set_idents <- function(req, name = "orig.ident") {
  e1$ident_idx <- which(colnames(e1$obj@meta.data) == name)
  return(levels(as.factor(
    e1$obj@meta.data[, e1$ident_idx]
  )))
}
