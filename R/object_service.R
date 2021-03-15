#' Get all Idents names from object
#'
#' @return array all idents names
#' @export
get_all_idents <- function() {
  all_idents <- list()
  for (i in seq_along(colnames(e1$obj@meta.data))) {
    this_ident <- colnames(e1$obj@meta.data)[i]
    this_levels <- levels(as.factor(e1$obj@meta.data[, i]))
    tmp_list <- list(
      ident = this_ident,
      levels = this_levels
    )
    all_idents <- rlist::list.append(all_idents, tmp_list)
  }
  return(all_idents)
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
  return(levels(as.factor(e1$obj@meta.data[, e1$ident_idx])))
}


#' Set subset obj or full_obj
#' @param req request payload
#' @param type string idents name
#' @return array levels of new idents
#' @export
#'
set_obj <- function(req, type = "full") {
  if (is.null(e1$full_obj)) {
    e1$full_obj <- e1$obj
  }
  if (is.null(e1$sub_obj)) {
    e1$sub_obj <- e1$obj
  }
  if (type == "full") {
    e1$obj <- e1$full_obj
  }
  if (type == "subset") {
    e1$obj <- e1$sub_obj
  }
  return(type)
}
