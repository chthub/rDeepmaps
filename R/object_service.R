#' Get all Idents names from object
#'
#' @return array all idents names
#' @export
get_all_idents <- function() {
  all_idents <- list()
  for (i in seq_along(colnames(e1$obj@meta.data))) {
    this_ident <- colnames(e1$obj@meta.data)[i]
    pattern <- c(
      "orig.ident",
      "nCount_*",
      "nFeature_*",
      "percent.*",
      "pct.*",
      "gex*",
      "atac*",
      "prediction*",
      "nucleosome*",
      "TSS.enrichment",
      "TSS.percentile",
      "high.tss",
      "SCT.weight",
      "blacklist*"
    )
    matches <- unique(grep(paste(pattern, collapse = "|"),
                           this_ident, value = F))
    if (length(matches) == 0) {
      this_levels <- levels(as.factor(e1$obj@meta.data[, i]))
      tmp_list <- list(ident = this_ident,
                       levels = this_levels)
      all_idents <- rlist::list.append(all_idents, tmp_list)
    }
  }
  all_idents
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

#' Get all assays names from object
#'
#' @return array all assays names
#' @export
get_all_assays <- function() {
  return(list(
    assay_idx = e1$assay_idx,
    all_assays = names(e1$obj@assays)
  ))
}

#' Set Seurat Assay by name
#' @param req request payload
#' @param name string idents name
#' @return array levels of new idents
#' @export
#'
set_assay <- function(req, name = "RNA") {
  e1$assay_idx <- which(names(e1$obj@assays) == name)
  return(list(
    assay_idx = e1$assay_idx,
    all_assays = names(e1$obj@assays)
  ))
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
