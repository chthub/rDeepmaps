

#' Get dimensional reduction information associated with an assay
#'
#' @param object A \code{Seurat} object
#' @param assay Name of assay that dimensional reduction objects should be
#' associated with
#' @param global Include global dimensional reductions
#'
#' @return A vector of dimensional reduction names
#'
#' @keywords internal
#' @author Paul Hoffman
#' @references https://github.com/satijalab/seurat-wrappers/wiki/Code-Guidelines
#'
AssociatedDimReducs <- function(object,
                                assay = DefaultAssay(object = object),
                                global = TRUE) {
  return(Filter(
    f = function(x) {
      check <- DefaultAssay(object = object[[x]]) == assay
      if (global) {
        check <- c(check, IsGlobal(object = object[[x]]))
      }
      return(any(check))
    },
    x = Reductions(object = object)
  ))
}

#' Find the default DimReduc
#'
#' Searches for DimReducs matching 'umap', 'tsne', or 'pca',
#' case-insensitive, and in that order. Priority given to
#' DimReducs matching the DefaultAssay or assay specified
#' (eg. 'pca' for the default assay weights higher than
#' 'umap' for a non-default assay)
#'
#' @param object A Seurat object
#' @param assay Name of assay to use;
#' defaults to the default assay of the object
#'
#' @return The default DimReduc, if possible
#' @author Paul Hoffman
#' @references https://github.com/satijalab/seurat-wrappers/wiki/Code-Guidelines
#'
DefaultDimReduc <- function(object, assay = NULL) {
  assay <- set_if_null(assay, DefaultAssay(object = object))
  drs.use <- c("umap", "tsne", "pca")
  dim.reducs <- Reductions(object = object)
  drs.assay <- Filter(
    f = function(x) {
      return(DefaultAssay(object = object[[x]]) == assay)
    },
    x = dim.reducs
  )
  if (length(x = drs.assay) > 0) {
    index <- lapply(
      X = drs.use,
      FUN = grep,
      x = drs.assay,
      ignore.case = TRUE
    )
    index <- Filter(f = length, x = index)
    if (length(x = index) > 0) {
      return(drs.assay[min(index[[1]])])
    }
  }
  index <- lapply(
    X = drs.use,
    FUN = grep,
    x = dim.reducs,
    ignore.case = TRUE
  )
  index <- Filter(f = length, x = index)
  if (length(x = index) < 1) {
    stop(
      "Unable to find a DimReduc matching one of '",
      paste(drs.use[1:(length(x = drs.use) - 1)], collapse = "', '"),
      "', or '",
      drs.use[length(x = drs.use)],
      "', please specify a dimensional reduction to use",
      call. = FALSE
    )
  }
  return(dim.reducs[min(index[[1]])])
}
