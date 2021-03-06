#' IRIS3
#'
#' @slot seurat
#' @slot upload_info UploadInfo
setClass(
  "IRIS3",
  slots = c(
    seurat = "Seurat",
    upload_info = "ANY"
  )
)


#' Create UploadInfo object
#'
#' @slot type ANY.
#' @slot expr_filename ANY.
#' @slot label_filename ANY.
setClass(
  "UploadInfo",
  slots = c(
    type = "character",
    expr_filename = "character",
    label_filename = "character"
  )
)


setMethod(
  f = "show",
  signature = "IRIS3",
  definition = function(object) {
    cat(format("Print IRIS3 object"), "\n")
    # invisible(x = NULL)
  }
)
