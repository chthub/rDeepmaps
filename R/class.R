#' Create UploadInfo object
#'
#' @slot type character
#' @slot expr_filename character
#' @slot label_filename character
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
  signature = "UploadInfo",
  definition = function(object) {
    cat(format("Print UploadInfo object"), "\n")
  }
)
