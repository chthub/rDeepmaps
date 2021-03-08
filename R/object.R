#' Create UploadInfo Object
#' @name create_uploadinfo_object
#' @title create_uploadinfo_object
#' @description Create IRIS3 UploadInfo object
#' @param upload_info upload data information
#' @return a S4 object.
#' @export
#'
create_uploadinfo_object <-
  function(upload_info = list()) {
    message("Creating IRIS3 object. \n")
    obj <-
      new(
        Class = "UploadInfo",
        type = "RNA",
        expr_filename = "expr.csv",
        label_filename = "label.csv"
      )
    return(obj)
  }
