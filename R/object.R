#' CreateIRIS3Object
#' @name CreateIRIS3Object
#' @title CreateIRIS3Object
#' @description Create IRIS3 object
#' @param x Seurat, a seurat object
#' @param upload_info upload data information
#' @return it should return a IRIS3 S4 object.
#' @export
#'
#' @examples
#' data("yan_2013")
#' seurat_obj <- Seurat::CreateSeuratObject(yan_2013$expr)
#' object <- CreateIRIS3Object(seurat_obj)
# globalVariables(c("input_matrix"))
CreateIRIS3Object <-
  function(x = seurat_object,
           upload_info = new(Class = "UploadInfo")) {
    message("Creating IRIS3 object. \n")
    obj <-
      new(Class = "IRIS3",
          seurat = x,
          upload_info = upload_info)
    # obj <- suppressMessages(AddMeta(obj))
    return(obj)
  }
