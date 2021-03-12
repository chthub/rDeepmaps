#' Run Seurat clustering
#'
#' @param req request payload
#' @param ident1_idx string
#' @param ident2_idx string
#' @param ident1 string
#' @param ident2 string
#' @param min_pct string
#' @param min_lfc string
#'
#' @return
#' @export
#'
calc_deg <-
  function(req,
           ident1_idx = 5,
           ident2_idx = 5,
           ident1 = 4,
           ident2 = 5,
           min_pct = 0.2,
           min_lfc = 0.5) {
    print("run deg")
    Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]
    this_markers <-
      FindMarkers(
        e1$obj,
        ident.1 = ident1,
        ident.2 = ident2,
        min.pct = min_pct,
        logfc.threshold = min_lfc
      )
    this_markers <- tibble::rownames_to_column(this_markers, "gene")
    return(list(this_markers))
  }
