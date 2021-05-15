#' Run Seurat clustering
#'
#' @param req request payload
#' @param ident1 string
#' @param ident2 string
#' @param min_pct string
#' @param min_lfc string
#' @param min_lfc string
#' @param assay string
#' @param pvalue string
#'
#' @return
#' @export
#'
calc_deg <-
  function(req,
           ident1 = 4,
           ident2 = 5,
           min_pct = 0.2,
           min_lfc = 0.5,
           assay = "RNA",
           pvalue = 0.05) {
    send_progress(paste0("Running differential gene expression analysis"))
    Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]
    this_markers <-
      FindMarkers(
        e1$obj,
        assay = assay,
        ident.1 = ident1,
        ident.2 = ident2,
        min.pct = min_pct,
        logfc.threshold = min_lfc
      ) %>%
      dplyr::filter(p_val_adj < pvalue)
    this_markers <- tibble::rownames_to_column(this_markers, "gene")
    e1$deg <- this_markers
    return(list(this_markers))
  }


#' Run Seurat clustering
#'
#' @param req request payload
#' @param ident1 string
#' @param ident2 string
#' @param min_pct string
#' @param min_lfc string
#' @param min_lfc string
#' @param assay string
#' @param pvalue string
#'
#' @return
#' @export
#'
calc_dr <-
  function(req,
           ident1 = 4,
           ident2 = 5
           ) {

    return(list(1))
  }

