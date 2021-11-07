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
#' @param direction string
#'
#' @return
#' @export
#'
calc_deg <-
  function(req,
           ident1 = c(0, 1),
           ident2 = c(2, 3),
           min_pct = 0.2,
           min_lfc = 0.5,
           assay = "RNA",
           pvalue = 0.05,
           direction = 'all') {
    send_progress(paste0("Running differential gene expression analysis"))
    Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]

    if(direction == 'up') {
      this_markers <-
        FindMarkers(
          e1$obj,
          ident.1 = ident1,
          ident.2 = ident2,
          min.pct = min_pct,
          assay = "RNA",
          only.pos = T,
          logfc.threshold = min_lfc
        ) %>%
        dplyr::filter(p_val_adj < pvalue) %>%
        dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
        tibble::rownames_to_column("gene")
    }
    if(direction == 'down') {
      this_markers <-
        FindMarkers(
          e1$obj,
          ident.1 = ident1,
          ident.2 = ident2,
          min.pct = min_pct,
          assay = "RNA",
          logfc.threshold = min_lfc
        ) %>%
        dplyr::filter(p_val_adj < pvalue) %>%
        dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::filter(avg_log2FC < 0)
    }
    if(direction == 'all') {
      this_markers <-
        FindMarkers(
          e1$obj,
          ident.1 = ident1,
          ident.2 = ident2,
          min.pct = min_pct,
          assay = "RNA",
          logfc.threshold = min_lfc
        ) %>%
        dplyr::filter(p_val_adj < pvalue) %>%
        dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
        tibble::rownames_to_column("gene")
    }
    this_markers <- this_markers %>%
      dplyr::rename(pct_1=pct.1, pct_2=pct.2)
    e1$deg <- this_markers
    return(this_markers)
  }

