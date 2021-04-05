#' Run GSEA enrichment
#'
#' @param req request payload
#' @param genes string
#' @param database string
#'
#' @return
#' @export
#'
calc_gsea_table <-
  function(req,
           genes,
           database) {
    print("run GSEA")
    Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]

    library(fgsea)
    library(msigdbr)
    m_df = msigdbr(species = "Mus musculus", category = "H")
    m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

    res <- as.numeric(seq_along(res))
    names(res) <- VariableFeatures(e1$obj)
    fgseaRes <- fgsea(pathways = m_list,
                      stats = res,
                      nperm = 1000)
    gseaTable <- fgseaRes %>%
      tibble::as_tibble() %>%
      dplyr::arrange(desc(NES)) %>%
      dplyr::select(-leadingEdge,-ES,-nMoreExtreme) %>%
      dplyr::arrange(padj) %>%
      tibble::as_data_frame()
    return(list(gseaTable))
  }

#' Generate GSEA plot
#'
#' @param req request payload
#' @param genes string
#' @param database string
#'
#' @return
#' @export
#'
calc_gsea_plot <-
  function(req,
           term) {
    print("run GSEA")
    library(fgsea)
    library(msigdbr)
    term <- examplePathways[["5991130_Programmed_Cell_Death"]]
    plot1 <-
      plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
                     exampleRanks) + ggplot2::labs(title = "Programmed Cell Death")
    return(plot1)
  }
