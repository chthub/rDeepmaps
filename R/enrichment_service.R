#' Run GSEA enrichment
#'
#' @param req request payload
#' @param genes array
#' @param database string
#'
#' @return
#' @export
#'
calc_gsea_table <-
  function(req,
           genes = c("CD74", "CD7"),
           database = "C2") {
    print("run GSEA")

    library(fgsea)
    library(msigdbr)
    if (e1$species == "Human") {
      this_species <- "Homo sapiens"
    } else {
      this_species <- "Mus musculus"
    }
    print(this_species)
    m_df <- msigdbr(species = this_species, category = database)
    m_list <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

    res <- e1$deg %>%
      dplyr::select(gene, avg_log2FC) %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      na.omit() %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene) %>%
      tibble::deframe() %>%
      sort(decreasing = T)

    fgseaRes <- fgsea(
      pathways = m_list,
      stats = res,
      nperm = 500
    )
    gseaTable <- fgseaRes %>%
      tibble::as_tibble() %>%
      dplyr::arrange(desc(NES)) %>%
      dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
      dplyr::arrange(padj) %>%
      tibble::as_data_frame()
    return(list(gseaTable))
  }

#' Run GSEA enrichment
#'
#' @param req request payload
#' @param genes array
#' @param database string
#'
#' @return
#' @export
#'
regulon_enrichment <-
  function(req,
           genes = c("CD74", "CD7"),
           database = "C2") {
    library(enrichR)
    if (length(e1$species == "Human")) {
      this_species <- "Homo sapiens"
    } else {
      this_species <- "Mus musculus"
    }
    message(this_species)
    if(database == "KEGG" && this_species == 'Homo sapiens') {
      dbs <- "KEGG_2019_Human"
    } else if (database == "KEGG" && this_species == 'Mus musculus'){
      dbs <- "KEGG_2019_Mouse"
    } else {
      dbs <- database
    }
    message(dbs)
    enriched_combined <- enrichr(genes,dbs)
    return(list(enriched_combined[[1]][,c(-3,-5,-6,-7)]))
  }

#' Generate GSEA plot
#'
#' @param req request payload
#' @param term string
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
      plotEnrichment(
        examplePathways[["5991130_Programmed_Cell_Death"]],
        exampleRanks
      ) + ggplot2::labs(title = "Programmed Cell Death")
    return(plot1)
  }
