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

    fgseaRes <- fgsea(pathways = m_list,
                      stats = res,
                      nperm = 500)
    gseaTable <- fgseaRes %>%
      tibble::as_tibble() %>%
      dplyr::arrange(desc(NES)) %>%
      dplyr::select(-leadingEdge,-ES,-nMoreExtreme) %>%
      dplyr::arrange(padj) %>%
      tibble::as_data_frame()
    return(gseaTable)
  }

#' Run Enrichr
#'
#' @param req request payload
#' @param genes array
#' @param database string
#'
#' @return
#' @export
#'
calc_enrichr_table <-
  function(req,
           genes = c("CD74", "CD7"),
           database = "KEGG") {
    library(enrichR)
    if (length(e1$species == "Human")) {
      this_species <- "Homo sapiens"
    } else {
      this_species <- "Mus musculus"
    }

    if (database == "KEGG" && this_species == 'Homo sapiens') {
      dbs <- "KEGG_2019_Human"
    } else if (database == "KEGG" &&
               this_species == 'Mus musculus') {
      dbs <- "KEGG_2019_Mouse"
    } else {
      dbs <- database
    }

    enriched_combined <- enrichr(genes, dbs)
    return(enriched_combined[[1]][, c(-3, -5, -6, -7)])
  }

#' Generate GSEA plot
#'
#' @param req request payload
#' @param term string
#' @return
#' @export
#'
plot_gsea <-
  function(req,
           term) {
    library(fgsea)
    library(msigdbr)
    this_idx <- sample.int(length(examplePathways), 1)
    term <- examplePathways[[this_idx]]
    plot1 <-
      plotEnrichment(term,
                     exampleRanks)
    return(print(plot1))
  }

#' Generate enrichment dot plot
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_size
#' @importFrom ggplot2 theme_bw ylab labs theme element_text scale_x_continuous
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_split str_replace_all str_remove_all
#' @param df
#' @param isPvalLog
#' @return

#' @export
#'
plot_enrichr_dot <-
  function(df=mtcars, isPvalLog = "false") {
    library(ggplot2)
    library(forcats)
    library(stringr)
    #write.csv(df,"df.csv")
    df <- read.csv("C:/Users/flyku/Documents/GitHub/iris3api/inst/endpoints/df.csv")
    new_df <- df %>%
      dplyr::mutate(
        Term = as.factor(stringr::str_replace_all(Term, " \\(GO.*", "")),
        len = lengths(str_split(Genes, ";")),
        pval = -log10(Adjusted.P.value)
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(gene_ratio = eval(parse(text = str_remove_all(Overlap, " ")))) %>%
      dplyr::select(Term, len, pval, gene_ratio, Adjusted.P.value) %>%
      dplyr::mutate(Term = fct_reorder(Term, Adjusted.P.value)) %>%
      dplyr::arrange(Adjusted.P.value)

    #new_df$Term <-
    #  factor(new_df$Term, levels = levels(factor(new_df$Term)))

    if (isPvalLog == "true") {
      new_df$Adjusted.P.value <- -1 * log10(new_df$Adjusted.P.value)
      pval_legend <- '-log10(adj. p-value)'
      low_color <- "blue"
      high_color <- "red"
      trans_color <- 'log10'
    } else {
      pval_legend <- 'adj. p-value'
      low_color <- "blue"
      high_color <- "red"
      trans_color <- 'reverse'
    }
    plot_dot <- ggplot(new_df,
                       aes(x = gene_ratio,
                           y = Term)) +
      geom_point(aes(color = Adjusted.P.value, size = len)) +
      scale_color_gradient(low = low_color,
                           high = high_color,
                           trans = trans_color) +
      theme_bw() +
      ylab("") +
      labs(size = "Overlapping count", color = pval_legend) +
      theme(
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)
      ) +
      scale_x_continuous(name = "Overlapping ratio") +
      scale_size(range = c(4, 8))

    return(print(plot_dot))
  }

#' Generate enrichment bar plot
#'

#' @importFrom ggplot2 ggplot aes geom_point scale_fill_gradient scale_size
#' @importFrom ggplot2 theme_bw ylab labs theme element_text scale_x_continuous
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_split str_replace_all str_remove_all
#' @param df
#' @param isPvalLog
#' @return
#' @export
#'
plot_enrichr_bar <-
  function(df, isPvalLog = "true") {
    library(ggplot2)
    library(forcats)
    library(stringr)
    new_df <- df %>%
      dplyr::mutate(
        Term = as.factor(stringr::str_replace_all(Term, " \\(GO.*", "")),
        len = lengths(stringr::str_split(Genes, ";")),
        pval = -log10(Adjusted.P.value)
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(gene_ratio = eval(parse(text = str_remove_all(Overlap, " ")))) %>%
      dplyr::select(Term, len, pval, gene_ratio, Adjusted.P.value) %>%
      dplyr::mutate(Term = fct_reorder(Term, Adjusted.P.value))

    new_df$Term <-
      factor(new_df$Term, levels = rev(levels(factor(new_df$Term))))

    if (isPvalLog == "true") {
      new_df$Adjusted.P.value <- -1 * log10(new_df$Adjusted.P.value)
      pval_legend <- '-log10(adj. p-value)'
      low_color <- "#EF9A9A"
      high_color <- "#F44336"
    } else {
      pval_legend <- 'adj. p-value'
      low_color <- "#F44336"
      high_color <- "#EF9A9A"
    }

    plot_bar <- ggplot(new_df,
                       aes(x = Adjusted.P.value, y = Term)) +
      geom_bar(stat = "identity", aes(fill = Adjusted.P.value)) +
      scale_fill_gradient(low = low_color, high = high_color, name = pval_legend) +
      theme_bw() +
      ylab("") +
      labs(size = "Overlapping count", color = '123') +
      theme(
        legend.title = element_text(size = 18, ),
        legend.text = element_text(size = 14, ),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)
      ) +
      scale_x_continuous(name = pval_legend) +
      scale_size(range = c(4, 8))

    return(print(plot_bar))
  }
