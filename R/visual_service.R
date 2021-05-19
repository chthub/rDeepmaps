




#' Gene-gene correlation static plot
#' @param req request payload
#' @param gene1 string gene1
#' @param gene2 string gene2
#' @return static plot
#' @export
#'
#'
gene_cor_plot <- function(req, gene1 = "Gad1", gene2 = "Gad2") {
  if ("ATAC" %in% names(e1$obj@assays)) {
    DefaultAssay(e1$obj) <- "SCT"
  }
  send_progress("Calculating gene-gene correlation")
  # gene1 <- VariableFeatures(e1$obj)[1]
  # gene2 <- VariableFeatures(e1$obj)[2]

  df <- FetchData(
    object = e1$obj,
    vars = c(gene1, gene2, group.by = NULL),
    cells = NULL,
    slot = "data"
  )

  plot1 <- ggpubr::ggscatter(
    df,
    x = gene1,
    y = gene2,
    color = "black",
    shape = 21,
    size = 3,
    # Points color, shape and size
    add = "reg.line",
    # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"),
    # Customize reg. line
    conf.int = TRUE,
    # Add confidence interval
    cor.coef = TRUE,
    # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.sep = "\n"
    )
  )
  # VlnPlot(e1$obj, features = c(gene1, gene2))
  return(print(plot1))
}


#' Return embeddings coordinates from a given embedding name
#' @param req request payload
#' @param categoryName string
#' @param embedding string
#'
#' @return json
#' @export
cluster_coords <-
  function(req,
           categoryName = "seurat_clusters",
           embedding = "umap") {
    this_ident_idx <-
      which(colnames(e1$obj@meta.data) == categoryName)[1]

    embedding <- names(e1$obj@reductions[e1$embedding_idx])

    Idents(e1$obj) <- e1$obj@meta.data[, this_ident_idx]

    res1 <- data.frame(
      id = rownames(Embeddings(e1$obj, reduction = embedding)),
      dim1 = Embeddings(e1$obj, reduction = embedding)[, 1],
      dim2 = Embeddings(e1$obj, reduction = embedding)[, 2],
      dim3 = Embeddings(e1$obj, reduction = embedding)[, 3],
      label = as.character(Idents(e1$obj)),
      index = as.integer(Idents(e1$obj)) - 1
    )
    legend <- levels(Idents(e1$obj))
    axis <- c(
      paste0(embedding, "_1"),
      paste0(embedding, "_2"),
      paste0(embedding, "_3")
    )
    dimension <- colnames(res1)
    coords <- as.matrix(res1)
    coords[, 2] <- as.numeric(coords[, 2])
    coords[, 3] <- as.numeric(coords[, 3])
    coords[, 4] <- as.numeric(coords[, 4])
    return(list(
      axis = axis,
      legend = legend,
      dimension = dimension,
      embedding = coords
    ))
  }

#' Return embeddings coordinates and feature expression
#' @param req request payload
#' @param categoryName string
#' @param embedding string
#'
#' @return json
#' @export
feature_coords <-
  function(req,
           gene = "Gad1", assay = "RNA") {
    embedding <- names(e1$obj@reductions[e1$embedding_idx])
    res1 <- data.frame(
      id = rownames(Embeddings(e1$obj, reduction = embedding)),
      dim1 = Embeddings(e1$obj, reduction = embedding)[, 1],
      dim2 = Embeddings(e1$obj, reduction = embedding)[, 2],
      dim3 = Embeddings(e1$obj, reduction = embedding)[, 3],
      expr = FetchData(
        object = e1$obj,
        vars = c(gene)
      )[, 1],
      index = 1
    )
    legend <- c(min(res1$expr), max(res1$expr))
    axis <- c(
      paste0(embedding, "_1"),
      paste0(embedding, "_2"),
      paste0(embedding, "_3")
    )
    dimension <- colnames(res1)
    coords <- as.matrix(res1)
    coords[, 2] <- as.numeric(coords[, 2])
    coords[, 3] <- as.numeric(coords[, 3])
    coords[, 4] <- as.numeric(coords[, 4])
    return(list(
      axis = axis,
      legend = legend,
      dimension = dimension,
      embedding = coords
    ))
  }

#' Plot static UMAP
#' @param req request payload
#' @param categoryName string
#'
#' @return static image
#' @export
umap_plot <- function(req, categoryName = "hgt_cluster") {
  send_progress(paste0("Plotting UMAP using cell category:", categoryName))
  this_ident_idx <-
    which(colnames(e1$obj@meta.data) == categoryName)[1]
  print(e1$embedding_idx)
  this_embedding <- names(e1$obj@reductions[e1$embedding_idx])

  Idents(e1$obj) <- e1$obj@meta.data[, this_ident_idx]
  plot <- DimPlot(e1$obj, reduction = this_embedding)

  # Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]
  return(print(plot))
}

#' Plot static gene UMAP
#' @param req request payload
#' @param gene string
#' @param assay string
#'
#' @return static image
#' @export
gene_umap_plot <- function(req, gene = "Gad1", assay = "RNA") {
  send_progress(paste0("Plotting UMAP gene:", gene))
  this_embedding <- names(e1$obj@reductions[e1$embedding_idx])

  current_assay <- e1$assay_idx
  e1$assay_idx <- which(names(e1$obj@assays) == assay)
  this_assay <- names(e1$obj@assays[e1$assay_idx])
  DefaultAssay(e1$obj) <- this_assay
  plot <-
    FeaturePlot(e1$obj, reduction = this_embedding, feature = gene)
  this_assay <- names(e1$obj@assays[current_assay])
  DefaultAssay(e1$obj) <- this_assay
  return(print(plot))
}


#' Plot violin gene UMAP
#' @param req request payload
#' @param gene string
#' @param split string
#' @param group string
#' @param assay string
#'
#' @return static image
#' @export
#'
violin_gene_plot <-
  function(req,
           gene = "Gad1",
           split = "sex",
           group = "cell_type",
           assay = "RNA") {
    # ident_idx=9
    send_progress(paste0("Plotting violin gene:", gene))
    if (split == "NULL") {
      Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]
      plot <-
        VlnPlot(e1$obj, gene, group.by = colnames(e1$obj@meta.data)[e1$ident_idx])
    } else {
      idx <- which(colnames(e1$obj@meta.data) == split)
      Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]

      current_assay <- e1$assay_idx
      e1$assay_idx <- which(names(e1$obj@assays) == assay)
      this_assay <- names(e1$obj@assays[e1$assay_idx])
      DefaultAssay(e1$obj) <- this_assay
      this_group_idx <-
        which(colnames(e1$obj@meta.data) == group)
      plot <-
        VlnPlot(
          e1$obj,
          gene,
          split.by = colnames(e1$obj@meta.data)[idx],
          group.by = colnames(e1$obj@meta.data)[this_group_idx]
        )
      this_assay <- names(e1$obj@assays[current_assay])
      DefaultAssay(e1$obj) <- this_assay
    }

    return(print(plot))
  }


#' Complex heatmap
#' @param req request payload
#' @param genes vector of string
#' @param meta metadata to annotate in columns
#' @param color color palette
#' @return static image
#' @export
#'
static_heatmap <-
  function(req,
           genes = c("CD74", "CD7"),
           meta = "cell_type",
           color = NULL) {
    # ident_idx=9
    library(ComplexHeatmap)
    send_progress(paste0("Plotting heatmap"))
    features <- VariableFeatures(e1$obj)[1:20]
    features <- genes[genes %in% rownames(e1$obj)]
    Idents(e1$obj) <- e1$obj$hgt_cluster

    cell_info <- Idents(e1$obj)
    cell_label <- cbind(colnames(e1$obj), as.character(cell_info))
    colnames(cell_label) <- c("cell_name", "label")
    cell_label <- cell_label[order(cell_label[, 1]), ]

    cell_label <- as.data.frame(cell_label)
    label_data <- cell_label[order(cell_label[, 2]), ]
    exp_data <- GetAssayData(e1$obj, slot = "data")
    cell_idx <- as.character(label_data[, 1])
    exp_data <- exp_data[, cell_idx]

    if (ncol(exp_data) > 500) {
      this_bin <- ncol(exp_data) %/% 500
      small_cell_idx <- seq(1, ncol(exp_data), by = this_bin)
      small_exp_data <<- t(apply(exp_data, 1, function(x) {
        BinMean(x, every = this_bin)
      }))
      small_cell_label <- label_data[small_cell_idx, ]
    }
    colnames(small_exp_data) <- small_cell_label[, 1]
    library(matrixStats)
    small_exp_data <- log1p(small_exp_data)
    small_exp_data <-
      (small_exp_data - rowMeans(small_exp_data)) / rowSds(small_exp_data)


    mat <-
      small_exp_data[match(features, rownames(small_exp_data)), ]

    library(circlize)
    col_fun <- as.character(Polychrome::palette36.colors(36)[-2])
    names(col_fun) <- seq_along(col_fun)
    col_fun <-
      as.character(Polychrome::palette36.colors(36)[-2])[seq_along(levels(Idents(e1$obj)))]
    names(col_fun) <- levels(Idents(e1$obj))

    ha <-
      rowAnnotation(foo = anno_mark(at = match(features, rownames(mat)), labels = features))

    ht <- Heatmap(
      mat,
      show_row_names = F,
      show_column_names = F,
      top_annotation = HeatmapAnnotation(
        Cluster = as.character(small_cell_label[, 2]),
        col = list(Cluster = col_fun)
      ),
      column_order = c(1:ncol(mat)),
      cluster_rows = F,
      name = "z-score",
      show_column_dend = F,
      right_annotation = ha
    )

    # plot <- draw(ht)

    return(draw(ht))
  }



#' Interactive complex heatmap
#' @param req request payload
#' @return html
#' @export
#'
interactive_heatmap <-
  function(req,
           genes = c("CD74", "CD7"),
           meta = "cell_type",
           color = NULL) {
    return(1)
  }


#' ATAC coverage static plot
#' @importFrom Signac CoveragePlot
#' @param req request payload
#' @param type string
#' @param gene string#'
#' @param flank string
#' @param chr string
#' @param start string
#' @param end string
#' @param is_annotation string
#' @param is_peak string
#' @return static plot
#' @export
#'
coverage_plot <-
  function(req,
           type = "gene",
           gene = "GAD2",
           flank = "10000",
           chr = "chr10",
           start = "26116307",
           end = "26316307",
           is_annotation = T,
           is_peak = F) {
    # Idents(e1$obj) <- e1$obj$cell_type
    if (type == "gene") {
      this_ranges <- get_gene_range(gene = gene, flank = flank)
    } else if (type == "region") {
      this_ranges <- paste(chr, start, end, sep = "-")
    }
    send_progress(paste0("Plotting coverage"))
    cov_plot <- Signac::CoveragePlot(
      object = e1$obj,
      assay = "ATAC",
      region = this_ranges,
      annotation = is_annotation,
      peaks = is_peak
    )
    send_progress(paste0("Plotting peaks"))
    peak_plot <- Signac::PeakPlot(
      object = e1$obj,
      assay = "ATAC",
      region = this_ranges
    )
    send_progress(paste0("Plotting tiles"))
    tile_plot <- Signac::TilePlot(
      object = e1$obj,
      assay = "ATAC",
      region = this_ranges
    )
    send_progress(paste0("Combining plots"))
    combine_plot <-
      Signac::CombineTracks(plotlist = list(cov_plot, tile_plot, peak_plot))
    return(print(combine_plot))
  }
