#' Load data to Seurat
#' @import Seurat
#' @param req request payload
#' @param filename string
#' @param min_cells number
#' @param min_genes number
#' @param nVariableFeatures number
#' @param percentMt number
#' @param removeRibosome boolean
#'
#' @return list of basic QC metrics
#' @export
#'
load_seurat <-
  function(req,
           filename,
           min_cells = 1,
           min_genes = 200,
           nVariableFeatures = 3000,
           percentMt = 5,
           removeRibosome = FALSE) {
    raw_expr_data <- iris3api::zeisel_2015$expr
    meta <- iris3api::zeisel_2015$meta
    raw_obj <- CreateSeuratObject(raw_expr_data)
    e1$obj <-
      CreateSeuratObject(raw_expr_data,
        min.cells = min_cells,
        min.features = min_genes
      )
    empty_ident <- as.factor(e1$obj$orig.ident)
    levels(empty_ident) <-
      rep("empty_ident", length(levels(empty_ident)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_ident, col.name = "empty_ident")

    e1$obj <-
      AddMetaData(e1$obj, meta$Label, col.name = "cell_type")
    e1$obj <- AddMetaData(e1$obj, meta$Sex, col.name = "sex")
    e1$obj <- AddMetaData(e1$obj, meta$Sample, col.name = "sample")
    e1$obj <-
      AddMetaData(e1$obj,
        PercentageFeatureSet(e1$obj, pattern = "^MT-"),
        col.name = "percent.mt"
      )

    Idents(e1$obj) <- e1$obj$orig.ident
    rb.genes <-
      rownames(e1$obj)[grep("^Rp[sl][[:digit:]]", rownames(e1$obj))]
    percent.ribo <-
      Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
    e1$obj <-
      AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")
    e1$obj <-
      subset(e1$obj, subset = `percent.mt` < as.numeric(percentMt))
    raw_percent_zero <-
      length(which((as.matrix(
        GetAssayData(raw_obj)
      ) > 0))) / length(GetAssayData(raw_obj))
    filter_percent_zero <-
      length(which((as.matrix(
        GetAssayData(e1$obj)
      ) > 0))) / length(GetAssayData(e1$obj))
    raw_mean_expr <- mean(as.matrix(GetAssayData(raw_obj)))
    filter_mean_expr <- mean(as.matrix(GetAssayData(e1$obj)))
    e1$obj <-
      FindVariableFeatures(
        e1$obj,
        selection.method = "vst",
        nfeatures = as.numeric(nVariableFeatures),
        verbose = F
      )
    return(
      list(
        raw_n_genes = dim(raw_obj)[1],
        raw_n_cells = dim(raw_obj)[2],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = dim(e1$obj)[1],
        filter_n_cells = dim(e1$obj)[2],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    )
  }
