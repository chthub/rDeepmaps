#' Load data to Seurat
#' @import Seurat
#' @param req request payload
#' @param idx sample index
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
load_single_rna <-
  function(req,
           idx = 1,
           jobid = "example",
           type = "eg",
           min_cells = 1,
           min_genes = 200,
           nVariableFeatures = 2000,
           percentMt = 5,
           removeRibosome = FALSE,
           expr = NULL,
           label = NULL,
           species = "Human") {
    # expr_type <- label_type <- 'application/vnd.ms-excel'
    # label_path <- "/mnt/c/Users/flyku/Documents/deepmaps-data/fbb36e99797c313276f8bca86539cb3c"
    # expr_path <- "/mnt/c/Users/flyku/Documents/deepmaps-data/b16af116a67fc89ec0991170f99b7503"

    # Yan: 1619284757781
    # Zeisel: 1619291336397
    TOTAL_STEPS <- 10
    send_progress(paste0("Start processing scRNA-seq dataset"))
    if (jobid != "example") {
      expr_type <- as.character(expr$mimetype[1])
      expr_path <- as.character(expr$filename[1])
      raw_expr_data <- read_deepmaps(expr_type, expr_path)

      if(length(label) > 0 ){
        label_type <- as.character(label$mimetype[1])
        label_path <- as.character(label$filename[1])
        e1$meta <- read_deepmaps(label_type, label_path)
      }
    } else {
      raw_expr_data <- iris3api::zeisel_2015$expr
      e1$meta <- iris3api::zeisel_2015$meta
      species <- "Mouse"
    }

    send_progress(paste0("Creating Seurat object"))
    raw_obj <- CreateSeuratObject(raw_expr_data)
    e1$obj <-
      CreateSeuratObject(
        raw_expr_data,
        min.cells = as.numeric(min_cells),
        min.features = as.numeric(min_genes)
      )
    e1$species <- species
    empty_category <- as.factor(e1$obj$orig.ident)
    levels(empty_category) <-
      rep("empty_category", length(levels(empty_category)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_category, col.name = "empty_category")

    for (idx in seq_along(colnames(e1$meta))) {
      this_meta_name <- colnames(e1$meta)[idx]
      e1$obj <-
        AddMetaData(e1$obj, e1$meta[, idx], col.name = this_meta_name)
    }

    e1$obj <-
      AddMetaData(e1$obj,
        PercentageFeatureSet(e1$obj, pattern = "^MT-"),
        col.name = "percent.mt"
      )
    send_progress("Calculating data summary statistics")
    Idents(e1$obj) <- e1$obj$orig.ident
    rb.genes <-
      rownames(e1$obj)[grep("^Rp[sl][[:digit:]]", rownames(e1$obj),
        ignore.case =
          TRUE
      )]
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
    send_progress("Finding highly variable features")
    e1$obj <-
      FindVariableFeatures(
        e1$obj,
        selection.method = "vst",
        nfeatures = as.numeric(nVariableFeatures),
        verbose = F
      )
    send_progress("Normalizing data")
    e1$obj <- NormalizeData(e1$obj, verbose = F)
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

#' Load multi scRNA-seq data to Seurat
#' @import Seurat
#' @param idx sample index
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
load_multi_rna <-
  function(req,
           idx = 1,
           jobid = "1642208223822",
           type = "multiRna",
           mode = "RNA",
           load = "load",
           filename,
           min_cells = 20000,
           min_genes = 0.001,
           nVariableFeatures = 2000,
           percentMt = 0.1,
           percentRb = 0.5,
           removeOutlier = FALSE,
           label = NULL,
           species = "Human") {
    send_progress("Start processing scRNA-seq dataset")
    expr <- data.frame(jobid = "1642208223822")
    expr$mimetype <- "application/vnd.ms-excel"
    expr$filename <- "c6c13e15a1f1dd3f0fdcba2ae3c65d20"
    if (jobid == "example1") {
      path <- gsub("/scratch/deepmaps","",as.character(path))
      expr_type <- as.character(expr$mimetype[idx])
      expr_path <- as.character(expr$filename[idx])
      raw_expr_data <- read_deepmaps(expr_type, expr_path)
      if(length(label) > 0 ){
        label_type <- as.character(label$mimetype[1])
        label_path <- as.character(label$filename[1])
        e1$meta <- read_deepmaps(label_type, label_path)
      }
    } else {
      raw_expr_data <- iris3api::ifnb_2800$expr
      e1$meta <- iris3api::ifnb_2800$meta
      # write.csv(raw_expr_data[,1:1400],"human_ifnb_sample1_expr.csv")
      # write.csv(raw_expr_data[,1401:2800],"human_ifnb_sample2_expr.csv")
      # write.csv(e1$meta,"human_ifnb_label.csv")
    }
    send_progress("Creating Seurat object")
    raw_obj <- CreateSeuratObject(raw_expr_data)
    e1$obj <-
      CreateSeuratObject(
        raw_expr_data,
        min.cells = as.numeric(min_genes) * ncol(raw_expr_data),
      )

    if(load == "Load" && file.exists(paste0(get_base_dir(), jobid, ".qsave"))) {
      e1$obj <- qs::qread(paste0(get_base_dir(), jobid, ".qsave"))
    }

    e1$species <- "Human"
    empty_category <- as.factor(e1$obj$orig.ident)
    levels(empty_category) <-
      rep("empty_category", length(levels(empty_category)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_category, col.name = "empty_category")

    for (index in seq_along(colnames(e1$meta))) {
      this_meta_name <- colnames(e1$meta)[index]
      e1$obj <-
        AddMetaData(e1$obj, e1$meta[, index], col.name = this_meta_name)
    }

    if (idx != 0) {
      this_idx <- as.character(levels(as.factor(e1$obj$sample))[idx])
      e1$obj <- subset(e1$obj, sample == this_idx)
    }

    e1$obj <-
      AddMetaData(e1$obj,
        PercentageFeatureSet(e1$obj, pattern = "^MT-"),
        col.name = "percent.mt"
      )

    Idents(e1$obj) <- e1$obj$orig.ident
    rb.genes <-
      rownames(e1$obj)[grep("^Rp[sl][[:digit:]]", rownames(e1$obj), ignore.case = T)]
    percent.ribo <-
      Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
    send_progress("Calculating data summary statistics")
    e1$obj <-
      AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")
    e1$obj <-
      subset(e1$obj, subset = `percent.mt` < as.numeric(percentMt) * 100)
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
    send_progress("Finding highly variable features")
    e1$obj <-
      FindVariableFeatures(
        e1$obj,
        selection.method = "vst",
        nfeatures = as.numeric(nVariableFeatures),
        verbose = F
      )
    send_progress("Normalizing data")
    e1$obj <- NormalizeData(e1$obj, verbose = F)
    if(load == "Calculate") {
      qs::qsave(e1$obj, paste0(get_base_dir(), jobid, ".qsave"))
    }
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


#' Load matched scRNA-seq and scATAC-seq data to Seurat / Signac
#' @import Seurat
#' @param req request payload
#' @param idx sample index
#' @param filename string
#' @param mode string
#' @param min_cells number
#' @param min_genes number
#' @param nVariableFeatures number
#' @param percentMt number
#' @param removeRibosome boolean
#'
#' @return list of basic QC metrics
#' @export
#'
load_multiome2 <-
  function(req,
           idx = 0,
           filename,
           jobid = "lymph",
           mode = "ATAC",
           min_cells = 1,
           min_genes = 200,
           nVariableFeatures = 2000,
           percentMt = 25,
           removeRibosome = FALSE,
           expr = NULL,
           label = NULL,
           species = "Human") {

    TOTAL_STEPS <- 6

    if (file.exists("/data")) {
      base_dir <- "/data/"
    } else {
      base_dir <- "F:/DeepMAPS/data/"
    }
    e1$obj <- qs::qread(paste0(base_dir, "lymph_obj.qsave"))
    raw_obj <- qs::qread(paste0(base_dir, "lymph_obj.qsave"))

    fragments <- CreateFragmentObject(
      path = paste0(base_dir, "lymph_node_lymphoma_14k_atac_fragments.tsv.gz"),
      cells = colnames(e1$obj),
      validate.fragments = FALSE
    )
    e1$obj@assays$ATAC@fragments[[1]] <- fragments
    Idents(e1$obj) <- e1$obj$orig.ident
    rb.genes <-
      rownames(e1$obj)[grep("^RP[SL][[:digit:]]", rownames(e1$obj),
                            ignore.case =
                              TRUE
      )]
    DefaultAssay(e1$obj) <- "RNA"
    percent.ribo <-
      Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
    e1$obj <-
      AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")


    e1$obj <-
      AddMetaData(e1$obj,
                  PercentageFeatureSet(e1$obj, pattern = "^MT-"),
                  col.name = "percent.mt"
      )
    e1$obj <-
      FindVariableFeatures(
        e1$obj,
        selection.method = "vst",
        nfeatures = as.numeric(2000),
        verbose = F
      )

    e1$obj <- NormalizeData(e1$obj, verbose = F)
    e1$obj <- ScaleData(e1$obj, verbose = F)


    # Add empety ident
    empty_category <- as.factor(e1$obj$orig.ident)
    levels(empty_category) <-
      rep("empty_category", length(levels(empty_category)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_category, col.name = "empty_category")
    raw_obj <- e1$obj
    send_progress(paste0("Calculating data summary statistics"))
    e1$species <- species
    if(mode == "RNA") {
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
      res <- list(
        raw_n_genes = dim(raw_obj)[1],
        raw_n_cells = dim(raw_obj)[2],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = dim(e1$obj)[1],
        filter_n_cells = dim(e1$obj)[2],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    } else {
      raw_percent_zero <-
        length(which((as.matrix(
          GetAssayData(raw_obj, assay="ATAC")
        ) > 0))) / length(GetAssayData(raw_obj, assay="ATAC"))
      filter_percent_zero <-
        length(which((as.matrix(
          GetAssayData(e1$obj, assay="ATAC")
        ) > 0))) / length(GetAssayData(e1$obj, assay="ATAC"))
      raw_mean_expr <- mean(as.matrix(GetAssayData(raw_obj)))
      filter_mean_expr <- mean(as.matrix(GetAssayData(e1$obj)))
      res <- list(
        raw_n_genes = nrow(GetAssayData(raw_obj, assay="ATAC"))[1],
        raw_n_cells = ncol(GetAssayData(raw_obj, assay="ATAC"))[1],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = nrow(GetAssayData(e1$obj, assay="ATAC"))[1] - 6000,
        filter_n_cells = ncol(GetAssayData(e1$obj, assay="ATAC"))[1],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    }
    return(
      res
    )
  }

#' Load matched scRNA-seq and scATAC-seq data to Seurat / Signac
#' @import Seurat
#' @param req request payload
#' @param idx sample index
#' @param filename string
#' @param mode string
#' @param min_cells number
#' @param min_genes number
#' @param nVariableFeatures number
#' @param percentMt number
#' @param removeRibosome boolean
#'
#' @return list of basic QC metrics
#' @export
#'
load_multiome3 <-
  function(req,
           idx = 0,
           filename,
           jobid = "pbmc_unsorted_10k",
           mode = "ATAC",
           min_cells = 1,
           min_genes = 200,
           nVariableFeatures = 2000,
           percentMt = 25,
           removeRibosome = FALSE,
           expr = NULL,
           label = NULL,
           species = "Human") {

    TOTAL_STEPS <- 6

    if (file.exists("/data")) {
      base_dir <- "/data/"
    } else {
      base_dir <- "C:/Users/flyku/Desktop/iris3/pbmc_match/pbmc/"
    }
    e1$obj <- qs::qread(paste0(base_dir, "pbmc_unsorted_10k_obj.qsave"))
    raw_obj <- qs::qread(paste0(base_dir, "pbmc_unsorted_10k_obj.qsave"))

    fragments <- CreateFragmentObject(
      path = paste0(base_dir, "pbmc_unsorted_10k_atac_fragments.tsv.gz"),
      cells = colnames(e1$obj),
      validate.fragments = FALSE
    )
    e1$obj@assays$ATAC@fragments[[1]] <- fragments
    Idents(e1$obj) <- e1$obj$orig.ident
    rb.genes <-
      rownames(e1$obj)[grep("^RP[SL][[:digit:]]", rownames(e1$obj),
                            ignore.case =
                              TRUE
      )]
    DefaultAssay(e1$obj) <- "RNA"
    percent.ribo <-
      Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
    e1$obj <-
      AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")


    e1$obj <-
      AddMetaData(e1$obj,
                  PercentageFeatureSet(e1$obj, pattern = "^MT-"),
                  col.name = "percent.mt"
      )
    e1$obj <-
      FindVariableFeatures(
        e1$obj,
        selection.method = "vst",
        nfeatures = as.numeric(2000),
        verbose = F
      )

    e1$obj <- NormalizeData(e1$obj, verbose = F)
    e1$obj <- ScaleData(e1$obj, verbose = F)


    # Add empety ident
    empty_category <- as.factor(e1$obj$orig.ident)
    levels(empty_category) <-
      rep("empty_category", length(levels(empty_category)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_category, col.name = "empty_category")
    raw_obj <- e1$obj
    send_progress(paste0("Calculating data summary statistics"))
    e1$species <- species
    if(mode == "RNA") {
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
      res <- list(
        raw_n_genes = dim(raw_obj)[1],
        raw_n_cells = dim(raw_obj)[2],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = dim(e1$obj)[1],
        filter_n_cells = dim(e1$obj)[2],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    } else {
      raw_percent_zero <-
        length(which((as.matrix(
          GetAssayData(raw_obj, assay="ATAC")
        ) > 0))) / length(GetAssayData(raw_obj, assay="ATAC"))
      filter_percent_zero <-
        length(which((as.matrix(
          GetAssayData(e1$obj, assay="ATAC")
        ) > 0))) / length(GetAssayData(e1$obj, assay="ATAC"))
      raw_mean_expr <- mean(as.matrix(GetAssayData(raw_obj)))
      filter_mean_expr <- mean(as.matrix(GetAssayData(e1$obj)))
      res <- list(
        raw_n_genes = nrow(GetAssayData(raw_obj, assay="ATAC"))[1],
        raw_n_cells = ncol(GetAssayData(raw_obj, assay="ATAC"))[1],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = nrow(GetAssayData(e1$obj, assay="ATAC"))[1] - 6000,
        filter_n_cells = ncol(GetAssayData(e1$obj, assay="ATAC"))[1],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    }
    return(
      res
    )
  }

#' Load matched scRNA-seq and scATAC-seq data to Seurat / Signac
#' @import Seurat
#' @param req request payload
#' @param idx sample index
#' @param filename string
#' @param mode string
#' @param min_cells number
#' @param min_genes number
#' @param nVariableFeatures number
#' @param percentMt number
#' @param removeRibosome boolean
#'
#' @return list of basic QC metrics
#' @export
#'
load_multiome4 <-
  function(req,
           idx = 0,
           filename,
           jobid = "lymph",
           mode = "ATAC",
           min_cells = 1,
           min_genes = 200,
           nVariableFeatures = 2000,
           percentMt = 25,
           removeRibosome = FALSE,
           expr = NULL,
           label = NULL,
           species = "Human") {

    TOTAL_STEPS <- 6

    if (file.exists("/data")) {
      base_dir <- "/data/"
    } else {
      base_dir <- "F:/DeepMAPS/data/"
    }

    raw_obj <- qs::qread(paste0(base_dir, "lymphoma_14k_raw_obj.qsave"))
    e1$obj <- qs::qread(paste0(base_dir, "lymphoma_14k_obj.qsave"))
    #tmp_ident <- e1$obj$cell_type
    #levels(tmp_ident)[1] <- "Exhausted CD4+ T cell"
    #levels(tmp_ident)[2] <- "Exhausted CD8+ T cell"
    #levels(tmp_ident)[3] <- "Macrophage"
    #levels(tmp_ident)[4] <- "Naive CD8+ T cell"
    #levels(tmp_ident)[5] <- "Effector-like CD4+ T cell_1"
    #levels(tmp_ident)[6] <- "DSLL state-2"
    #levels(tmp_ident)[7] <- "Effector-like CD4+ T cell_2"
    #levels(tmp_ident)[8] <- "DSLL state-1"
    #levels(tmp_ident)[9] <- "Normal B cell"
    #levels(tmp_ident)[10] <- "Effector-like CD8+ T cell"
    #levels(tmp_ident)[11] <- "DC"
    #e1$obj <- AddMetaData(e1$obj, tmp_ident, col.name = "cell_type")
    #DimPlot(e1$obj, group.by = "cell_type", label = T)
    #qs::qsave(e1$obj, "lymphoma_14k_obj.qsave")

    fragments <- CreateFragmentObject(
      path = paste0(base_dir, "lymph_node_lymphoma_14k_atac_fragments.tsv.gz"),
      cells = colnames(e1$obj),
      validate.fragments = FALSE
    )
    e1$obj@assays$ATAC@fragments[[1]] <- fragments
    raw_obj@assays$ATAC@fragments[[1]] <- fragments
    Idents(e1$obj) <- e1$obj$orig.ident
    rb.genes <-
      rownames(e1$obj)[grep("^RP[SL][[:digit:]]", rownames(e1$obj),
                            ignore.case =
                              TRUE
      )]
    DefaultAssay(e1$obj) <- "RNA"
    percent.ribo <-
      Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
    e1$obj <-
      AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")


    e1$obj <-
      AddMetaData(e1$obj,
                  PercentageFeatureSet(e1$obj, pattern = "^MT-"),
                  col.name = "percent.mt"
      )
    e1$obj <-
      FindVariableFeatures(
        e1$obj,
        selection.method = "vst",
        nfeatures = as.numeric(2000),
        verbose = F
      )

    e1$obj <- NormalizeData(e1$obj, verbose = F)
    e1$obj <- ScaleData(e1$obj, verbose = F)


    # Add empety ident
    empty_category <- as.factor(e1$obj$orig.ident)
    levels(empty_category) <-
      rep("empty_category", length(levels(empty_category)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_category, col.name = "empty_category")

    e1$species <- species
    if(mode == "RNA") {
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
      res <- list(
        raw_n_genes = dim(raw_obj)[1],
        raw_n_cells = dim(raw_obj)[2],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = dim(e1$obj)[1],
        filter_n_cells = dim(e1$obj)[2],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    } else {
      raw_percent_zero <-
        length(which((as.matrix(
          GetAssayData(raw_obj, assay="ATAC")
        ) > 0))) / length(GetAssayData(raw_obj, assay="ATAC"))
      filter_percent_zero <-
        length(which((as.matrix(
          GetAssayData(e1$obj, assay="ATAC")
        ) > 0))) / length(GetAssayData(e1$obj, assay="ATAC"))
      raw_mean_expr <- mean(as.matrix(GetAssayData(raw_obj, assay="ATAC")))
      filter_mean_expr <- mean(as.matrix(GetAssayData(e1$obj, assay="ATAC")))
      res <- list(
        raw_n_genes = nrow(GetAssayData(raw_obj, assay="ATAC"))[1],
        raw_n_cells = ncol(GetAssayData(raw_obj, assay="ATAC"))[1],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = nrow(GetAssayData(e1$obj, assay="ATAC"))[1],
        filter_n_cells = ncol(GetAssayData(e1$obj, assay="ATAC"))[1],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    }
    return(
      res
    )
  }


#' Load matched scRNA-seq and scATAC-seq data to Seurat / Signac
#' @import Seurat
#' @param req request payload
#' @param idx sample index
#' @param filename string
#' @param mode string
#' @param min_cells number
#' @param min_counts number
#' @param nVariableFeatures number
#' @param percentMt number
#' @param percentMt number
#' @param removeOutlier boolean
#'
#' @return list of basic QC metrics
#' @export
#'
load_multiome <-
  function(req,
           idx = 0,
           filename,
           jobid = "example",
           mode = "RNA",
           min_cells = 0.001,
           min_counts = 20000,
           removeOutlier = T,
           nVariableFeatures = 2000,
           percentMt = 0.25,
           percentRb = 0.5,
           expr = NULL,
           label = NULL,
           species = "Human",
           destination = NULL
           ) {
    # expr_type <- label_type <- 'application/octet-stream'
    # label_path <- ""
    # expr_path <- "9a9841a85c48b692e70bc03db811ccdc"

    # Brain: 1619298241128
    # PBMC 3K: 1619297987450
    # PBMC : 1619311996943

    TOTAL_STEPS <- 6
    if (jobid == "example") {
      send_progress(paste0("Loading example data"))
      if (file.exists("/data")) {
        sleep(0.1)
        e1$obj <- qs::qread("/data/pbmc_sorted_3k.qsave")
        raw_obj <- qs::qread("/data/pbmc_match_3k.qsave")
        e1$obj@assays$ATAC@fragments[[1]]@path <-
          "/data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
      } else {
        e1$obj <-
          qs::qread(
            "C:/Users/flyku/Documents/GitHub/iris3api/inst/extdata/tmp1.qsave"
          )
        raw_obj <- e1$obj

        e1$obj@assays$ATAC@fragments[[1]]@path <-
          "C:/Users/flyku/Desktop/iris3/eg2/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
        iris3api::set_embedding(name = "umap.rna")

      }
      Idents(e1$obj) <- e1$obj$orig.ident
      rb.genes <-
        rownames(e1$obj)[grep("^R[Pp][slSL][[:digit:]]", rownames(e1$obj),
                              ignore.case =
                                TRUE)]
      percent.ribo <-
        Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
      e1$obj <-
        AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")
    } else {
      #path <- "/mnt/c/Users/flyku/Documents/deepmaps-data/0cc4610f4c9b4c0f97ec5a84c2e19e30"
      #0cc4610f4c9b4c0f97ec5a84c2e19e30
      #
      #path <- gsub("/mnt/c","c:/",as.character(path))
      path <- gsub("/mnt/c","c:/",as.character(expr$path[1]))
      path <- gsub("/scratch/deepmaps","",as.character(path))
      print(path)
      raw_expr_data <- Read10X_h5(paste0(path))
      raw_obj <- CreateSeuratObject(raw_expr_data$`Gene Expression`)
      e1$obj <-
        CreateSeuratObject(
          raw_expr_data$`Gene Expression`,
          min.cells = as.numeric(0.001) *  ncol(raw_expr_data$`Gene Expression`)
        )
      atac_obj <- CreateChromatinAssay(counts = raw_expr_data$Peaks[,colnames(e1$obj)],
                                       sep = c(":", "-"))

      e1$obj[["ATAC"]] <- atac_obj
      e1$obj <-
        AddMetaData(e1$obj,
                    PercentageFeatureSet(e1$obj, pattern = "^MT-"),
                    col.name = "percent.mt"
        )

      Idents(e1$obj) <- e1$obj$orig.ident
      rb.genes <-
        rownames(e1$obj)[grep("^R[p][sl][[:digit:]]", rownames(e1$obj),
                              ignore.case =
                                TRUE)]
      percent.ribo <-
        Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
      e1$obj <-
        AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")

    }
    e1$obj <-
      subset(e1$obj, subset = `percent.mt` < as.numeric(percentMt) * 100)
    #e1$obj <-
    #  subset(e1$obj, subset = `percent.ribo` < as.numeric(percentRb) * 100)
    # Add empety ident
    empty_category <- as.factor(e1$obj$orig.ident)
    levels(empty_category) <-
      rep("empty_category", length(levels(empty_category)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_category, col.name = "empty_category")
    raw_obj <- e1$obj
    e1$species <- species
    if(mode == "RNA") {
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
      res <- list(
        raw_n_genes = dim(raw_obj)[1],
        raw_n_cells = dim(raw_obj)[2],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = dim(e1$obj)[1],
        filter_n_cells = dim(e1$obj)[2],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    } else {
      raw_percent_zero <-
        length(which((as.matrix(
          GetAssayData(raw_obj, assay="ATAC")
        ) > 0))) / length(GetAssayData(raw_obj, assay="ATAC"))
      filter_percent_zero <-
        length(which((as.matrix(
          GetAssayData(e1$obj, assay="ATAC")
        ) > 0))) / length(GetAssayData(e1$obj, assay="ATAC"))
      raw_mean_expr <- mean(as.matrix(GetAssayData(raw_obj)))
      filter_mean_expr <- mean(as.matrix(GetAssayData(e1$obj)))
      res <- list(
        raw_n_genes = nrow(GetAssayData(raw_obj, assay="ATAC"))[1],
        raw_n_cells = ncol(GetAssayData(raw_obj, assay="ATAC"))[1],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = nrow(GetAssayData(e1$obj, assay="ATAC"))[1] * 0.8,
        filter_n_cells = ncol(GetAssayData(e1$obj, assay="ATAC"))[1],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    }
    return(
      res
    )
  }

#' Load CITE-seq data
#' @import Seurat
#' @param req request payload
#' @param idx sample index
#' @param filename string
#' @param mode string
#' @param min_cells number
#' @param min_counts number
#' @param nVariableFeatures number
#' @param percentMt number
#' @param percentMt number
#' @param removeOutlier boolean
#'
#' @return list of basic QC metrics
#' @export
#'
load_citeseq <-
  function(req,
           idx = 0,
           filename,
           jobid = "example",
           mode = "RNA",
           min_cells = 0.001,
           min_counts = 20000,
           removeOutlier = T,
           nVariableFeatures = 2000,
           percentMt = 0.25,
           percentRb = 0.5,
           expr = NULL,
           label = NULL,
           species = "Human",
           destination = NULL
  ) {
    if (jobid !="example1") {
      if (file.exists("/data")) {
        e1$obj <- qs::qread("/data/PBMCandLung_obj.qsave")
        raw_obj <- qs::qread("/data/PBMCandLung_obj.qsave")
      } else {
        e1$obj <-
          qs::qread(
            "F:/DeepMAPS/data/PBMCandLung_obj.qsave"
          )
        raw_obj <- e1$obj
        iris3api::set_embedding(name = "umap.hgt")

      }
      Idents(e1$obj) <- e1$obj$orig.ident
      rb.genes <-
        rownames(e1$obj)[grep("^R[Pp][slSL][[:digit:]]", rownames(e1$obj),
                              ignore.case =
                                TRUE)]
      percent.ribo <-
        Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
      e1$obj <-
        AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")
    } else {

      path <- gsub("/mnt/c","c:/",as.character(expr$path[1]))
      path <- gsub("/scratch/deepmaps","",as.character(path))
      print(path)
      raw_expr_data <- Read10X_h5(paste0(path))
      raw_obj <- CreateSeuratObject(raw_expr_data$`Gene Expression`)
      e1$obj <-
        CreateSeuratObject(
          raw_expr_data$`Gene Expression`,
          min.cells = as.numeric(0.001) *  ncol(raw_expr_data$`Gene Expression`)
        )
      atac_obj <- CreateChromatinAssay(counts = raw_expr_data$Peaks[,colnames(e1$obj)],
                                       sep = c(":", "-"))

      e1$obj[["ATAC"]] <- atac_obj
      e1$obj <-
        AddMetaData(e1$obj,
                    PercentageFeatureSet(e1$obj, pattern = "^MT-"),
                    col.name = "percent.mt"
        )

      Idents(e1$obj) <- e1$obj$orig.ident
      rb.genes <-
        rownames(e1$obj)[grep("^R[p][sl][[:digit:]]", rownames(e1$obj),
                              ignore.case =
                                TRUE)]
      percent.ribo <-
        Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
      e1$obj <-
        AddMetaData(e1$obj, percent.ribo, col.name = "percent.ribo")

    }
    #e1$obj <-
    #  subset(e1$obj, subset = `percent.mt` < as.numeric(percentMt) * 100)
    #e1$obj <-
    #  subset(e1$obj, subset = `percent.ribo` < as.numeric(percentRb) * 100)
    # Add empety ident
    empty_category <- as.factor(e1$obj$orig.ident)
    levels(empty_category) <-
      rep("empty_category", length(levels(empty_category)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_category, col.name = "empty_category")
    raw_obj <- e1$obj
    e1$species <- species
    if(mode == "RNA") {
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
      res <- list(
        raw_n_genes = dim(raw_obj)[1],
        raw_n_cells = dim(raw_obj)[2],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = dim(e1$obj)[1],
        filter_n_cells = dim(e1$obj)[2],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    } else {
      raw_percent_zero <-
        length(which((as.matrix(
          GetAssayData(raw_obj, assay="ADT")
        ) > 0))) / length(GetAssayData(raw_obj, assay="ADT"))
      filter_percent_zero <-
        length(which((as.matrix(
          GetAssayData(e1$obj, assay="ADT")
        ) > 0))) / length(GetAssayData(e1$obj, assay="ADT"))
      raw_mean_expr <- mean(as.matrix(GetAssayData(raw_obj)))
      filter_mean_expr <- mean(as.matrix(GetAssayData(e1$obj)))
      res <- list(
        raw_n_genes = nrow(GetAssayData(raw_obj, assay="ADT"))[1],
        raw_n_cells = ncol(GetAssayData(raw_obj, assay="ADT"))[1],
        raw_percent_zero = raw_percent_zero,
        raw_mean_expr = raw_mean_expr,
        filter_n_genes = nrow(GetAssayData(e1$obj, assay="ADT"))[1],
        filter_n_cells = ncol(GetAssayData(e1$obj, assay="ADT"))[1],
        filter_percent_zero = filter_percent_zero,
        filter_mean_expr = filter_mean_expr
      )
    }
    return(
      res
    )
  }
