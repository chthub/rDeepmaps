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
    # expr_path <- "/mnt/c/Users/flyku/Documents/deepmaps-data/3c0955e700b4bf2d177442236d0c2395"

    # Yan: 1619284757781
    # Zeisel: 1619291336397
    TOTAL_STEPS <- 10
    send_progress(paste0("Start processing scRNA-seq dataset"))
    if (jobid != "example") {
      expr_type <- as.character(expr$mimetype[1])
      expr_path <- as.character(expr$filename[1])
      raw_expr_data <- read_deepmaps(expr_type, expr_path)

      label_type <- as.character(label$mimetype[1])
      label_path <- as.character(label$filename[1])
      e1$meta <- read_deepmaps(label_type, label_path)
    } else {
      raw_expr_data <- iris3api::zeisel_2015$expr
      e1$meta <- iris3api::zeisel_2015$meta
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
    empty_ident <- as.factor(e1$obj$orig.ident)
    levels(empty_ident) <-
      rep("empty_category", length(levels(empty_ident)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_ident, col.name = "empty_category")

    for (idx in seq_along(colnames(e1$meta))) {
      this_meta_name <- colnames(e1$meta)[idx]
      e1$obj <-
        AddMetaData(e1$obj, e1$meta[, idx], col.name = this_meta_name)
    }

    e1$obj <-
      AddMetaData(e1$obj,
        PercentageFeatureSet(e1$obj, pattern = "^MT-") + 0.001,
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
           idx = 2,
           jobid = "example",
           type = "example",
           filename,
           min_cells = 1,
           min_genes = 200,
           nVariableFeatures = 3000,
           percentMt = 5,
           removeRibosome = FALSE,
           label = NULL,
           species = "Human") {
    send_progress("Start processing scRNA-seq dataset")
    if (jobid == "example1") {
      expr_type <- as.character(expr$mimetype[1])
      expr_path <- as.character(expr$filename[1])
      raw_expr_data <- read_deepmaps(expr_type, expr_path)

      label_type <- as.character(label$mimetype[1])
      label_path <- as.character(label$filename[1])
      e1$meta <- read_deepmaps(label_type, label_path)
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
        min.cells = as.numeric(min_cells),
        min.features = as.numeric(min_genes)
      )
    e1$species <- "Human"
    empty_ident <- as.factor(e1$obj$orig.ident)
    levels(empty_ident) <-
      rep("empty_category", length(levels(empty_ident)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_ident, col.name = "empty_category")

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
        PercentageFeatureSet(e1$obj, pattern = "^MT-") + 0.001,
        col.name = "percent.mt"
      )

    Idents(e1$obj) <- e1$obj$orig.ident
    rb.genes <-
      rownames(e1$obj)[grep("^Rp[sl][[:digit:]]", rownames(e1$obj))]
    percent.ribo <-
      Matrix::colSums(e1$obj[rb.genes, ]) / Matrix::colSums(e1$obj) * 100
    send_progress("Calculating data summary statistics")
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


#' Load matched scRNA-seq and scATAC-seq data to Seurat / Signac
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
load_multiome <-
  function(req,
           idx = 0,
           filename,
           jobid = "example",
           min_cells = 1,
           min_genes = 200,
           nVariableFeatures = 2000,
           percentMt = 25,
           removeRibosome = FALSE,
           expr = NULL,
           label = NULL,
           species = "Human") {
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
        e1$obj <- qs::qread("/data/pbmc_sorted_3k.qsave")
        raw_obj <- qs::qread("/data/pbmc_match_3k.qsave")
        e1$obj@assays$ATAC@fragments[[1]]@path <-
          "/data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
      } else {
        e1$obj <-
          qs::qread(
            "C:/Users/flyku/Documents/GitHub/iris3api/inst/extdata/pbmc_sorted_3k.qsave"
          )
        raw_obj <- qs::qread("C:/Users/flyku/Documents/GitHub/iris3api/inst/extdata/pbmc_match_3k.qsave")
        #qs::qsave(e1$obj,
        #  "C:/Users/flyku/Documents/GitHub/iris3api/inst/extdata/pbmc_match_3k.qsave"
        #)
        e1$obj@assays$ATAC@fragments[[1]]@path <-
          "C:/Users/flyku/Desktop/iris3/eg2/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
        iris3api::set_embedding(name = "umap.rna")


        #e1$meta <- e1$obj@meta.data[, c("cell_type", "sex")]
        #e1$meta$disease <-
        #  rep(c("disease", "control"), nrow(e1$meta) / 2)
        # e1$embedding_idx <- which(names(e1$obj@reductions) == 'HGT')
      }
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
    } else {
      send_progress(paste0("Start processing single-cell multiome dataset"))
      expr_type <- as.character(expr$mimetype[1])
      expr_path <- as.character(expr$filename[1])
      raw_expr_data <- read_deepmaps(expr_type, expr_path)

      label <- list()
      label_type <- as.character(label$mimetype[1])
      label_path <- as.character(label$filename[1])
      e1$meta <- NULL
      if (length(label_type) > 0) {
        e1$meta <- read_deepmaps(label_type, label_path)
      }
      rna_counts <- raw_expr_data$`Gene Expression`
      atac_counts <- raw_expr_data$Peaks


      # Use peaks in standard chromosomes
      grange.counts <-
        Signac::StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
      grange.use <-
        as.character(GenomicRanges::seqnames(grange.counts)) %in% GenomeInfoDb::standardChromosomes(grange.counts)
      atac_counts <- atac_counts[as.vector(grange.use), ]

      # library(EnsDb.Hsapiens.v86)
      # annotations <-
      #  Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
      # seqlevelsStyle(annotations) <- 'UCSC'
      # genome(annotations) <- "hg38"
      # qs::qsave(annotations, "hg38_annotations.qsave")

      send_progress(paste0("Creating Seurat object"))

      if (file.exists("/data")) {
        annotations <- qs::qread("/data/hg38_annotations.qsave")
      } else {
        annotations <-
          qs::qread("C:/Users/flyku/Desktop/iris3/pbmc_match/db/hg38_annotations.qsave")
      }
      chrom_assay <- Signac::CreateChromatinAssay(
        counts = atac_counts,
        sep = c(":", "-"),
        genome = "hg38",
        min.cells = 5,
        # min.feature = 200,
        annotation = annotations,
      )

      e1$obj <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "ATAC"
      )
      exp_assay <- CreateAssayObject(counts = rna_counts)
      e1$obj[["RNA"]] <- exp_assay

      # DefaultAssay(e1$obj) <- "ATAC"
      # Need fragments
      # e1$obj <- Signac::NucleosomeSignal(object = e1$obj)
      # e1$obj <- Signac::TSSEnrichment(object = e1$obj, fast = FALSE)
      # downstream plotting of the TSS enrichment signal for different groups of cells.
      # e1$obj$high.tss <- ifelse(e1$obj$TSS.enrichment > 2, 'High', 'Low')
      # e1$obj$pct_reads_in_peaks <-  e1$obj$atac_peak_region_fragments/ e1$obj$atac_fragments*100
      # e1$obj$blacklist_ratio <-  e1$obj$blacklist_region_fragments / ( e1$obj$atac_peak_region_fragments+0.01)

      DefaultAssay(e1$obj) <- "RNA"
      e1$obj <-
        AddMetaData(e1$obj,
          PercentageFeatureSet(e1$obj, pattern = "^MT-") + 0.001,
          col.name = "percent.mt"
        )

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
      send_progress(paste0("Finding variable features"))

      e1$obj <-
        FindVariableFeatures(
          e1$obj,
          selection.method = "vst",
          nfeatures = as.numeric(2000),
          verbose = F
        )
      send_progress(paste0("Normalizing data"))
      e1$obj <- NormalizeData(e1$obj, verbose = F)
      e1$obj <- SCTransform(e1$obj, verbose = F)
      e1$obj <- ScaleData(e1$obj, verbose = F)
      # END custom
    }

    # Add empety ident
    empty_ident <- as.factor(e1$obj$orig.ident)
    levels(empty_ident) <-
      rep("empty_category", length(levels(empty_ident)))
    e1$obj <-
      AddMetaData(e1$obj, metadata = empty_ident, col.name = "empty_category")

    send_progress(paste0("Calculating data summary statistics"))
    e1$species <- species
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
