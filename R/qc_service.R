#' Get QC metrics list for single scRNA-seq dataset
#'
#' @return list of QC metrics for plotting
#' @export
#'
rna_qc_list <- function() {
  TOTAL_STEPS <- 6
  send_progress(paste0("Calculating scRNA-seq QC metrics"))
  if (length(names(e1$obj@assays)) > 2) {
    DefaultAssay(e1$obj) <- "SCT"
  }

  vargenes <- VariableFeatures(e1$obj)
  Idents(e1$obj) <- e1$obj@meta.data$empty_ident
  this_obj <-
    as.matrix(GetAssayData(subset(e1$obj, features = vargenes), slot = "data"))
  row_min <- apply(this_obj, 1, min)
  row_mean <- HVFInfo(e1$obj)[VariableFeatures(e1$obj), ][, 1]
  row_sd <- HVFInfo(e1$obj)[VariableFeatures(e1$obj), ][, 2]
  row_residual_variance <-
    HVFInfo(e1$obj)[VariableFeatures(e1$obj), ][, 3]
  row_max <- apply(this_obj, 1, max)
  n_features_per_cell <- e1$obj$nFeature_RNA
  n_reads_per_cell <- e1$obj$nCount_RNA
  n_cells_per_gene <- apply(this_obj, 1, function(x) {
    sum(x > 0)
  })
  pct_ribo_per_gene <- e1$obj$percent.ribo
  pct_mito_per_gene <- e1$obj$percent.mt
  HVFInfo(e1$obj)[VariableFeatures(e1$obj), ]
  gene_result <-
    data.frame(
      gene = rownames(this_obj),
      min = format(round(row_min, 2), nsmall = 2),
      max = format(round(row_max, 2), nsmall = 2),
      mean = format(round(rowMeans(this_obj), 2), nsmall = 2),
      std = format(round(row_sd, 2), nsmall = 2),
      residual_variance = format(round(row_residual_variance, 2), nsmall = 2),
      n_cells_per_gene = n_cells_per_gene
    )
  cell_result <-
    data.frame(
      n_reads_per_cell = n_reads_per_cell,
      n_features_per_cell = n_features_per_cell,
      pct_ribo_per_gene = pct_ribo_per_gene,
      pct_mito_per_gene = pct_mito_per_gene
    )
  hist_features_per_cell <- data.frame(
    breaks = hist(n_features_per_cell, plot = F)$breaks[-1],
    counts = hist(n_features_per_cell, plot = F)$counts
  )
  hist_reads_per_cell <- data.frame(
    breaks = hist(n_reads_per_cell, plot = F)$breaks[-1],
    counts = hist(n_reads_per_cell, plot = F)$counts
  )
  hist_cells_per_gene <- data.frame(
    breaks = hist(n_cells_per_gene, plot = F)$breaks[-1],
    counts = hist(n_cells_per_gene, plot = F)$counts
  )


  meta_list <- list()
  for (idx in seq_along(colnames(e1$meta))) {
    this_title <- colnames(e1$meta)[idx]
    this_name <- names(table(e1$meta[, idx]))
    this_val <- as.vector(table(e1$meta[, idx]))
    tmp <-
      list(
        title = this_title,
        name = this_name,
        val = data.frame(name = this_name, value = this_val)
      )
    meta_list <- append(meta_list, list(tmp))
  }

  # jsonlite::toJSON(meta_list)

  result <- list(
    gene_result = gene_result,
    cell_result = cell_result,
    hist_features_per_cell = hist_features_per_cell,
    hist_reads_per_cell = hist_reads_per_cell,
    hist_cells_per_gene = hist_cells_per_gene,
    meta_list = meta_list
  )
  return(result)
}

#' Variable genes scatter plot
#'
#' @return static plot
#' @export
#'
rna_qc_plot <- function() {
  if ("ATAC" %in% names(e1$obj@assays)) {
    DefaultAssay(e1$obj) <- "SCT"
  }

  top10 <- head(VariableFeatures(e1$obj), 10)
  plot1 <- VariableFeaturePlot(e1$obj)
  plot2 <- LabelPoints(
    plot = plot1,
    points = top10,
    repel = TRUE
  )
  return(print(plot2))
}


#' Get QC metrics list for single scRNA-seq dataset
#'
#' @return list of QC metrics for plotting
#' @export
#'
atac_qc_list <- function() {
  TOTAL_STEPS <- 6
  send_progress(paste0("Calculating scATAC-seq QC metrics"))
  if (length(names(e1$obj@assays)) > 0) {
    DefaultAssay(e1$obj) <- "ATAC"
  }

  n_features_per_cell <- e1$obj$nFeature_ATAC
  n_reads_per_cell <- e1$obj$nCount_ATAC

  # pct_reads_in_peaks <- e1$obj$pct_reads_in_peaks
  # atac_peak_region_fragments <- e1$obj$atac_peak_region_fragments
  # blacklist_ratio <- e1$obj$blacklist_ratio
  # nucleosome_signal <- e1$obj$nucleosome_signal
  # tss_enrichment <- e1$obj$TSS.enrichment
  pct_reads_in_peaks <- n_features_per_cell
  atac_peak_region_fragments <- n_features_per_cell
  blacklist_ratio <- n_features_per_cell
  nucleosome_signal <- n_features_per_cell
  tss_enrichment <- n_features_per_cell

  cell_result <-
    data.frame(
      n_reads_per_cell = n_reads_per_cell,
      n_features_per_cell = n_features_per_cell,
      pct_reads_in_peaks = pct_reads_in_peaks,
      atac_peak_region_fragments = atac_peak_region_fragments,
      blacklist_ratio = blacklist_ratio,
      nucleosome_signal = nucleosome_signal
    )
  hist_features_per_cell <- data.frame(
    breaks = hist(n_features_per_cell, plot = F)$breaks[-1],
    counts = hist(n_features_per_cell, plot = F)$counts
  )
  hist_reads_per_cell <- data.frame(
    breaks = hist(n_reads_per_cell, plot = F)$breaks[-1],
    counts = hist(n_reads_per_cell, plot = F)$counts
  )
  hist_atac_fragments <- data.frame(
    breaks = hist(atac_peak_region_fragments, plot = F)$breaks[-1],
    counts = hist(atac_peak_region_fragments, plot = F)$counts
  )

  hist_tss_enrichment <- data.frame(
    breaks = hist(tss_enrichment, plot = F)$breaks[-1],
    counts = hist(tss_enrichment, plot = F)$counts
  )

  meta_list <- list()
  for (idx in seq_along(colnames(e1$meta))) {
    this_title <- colnames(e1$meta)[idx]
    this_name <- names(table(e1$meta[, idx]))
    this_val <- as.vector(table(e1$meta[, idx]))
    tmp <-
      list(
        title = this_title,
        name = this_name,
        val = data.frame(name = this_name, value = this_val)
      )
    meta_list <- append(meta_list, list(tmp))
  }

  result <- list(
    cell_result = cell_result,
    hist_features_per_cell = hist_features_per_cell,
    hist_reads_per_cell = hist_reads_per_cell,
    hist_atac_peak_region_fragments = hist_atac_fragments,
    hist_tss_enrichment = hist_tss_enrichment,
    meta_list = meta_list
  )
  return(result)
}
