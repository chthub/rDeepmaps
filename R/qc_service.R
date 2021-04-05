#' Get QC metrics list for single scRNA-seq dataset
#'
#' @return list of QC metrics for plotting
#' @export
#'
rna_qc_list <- function() {
  if (length(names(e1$obj@assays)) > 2) {
    DefaultAssay(e1$obj) <- "SCT"
  }

  vargenes <- VariableFeatures(e1$obj)
  Idents(e1$obj) <- e1$obj@meta.data$empty_ident
  this_obj <-
    as.matrix(GetAssayData(subset(e1$obj, features = vargenes), slot = "data"))
  row_min <- apply(this_obj, 1, min)
  row_sd <- apply(this_obj, 1, stats::sd)
  row_max <- apply(this_obj, 1, max)
  n_features_per_cell <- e1$obj$nFeature_RNA
  n_reads_per_cell <- e1$obj$nCount_RNA
  n_cells_per_gene <- apply(this_obj, 1, function(x) {
    sum(x > 0)
  })
  pct_ribo_per_gene <- e1$obj$percent.ribo
  pct_mito_per_gene <- e1$obj$percent.mt
  gene_result <-
    data.frame(
      gene = rownames(this_obj),
      mean = rowMeans(this_obj),
      std = row_sd,
      min = row_min,
      max = row_max,
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

  meta1_title <- "Cell type"
  meta1_name <- names(table(e1$obj$cell_type))
  meta1_val <- as.vector(table(e1$obj$cell_type))

  meta2_title <- "Sex"
  meta2_name <- names(table(e1$obj$sex))
  meta2_val <- as.vector(table(e1$obj$sex))

  meta3_title <- "Sample"
  meta3_name <- names(table(e1$obj$sample))
  meta3_val <- as.vector(table(e1$obj$sample))

  return(
    list(
      gene_result = gene_result,
      cell_result = cell_result,
      hist_features_per_cell = hist_features_per_cell,
      hist_reads_per_cell = hist_reads_per_cell,
      hist_cells_per_gene = hist_cells_per_gene,
      meta1_title = meta1_title,
      meta1_name = meta1_name,
      meta1_val = data.frame(name = meta1_name, value = meta1_val),
      meta2_title = meta2_title,
      meta2_name = meta2_name,
      meta2_val = data.frame(name = meta2_name, value = meta2_val),
      meta3_title = meta3_title,
      meta3_name = meta3_name,
      meta3_val = data.frame(name = meta3_name, value = meta3_val)
    )
  )
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
  if (length(names(e1$obj@assays)) > 0) {
    DefaultAssay(e1$obj) <- "ATAC"
  }

  n_features_per_cell <- e1$obj$nFeature_ATAC
  n_reads_per_cell <- e1$obj$nCount_ATAC
  pct_reads_in_peaks <- e1$obj$pct_reads_in_peaks
  atac_peak_region_fragments <- e1$obj$atac_peak_region_fragments
  blacklist_ratio <- e1$obj$blacklist_ratio
  nucleosome_signal <- e1$obj$nucleosome_signal

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
    breaks = hist(e1$obj$TSS.enrichment, plot = F)$breaks[-1],
    counts = hist(e1$obj$TSS.enrichment, plot = F)$counts
  )

  meta1_title <- "Cell type"
  meta1_name <- names(table(e1$obj$cell_type))
  meta1_val <- as.vector(table(e1$obj$cell_type))

  meta2_title <- "Sex"
  meta2_name <- names(table(e1$obj$sex))
  meta2_val <- as.vector(table(e1$obj$sex))

  meta3_title <- "Sample"
  meta3_name <- names(table(e1$obj$sample))
  meta3_val <- as.vector(table(e1$obj$sample))

  return(
    list(
      cell_result = cell_result,
      hist_features_per_cell = hist_features_per_cell,
      hist_reads_per_cell = hist_reads_per_cell,
      hist_atac_peak_region_fragments = hist_atac_fragments,
      hist_tss_enrichment = hist_tss_enrichment,
      meta1_title = meta1_title,
      meta1_name = meta1_name,
      meta1_val = data.frame(name = meta1_name, value = meta1_val),
      meta2_title = meta2_title,
      meta2_name = meta2_name,
      meta2_val = data.frame(name = meta2_name, value = meta2_val),
      meta3_title = meta3_title,
      meta3_name = meta3_name,
      meta3_val = data.frame(name = meta3_name, value = meta3_val)
    )
  )
}

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
  plot1 <- FeatureScatter(object = e1$obj, feature1 = gene1, feature2 = gene2, slot = "data")
  return(print(plot1))
}
