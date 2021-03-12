#' Get QC metrics list
#'
#' @return list of QC metrics for plotting
#' @export
#'
qc_list <- function() {
  vargenes <- VariableFeatures(e1$obj)
  Idents(e1$obj) <- e1$obj@meta.data$empty_ident
  this_obj <-
    as.matrix(GetAssayData(subset(e1$obj, features = vargenes), slot = "data"))
  row_min <- apply(this_obj, 1, min)
  row_sd <- apply(this_obj, 1, stats::sd)
  row_max <- apply(this_obj, 1, max)
  n_genes_per_cell <- e1$obj$nFeature_RNA
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
      n_genes_per_cell = n_genes_per_cell,
      pct_ribo_per_gene = pct_ribo_per_gene,
      pct_mito_per_gene = pct_mito_per_gene
    )
  hist_genes_per_cell <- data.frame(
    breaks = hist(n_genes_per_cell, plot = F)$breaks[-1],
    counts = hist(n_genes_per_cell, plot = F)$counts
  )
  hist_reads_per_cell <- data.frame(
    breaks = hist(n_reads_per_cell, plot = F)$breaks[-1],
    counts = hist(n_reads_per_cell, plot = F)$counts
  )
  hist_cells_per_gene <- data.frame(
    breaks = hist(n_cells_per_gene, plot = F)$breaks[-1],
    counts = hist(n_cells_per_gene, plot = F)$counts
  )
  return(
    list(
      gene_result,
      cell_result,
      hist_genes_per_cell,
      hist_reads_per_cell,
      hist_cells_per_gene
    )
  )
}
