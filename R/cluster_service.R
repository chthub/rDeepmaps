#' Run Seurat clustering
#'
#' @param req request payload
#' @param nPCs string
#' @param resolution string
#' @param neighbor string
#'
#' @return
#' @export
#'
cluster_seurat <- function(req,
                           nPCs = 20,
                           resolution = 0.5,
                           neighbor = 25) {
  message(glue(
    "Run clustering. nPC={nPCs}, resolution={resolution}, neighbor={neighbor}"
  ))
  nPCs <- as.numeric(nPCs)
  resolution <- as.numeric(resolution)
  neighbor <- as.numeric(neighbor)
  e1$obj <- NormalizeData(e1$obj, verbose = F)
  e1$obj <-
    ScaleData(e1$obj, features = rownames(e1$obj), verbose = F)
  variable_genes <- VariableFeatures(e1$obj)

  e1$obj <-
    RunPCA(e1$obj,
      features = variable_genes,
      npcs = nPCs,
      verbose = F
    )
  e1$obj <-
    FindNeighbors(e1$obj,
      dims = 1:nPCs,
      k.param = neighbor,
      verbose = F
    )
  e1$obj <-
    FindClusters(e1$obj, resolution = resolution, verbose = F)
  e1$obj <- RunUMAP(
    e1$obj,
    dims = 1:nPCs,
    n.neighbors = 15,
    verbose = F
  )
  return(list(
    n_seurat_clusters = length(levels(e1$obj$seurat_clusters)),
    umap_pts = data.frame(
      umap1 = as.vector(Embeddings(e1$obj, reduction = "umap")[, 1]),
      umap2 = as.vector(Embeddings(e1$obj, reduction = "umap")[, 2])
    )
  ))
}

#' Return active cell idents labels
#'
#' @return
#' @export
#'
active_label <- function() {
  return(1)
}

#' Merge active clusters
#' @param req request payload
#' @param newClusterIds array
#' @return array levels of renamed new idents
#' @export
#'
merge_idents <- function(req, newClusterIds) {
  print(newClusterIds)
  message(glue("Renaming idents: {e1$ident_idx} at ID: {newClusterIds} "))
  this_meta_name <- glue("new_ident_{e1$meta_counter}")
  this_idents <- as.factor(e1$obj@meta.data[, e1$ident_idx])
  this_idents_levels <- levels(this_idents)
  for (i in seq_along(this_idents_levels)) {
    if (this_idents_levels[i] %in% newClusterIds) {
      this_idents_levels[i] <- newClusterIds[1]
    }
  }
  levels(this_idents) <- this_idents_levels
  e1$obj <-
    AddMetaData(e1$obj, this_idents, col.name = this_meta_name)
  e1$meta_counter <- e1$meta_counter + 1

  return(list(
    new_ident = this_meta_name,
    new_levels = levels(this_idents)
  ))
}
