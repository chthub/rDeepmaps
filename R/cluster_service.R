#' Run Seurat clustering
#'
#' @param req test
#' @param nPCs test
#' @param resolution test
#'
#' @return
#' @export
#'
cluster_seurat <- function(req,
                           nPCs = 20,
                           resolution = 0.5) {
  message(glue("Run clustering. nPC={nPCs}, resolution={resolution}"))
  nPCs <- as.numeric(nPCs)
  resolution <- as.numeric(resolution)
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
  e1$obj <- FindNeighbors(e1$obj, dims = 1:nPCs, verbose = F)
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
