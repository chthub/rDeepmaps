#' Run Seurat clustering
#'
#' @param cluster string
#'
#' @return json
#' @export
#'
example_regulon_network <- function(cluster = "0") {
  data(dt)
  set.seed(42)
  send_progress("Start calculation")
  Sys.sleep(1)
  send_progress("Calculating regulons")
  tmp_regulon <- dt$ct_regulon
  rand_num <- sample.int(8,1) - 1
  if(!cluster %in% 0:7) {
    cluster <- rand_num
  }
  if(e1$regulon_ident == 'other') {
    e1$regulon_ident <- 'hgt_cluster'
  }
  #e1$regulon_ident <- 'seurat_clusters'
  active_idents <-
    as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == e1$regulon_ident)])
  if (length(active_idents) == 0) {
    active_idents <-
      as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == 'hgt_cluster')])
  }

  if (length(active_idents) == 0) {
    active_idents <-
      as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == 'seurat_clusters')])
  }

  all_network <- tibble::tibble()
  i=1
  for (i in seq_along(tmp_regulon)) {
    this_tf <- unlist(strsplit(names(tmp_regulon[i]), split = "_"))[1]
    this_ct <- unlist(strsplit(names(tmp_regulon[i]), split = "_"))[2]
    this_ct <- gsub("ct","", this_ct)
    if(this_ct == as.numeric(this_ct)) {
      this_ct <- as.numeric(this_ct) - 1
    }
    if(length(tmp_regulon[[i]]) > 300) {
      max_int  <- 300
    } else {
      max_int  <- length(tmp_regulon[[i]])
    }

    this_target <- as.character(tmp_regulon[[i]])[sample.int(length(tmp_regulon[[i]]), max_int)]
    this_network <-
      tibble::tibble(tf = this_tf, target = this_target, ct = this_ct)
    all_network <- dplyr::bind_rows(all_network, this_network)
  }


  Sys.sleep(0)
  send_progress("Calculating regulon intensity")
  all_network <- all_network %>%
    dplyr::mutate(id = dplyr::group_indices(., tf)) %>%
    dplyr::group_by(tf) %>%
    dplyr::mutate(idx = seq_along(tf)) %>%
    dplyr::filter(ct == cluster)

  g <- igraph::graph.data.frame(all_network)

  #coords <- layout.norm(layout.auto(g))
  this_tf <- as.character(unique(all_network$tf))
  this_edges <- as.data.frame(igraph::get.edgelist(g))


  Sys.sleep(0)
  send_progress("Calculating regulon networks")
  Sys.sleep(0)
  nodes <-
    tibble(
      name = as.character(igraph::V(g)$name),
      centrality = as.numeric(scales::rescale(igraph::eigen_centrality(g)$vector, to = c(0.1, 1)))
    ) %>%
    dplyr::mutate(
      index = seq_along(name),
      color_index = index %% 34 + 1,
      id = name,
      category = dplyr::case_when(
        name %in% this_tf ~ paste0('Regulon-', name),
        !name %in% this_tf ~ "Gene"
      ),
      color = dplyr::case_when(
        name %in% this_tf ~ as.character(Polychrome::palette36.colors(36)[-2])[color_index],!name %in% this_tf ~ 'grey',
        TRUE ~ 'grey'
      ),
      geneSymbol = dplyr::case_when(name %in% this_tf ~ paste0(name), !name %in% this_tf ~ "Gene"),
      type = dplyr::case_when(name %in% this_tf ~ "tf", !name %in% this_tf ~ "gene")
    )

  colnames(this_edges) <- c("source", "target")

  this_edges <- this_edges %>%
    dplyr::mutate(id = paste0(source, "-", target)) %>%
    unique()

  this_genes <- all_network %>%
    dplyr::group_by(tf) %>%
    dplyr::mutate(genes = paste0(target, collapse = ",")) %>%
    dplyr::select(tf, genes, ct) %>%
    unique()

  this_centrality <- nodes %>%
    dplyr::select(name, centrality) %>%
    dplyr::rename(tf=name)


  this_regulon <- all_network %>%
    dplyr::select(tf) %>%
    dplyr::group_by(tf) %>%
    dplyr::count() %>%
    dplyr::left_join(this_genes, by = "tf") %>%
    dplyr::left_join(this_centrality, by = "tf") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(index = dplyr::row_number())

  this_dr <- dt$dr[which(dt$dr$gene %in% this_regulon$tf),] %>%
    dplyr::rename(tf=gene) %>%
    dplyr::group_by(tf) %>%
    dplyr::arrange(tf) %>%
    dplyr::filter(dplyr::row_number()==1)

  this_vr <- dt$VR %>%
    as.data.frame() %>%
    tibble::rownames_to_column('tf')

  this_regulon_score <- this_regulon %>%
    dplyr::left_join(this_dr, by="tf") %>%
    dplyr::left_join(this_vr, by="tf")


  result <- list()
  result$idents <- levels(active_idents)
  result$nodes <- nodes
  result$edges <- this_edges
  result$regulons <- this_regulon_score

  #library(jsonlite)
  #write(toJSON(this_regulon), paste0("C:/Users/flyku/Documents/GitHub/iris3-frontend/static/json/regulon/example_regulon.json"))
  #write(toJSON(nodes), "C:/Users/flyku/Documents/GitHub/iris3-frontend/static/json/regulon/example_cyto_nodes.json")
  #write(toJSON(this_edges), "C:/Users/flyku/Documents/GitHub/iris3-frontend/static/json/regulon/example_cyto_edges.json")
  return(result)
}

#' RAS
#'
#' @param gene string
#'
#' @return
#' @export
#'
example_cluster_coords <- function() {
  if(e1$regulon_ident == 'other') {
    e1$regulon_ident <- 'hgt_cluster'
  }
  active_idents <-
    as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == e1$regulon_ident)])
  if (length(active_idents) == 0) {
    active_idents <-
      as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == 'hgt_cluster')])
  }
  if (length(active_idents) == 0) {
    active_idents <-
      as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == 'seurat_clusters')])
  }

  embedding <- names(e1$obj@reductions[e1$embedding_idx])

  Idents(e1$obj) <- active_idents

  res1 <- data.frame(
    id = rownames(Embeddings(e1$obj, reduction = embedding)),
    dim1 = Embeddings(e1$obj, reduction = embedding)[, 1],
    dim2 = Embeddings(e1$obj, reduction = embedding)[, 2],
    label = as.character(Idents(e1$obj)),
    index = as.integer(Idents(e1$obj)) - 1
  )
  legend <- levels(Idents(e1$obj))
  axis <- c(
    paste0(embedding, "_1"),
    paste0(embedding, "_2")
  )
  dimension <- colnames(res1)
  coords <- as.matrix(res1)
  coords[, 2] <- as.numeric(coords[, 2])
  coords[, 3] <- as.numeric(coords[, 3])
  return(list(
    axis = axis,
    legend = legend,
    dimension = dimension,
    embedding = coords
  ))

}
#' RAS
#'
#' @param gene string
#'
#' @return
#' @export
#'
example_ras <- function(gene = "MAFF", assay = "RNA") {
  send_progress(paste0("Loading regulon: ", gene))
  if('Gad1' %in% rownames(e1$obj)) {
    gene <- stringr::str_to_title(gene)
  }
  if(!gene %in% rownames(e1$obj)) {
    gene <- rownames(e1$obj)[1001]
  }
  embedding <- names(e1$obj@reductions[e1$embedding_idx])
  res1 <- data.frame(
    id = rownames(Embeddings(e1$obj, reduction = embedding)),
    dim1 = Embeddings(e1$obj, reduction = embedding)[, 1],
    dim2 = Embeddings(e1$obj, reduction = embedding)[, 2],
    dim3 = Embeddings(e1$obj, reduction = embedding)[, 2],
    expr = log1p(FetchData(
      object = e1$obj,
      vars = c(gene)
    )[, 1])*2,
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

#' GAS
#'
#' @param gene string
#'
#' @return
#' @export
#'
example_gas <- function(gene = "Gad1", assay = "RNA") {
  if('Gad1' %in% rownames(e1$obj)) {
    gene <- stringr::str_to_title(gene)
  }
  if(!gene %in% rownames(e1$obj)) {
    gene <- rownames(e1$obj)[1001]
  }
  embedding <- names(e1$obj@reductions[e1$embedding_idx])
  res1 <- data.frame(
    id = rownames(Embeddings(e1$obj, reduction = embedding)),
    dim1 = Embeddings(e1$obj, reduction = embedding)[, 1],
    dim2 = Embeddings(e1$obj, reduction = embedding)[, 2],
    dim3 = Embeddings(e1$obj, reduction = embedding)[, 2],
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


#' Run DR
#' @param tf
#' @param ct1
#' @param ct2
#' @return
#' @export
#'
example_dr <- function(tf='CTCF', ct1=1, ct2='all other cell types') {
  data(dt)
  score<-dt$Dregulon[[2]]
  pvale<-dt$Dregulon[[1]]
  set.seed(42)
  rand_num <- sample.int(7,2)
  if(!ct1 %in% colnames(score[[1]])) {
    ct1 <- colnames(score[[1]])[rand_num[1]]
  }
  if (ct2 == 'all other cell types')
  {
    result <- data.frame()
    for(i in tf) {
      tmp <- data.frame(tf=i, logfc=score[[i]][ct1,ct2],pvalue=pvale[[i]][ct1,ct2])
      result <- rbind(result, tmp)
    }
  } else if (!ct2 %in% colnames(score[[1]])) {
    ct2 <- colnames(score[[1]])[rand_num[2]]
    result <- data.frame()
    for(i in tf) {
      tmp <- data.frame(tf=i, logfc=score[[i]][ct1,ct2],pvalue=pvale[[i]][ct1,ct2])
      result <- rbind(result, tmp)
    }
  }



  return (result)
}

#' Run RI heatmap
#' @param tf
#' @param genes
#' @return
#' @export
#'
example_ri_heatmap <- function(tf='CTCF', genes) {

  #genes <- dt$ct_regulon[[1]][1:20]
  #heatmap_rownames <- paste0(tf, "_", genes)
  #
  #heat_idx <- which(rownames(dt$RI_CT) %in% heatmap_rownames)
  #heatmap_mat <-dt$RI_CT[heat_idx, ]

  send_progress(paste0("Loading regulon itensity heatmap: ", tf))
  if('Gad1' %in% rownames(e1$obj)) {
    genes <- stringr::str_to_title(genes)
  }

  if(e1$regulon_ident == 'other') {
    e1$regulon_ident <- 'hgt_cluster'
  }
  active_idents <-
    as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == e1$regulon_ident)])
  if (length(active_idents) == 0) {
    active_idents <-
      as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == 'hgt_cluster')])
  }
  if (length(active_idents) == 0) {
    active_idents <-
      as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == 'seurat_clusters')])
  }
  embedding <- names(e1$obj@reductions[e1$embedding_idx])
  Idents(e1$obj) <- active_idents

  res1 <- as.matrix(log1p(log1p(AverageExpression(e1$obj, features = genes)$RNA)),to = c(0.01, 1))
  res1 <- (res1 - rowMeans(res1)) / matrixStats::rowSds(res1)
  res1 <- as.data.frame(res1)

  res2 <- data.frame()
  for(i in 1:ncol(res1)) {
    for(j in 1:nrow(res1)) {
      tmp <- data.frame(i-1, j, res1[j,i])
      res2 <- rbind(res2, tmp)
    }
  }
  res2 <- as.matrix(res2)
  legend <- c(min(res1), max(res1))

  result <- list(
    column = genes,
    row = levels(Idents(e1$obj)),
    data = res2,
    legend = legend
  )
  #jsonlite::toJSON(result)
  return(result)
}


