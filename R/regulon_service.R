#' Run Seurat clustering
#'
#' @param cluster string
#'
#' @return json
#' @export
#'
calc_regulon_network <- function(dat = "lymph", clust = "2") {
  set.seed(42)
  dt <- get(dat)
  this_ct_name <- paste0('ct', clust)
  tmp_regulon <- dt$ct_regulon

  if (e1$regulon_ident == 'other') {
    e1$regulon_ident <- 'hgt_cluster'
  }
  #e1$regulon_ident <- 'seurat_clusters'
  active_idents <-
    droplevels(as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == e1$regulon_ident)]))
  if (length(active_idents) == 0) {
    active_idents <-
      droplevels(as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == 'hgt_cluster')]))
  }

  if (length(active_idents) == 0) {
    active_idents <-
      droplevels(as.factor(e1$obj@meta.data[, which(colnames(e1$obj@meta.data) == 'seurat_clusters')]))
  }
  this_ct_idx <- which(levels(active_idents) == clust)
  all_network <- tibble::tibble()
  i = 1
  for (i in seq_along(tmp_regulon)) {
    this_tf <- unlist(strsplit(names(tmp_regulon[i]), split = "_"))[1]
    this_ct <-
      unlist(strsplit(names(tmp_regulon[i]), split = "_"))[2]
    this_ct <- gsub("ct", "", this_ct)
    if (this_ct == as.numeric(this_ct)) {
      this_ct <- as.numeric(this_ct)
    }
    if (length(tmp_regulon[[i]]) > 500) {
      max_int  <- 500
    } else {
      max_int  <- length(tmp_regulon[[i]])
    }

    this_target <-
      as.character(tmp_regulon[[i]])[sample.int(length(tmp_regulon[[i]]), max_int)]
    this_network <-
      tibble::tibble(tf = this_tf,
                     target = this_target,
                     ct = this_ct)
    all_network <- dplyr::bind_rows(all_network, this_network)
  }

  Sys.sleep(0)
  all_network <- all_network %>%
    dplyr::mutate(id = dplyr::group_indices(., tf)) %>%
    dplyr::group_by(tf) %>%
    dplyr::mutate(idx = seq_along(tf)) %>%
    dplyr::filter(ct == clust)

  g <- igraph::graph.data.frame(all_network)
  #coords <- layout.norm(layout.auto(g))
  this_tf <- as.character(unique(all_network$tf))
  this_edges <- as.data.frame(igraph::get.edgelist(g))

  this_tf_cen <- dt$TF_cen[[this_ct_name]]
  this_gene_cen <- dt$gene_cen[[this_ct_name]]
  this_node_cen <- c(this_tf_cen, this_gene_cen)
  this_node_cen <- this_node_cen[order(names(this_node_cen))]

  nodes <-
    tibble(name = names(this_node_cen),
           centrality = this_node_cen) %>%
    dplyr::mutate(
      index = seq_along(name),
      color_index = index %% 34 + 1,
      id = name,
      category = dplyr::case_when(
        name %in% this_tf ~ paste0('Regulon-', name),!name %in% this_tf ~ "Gene"
      ),
      color = dplyr::case_when(
        name %in% this_tf ~ as.character(Polychrome::palette36.colors(36)[-2])[color_index],
        !name %in% this_tf ~ 'grey',
        TRUE ~ 'grey'
      ),
      geneSymbol = dplyr::case_when(name %in% this_tf ~ paste0(name),!name %in% this_tf ~ "Gene"),
      type = dplyr::case_when(name %in% this_tf ~ "tf",!name %in% this_tf ~ "gene")
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
    dplyr::rename(tf = name)


  this_regulon <- all_network %>%
    dplyr::select(tf) %>%
    dplyr::group_by(tf) %>%
    dplyr::count() %>%
    dplyr::left_join(this_genes, by = "tf") %>%
    dplyr::left_join(this_centrality, by = "tf") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(index = dplyr::row_number())

  this_dr <- dt$dr %>%
    dplyr::filter(cluster == as.numeric(clust)) %>%
    dplyr::rename(tf = gene) %>%
    dplyr::group_by(tf) %>%
    dplyr::arrange(tf) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    tidyr::separate(tf, c("tf",'ct2'), "-") %>%
    dplyr::select(-ct2)


  this_ras <- dt$RAS[, this_ct_idx] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("tf") %>%
    dplyr::rename(ras = 2)

  this_regulon_score <- this_regulon %>%
    dplyr::left_join(this_dr, by = "tf") %>%
    dplyr::left_join(this_ras, by = "tf") %>%
    dplyr::mutate(avg_log2FC = tidyr::replace_na(avg_log2FC, NA_real_)) %>%
    dplyr::mutate(p_val_adj  = tidyr::replace_na(p_val_adj , NA_real_)) %>%
    dplyr::mutate(isCtsr = dplyr::case_when(avg_log2FC > 0.25 &
                                       p_val_adj < 0.05 ~ 'yes',
                                     T ~ 'no'))

  result <- list()
  result$idents <- levels(active_idents)
  result$nodes <- nodes
  result$edges <- this_edges
  result$regulons <- this_regulon_score

  return(result)
}

#' Run Seurat clustering
#'
#'
#' @return json
#' @export
#'
list_regulon_network <- function(dat = "lymph") {
  set.seed(42)
  dt <- get(dat)
  tmp_regulon <- dt$ct_regulon

  all_network <- tibble::tibble()
  i = 1
  for (i in seq_along(tmp_regulon)) {
    this_tf <- unlist(strsplit(names(tmp_regulon[i]), split = "_"))[1]
    this_ct <-
      unlist(strsplit(names(tmp_regulon[i]), split = "_"))[2]
    this_ct <- gsub("ct", "", this_ct)
    if (this_ct == as.numeric(this_ct)) {
      this_ct <- as.numeric(this_ct)
    }
    if (length(tmp_regulon[[i]]) > 500) {
      max_int  <- 500
    } else {
      max_int  <- length(tmp_regulon[[i]])
    }

    this_target <-
      as.character(tmp_regulon[[i]])[sample.int(length(tmp_regulon[[i]]), max_int)]
    this_network <-
      tibble::tibble(tf = this_tf,
                     target = this_target,
                     ct = this_ct)
    all_network <- dplyr::bind_rows(all_network, this_network)
  }
  all_network <-  all_network %>% tibble::rowid_to_column("id")
  return(all_network)
}

#' RAS
#'
#' @param gene string
#'
#' @return
#' @export
#'
example_cluster_coords <- function() {
  if (e1$regulon_ident == 'other') {
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
  axis <- c(paste0(embedding, "_1"),
            paste0(embedding, "_2"))
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
example_ras <- function(gene = "MAFF", assay = "RNA", clust = "0") {
  send_progress(paste0("Loading regulon: ", gene))
  if ('Gad1' %in% rownames(e1$obj)) {
    gene <- stringr::str_to_title(gene)
  }
  if (!gene %in% rownames(e1$obj)) {
    gene <- rownames(e1$obj)[1001]
  }
  this_ct_name <- paste0("ct",clust)
  this_ras_rowid <- data.frame(rowname = rownames(dt$RAS_C)) %>%
    tidyr::separate(rowname, c("tf",'ct'), "_") %>%
    tibble::rowid_to_column() %>%
    dplyr::filter(ct %in% this_ct_name & tf == gene) %>%
    dplyr::pull(rowid)

  embedding <- names(e1$obj@reductions[e1$embedding_idx])
  res1 <- data.frame(
    id = rownames(Embeddings(e1$obj, reduction = embedding)),
    dim1 = Embeddings(e1$obj, reduction = embedding)[, 1],
    dim2 = Embeddings(e1$obj, reduction = embedding)[, 2],
    dim3 = Embeddings(e1$obj, reduction = embedding)[, 2],
    expr = as.numeric(dt$RAS_C[this_ras_rowid, ]),
    index = 1
  )
  legend <- c(min(res1$expr), max(res1$expr))
  axis <- c(paste0(embedding, "_1"),
            paste0(embedding, "_2"),
            paste0(embedding, "_3"))
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
  if ('Gad1' %in% rownames(e1$obj)) {
    gene <- stringr::str_to_title(gene)
  }
  if (!gene %in% rownames(e1$obj)) {
    gene <- rownames(e1$obj)[1001]
  }
  embedding <- names(e1$obj@reductions[e1$embedding_idx])
  res1 <- data.frame(
    id = rownames(Embeddings(e1$obj, reduction = embedding)),
    dim1 = Embeddings(e1$obj, reduction = embedding)[, 1],
    dim2 = Embeddings(e1$obj, reduction = embedding)[, 2],
    dim3 = Embeddings(e1$obj, reduction = embedding)[, 2],
    expr = FetchData(object = e1$obj,
                     vars = c(gene))[, 1],
    index = 1
  )
  legend <- c(min(res1$expr), max(res1$expr))
  axis <- c(paste0(embedding, "_1"),
            paste0(embedding, "_2"),
            paste0(embedding, "_3"))
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
example_dr1 <- function(tf = c('CTCF', 'DEAF1'),
                        ct1 = c(0, 1),
                        ct2 = c(2, 3)) {
  data(dt)
  score <- dt$Dregulon[[2]]
  pvale <- dt$Dregulon[[1]]
  set.seed(42)
  rand_num <- sample.int(7, 2)
  example_ct <-
    as.numeric(stringr::str_extract(colnames(score[[1]]), "\\d+"))

  if (!ct1 %in% colnames(score[[1]])) {
    ct1 <- colnames(score[[1]])[rand_num[1]]
  }

  if (ct2 == 'all other cell types')
  {
    result <- data.frame()
    for (i in tf) {
      tmp <-
        data.frame(tf = i,
                   logfc = score[[i]][ct1, ct2],
                   pvalue = pvale[[i]][ct1, ct2])
      result <- rbind(result, tmp)
    }
  } else if (!ct2 %in% colnames(score[[1]])) {
    ct2 <- colnames(score[[1]])[rand_num[2]]
    result <- data.frame()
    for (i in tf) {
      tmp <-
        data.frame(tf = i,
                   logfc = score[[i]][ct1, ct2],
                   pvalue = pvale[[i]][ct1, ct2])
      result <- rbind(result, tmp)
    }
  }

  return (result)
}

#' Run DR
#' @param tf
#' @param ct1
#' @param ct2
#' @return
#' @export
#'
calc_dr <- function(dat = "lymph",
                       tf = c('CTCF', 'ELF1', 'MEF2C', 'E2F6', 'EGR1'),
                       ct1 = c(4),
                       ct2 = c(1)) {

  dt <- get(dat)

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


  this_tf_names <- paste0("ct", ct1)
  this_ras_rowid <- data.frame(rowname = rownames(dt$RAS_C)) %>%
    tidyr::separate(rowname, c("tf",'ct'), "_") %>%
    tibble::rowid_to_column() %>%
    dplyr::filter(ct %in% this_tf_names) %>%
    dplyr::pull(rowid)

  this_ras <- dt$RAS_C[this_ras_rowid, ]
  rownames(this_ras) <- stringr::str_remove(rownames(this_ras), "_.*")

  ras_obj <- CreateSeuratObject(dt$RAS_C)
  ras_obj <-
    AddMetaData(ras_obj, active_idents, col.name = "hgt_cluster")
  Idents(ras_obj) <- 'hgt_cluster'
  dr <-
    FindMarkers(
      ras_obj,
      ident.1 = ct1,
      ident.2 = ct2,
      logfc.threshold = 0.25,
      min.pct = 0.25,
      only.pos = T
    )
  dr <- tibble::rownames_to_column(dr, "tf") %>%
    tidyr::separate(tf, c("tf", "ct"), "-") %>%
    dplyr::filter(ct %in% this_tf_names)
  return (dr)
}

#' Run DR peak-gene figure
#' @param genes
#' @param ct1
#' @param ct2
#' @return
#' @export
#'
example_dr_figure <- function(genes = c('CTCF', 'ELF1', 'MEF2C','E2F6','EGR1'),
                       ct1 = c(4),
                       ct2 = c(1)) {
  library(Signac)
  library(tidyverse)
  library(ggpubr)
  library(GenomicRanges)
  #ct1 = c(4)
  #ct2 = c(1)
  obj2<- qs::qread("C:/Users/flyku/Documents/GitHub/iris3api/example_peak_fig.qsave")
  #genes <- dt$ct_regulon[tmp1][[4]]
  this_tf <- strsplit(names(dt$ct_regulon), split = "_")
  tmp1 <- sapply(this_tf, function(x){
    if(x[1] =="PRDM1") {
      return(T)
    } else {
      return(F)
    }
  })



  this_ct <-
    unlist(strsplit(names(dt$ct_regulon), split = "_"))
  this_ct <- gsub("ct", "", this_ct)



  DefaultAssay(obj2) <- "RNA"
  this_obj <- subset(obj2, features = genes)

  DefaultAssay(obj2) <- "ATAC"
  gene.ranges <- GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  this_ranges <- gene.ranges[which(gene.ranges$symbol %in% genes)]
  this_start <- IRanges::start(this_ranges) - 2000
  this_end <- IRanges::end(this_ranges)
  #result <- paste0("chr", this_ranges@seqnames, "-", this_start, "-", this_end)
  levels(this_ranges@seqnames) <- paste0('chr',levels(this_ranges@seqnames))
  gr2 <- GRanges(seqnames=this_ranges@seqnames, ranges=IRanges(start=this_start, end = this_end), strand="+")
  overlap1 <- findOverlaps(granges(obj2), gr2)

  obj1 <- Seurat::GetAssayData(
    object = obj2, assay = 'ATAC', slot = 'data'
  )[overlap1@from,, drop=F]

  obj1 <- CreateChromatinAssay(counts = obj1, genome = "hg38")
  obj1 <- CreateSeuratObject(counts = obj1)
  obj1 <- AddMetaData(obj1, obj2$hgt_cluster, col.name="hgt_cluster")
  Idents(obj1) <- obj1$hgt_cluster

  deg <-
    FindMarkers(
      this_obj,
      ident.1 = ct1,
      ident.2 = ct2,
      min.pct = 0,
      logfc.threshold = 0,
      assay = "RNA"
    ) %>%
  tibble::rownames_to_column("gene") %>%
    dplyr::mutate(type = 'deg',
                  isSignificant = dplyr::case_when(
                    p_val_adj < 0.05 ~ T,
                    T ~ F
                  ))%>%
    dplyr::mutate(
      isSignificant = as_factor(isSignificant))


  DefaultAssay(obj2) <- "ATAC"
  #gene.activities <- GeneActivity(obj2)
  #obj2[['gas']] <- CreateAssayObject(counts = gene.activities)
  #qs::qsave(obj2, "example_peak_fig.qsave")
  DefaultAssay(obj2) <- 'gas'
  this_gas_obj <- subset(obj2, features = genes)


  dp <-
    FindMarkers(
      this_gas_obj,
      ident.1 = ct1,
      ident.2 = ct2,
      min.pct = 0,
      logfc.threshold = 0,
      assay = "gas"
    ) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(type = 'gas',
                  isSignificant = dplyr::case_when(p_val_adj < 0.05 ~ T,
                                            T ~ F)) %>%
    dplyr::mutate(
                  isSignificant = as_factor(isSignificant))


  dp2 <-
    FindMarkers(
      obj1,
      ident.1 = ct1,
      ident.2 = ct2,
      min.pct = 0,
      logfc.threshold = 0,
      assay = "RNA"
    ) %>%
    tibble::rownames_to_column("peak") %>%
    dplyr::mutate(type = 'peak',
                  gene = genes[overlap1@to],
                  isSignificant = dplyr::case_when(p_val_adj < 0.05 ~ T,
                                            T ~ F)) %>%
    dplyr::mutate(
      isSignificant = as_factor(isSignificant))




  df <- deg %>%
    dplyr::bind_rows(dp)
  library(ggplot2)

  p1 <- ggplot(dp,
         aes(x = avg_log2FC,
             y = gene)) +
    geom_point(aes(color = isSignificant), size = 4) +
    scale_colour_manual(values = c("grey","red")) +
    theme_bw() +
    ylab("") +
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.position = "none"
    ) +
    scale_x_continuous(name = "Peak: log2FC") +
    scale_size(range = c(6, 14))

  p2 <- ggplot(deg,
         aes(x = avg_log2FC,
             y = gene)) +
    geom_point(aes(color = isSignificant), size = 4) +
    scale_colour_manual(values = c("grey","red")) +
    theme_bw() +
    ylab("") +
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text.y = element_blank(),
      legend.position = "none"
    ) +
    scale_x_continuous(name = "Gene: log2FC") +
    scale_size(range = c(6, 14))

  p3 <- ggplot(dp2,
               aes(x = avg_log2FC,
                   y = gene)) +
    geom_point(aes(color = isSignificant), size = 4) +
    scale_colour_manual(values = c("grey","red")) +
    theme_bw() +
    ylab("") +
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.position = "none"
    ) +
    scale_x_continuous(name = "Peak: log2FC") +
    scale_size(range = c(6, 14))



  figure <- ggarrange(p3, p2,
                      ncol = 2, nrow = 1)
  return (print(figure))
}


#' Run RI heatmap
#' @param tf
#' @param genes
#' @return
#' @export
#'
example_ri_heatmap2 <-
  function(tf = 'CTCF',
           genes = c("Arhgap27", "Zbtb22", "Fam126a", "Mfsd14b")) {
    #genes <- dt$ct_regulon[[1]][1:20]
    #heatmap_rownames <- paste0(tf, "_", genes)
    #
    #heat_idx <- which(rownames(dt$RI_CT) %in% heatmap_rownames)
    #heatmap_mat <-dt$RI_CT[heat_idx, ]

    send_progress(paste0("Loading regulon itensity heatmap: ", tf))
    if ('Gad1' %in% rownames(e1$obj)) {
      genes <- stringr::str_to_title(genes)
    }

    if (e1$regulon_ident == 'other') {
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

    res1 <-
      as.matrix(log1p(log1p(
        AverageExpression(e1$obj, features = genes)$RNA
      )), to = c(0.01, 1))
    #res1 <- (res1 - rowMeans(res1)) / matrixStats::rowSds(res1)
    res1 <- as.matrix((res1 > 0) + 0)
    res1 <- as.data.frame(res1)
    library(tidyverse)
    js <-
      jsonlite::fromJSON('C:/Users/flyku/Desktop/iris3/tmp/NRF1.json')
    js$cat_colors$row$`cat-0` <- "regulon"
    js$row_nodes$`cat-0` <- "regulon"
    js$col_nodes <- js$col_nodes[1:18, ]
    js$col_nodes$name <- 1:18
    js$col_nodes$`cat-0` <- 1:18
    jsonlite::write_json(js, "CTCF.json", simplifyVector = T)

    writeLines(toJSON(js, auto_unbox = T), "CTCF.json")
    library(jsonlite)

    js$mat <- js$mat[, 1:18]

    #result <- list(
    #  mat = as.matrix(res1),
    #  links = list(),
    #  views = list(
    #    N_row_sum = c("all", NA),
    #    dist = c("cos","cos")
    #    nodes = list(
    #      row_nodes = list(),
    #      col_nodes = list
    #    )
    #  )
    #)

    #jsonlite::toJSON(result)
    return(1)
  }

#' Run RI heatmap
#' @param tf
#' @param genes
#' @return
#' @export
#'
example_ri_heatmap <- function(tf = 'CTCF', genes) {
  #genes <- dt$ct_regulon[[1]][1:20]
  #heatmap_rownames <- paste0(tf, "_", genes)
  #
  #heat_idx <- which(rownames(dt$RI_CT) %in% heatmap_rownames)
  #heatmap_mat <-dt$RI_CT[heat_idx, ]

  send_progress(paste0("Loading regulon itensity heatmap: ", tf))
  if ('Gad1' %in% rownames(e1$obj)) {
    genes <- stringr::str_to_title(genes)
  }

  if (e1$regulon_ident == 'other') {
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

  res1 <-
    as.matrix(log1p(log1p(
      AverageExpression(e1$obj, features = genes)$RNA
    )), to = c(0.01, 1))
  #res1 <- (res1 - rowMeans(res1)) / matrixStats::rowSds(res1)
  res1 <- as.matrix((res1 > 0) + 0)
  res1 <- as.data.frame(res1)

  res2 <- data.frame()
  for (i in 1:ncol(res1)) {
    for (j in 1:nrow(res1)) {
      tmp <- data.frame(i - 1, j, res1[j, i])
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
