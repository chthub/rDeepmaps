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
                           neighbor = 20) {
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

  e1$ident_idx <-
    which(colnames(e1$obj@meta.data) == "seurat_clusters")
  Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]

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
  this_meta_name <- glue("new_merge_{e1$new_meta_counter}")
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
  e1$new_meta_counter <- e1$new_meta_counter + 1

  return(list(
    new_ident = this_meta_name,
    new_levels = levels(this_idents)
  ))
}


#' Set or select category to assign labels
#'
#' @param req request payload
#' @param categoryName string

#' @return active category name and levels, all available categories,
#' @export
#'
select_category <- function(req, categoryName = "name3") {
  if (categoryName == "") {
    this_idents <- as.factor("")
  } else if (categoryName %in% colnames(e1$obj@meta.data)) {
    this_category_idx <-
      which(colnames(e1$obj@meta.data) == categoryName)
    this_idents <- e1$obj@meta.data[, this_category_idx]
  } else {
    this_idents <- as.factor(e1$obj$empty_ident)
    levels(this_idents) <- "other"
    e1$obj <-
      AddMetaData(e1$obj, this_idents, col.name = categoryName)
    e1$obj@tools$active_category <- categoryName
    e1$obj@tools$available_category <-
      unique(c(e1$obj@tools$available_category, categoryName))
  }

  return(
    list(
      active_category = categoryName,
      active_category_levels = levels(this_idents),
      available_category = e1$obj@tools$available_category
    )
  )
}

#' Select cells based on criteria
#'
#' @param req request payload
#' @param newLevelName array
#' @param filterPayload array
#' @return filtered cell labels, obj will store current cells
#' @export
#' @examples
#' \dontrun{
#' filterPayload <- data.frame(
#'   type = c("gene", "cluster"),
#'   name = c("Gad1", NA),
#'   direction = c(">", "in"),
#'   thres = c("1", NA),
#'   category = c(NULL, "cell_type"),
#'   level = c(NA, "1_oligodendrocytes")
#' )
#' newLevelName <- "label1"
#'
#' this_filter <- list(
#'   type = "cluster",
#'   direction = "in",
#'   category = "cell_type",
#'   level = "1_oligodendrocytes"
#' )
#' }
#'
select_cells <- function(req, newLevelName = "ct1", filterPayload) {
  message(glue("Select cells..."))
  list_cells <- list()
  print(filterPayload)
  for (i in seq_len(nrow(filterPayload))) {
    this_filter <- filterPayload[i, ]
    if (this_filter$type == "gene") {
      this_cells <-
        eval(parse(
          text = paste0(
            "WhichCells(e1$obj, expression=",
            this_filter$name,
            this_filter$direction,
            this_filter$thres,
            ")"
          )
        ))
    } else if (this_filter$type == "cluster") {
      if (this_filter$direction == "in") {
        tmp_direction <- "ident"
      } else {
        tmp_direction <- "ident.remove"
      }

      # Set temp ident
      Idents(e1$obj) <-
        eval(parse(text = paste0("e1$obj$", this_filter$category)))

      this_cells <-
        eval(parse(
          text = paste0(
            "WhichCells(e1$obj, ",
            tmp_direction,
            " = '",
            this_filter$level,
            "')"
          )
        ))

      # Set back to active ident
      Idents(e1$obj) <- as.factor(e1$obj@meta.data[, e1$ident_idx])
    }
    list_cells[[i]] <- this_cells
  }

  # Intersect all filtered cells
  result_cells <- Reduce(union, list_cells)

  # Update active category
  active_cate <- e1$obj@tools$active_category
  active_cate_idx <-
    which(colnames(e1$obj@meta.data) == active_cate)

  # Apply new levels to filtered cells
  filtered_cells <- colnames(e1$obj) %in% result_cells
  tmp_ident <- as.character(e1$obj@meta.data[, active_cate_idx])
  tmp_ident[which(filtered_cells)] <- newLevelName
  e1$obj@meta.data[, active_cate_idx] <- as.factor(tmp_ident)
  active_cate_levels <- levels(e1$obj@meta.data[, active_cate_idx])
  return(
    list(
      active_category = active_cate,
      active_category_levels = active_cate_levels,
      active_category_idents <- e1$obj@meta.data[, active_cate_idx]
    )
  )
}

#' Select cells based on criteria
#'
#' @param req request payload
#' @param selectionPayload array
#' @return select usbset cell labels, obj will store current cells
#' @export

subset_cells <- function(req, selectionPayload) {
  message(glue("Subset cells..."))
  list_cells <- list()
  print(selectionPayload)

  for (i in seq_len(nrow(selectionPayload))) {
    this_filter <- selectionPayload[i, ]
    if (this_filter$type == "gene") {
      this_cells <-
        eval(parse(
          text = paste0(
            "WhichCells(e1$obj, expression=",
            this_filter$name,
            this_filter$direction,
            this_filter$thres,
            ")"
          )
        ))
    } else if (this_filter$type == "cluster") {
      tmp_direction <- "ident"

      # Set temp ident
      Idents(e1$obj) <-
        eval(parse(text = paste0("e1$obj$", this_filter$category)))

      this_cells <-
        eval(parse(
          text = paste0(
            "WhichCells(e1$obj, ",
            tmp_direction,
            " = '",
            this_filter$level,
            "')"
          )
        ))

      # Set back to active ident
      Idents(e1$obj) <- as.factor(e1$obj@meta.data[, e1$ident_idx])
    }
    list_cells[[i]] <- this_cells
  }

  # Intersect all filtered cells
  result_cells <- Reduce(union, list_cells)

  # Perform backup, subset, reset
  e1$full_obj <- e1$obj
  e1$sub_obj <- subset(e1$obj, cells = result_cells)
  e1$obj <- e1$sub_obj

  return(
    list(
      active_category = e1$obj@tools$active_category,
      n_cells = ncol(e1$obj)
    )
  )
}


#' Plot static UMAP
#' @param req request payload
#' @param categoryName string
#'
#' @return
#' @export

umap_plot <- function(req, categoryName = "seurat_clusters") {
  print(categoryName)
  this_ident_idx <-
    which(colnames(e1$obj@meta.data) == categoryName)

  Idents(e1$obj) <- e1$obj@meta.data[, this_ident_idx]
  plot <- DimPlot(e1$obj, reduction = "umap")

  Idents(e1$obj) <- e1$obj@meta.data[, e1$ident_idx]
  return(print(plot))
}
