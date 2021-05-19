## code to prepare `DATASET` dataset goes here

#' load_zeisel_2015
#' @description Zeisel 2015: Mouse Brain, 7 cell types, 3005 cells
#' @return NULL
#'

load_zeisel_2015 <- function() {
  expr <- fst::read_fst("data-raw/Zeisel_expression.fst")
  meta <- readr::read_csv("data-raw/Zeisel_index_label.csv")
  rownames(expr) <- NULL
  rownames(meta) <- NULL
  expr <- tibble::column_to_rownames(expr, "X1")
  meta <- tibble::column_to_rownames(meta, "Cell")
  zeisel_2015 <- list(expr = expr, meta = meta)
  usethis::use_data(zeisel_2015, overwrite = TRUE)
}


#' load_yan_2013
#' @description Yan 2013: Human embryo, 7 cell types, 90 cells
#' @return null
#'
load_yan_2013 <- function() {
  expr <- readr::read_csv("data-raw/Yan_2013_expression.csv")
  meta <- readr::read_csv("data-raw/Yan_2013_label.csv")
  rownames(expr) <- NULL
  rownames(meta) <- NULL
  expr <- tibble::column_to_rownames(expr, "Gene_ID")
  meta <- tibble::column_to_rownames(meta, "Cell_type")
  yan_2013 <- list(expr = expr, meta = meta)
  usethis::use_data(yan_2013, overwrite = TRUE)
}


#' load_ifnb_2800
#' @description Yan 2013: Human embryo, 7 cell types, 90 cells
#' @return null
#'
load_ifnb <- function() {
  library(Seurat)
  library(SeuratData)
  library(patchwork)
  # load dataset
  LoadData("ifnb")


  # split the dataset into a list of two seurat objects (stim and CTRL)
  ifnb.list <- SplitObject(ifnb, split.by = "stim")
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = ifnb.list)
  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
  immune.combined <- IntegrateData(anchorset = immune.anchors)
  DefaultAssay(immune.combined) <- "integrated"

  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindClusters(immune.combined, resolution = 0.4)

  p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
  p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
  p1 + p2

  immune.combined <- RenameIdents(immune.combined,
    `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
    `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
    `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC"
  )
  DimPlot(immune.combined, label = TRUE)

  random_cells <- seq(1, ncol(ifnb), 5)
  ifnb_subset <- ifnb[, random_cells]

  expr <- as.matrix(GetAssayData(ifnb_subset))
  meta <- ifnb_subset@meta.data[, c("stim", "seurat_annotations")]
  colnames(meta) <- c("sample", "cell_type")
  meta$sex <- rep(c("male", "female"), ncol(ifnb_subset) / 2)

  qs::qsave(ifnb, "ifnb.qsave")
  ifnb_2800 <- list(expr = expr, meta = meta)
  usethis::use_data(ifnb_2800, overwrite = TRUE)
}


#' load_pbmc_match_3k
#' @description pbmc_match_3k
#' @return null
#'
load_pbmc_match_3k <- function() {
  pbmc <- qs::qread("./data/pbmc_match_3k.qsave")

  empty_ident <- as.factor(pbmc$orig.ident)
  levels(empty_ident) <-
    rep("empty_ident", length(levels(empty_ident)))
  pbmc <-
    AddMetaData(pbmc, metadata = empty_ident, col.name = "empty_ident")

  pbmc@meta.data$sex <- rep(c("male", "female"), ncol(pbmc) / 2)
  pbmc@meta.data$sample <- rep(c("sample1"), ncol(pbmc) / 1)
  pbmc@meta.data$cell_type <- pbmc$predicted.id


  pbmc <-
    AddMetaData(pbmc,
      PercentageFeatureSet(pbmc, pattern = "^MT-"),
      col.name = "percent.mt"
    )

  Idents(pbmc) <- pbmc$orig.ident
  rb.genes <-
    rownames(pbmc)[grep("^R[P][[:digit:]]", rownames(pbmc))]
  percent.ribo <-
    Matrix::colSums(pbmc[rb.genes, ]) / Matrix::colSums(pbmc) * 100
  pbmc <-
    AddMetaData(pbmc, percent.ribo, col.name = "percent.ribo")


  pbmc@meta.data$sex <- rep(c("male", "female"), ncol(pbmc) / 2)
  pbmc@meta.data$sample <- rep(c("sample1"), ncol(pbmc) / 1)
  pbmc@meta.data$cell_type <- pbmc$predicted.id

  qs::qsave(pbmc, "pbmc_match_3k.qsave")
  usethis::use_data(pbmc, overwrite = TRUE)
}


#' example regulon data
#' @description pbmc_match_3k
#' @return null
#'
load_pbmc_match_3k <- function() {
  PATH <- 'C:/Users/flyku/Desktop/iris3/pbmc_match/'
  dt <- list()
  dt$RAS <- as.matrix(readRDS(paste0(PATH, "RAS.rds")))
  GAS <- as.matrix(readRDS(paste0(PATH, "GAS.rds")))
  dt$RI_CT <- as.matrix(readRDS(paste0(PATH, "RI_CT.rds")))
  dt$Dregulon <- readRDS(paste0(PATH, "Dregulon.rds"))
  dt$ct_regulon <- readRDS(paste0(PATH, "ct_regulon.rds"))
  usethis::use_data(dt, overwrite = TRUE)
}
