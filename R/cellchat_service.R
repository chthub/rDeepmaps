#' Run CellChat
#'
#' @param active_idents string
#'
#' @return json
#' @export
#'
run_cellchat <- function(job = "single_rna_example", active_idents = "hgt_cluster") {
  if(job == "single_rna_example" &
     active_idents =="seurat_clusters") {
    cellchat <- iris3api::cellchat_multiome_hgt
    # cellchat_single_cluster <- cellchat
    # usethis::use_data(cellchat_single_cluster)
  } else if(job == "single_rna_example" &
       active_idents =="Label") {
    cellchat <- iris3api::cellchat_single_Label
    # cellchat_single_Label <- cellchat
    # usethis::use_data(cellchat_single_Label)

    # cellchat_multiome_hgt <- cellchat
    # usethis::use_data(cellchat_multiome_hgt)
  } else{
    tictoc::tic()
    cellchat <-
      CellChat::createCellChat(
        object = e1$obj,
        meta = e1$obj@meta.data,
        group.by = active_idents
      )

    if (e1$species == "Human") {
      CellChatDB <- CellChat::CellChatDB.human
      ppi <- CellChat::PPI.human
    } else {
      CellChatDB <- CellChat::CellChatDB.mouse
      ppi <- CellChat::PPI.mouse
    }
    all_dbs <- dplyr::glimpse(CellChatDB$interaction)
    all_dbs_names <- unique(all_dbs$annotation)

    CellChatDB.use <-
      CellChat::subsetDB(CellChatDB, search = all_dbs_names)
    cellchat@DB <- CellChatDB.use
    cellchat <- CellChat::subsetData(cellchat)
    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
    cellchat <-
      CellChat::identifyOverExpressedInteractions(cellchat)
    cellchat <- CellChat::projectData(cellchat, ppi)
    cellchat <- CellChat::computeCommunProb(cellchat)
    cellchat <-
      CellChat::filterCommunication(cellchat, min.cells = 10)
    cellchat <- CellChat::computeCommunProbPathway(cellchat)
    cellchat <- CellChat::aggregateNet(cellchat)
    tictoc::toc()
  }

  all_idents <- levels(as.factor(cellchat@idents))
  pathway_list <- cellchat@netP$pathways
  enrichedLR <-
    CellChat::extractEnrichedLR(cellchat, signaling = pathway_list, geneLR.return = FALSE)[,1]
  df.netp <-
    CellChat::subsetCommunication(cellchat, slot.name = "netP")
  df.net <-
    CellChat::subsetCommunication(cellchat, slot.name = "net")
  res <- list(
    idents = all_idents,
    enrichedLR = enrichedLR,
    pathway = pathway_list,
    netp = df.netp,
    net = df.net
  )

  e1$cellchat <- cellchat
  return(res)
}


#' Run CellChat
#'
#' @param db string
#'
#' @return json
#' @export
#'
plot_cellchat <-
  function(mode = "single_pathway",
           slot = "net",
           pathway_show = c('AGRN'),
           lr_show = "WNT4_FZD3_LRP6",
           ident1 = c(2),
           ident2 = c(3)) {
    # slot: L-R: net; pathway: netP
    # mode: single_pathway
    # circle type: chord, circle
    par(mfrow = c(1, 1))

    if (mode == "single_pathway") {
      result <-
        CellChat::netVisual_aggregate(e1$cellchat, signaling = pathway_show, layout = "circle")
    } else if (mode == "all_network") {
      result <- CellChat::netVisual_chord_gene(
        e1$cellchat,
        sources.use = ident1,
        targets.use = ident2,
        slot.name = slot
      )
    } else if (mode == "some_network") {
      result <- CellChat::netVisual_chord_gene(
        e1$cellchat,
        sources.use = ident1,
        targets.use = ident2,
        signaling = pathway_show,
        slot.name = slot
      )
    } else if (mode == "single_lr_pathway") {
      result <-
        CellChat::netVisual_individual(
          e1$cellchat,
          signaling = pathway_show,
          pairLR.use = lr_show,
          layout = "circle"
        )
    } else if (mode == "aggregated_network_weight") {
      mat <- e1$cellchat@net$weight
      mat2 <-
        matrix(
          0,
          nrow = nrow(mat),
          ncol = ncol(mat),
          dimnames = dimnames(mat)
        )
      mat2[ident1,] <- mat[ident1,]
      result <- CellChat::netVisual_circle(
        mat2,
        vertex.weight = as.numeric(table(as.factor(
          e1$cellchat@idents
        ))),
        weight.scale = T,
        label.edge = F,
        title.name = paste0(rownames(mat)[ident1], ": number of interactions")
      )
    } else if (mode == "aggregated_network_count") {
      mat <- e1$cellchat@net$count
      mat2 <-
        matrix(
          0,
          nrow = nrow(mat),
          ncol = ncol(mat),
          dimnames = dimnames(mat)
        )
      mat2[ident1,] <- mat[ident1,]
      result <- CellChat::netVisual_circle(
        mat2,
        vertex.weight = as.numeric(table(as.factor(
          e1$cellchat@idents
        ))),
        weight.scale = T,
        label.edge = F,
        title.name = paste0(rownames(mat)[ident1], ": interaction weights/strength")
      )
    }

    return(print(result))
  }
