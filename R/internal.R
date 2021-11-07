#' Get genomic range of a gene and extend a flank region (Multiome project)
#' @importFrom GenomicFeatures genes
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom EnsDb.Mmusculus.v79 EnsDb.Mmusculus.v79
#' @importFrom IRanges start
#' @importFrom IRanges end
#' @param gene string
#' @param flank string
#' @return
#'
get_gene_range <- function(gene = "GAD1", flank = "10000") {
  if (e1$species == "Human") {
    gene.ranges <- GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  } else {
    gene.ranges <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
  }
  this_ranges <- gene.ranges[which(gene.ranges$symbol == gene)]
  this_start <- IRanges::start(this_ranges) - as.numeric(flank)
  this_end <- IRanges::end(this_ranges) + as.numeric(flank)
  result <- paste0("chr", this_ranges@seqnames, "-", this_start, "-", this_end)
  return(result[1])
}
