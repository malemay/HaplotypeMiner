#' Title
#'
#' Description
#'
#' Details
#'
#' @param snp_original To be completed
#' @param snp_target To be completed
#'
#' @return To be completed
#'
#' @examples
#' NULL
#'
restore_markers <- function(snp_original, snp_target) {

  kept_clusters <- unique(snp_target$Clusters[as.character(snp_target$Markers$rs)])
  kept_markers <- names(snp_target$Clusters)[snp_target$Clusters %in% kept_clusters]

  return(marker_subset(snp_original, snp_original$Markers$rs %in% kept_markers))
}
