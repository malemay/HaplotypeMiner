#' Subsets a snp_data by keeping only certain markers.
#'
#' Description
#'
#' Details
#'
#' @param snp_data The snp_data to be subset.
#' @param indices Indices of the markers to be kept.
#'
#' @return To be completed
#'
#' @examples
#' NULL
#'
marker_subset <- function(snp_data, indices) {

  if(is.null(indices) || sum(indices) == 0) {
    return(NULL)
  }

  # All snp_data at least have a matrix which we must subset.
  snp_data$Genotypes <- snp_data$Genotypes[ ,indices]

  # If metadata are present, subset those.
  if(!is.null(snp_data$Markers)) {
    snp_data$Markers <- snp_data$Markers[indices, ]
  }

  # If an LD matrix is present, subset it.
  if(!is.null(snp_data$LD)) {
    snp_data$LD <- snp_data$LD[indices, indices]
  }

  if(!is.null(snp_data$VCF)) {
    snp_data$VCF <- snp_data$VCF[indices, ]
  }

  # Return the subset object.
  return(snp_data)
}
