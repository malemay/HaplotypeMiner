#' Generate clusters of markers in high linkage disequilibrium
#'
#' This function clusters groups of markers that are in high linkage
#' disequilibrium on a single side of the central gene position considered.
#' The objective of the clustering step is to reduce the number of markers
#' considered and thus allow for a simplification of haplotype definition
#' problem.
#'
#' This function is used internally by the package and is not normally
#' available for end users.
#'
#' @param snp_data A list: the snp data to be clustered. Must contain at least
#'   \code{Genotypes}, \code{LD} and \code{Markers} fields for clustering to
#'   occur properly.
#' @param center_pos The central position on either side of which clustering
#'   should be performed.
#' @param block_threshold The minimum R^2 for two adjacent markers to be
#'   eligible for blocking (defaults to 1, i.e. perfect linkage disequilibrium).
#'
#' @return A subset of the original \code{snp_data} containing only one
#'  representative marker per LD block.
#'
#' @examples
#' NULL
#'
cluster <- function(snp_data, center_pos, block_threshold = 1) {

  # Split the markers into 5' and 3' groups
  split_snp <- list(five_prime  = marker_subset(snp_data, snp_data$Markers$pos < center_pos),
                    three_prime = marker_subset(snp_data, snp_data$Markers$pos >= center_pos))

  # Generate a 5' block of markers in high linkage disequilibrium
  five_block  <- subcluster(snp_data = split_snp$five_prime,
                            block_threshold = block_threshold,
                            reverse = TRUE)

  # Doing the same thing with 3' markers
  three_block <- subcluster(snp_data = split_snp$three_prime,
                            block_threshold = block_threshold,
                            reverse = FALSE)

  # "Genotypes" and "Markers" elements of the two lists are combined
  clustered_markers <-
    list(Genotypes = BiocGenerics::cbind(five_block$Genotypes, three_block$Genotypes),
         Markers = rbind(five_block$Markers, three_block$Markers))

  # The LD matrix is retrieved from the unclustered data and subset
  kept_markers <- BiocGenerics::colnames(clustered_markers$Genotypes)
  clustered_markers$LD <- snp_data$LD[rownames(snp_data$LD) %in% kept_markers,
                                      colnames(snp_data$LD) %in% kept_markers]

  # Also combine VCF data, if found
  if(!is.null(five_block$VCF) && !is.null(three_block$VCF)) {
    clustered_markers$VCF <- BiocGenerics::rbind(five_block$VCF, three_block$VCF)
  }

  # Also combining clusters
  clustered_markers$Clusters <- c(five_block$Clusters, three_block$Clusters + max(five_block$Clusters) + 1)

  return(clustered_markers)
}
