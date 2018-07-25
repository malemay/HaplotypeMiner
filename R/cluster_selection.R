#' Selecting marker clusters according to filtering criteria
#'
#' This function selects marker clusters to be used for haplotype generation.
#' Clusters are selecte based on three different analysis criteria : high enough
#' linkage disequilibrium between marker pairs (parameter \code{ld_threshold}),
#' proximity between markers in a pair (parameter \code{max_pair_distance}),
#' and proximity between markers and gene center position (parameter
#' \code{max_marker_distance}).
#'
#' This function is meant for internal package use only. It is called by function
#' \code{\link{haplo_selection}} in the final stage of marker selection.
#'
#' @param snp_data A list describing snp_data and containing at least two items :
#'   \itemize{
#'   \item{Markers} {A data.frame containing metadata regarding the markers to be
#'   selected, as output by functions \code{\link{read_hapmap}} and \code{\link{read_vcf}}
#'   }.
#'   \item{LD} {A square matrix containing linkage disequilibrum measures (could be either
#'   R2, R2v, R2s or R2vs) between all gene pairs.}
#'   }
#' @param center_pos A numeric of length one. The central position of the gene of interest.
#' @param ld_threshold A numeric between 0 and 1. The minimum LD value to pass the filter. Defaults to 0.5.
#' @param max_pair_distance A numeric of length one. The maximum distance (in base pairs)
#'   between two markers in a pair for markers in that pair to be selected.
#' @param max_marker_distance A numeric of length one. The maximum distance (in
#'   base pairs) between the marker and the gene center for a marker to be selected.
#'
#' @return A list similar to that provided as input, but with a subset of markers
#'   that passed the filters.
#'
#' @examples
#' NULL
#'
cluster_selection <- function(snp_data, center_pos, ld_threshold = 0.5,
                              max_pair_distance = 10^10, max_marker_distance = 10^10) {

  # The function will need to access marker metadata and LD data
  if(is.null(snp_data$Markers) || is.null(snp_data$LD)) {
    stop("Unsuitable snp_data object passed to function cluster_selection.")
  }

  # Checking parameter center_pos
  stopifnot(is.numeric(center_pos) || is.integer(center_pos))
  stopifnot(length(center_pos) == 1)

  # Checking parameter ld_threshold
  if(ld_threshold < 0 || ld_threshold > 1) {
    stop("Linkage disequilibrium threshold for cluster selection must be between 0 and 1.")
  }

  # Logical vectors indicating which markers are upstream and
  #  downstream of the central position
  five_prime  <- snp_data$Markers$pos <= center_pos
  three_prime <- snp_data$Markers$pos > center_pos

  # Generate a data.frame of the LD for each pair
  LD_df <- reshape2::melt(as.matrix(snp_data$LD),
                          varnames = c("SNP_1", "SNP_2"),
                          value.name = "R2")

  # Adding the position (in base pairs) of SNP_1 and SNP_2
  LD_df[["pos1"]] <- snp_data$Markers$pos[match(LD_df$SNP_1, snp_data$Markers$rs)]
  LD_df[["pos2"]] <- snp_data$Markers$pos[match(LD_df$SNP_2, snp_data$Markers$rs)]
  # Adding the distance between the two markers
  LD_df[["pair_distance"]] <- abs(LD_df$pos2 - LD_df$pos1)
  # Determining which markers (logical) are close enough to central gene position
  LD_df[["gene_distance"]] <- abs(LD_df$pos1 - center_pos) < max_marker_distance &
                              abs(LD_df$pos2 - center_pos) < max_marker_distance

  # Three logical matrices are to be generated to select the final markers :
  #  1- A matrix indicating which marker paris are in high enough linkage
  #     disequilibrium (based on parameter ld_threshold).
  #  2- A matrix indicating which markers are close enough to the gene center
  #     (based on parameter max_marker_distance).
  #  3- A matrix indicating which marker pairs are close enough one the other
  #     (based on parameter max_pair_distance)
  #  All these matrices will have as many rows as there are markers upstream
  #   of the gene (five_prime) and as many columns as there are markers
  #   downstream of the gene (three_prime).

  # Generating the matrix for the first criterion (ld_threshold)
  ld_filter <- Matrix::forceSymmetric(snp_data$LD)[five_prime, three_prime, drop = FALSE] >= ld_threshold
  # Generating the matrix for the second criterion (max_marker_distance)
  gene_distance_filter <- reshape2::acast(LD_df, SNP_1 ~ SNP_2, value.var = "gene_distance")
  gene_distance_filter <- gene_distance_filter[five_prime, three_prime, drop = FALSE]
  # Generating the matrix for the third criterion (max_pair_distance)
  pair_distance_filter <- reshape2::acast(LD_df, SNP_1 ~ SNP_2, value.var = "pair_distance")
  pair_distance_filter <- pair_distance_filter[five_prime, three_prime, drop = FALSE] <= max_pair_distance
  # All are combined into a single logical matrix indicating marker pairs for which all
  #  conditions are satisfied
  filter_matrix <- ld_filter & gene_distance_filter & pair_distance_filter

  # Determine which upstream and downstream markers satisfy the conditions
  #  with at least one marker of the opposite group
  upstream_markers  <- apply(filter_matrix, 1, any, na.rm = TRUE)
  downstream_markers <- apply(filter_matrix, 2, any, na.rm = TRUE)
  kept_markers <- c(upstream_markers, downstream_markers)

  # Return a subset of markers
  return(marker_subset(snp_data, kept_markers))
}
