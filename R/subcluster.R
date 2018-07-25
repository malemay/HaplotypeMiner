#' Drop markers in perfect LD with adjacent markers, keeping only a single
#' representative marker.
#'
#' Description
#'
#' Details
#'
#' @param snp_data The snp_data to be filtered.
#' @param block_threshold The minimum R^2 for two adjacent markers to be
#'   eligible for blocking (default: 1).
#' @param reverse If true, blocks start at the most 3' position.
#'
#' @return A subset of the original \code{snp_data} containing only one
#'      representative marker per LD block.
#'
#' @examples
#' NULL
#'
subcluster <- function(snp_data, block_threshold = 1, reverse = FALSE) {

  # If there is only one marker, it is alone in its cluster and it is returned
  # immediately
  if(nrow(snp_data$Markers) == 1) {
    results <- snp_data
    clusters <- 1
    names(clusters) <- snp_data$Markers$rs
    results$Clusters <- clusters
    return(results)
  }

  # Make a copy of the original snp_data.
  work_snp <- snp_data

  # If we're operating in reverse, reverse both the SnpMatrix and the LD matrix.
  if(reverse) {
    work_snp$Genotypes <- snp_data$Genotypes[ , ncol(snp_data$Genotypes):1]
    work_snp$LD <- Matrix::forceSymmetric(snp_data$LD)[nrow(snp_data$LD):1, ncol(snp_data$LD):1]
  }

  # By default, we won't keep any marker but the very first one.
  kept_indices <- rep(TRUE, ncol(work_snp$Genotypes))
  kept_indices[1] <- TRUE

  # Vector for keeping track of which cluster each marker is assigned to.
  clusters <- 1:ncol(work_snp$Genotypes)

  # Loop over all markers
  for(i in 1:(nrow(work_snp$LD) - 1)) {

    if(kept_indices[i]) {

      # work_snp$LD is a sparse matrix, and accessing an item by row and column index
      # is very slow. We speed things up by a factor of 200 by making a copy of the
      # ith row then accessing that vector instead of accessing the sparse matrix
      # in each iteration of the loop.
      rowTemp <- work_snp$LD[i, ]

      for(j in (i + 1):ncol(work_snp$Genotypes)) {

        if(!is.na(rowTemp[j]) && round(rowTemp[j], 5) >= block_threshold) {
          clusters[j] <- i
          kept_indices[j] <- FALSE
        }

      }
    }
  }

  # If we were operating in reverse, reverse the kept indices vector before
  # applying to the original snp_data.
  if(reverse) {
    kept_indices <- rev(kept_indices)
    clusters <- rev(clusters)
  }

  names(clusters) <- snp_data$Markers$rs

  # Subset the markers, keeping only the block heads.
  results <- marker_subset(snp_data, kept_indices)
  results$Clusters <- clusters

  return(results)
}
