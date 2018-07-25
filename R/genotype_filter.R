#' Filter SNP data
#'
#' This function applies different filters on SNP data so as to generate as
#' set of markers suitable for haplotype analysis. Although the function can
#' be called separately for filtering purposes, it was thought and designed
#' specifically for the needs of package \code{autohaplo} and might not suit
#' the needs of the general. See section \code{Details} for a discussion of
#' the different mandatory and optional filters applied by this function.
#'
#' \code{genotype_filter} applies automaticaly three filters to the
#' \code{snp_data} object provided as an argument. These filters are applied
#' both when the function is used separately and when used internally inside
#' a call to \code{\link{haplo_selection}}. The three filters are :
#'
#' \itemize{
#' \item{chrom} {Only markers lying on the chromosome of interest are kept
#'   for further analysis.}
#' \item{distance} {Only markers less than \code{max_distance_to_gene} base pairs
#'   away from the central position of interest (usually the center of the gene
#'   for which haplotypes are beign generated) are kept for further analysis. The
#'   value of this distance is 1 Gb by default, which should mean that all
#'   markers on a given chromosome are kept by default (a more sensible
#'   default value could be used).}
#' \item{multiallelism} {Markers that are not biallelic (i.e. either triallelic
#'   or tetraallelic) are automatically removed from the dataset as package
#'   \code{autohaplo} does not know how to handle these markers yet.}
#' }
#'
#' Other filters are optional and are not applied by default, although is is
#' recommended that users do apply these filters either prior to the analysis,
#' externally to package \code{autohaplo}, or as part of the analysis pipeline
#' implemented by function \code{\link{haplo_selection}}. These four filters
#' are :
#'
#' \itemize{
#'   \item{Missing data} {Markers harbouring a missing data rate higher than
#'     \code{max_missing_threshold} can be selectively removed during the
#'     analysis.}
#'   \item{Heterozygosity} {Markers harbouring a heterozygosity rate higher
#'     than \code{max_het_threshold} can be selectively removed during the
#'     analysis. This may not be relevant for species found in the while, but
#'     is relevant e.g. for crop species which are expected to by homozygous
#'     at all loci.}
#'   \item{Minor allele frequency} {Markers harbouring a minor
#'     allele frequency (MAF) lower than \code{min_alt_treshold} can be selectively
#'     removed from the analysis.}
#'   \item{Minor allele count} {Markers harbouring a minor allele count (MAC)
#'     lower than \code{min_allele_count} can be selectively removed from the
#'     analysis.}
#' }
#'
#' @param snp_data A list of data pertaining so SNP markers and having at least
#'   elements \code{Markers} and \code{Genotypes}. A more detailed account of
#'   the contents of this object can be found in the documentation for functions
#'   \code{\link{read_hapmap}} and \code{\link{read_vcf}}.
#' @param chrom A character of length one. The name of the chromosome for which
#'   markers should be kept.
#' @param center_pos A numeric of length one. The central position (in base pairs)
#'   of the gene of interest.
#' @param max_distance_to_gene A numeric of length one. The maximum distance
#'   (in base pairs) between \code{center_pos} and the position of a marker
#'   for this marker to be kept in the analysis.
#' @param max_missing_threshold A numeric of length one between 0 and 1, or
#'   \code{NULL}. If \code{NULL} (default), no such filter is applied. Otherwise,
#'   all markers with a missing data rate over this value are removed.
#' @param max_het_threshold A numeric of length one between 0 and 1, or
#'   \code{NULL}. If \code{NULL} (default), no such filter is applied. Otherwise,
#'   all markers with a heterozygosity rate over this value are
#'   removed.
#' @param min_alt_threshold A numeric of length one between 0 and 1, or
#'   \code{NULL}. If \code{NULL} (default), no such filter is applied. Otherwise,
#'   all markers with a minor allele frequency below this value are removed.
#' @param min_allele_count A positive numeric value, or \code{NULL}. If
#'   \code{NULL} (default), no such filter is applied. Otherwise, all markers
#'   with a minor allele count below this threshold are removed.
#' @param verbose Logical. Should information regarding the filtering process
#'   be printed to screen? Defaults to TRUE.
#'
#' @return A list containing 3 or 4 elements depending on the \code{snp_data}
#'   object used as input :
#'
#'   \itemize{
#'     \item{Genotypes}{An object of class \code{snpMatrix} containing the
#'       genotypes corresponding to the various markers for every individual.
#'       This is essentially a subset of \code{snp_data$Genotypes} that
#'       contains only markers that have been selected.}
#'     \item{Markers}{A \code{data.frame} containing metadata relative to the
#'       genotyped markers. This is essentially a subset of
#'       \code{snp_data$Markers} that contains only markers that have been
#'       selected.}
#'     \item{Filters}{A list of eight integer vectors indicating how many
#'       markers remained following different filtering steps : (1) the total
#'       number of markers, (2) the number of markers located on the chromosome
#'       of interest, (3) the number of markers located close enough to the
#'       central gene position, (4) the number of biallelic markers, (5) the
#'       number of markers passing the missing data filter, (6) the number of
#'       markers passing the heterozygosity filter, (7) the number of markers
#'       passing the MAF filter, and (8) the number of markers passing the MAC
#'       filter. All these numbers are the number of markers remaining after
#'       every preceding step and not the absolute number of markers passing
#'       this filter.}
#'     \item{VCF}{If a VCF element was present in the initial \code{snp_data}
#'       object, this element is a subset of it containing only the markers
#'       remaining following filtering.}
#'   }
#' @export
#'
#' @examples
#' NULL
#'
genotype_filter  <- function(snp_data, chrom, center_pos,
                             max_distance_to_gene = 10^9,
                             max_missing_threshold = NULL,
                             max_het_threshold = NULL,
                             min_alt_threshold = NULL,
                             min_allele_count = NULL,
                             verbose = TRUE) {

  # Input checking (snp_data)
  if(is.null(snp_data$Markers) || is.null(snp_data$Genotypes)) {
    stop("Improper snp_data argument to function genotype_filter.")
  }

  # Input checking (chrom)
  if(!is.character(chrom) || length(chrom) != 1) {
    stop("Argument chrom to function genotype_filter must be a character of length 1.")
  }

  # Input checking (center_pos)
  if(!is.numeric(center_pos) || center_pos < 0 || length(center_pos) != 1) {
    stop("Argument center_pos to function genotype_filter must be a positive integer of length 1.")
  }

  # Input checking (max_distance_to_gene)
  if(!is.numeric(max_distance_to_gene) || length(max_distance_to_gene) != 1) {
    stop("Argument max_distance_to_gene to function genotype_filter must be a numeric of length 1.")
  }

  # Input checking (max_missing_threshold)
  if(!is.null(max_missing_threshold)) {
    if(max_missing_threshold < 0 || max_missing_threshold > 1) {
      stop("Argument max_missing_threshold to function genotype_filter must take values between 0 and 1.")
    }
  }

  # Input checking (max_het_threshold)
  if(!is.null(max_het_threshold)) {
    if(max_het_threshold < 0 || max_het_threshold > 1) {
      stop("Argument max_het_threshold to function genotype_filter must take values between 0 and 1.")
    }
  }

  # Input checking (min_alt_threshold)
  if(!is.null(min_alt_threshold)) {
    if(min_alt_threshold < 0 || min_alt_threshold > 1) {
      stop("Argument min_alt_threshold to function genotype_filter must take values between 0 and 1.")
    }
  }

  # Input checking (min_allele_count)
  if(!is.null(min_allele_count)) {
    if(min_allele_count < 0) {
      stop("Argument min_allele_count to function genotype_filter must take values between 0 and 1.")
    }
  }

  # Informing the user (there should be a "verbose" option)
  total_number <- nrow(snp_data$Markers)
  if(verbose) cat("Total number of markers :", nrow(snp_data$Markers), "\n")

  # We make sure index numbers are stored as rownames (we will need them later on)
  rownames(snp_data$Markers) <- as.character(1:nrow(snp_data$Markers))

  # We filter first according to chromosome (mandatory step)
  kept_markers <- snp_data$Markers[snp_data$Markers$chrom == chrom, ]
  chrom_filter <- nrow(kept_markers)
  if(verbose) cat("Markers located on chromosome", as.character(chrom), ":", chrom_filter, "\n")

  # Remove markers too far from the gene position (mandatory step)
  kept_markers <- kept_markers[abs(kept_markers$pos - center_pos) <= max_distance_to_gene, ]
  distance_filter <- nrow(kept_markers)
  if(verbose) cat("Markers less than", max_distance_to_gene, "bp from gene center :", distance_filter, "\n")

  # Remove markers that have more than two different alleles (mandatory step)
  kept_markers <- kept_markers[stringr::str_count(kept_markers$alleles, "/") == 1, ]
  biallelic_filter <- nrow(kept_markers)
  if(verbose) cat("Number of biallelic markers :", biallelic_filter, "\n")

  # Get a summary of data about marker genotypes; used for the next three filtering steps
  genotype_subset <- snp_data$Genotypes[ , as.numeric(rownames(kept_markers))]
  snp_summary <- snpStats::col.summary(genotype_subset)

  # Remove markers with missing data (optional)
  if(!is.null(max_missing_threshold)) {
    missing_rate <- (1 - snp_summary$Call.rate)
    missing_pass <- missing_rate <= max_missing_threshold
    missing_filter <- sum(missing_pass)
    if(verbose) cat("Markers passing missing data filter :", missing_filter, "\n")
  } else {
    # All markers are otherwise kept
    missing_pass <- rep(TRUE, nrow(kept_markers))
    missing_filter <- sum(missing_pass)
    if(verbose) cat("No missing data filter applied.\n")
  }

  # Remove markers with heterozygous alleles (optional)
  if(!is.null(max_het_threshold)) {
    heterozygosity_pass <- snp_summary$P.AB <= max_het_threshold
    het_filter <- sum(missing_pass & heterozygosity_pass)
    if(verbose) cat("Markers passing heterozygosity filter :", het_filter, "\n")
  } else {
    # All markers are otherwise kept
    heterozygosity_pass <- rep(TRUE, nrow(kept_markers))
    het_filter <- sum(missing_pass & heterozygosity_pass)
    if(verbose) cat("No heterozygosity filter applied")
  }

  # Remove markers with an alternative allele frequency below threshold (optional)
  if(!is.null(min_alt_threshold)) {
    MAF_pass <- !(snp_summary$P.AA <= min_alt_threshold | snp_summary$P.AA >=  (1 - min_alt_threshold))
    MAF_filter <- sum(missing_pass & heterozygosity_pass & MAF_pass)
    if(verbose) cat("Markers passing MAF filter :", MAF_filter, "\n")
  } else {
    # All markers are otherwise kept
    MAF_pass <- rep(TRUE, nrow(kept_markers))
    MAF_filter <- sum(missing_pass & heterozygosity_pass & MAF_pass)
    if(verbose) cat("No MAF filter applied.\n")
  }

  # Subset the markers according to the last three statistics before applying MAC filter
  kept_markers <- kept_markers[missing_pass & heterozygosity_pass & MAF_pass, ]
  genotype_subset <- snp_data$Genotypes[ , as.numeric(rownames(kept_markers))]
  # Remove markers with a minor allele count below the minimum threshold (optional)
  if(!is.null(min_allele_count)) {
    kept_markers <- kept_markers[mac(genotype_subset) >= min_allele_count, ]
    MAC_filter <- nrow(kept_markers)
    if(verbose) cat("Markers passing MAC filter :", MAC_filter, "\n")
  } else {
    MAC_filter <- nrow(kept_markers)
    if(verbose) cat("No MAC filter applied.\n")
  }

  # Generating the indices to be passed to function marker_subset
  kept_indices <- logical(nrow(snp_data$Markers))
  kept_indices[as.numeric(rownames(kept_markers))] <- TRUE

  # Subsetting the markers
  result <- marker_subset(snp_data = snp_data, indices = kept_indices)

  # These number will be needed when writing the logfile
  filters <- list(Total_number = total_number,
                  Chromosome_filter = chrom_filter,
                  Distance_filter  = distance_filter,
                  Biallelic_filter =  biallelic_filter,
                  Missing_filter = missing_filter,
                  Heterozygous_filter = het_filter,
                  Alternate_allele_filter = MAF_filter,
                  Allele_count_filter = MAC_filter)

  # Output depends on whether there is a VCF element
  if(!is.null(result$VCF)) {
    result <- list(Genotypes = result$Genotypes,
                   Markers = result$Markers,
                   Filters = filters,
                   VCF = result$VCF)
  } else {
    result <- list(Genotypes = result$Genotypes,
                   Markers = result$Markers,
                   Filters = filters)
  }

  return(result)
}
