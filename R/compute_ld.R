#' Linkage disequilibrium computation
#'
#' This function computes linkage disequilibrium (LD) for a set of selected loci.
#' The measure can be either r2, r2v, r2s, or r2vs. Function \code{ld} from
#' package \code{snpStats} is used when computing r2, whereas r2 values corrected
#' with kinship data (r2v), population structure (r2s), or both (r2vs) make use
#' of function \code{LD.Measures} from package \code{LDcorSV}.
#'
#' An error will be thrown if the function is asked to compute an r2 measure
#' without all required inputs (kinship and/or structure data). If kinship
#' and/or structure data are given as input but the computed measure makes no
#' use of these data, they will be ignored with a warning.
#'
#' @param snp_data A list of SNP data comprising at least an element called
#'   \code{Genotypes}, a \code{SnpMatrix} object containing the genotype of
#'   individuals included in the analysis.
#' @param measure {A character of length one. The r2 measure that was used
#'   for computing the linkage disequilibrium values in \code{LD}. Can be either
#'   \code{"r2"}, \code{"r2v"}, \code{"r2s"}, \code{"r2vs"}.}
#' @param kinship A square numeric matrix of kinship values with as many rows and
#'   columns as there are individuals included in the analysis. More specifically,
#'   \code{kinship} must have as many rows and columns as there are rows in
#'   \code{snp_data$Genotypes}.
#' @param structure A numeric matrix with as many rows as there are individuals
#'   in the analysis and as many columns as there are putative subpopulations in
#'   the dataset.
#'
#' @return A list with two components :
#'   \itemize{
#'   \item{LD} {A symmetric Matrix of class \code{dsCMatrix} containing linkage
#'   disequilibrium values for all marker pairs.}
#'   \item{LD_measure} {A character of length one. The r2 measure that was used
#'   for computing the linkage disequilibrium values in \code{LD}. Can be either
#'   \code{"r2"}, \code{"r2v"}, \code{"r2s"}, \code{"r2vs"}.}
#'   }
#'
#' @export
#' @examples
#' NULL
#'
compute_ld <- function(snp_data, measure, kinship = NULL, structure = NULL) {

  # Assessing whether kinship and/or structure are included in the analysis
  has_kinship   <- !is.null(kinship)
  if(!has_kinship) kinship <- NA

  has_structure <- !is.null(structure)
  if(!has_structure) structure <- NA

  # Input checking
  if(measure %in% c("r2s", "r2vs") && !has_structure) {
    stop("Structure data must be provided in order to compute ", measure, ".")
  }

  if(measure %in% c("r2v", "r2vs") && !has_kinship) {
    stop("Kinship data must be provided in order to compute ", measure, ".")
  }

  # Warnings given to the user if the r2 measure computed might not be that expected
  if(has_kinship && measure %in% c("r2", "r2s")) {
    warning("Kinship data is provided but is ignored when computing ", measure, ".")
  }

  if(has_structure && measure %in% c("r2", "r2v")) {
    warning("Structure data is provided but is ignored when computing ", measure, ".")
  }

  # Computing simple r2 is done with package snpStats
  if(measure == "r2") {

    LD <- snpStats::ld(snp_data$Genotypes,
                       depth = ncol(snp_data$Genotypes) - 1,
                       stats = "R.squared",
                       symmetric = TRUE)

  } else if (measure %in% c("r2s", "r2v", "r2vs")) {
    # Computing the ld using package LDcorSV with the appropriate method
    ld_df <- LDcorSV::LD.Measures(as(snp_data$Genotypes, "numeric"),
                                  V = kinship, S = structure,
                                  data = "G", na.presence = TRUE)

    # Ensuring that the proper order will be kept
    ld_df$loc1 <- as.character(ld_df$loc1)
    ld_df$loc2 <- as.character(ld_df$loc2)
    ld_df$loc1 <- factor(ld_df$loc1, levels = colnames(snp_data$Genotypes@.Data))
    ld_df$loc2 <- factor(ld_df$loc2, levels = colnames(snp_data$Genotypes@.Data))

    # Turn it into a matrix.
    ld_matrix <- reshape2::acast(ld_df, loc1 ~ loc2,
                                 value.var = measure,
                                 fun.aggregate = mean)

    # NaN values are set to 0
    ld_matrix[is.na(ld_matrix)] <- 0
    # Adding the first column filled with 0 (has the name of the first row)
    ld_matrix <- cbind(rep(0, nrow(ld_matrix)), ld_matrix)
    colnames(ld_matrix)[1] <- rownames(ld_matrix)[1]
    # Adding the last row filled with 0 (has the name of the last column)
    ld_matrix <- rbind(ld_matrix, rep(0, ncol(ld_matrix)))
    rownames(ld_matrix)[nrow(ld_matrix)] <- colnames(ld_matrix)[ncol(ld_matrix)]
    # Coercing the result to a symmetric Matrix
    LD <- suppressMessages(Matrix::forceSymmetric(Matrix::Matrix(ld_matrix)))

  } else {
    # Only r2, r2v, r2s and r2vs are supported
    stop(measure, " is not an appropriate R-squared measure.")
  }

  # Return both the LD matrix and the computed measure
  return(list(LD = LD, LD_measure = measure))
}
