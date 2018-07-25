#' Plot of linkage disequilibrium as a function of the distance between markers
#'
#' This is a high-level plotting function used to plot the linkage disequilibrium
#' between pairs of markers as a function of the distance between them. One
#' should expect a pattern of decreasing linkage disequilibrium with increasing
#' distance. It is important to note that only pairs of markers that are
#' across the gene are represented, as these are the only ones relevant to
#' haplotype generation.
#'
#' As the output of this function is a graphical object from package \code{ggplot2},
#' it is possible to customize the appearance of the graph and add layers to
#' the graph using the syntax of this package.
#'
#' @param snp_data A list describing SNP data and containing at least the
#'   elements \code{LD}, \code{\link{read_hapmap}} for more information on
#'   the contents of this object.
#' @param center_pos The central position (in base pairs) of the gene for
#'   which haplotypes are being defined.
#' @param r2_threshold A numeric of length one. The r2 value that is considered
#'   as a treshold for selecting marker pairs in a given analysis. Results
#'   in a red horizontal line being added to the graph. \code{NULL} by default,
#'   in which case no horizontal line is added to the graph.
#'
#' @return A graphical object inheriting from classes \code{gg} and \code{ggplot}
#'   from package \code{ggplot2}, which can be either assigned to an object
#'   or printed to a graphical device.
#'
#' @export
#' @examples
#' NULL
#'
distance_plot <- function(snp_data, center_pos, r2_threshold = NULL) {

  # Transforming the LD matrix into a data.frame
  #  First column : name of the first marker
  #  Second column : name of the second marker
  #  Third column : Rsquared
  ld_df <- reshape2::melt(as.matrix(snp_data$LD),
                      varnames = c("SNP_1", "SNP_2"),
                      value.name = "Rsquared")

  # Associate positions of SNP_1 and SNP_2
  ld_df$pos1 <- snp_data$Markers$pos[match(ld_df$SNP_1, snp_data$Markers$rs)]
  ld_df$pos2 <- snp_data$Markers$pos[match(ld_df$SNP_2, snp_data$Markers$rs)]
  # Only markers that are across the gene center are kept
  # It also makes sure that a given marker pair is not considered twice
  ld_df <- ld_df[ld_df$pos1 < center_pos & ld_df$pos2 >= center_pos, ]

  # Rescaling the Distance variable
  ld_df$Distance <- abs(ld_df$pos2 - ld_df$pos1) / 10^6

  dplot <-
    # Initializing the plot with LD data
    ggplot2::ggplot(ld_df, ggplot2::aes_string(x = "Distance", y = "Rsquared")) +
    # Adding points corresponding to the given aesthetics
    ggplot2::geom_point() +
    # Not clear what this function call does
    ggplot2::expand_limits(x = 0) +
    # Setting y-axis limits (allowing values y == 1)
    ggplot2::ylim(c(0, 1.001)) +
    # Setting the axis labels
    ggplot2::labs(x = "Distance (Mb)", y = r2expr(snp_data[["LD_measure"]])) +
    # Setting the default theme
    ggplot2::theme_bw()

  # Adding a red horizontal line if r2_threshold is provided
  if(!is.null(r2_threshold)) dplot <- dplot + ggplot2::geom_hline(yintercept = r2_threshold, color = "red")

  return(dplot)
}
