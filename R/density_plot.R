#' Plot of marker density along chromosome
#'
#' This is a high-level function to plot marker density along the chromosome.
#' This may help in identifying regions that are poorly covered by markers and
#' is thus helpful for troubleshooting.
#'
#' As the output of this function is a graphical object from package \code{ggplot2},
#' it is possible to customize the appearance of the graph and add layers to
#' the graph using the syntax of this package.
#'
#' @param snp_data A list comprising at least an element called \code{Markers},
#'   a \code{data.frame} giving information about a set of markers. See
#'   \code{\link{read_hapmap}} and \code{\link{read_vcf}} for more information
#'   about the structure of \code{snp_data}.
#' @param chr_length A numeric of length one. The length of the chromosome
#'   (in base pairs). This is used to set the plotting window of the x-axis
#'   from 0 to the length of the chromosome. If \code{NULL} (default), the
#'   plotting window will only cover those markers present in the dataset.
#' @param center_pos A numeric of length one. The position of the center of
#'   the gene of interest (in base pairs). This is used to add a red vertical
#'   line at gene position. If \code{NULL}, no vertical line will be added.
#'
#' @return A graphical object inheriting from classes \code{gg} and \code{ggplot}
#'   from package \code{ggplot2}, which can be either assigned to an object
#'   or printed to a graphical device.
#'
#' @export
#' @examples
#' NULL
#'
density_plot <- function(snp_data, chr_length = NULL, center_pos = NULL) {

  # Rescaling values on the x-axis by 10^-6
  snp_data$Markers$pos <- snp_data$Markers$pos / 10^6
  if(!is.null(chr_length)) chr_length <- chr_length / 10^6
  if(!is.null(center_pos)) center_pos <- center_pos / 10^6

  dplot <-
    # Initializing the plot with position on the x-axis
    ggplot2::ggplot(snp_data$Markers, ggplot2::aes_string(x = "pos")) +
    # Adding the density geom
    ggplot2::geom_density(fill = "blue3", alpha = 0.6) +
    # Setting appropriate axis labels
    ggplot2::labs(x = "Position on chromosome (Mb)", y = "Marker density") +
    # Setting the theme
    ggplot2::theme_bw()

  # Setting the axis limits if chr_length is provided
  if(!is.null(chr_length)) dplot <- dplot + ggplot2::xlim(c(0, chr_length))
  # Adding a vertical line at gene central position if provided
  if(!is.null(center_pos)) dplot <- dplot + ggplot2::geom_vline(xintercept = center_pos, color = "red")

  return(dplot)
}
