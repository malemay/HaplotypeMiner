#' Generate graphical output of a haplotype analysis
#'
#' This is a high-level function for taking an object obtained as output from
#' the function \code{\link{haplo_selection}} and generating graphical output
#' written to files. The \code{graphs} argument allows a high degree of
#' coustomization as to the graphs that will be generated. Some graphs input/output
#' combinations are not allowed (see \code{Details}).
#'
#' The default \code{graphs} option currently results in no graphs being generated.
#' Sensible defaults still have to be defined.
#'
#' The output is directed by default to the current directory. An error is thrown
#' if the suggested output directory already exists, else a new directory with this
#' name is created. This is meant to avoid name conflicts and file overwriting,
#' but this could be changed in the future.
#'
#' This function calls high-level plotting functions defined in the package itself :
#' \code{\link{density_plot}}, \code{\link{ld_plot}}, \code{\link{distance_plot}},
#' \code{\link{genotype_plot}}. These functions are called in the context of
#' \code{haplo_graphs} with their default values, but finer control could be
#' obtained by using these functions directly.
#'
#' Arguments can be passed to function \code{ggplot::ggsave}, which directs the
#' output (\code{width} and \code{height} default to 7). PDF files are generated
#' by default.
#'
#' @param results A list generated as output from the function  \code{\link{haplo_selection}}
#'   and thus containing the results of a haplotype analysis.
#' @param output_dir A character. The name of the directory to which the files
#'   should be written.
#' @param graphs A named list containing at most 6 elements corresponding to the
#'   names of elements in a list output by function \code{\link{haplo_selection}} :
#'   \code{All_markers}, \code{Filtered_markers}, \code{Clustered_markers},
#'   \code{Selected_clusters}, \code{Selected_markers}, \code{Haplotypes}. Each
#'   of the elements of these lists must contain a character vector cotaining at
#'   most 4 elements corresponding to the 4 possible graphs that can be output :
#'   \code{"density"}, \code{"matrix"}, \code{"distance"}, \code{genotypes}. The
#'   list provided will determine which graphs will be output from which objects.
#' @param ... Arguments passed to function \code{ggplot::ggsave} which writes
#'   the image file to disc.
#'
#' @return \code{NULL}, invisibly. The function is called for its side effect of
#'   writing graphs to files.
#' @export
#'
#' @examples
#' NULL
#'
haplo_graphs <- function(results, output_dir = ".", graphs = "default", ...) {

  # Creating a directory with a name that already exists is not supported.
  #  Writing the results in the current directory is allowed.
  if(dir.exists(output_dir) && output_dir != ".") {
    stop("Directory ", output_dir, " already exists. Overwriting not allowed.")
  } else if (output_dir != ".") {
    dir.create(output_dir, recursive = TRUE)
  }

  # Register on.exit handler to restore previous directory
  oldwd <- setwd(output_dir)
  on.exit(setwd(oldwd))

  # No graphs are written if graphs = NULL
  if(is.null(graphs)) {
    # .null_graphs is an object defined in the package
    graphs <- .null_graphs
    # .default_graphs is an object defined in the package
    # at the moment, .default_graphs == .null_graphs
  } else if (identical(graphs, "default")) {
    graphs <- .default_graphs
  }

  # It is an error to ask for input data that is not among the possible options
  #  .input_options is an object defined in the package
  if(!all(names(graphs) %in% .input_options)) {
    which_error <- names(graphs)[!(names(graphs) %in% .input_options)]
    stop(which_error[1], " is not a valid graph input option.")
  }

  # The foor loop is applied over all possible inputs.
  for(i in names(graphs)) {
    # Assigning data and output options for this iteration to local objects
    i_data    <- results[[i]]
    i_options <- graphs[[i]]

    # No graphs are generated if no options are specified
    if(!length(i_options)) next

    # It is an error to ask for an output that is not among the possible options
    #  output_options is an object defined in the package
    if(!all(i_options %in% .output_options)) {
      which_error <- i_options[!(i_options %in% .output_options)]
      stop(which_error[1], " is not a valid graph output option.")
    }

    # Generating the density plot if it is specified as an option
    if("density" %in% i_options) {
      p <- density_plot(snp_data   = i_data,
                        center_pos = results$Parameters$Gene_center,
                        chr_length = results$Parameters$Chromosome_length)

      # Naming and writing the file
      filename <- paste0("Density plot of ", tolower(sub("_", " ", i)), ".pdf")
      ggplot2::ggsave(filename = filename, plot = p, width = 7, height = 7, ...)
    }

    # Generating the LD matrix if it is specified as an option
    if("matrix" %in% i_options) {

      if(nrow(i_data$Markers) < 100) {
        p <- ld_plot(snp_data   = i_data,
                     center_pos = results$Parameters$Gene_center,
                     kept_markers = results$Haplotypes$Markers)

        # Naming and writing the file
        filename <- paste0("LD plot of ", tolower(sub("_", " ", i)), ".pdf")
        ggplot2::ggsave(filename = filename, plot = p, width = 7, height = 7, ...)
      } else {
        # The plot would be too large in this case.
        #  Warnings like this should be given for other plots.
        #  Maximum dimensions should be defined for other plots.
        warning("LD plot for ", i, " could not be computed. The number of markers is ", nrow(i_data$Markers), " whereas the maximum allowed is 100.")
      }
    }

    # Generating the LD-distance plot if it is specified as an option
    if("distance" %in% i_options) {
      p <- distance_plot(snp_data   = i_data,
                         center_pos = results$Parameters$Gene_center,
                         r2_threshold = results$Parameters$Marker_independence_threshold)

      # Naming and writing the file
      filename <- paste0("LD vs distance plot of ", tolower(sub("_", " ", i)), ".pdf")
      ggplot2::ggsave(filename = filename, plot = p, width = 7, height = 7, ...)
    }

    # Generating the genotype plot if it is specified as an option
    if("genotypes" %in% i_options) {

      # Plotting the haplotypes requires different options from the other
      #  input data types
      if(i == "Haplotypes") {
        p <- genotype_plot(snp_data = i_data,
                           gene_pos = results$Parameters$Gene_center,
                           kept_markers = results$Haplotypes$Markers,
                           assignment = NULL,
                           name_order = TRUE)
      } else {
        p <- genotype_plot(snp_data = i_data,
                           gene_pos = results$Parameters$Gene_center,
                           kept_markers = results$Haplotypes$Markers,
                           assignment = results$Haplotypes$Assignment,
                           name_order = FALSE)
      }

      # Naming and writing the file
      filename <- paste0("Genotype plot of ", tolower(sub("_", " ", i)), ".pdf")
      ggplot2::ggsave(filename = filename, plot = p, width = 7, height = 7, ...)
    }

  }

  return(invisible(NULL))
}
