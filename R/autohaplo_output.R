#' Output the results of a haplotype analysis to disk
#'
#' This function takes an object returned by function \code{\link{haplo_selection}}
#' as input and extracts information from this object so as the make the
#' information easily interpretable by the user. The function directs the
#' output of essentially three elements : the analysis logfile, the variant
#' files (vcf or hapmap) containing the marker set at different steps of the
#' analysis, and various kinds of graphs
#'
#' More control over the files output by this function can be obtained by
#' using the underlying functions : \code{\link{autohaplo_logfile}},
#' \code{\link{autohaplo_graphs}}, \code{\link{write_hapmap}}
#'
#' @param results A list generated as output from the function  \code{\link{haplo_selection}}
#'   and thus containing the results of a haplotype analysis.
#' @param output_dir A list generated as output from the function  \code{\link{haplo_selection}}
#'   and thus containing the results of a haplotype analysis.
#' @param variant_files A character vector. The name of the list elements of
#'   \code{results} that should used to output the variant data to disk. Can
#'   be any of \code{All_markers}, \code{Filtered_markers}, \code{Clustered_markers},
#'   \code{Selected_clusters}, \code{Selected_markers}, \code{Haplotypes}.
#'   Defaults to the last three of these elements. If null, no variant file
#'   will be written to disk.
#' @param variant_format A character. The format to be used for writing
#'   variant files. Must be either "vcf" or "hapmap". Other formats are not
#'   supported and will result in an error.
#' @param graphs A list describing which graphical elements should be output.
#'   See \code{\link{autohaplo_graphs}} for more details.
#' @param logfile A character. The name of the file to which the logfile
#'   should be written. Defaults to "Log.txt". If \code{NULL}, no logfile
#'   is written.
#'
#' @return \code{NULL}, invisibly. This function is called for its side
#'   effect of writing output.
#' @export
#'
#' @examples
#' NULL
#'
autohaplo_output <- function(results, output_dir = NULL,
                             logfile = "Log.txt",
                             variant_files = "default",
                             variant_format = results$Parameters$File_format,
                             graphs = "default") {

  # Create output directory and move working directory there to output files in it
  if(is.null(output_dir)) {
    output_dir <- results$Parameters$Gene_name
  }

  # An error is thrown if the name of the directory specified already exists.
  #  This behaviour could be changed in the future.
  if(dir.exists(output_dir)) {
    stop("Directory ", output_dir, " already exists. Overwriting not allowed.")
  } else {
    dir.create(output_dir, recursive = TRUE)
  }

  # Register on.exit handler to restore previous directory
  oldwd <- setwd(output_dir)
  on.exit(setwd(oldwd))

  # The logfile is written with name logfile unless logfile = NULL
  if(!is.null(logfile)) {
    autohaplo_logfile(results = results, filename = logfile, to_file = TRUE)
  }

  # If not NULL, variant files will be written to disk.
  #  Default options write variant files for selecte markers/clusters and haplotypes
  if(is.null(variant_files)) {
    variant_files <- character()
  } else if (identical(variant_files, "default")) {
    variant_files <- c("Selected_clusters", "Selected_markers", "Haplotypes")
  }

  # Writing the files to disk. Output format is supplied by the user.
  for(i in variant_files) {
    if(variant_format == "vcf") {
      write_vcf(results[[i]], paste0(i, ".vcf"))
    } else if (variant_format == "hapmap") {
      write_hapmap(results[[i]], paste0(i, ".hmp.txt"))
    } else {
      # An error is thrown if the format is neither vcf nor hapmap
      stop(variant_format, " is not a valid variant format supported by package autohaplo.")
    }
  }

 # Calling function autohaplo_graphs to write the graphs to files.
  autohaplo_graphs(results, output_dir = ".", graphs = graphs)

  return(invisible(NULL))
}
