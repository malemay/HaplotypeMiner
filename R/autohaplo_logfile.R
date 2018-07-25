#' Generate the logfile of a haplotype analysis
#'
#' This function inspects the output of the \code{\link{haplo_selection}} function
#' in order to generate a logfile summarizing details of the analysis. The logfile
#' can be output to a file (defaults to "Log.txt") or to optionally to standard
#' output.
#'
#' The logfile is divided into 4 parts :
#' \itemize{
#' \item{Analysis parameters} {A description of the input files and parameters that
#'   were used for the analysis.}
#' \item{Marker clustering and filtering results} {A description of the number of
#'   markers that were kept following different filtering steps.}
#' \item{Selecting marker pairs in LD across gene} {Data regarding the final set
#'   of markers kept for haplotype generation.}
#' \item{Haplotype assignment} {A summary of the number of different haplotypes
#'   generated and the number of individuals unambiguously assigned to each
#'   haplotype, plus the number of individuals that could not be unambiguously
#'   assigned to a haplotype due to missing data.}
#' }
#'
#' @param results A list containing a variable number of items among the following :
#'   \code{Parameters}, \code{Kinship}, \code{Structure}, \code{All_markers},
#'   \code{Filtered_markers}, \code{Clustered_markers}, \code{Selected_clusters},
#'   \code{Selected_markers}, \code{Haplotypes}. The precise output of the function
#'   will depend on the number of these elements given as input. The function
#'   will throw an error if any element other than those aforementioned is present
#'   in the input list.
#' @param filename A character. The name of the file to which the logfile should
#'   be written. Ignored if \code{to_file} is FALSE.
#' @param to_file A logical. Should the output be written to file? Defaults to
#'   TRUE, otherwise output is directed to standard output.
#'
#' @return NULL, invisibly. This function is called for its side effect of writing
#'   to file or to standard output.
#' @export
#' @examples
#' NULL
#'
autohaplo_logfile <- function(results, filename = "Log.txt", to_file = TRUE) {

  # INPUT CHECKING
  # Checking that all input list elements are standard
  if(!all(names(results) %in% c("Parameters", "Kinship", "Structure", "All_markers",
                                "Filtered_markers", "Clustered_markers", "Selected_clusters",
                                "Selected_markers", "Haplotypes"))) {
    stop("The results object is not suitable as input to autohaplo_logfile.")
  }

  # Redirect output to the logfile
  if(to_file) {
    if(file.exists(filename)) {
      stop("File ", filename, " already exists. Overwriting not allowed.")
    } else {
      sink(filename)
      on.exit(sink(NULL), add = TRUE)
    }
  }

  # Writing analysis parameters to file =========================================
  cat("Analysis parameters:\n")

  # Writing of the log file stops here if no parameters are found
  if(is.null(results$Parameters)) {
    cat("\tNo parameters found in results object. Analysis output stops here.")
    return(invisible(NULL))
  }

  for(parameter in names(results$Parameters)) {
    if(!(parameter %in% c("Gene_database", "Chromosome_database", "Chromosome_length"))) {
      # Removing the dots from the parameter names
      display_name <- gsub("\\.", " ", parameter, fixed = TRUE)
      # Assigning the value corresponding to this parameter
      value <- results$Parameters[[parameter]]
      # If NULL, will print NA (for example if structure or kinship are NULL)
      if(is.null(value)) value <- NA
      if(parameter == "Input_file") cat("\t", "INPUT FILES", "\n")
      if(parameter == "Gene_name") cat("\t", "GENE CHARACTERISTICS", "\n")
      if(parameter == "Minimum_alternate_allele_frequency") cat("\t", "FILTERING PARAMETERS", "\n")
      if(parameter == "Marker_cluster_threshold") cat("\t", "PAIR SELECTION PARAMETERS", "\n")
      # Writing to file
      cat("\t\t", display_name, ":", value, "\n")
      }
    }


  # Writing filtering results to file ===========================================
  cat("Marker clustering and filtering results:\n")

  # Writing of the log file stops here if no markers were kept
  if(is.null(results$Filtered_markers)) {
    cat("No markers remained following filtering. Analysis output stops here.")
    return(invisible(NULL))
  }

  # Shortening the access to some objects
  filters <- results$Filtered_markers$Filters

  cat("\tTotal number of markers:", filters$Total_number, "\n")
  cat("\tNumber of markers located on chromosome", as.character(results$Parameters[["Gene_chromosome"]]), ":", filters$Chromosome_filter, "\n")
  cat("\tNumber of markers less than", results$Parameters[["Maximum_marker_to_gene_distance"]], "bp from gene center :", filters$Distance_filter, "\n")
  cat("\tNumber of biallelic markers :", filters$Biallelic_filter, "\n")
  cat("\tNumber of markers passing missing data filter :", filters$Missing_filter, "\n")
  cat("\tNumber of markers passing heterozygosity filter :", filters$Heterozygous_filter, "\n")
  cat("\tNumber of markers passing MAF filter :", filters$Alternate_allele_filter, "\n")
  cat("\tNumber of markers passing MAC filter (final number kept for clustering):", filters$Allele_count_filter, "\n")

  # Writing the results of the clustering process to file =====================
  cat("\tTotal markers kept following clustering:", nrow(results$Clustered_markers$Markers), "\n")

  # Writing of the log file stops here if no markers were kept following clustering
  #  (This should not happen if filtered markers remained)
  if(is.null(results$Clustered_markers)) {
    cat("No markers remained following clustering. Analysis output stops here.")
    return(invisible(NULL))
  }

  # Writing the results pertaining to markers across the gene to file ===========
  if(is.null(results$Selected_clusters)) {
    cat("No markers in linkage disequilibrium across gene. Analysis output stops here.")
    return(invisible(NULL))
  }

  cpos <- results$Selected_clusters$Markers$pos
  n_across <- nrow(results$Selected_clusters$Markers)
  haplotype_size <- max(cpos, na.rm = TRUE) - min(cpos, na.rm = TRUE)

  cat("Selecting marker pairs in LD across gene:\n")
  cat("\tTotal markers kept:", n_across, "\n")
  cat("\tHaplotype size (distance between two farthest markers):", haplotype_size, "\n")

  # Haplotype assignment
  cat("Haplotype assignment :\n")
  cat("\tNumber of distinct haplotypes :", nrow(results$Haplotypes$Genotypes), "\n")

  for(i in sub("haplotype_", "", rownames(results$Haplotypes$Genotypes@.Data))) {
    cat("\tNumber of individuals assigned to haplotype", i, ":", sum(results$Haplotypes$Assignment == i, na.rm = TRUE), "\n")
  }

  cat("\tNumber of individuals not unambiguously assigned a haplotype:", sum(is.na(results$Haplotypes$Assignment)), "\n")

  # The function returns NULL, invisibly.
  return(invisible(NULL))
}
