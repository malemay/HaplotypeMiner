#' Title
#'
#' Description
#'
#' Details
#'
#' @param input_file To be completed
#' @param kinship_file To be completed
#' @param structure_file To be completed
#' @param gene_db_file To be completed
#' @param chr_db_file To be completed
#' @param gene_name To be completed
#' @param gene_chrom To be completed
#' @param gene_start To be completed
#' @param gene_end To be completed
#' @param chr_length To be completed
#' @param file_format To be completed
#' @param R2_measure To be completed
#' @param min_alt_threshold To be completed
#' @param max_het_threshold To be completed
#' @param max_missing_threshold To be completed
#' @param min_allele_count To be completed
#' @param cluster_threshold To be completed
#' @param cluster_R2 To be completed
#' @param marker_independence_threshold To be completed
#' @param max_flanking_pair_distance To be completed
#' @param max_marker_to_gene_distance To be completed
#'
#' @return To be completed
#' @export
#'
#' @examples
#' NULL
#'
haplo_params <- function(input_file = NULL, kinship_file = NULL,
                         structure_file = NULL, gene_db_file = NULL,
                         chr_db_file = NULL, gene_name = NULL,
                         gene_chrom = NULL, gene_start = NULL,
                         gene_end = NULL, chr_length = NULL,
                         file_format = NULL, R2_measure = "r2",
                         min_alt_threshold = NULL, max_het_threshold = 0,
                         max_missing_threshold = 0, min_allele_count = 1,
                         cluster_threshold = 0.95, cluster_R2 = R2_measure,
                         marker_independence_threshold = 0.5,
                         max_flanking_pair_distance = 3 * 10^6,
                         max_marker_to_gene_distance = 10^12) {

  # Load gene database, and infer gene position if gene_name was given.
  if(!is.null(gene_db_file)) {

    gene_db <- read.table(gene_db_file, header = TRUE, row.names = 1)
    rownames(gene_db) <- toupper(rownames(gene_db))

    if(is.null(gene_chrom) && !is.null(gene_name)) {
      gene_chrom <- gene_db[toupper(gene_name), "chr"]
      gene_start <- gene_db[toupper(gene_name), "start"]
      gene_end   <- gene_db[toupper(gene_name), "end"]
    }

  }

  # Load chromosome database, and infer chromosome length.
  if(!is.null(chr_db_file)) {

    chr_db <- read.table(chr_db_file, header = TRUE, row.names = 1)
    rownames(chr_db) <- toupper(rownames(chr_db))

    if(is.null(chr_length) && !is.null(gene_chrom)) {
      chr_length <- chr_db[toupper(as.character(gene_chrom)), "length"]
    }

  }

  # Compute gene center position.
  if(!is.null(gene_start) && !is.null(gene_end)) gene_pos <- gene_start + ((gene_end - gene_start) / 2)
  # Determine input file format.
  if(is.null(file_format) && !is.null(input_file)) file_format <- ifelse(grepl("\\.vcf$", input_file), "vcf", "hapmap")
  if(is.null(gene_chrom)) stop("No chromosome provided")
  if(is.null(gene_pos)) stop("Gene position could not be inferred.")
  if(is.null(chr_length)) stop("Chromosome length could not be inferred.")

  list(
    # Input files
    Input_file = input_file,
    File_format = file_format,
    Structure_file = structure_file,
    Kinship_file = kinship_file,
    Gene_database_file = gene_db_file,
    Chromosome_database_file = chr_db_file,
    # Databases
    Gene_database = gene_db,
    Chromosome_database = chr_db,
    # Gene and chromosome information.
    Gene_name = gene_name,
    Gene_chromosome = as.character(gene_chrom),
    Gene_start = gene_start,
    Gene_end = gene_end,
    Gene_center = gene_pos,
    Chromosome_length = chr_length,
    # Analysis parameters
    Cluster_R2_measure = cluster_R2,
    R2_measure = R2_measure,
    Minimum_alternate_allele_frequency = min_alt_threshold,
    Maximum_heterozygous_frequency = max_het_threshold,
    Maximum_missing_allele_frequency = max_missing_threshold,
    Minimum_allele_count = min_allele_count,
    Marker_cluster_threshold = cluster_threshold,
    Marker_independence_threshold = marker_independence_threshold,
    Maximum_flanking_pair_distance = max_flanking_pair_distance,
    Maximum_marker_to_gene_distance = max_marker_to_gene_distance
    )
}
