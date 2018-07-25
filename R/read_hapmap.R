#' Parse a hapmap file.
#'
#' This function takes a hapmap file as input and converts it to a \code{SnpMatrix}
#' object.
#'
#' Details
#'
#' @param filename A character. The name of the file to be parsed.
#'
#' @return To be completed.
#'
#' @export
#' @examples
#' NULL
#'
read_hapmap <- function(filename) {

  # Read the file
  raw_data <- read.table(filename, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE)
  raw_data$chrom <- as.character(raw_data$chrom)

  # There should be some input checks here
  # stopifnot(condition)

  # Remove the metadata, keeping only genotypes
  genotype_data <- raw_data[ ,-(1:11)]

  # Determine what the reference and alternate allele are for every variant
  ref_alleles <- gsub("/.*", "", raw_data$alleles) #keeps the allele prior to the /
  other_alleles <- gsub(".*/", "", raw_data$alleles) # keeps allele after the /

  # Preparing vectors necessary for recoding the genotype matrix
  ref <- paste0(ref_alleles, ref_alleles) #homozygous for the reference allele
  het_1 <- paste0(ref_alleles, other_alleles) # heterozygoys A/B
  het_2 <- paste0(other_alleles, ref_alleles) # heterozygous B/A
  alt <- paste0(other_alleles, other_alleles) #homozygous for the alternate allele

  # Recode the matrix with integer values:
  #   1 = Homozygous for the reference allele
  #   2 = Heterozygous
  #   3 = Homozygous for the alternate allele.

  # Initializing the recoded matrix
  int_genotypes <- matrix(0, ncol = ncol(genotype_data), nrow = nrow(genotype_data))
  # Setting value to 1 for reference allele homozygosity
  int_genotypes[apply(genotype_data, 2, "==", ref)] <- 1
  # Setting value to 2 for heterozygosity (two possibilities)
  int_genotypes[apply(genotype_data, 2, "==", het_1) | apply(genotype_data, 2, "==", het_2)] <- 2
  # Setting value to 3 for alternate allele homozygosity
  int_genotypes[apply(genotype_data, 2, "==", alt)] <- 3

  # Row names of the recoded matrix are marker names
  rownames(int_genotypes) <- raw_data$rs
  # Colmumn names of the recoded matrix are the names of the cultivars
  colnames(int_genotypes) <- colnames(raw_data)[-(1:11)]

  # The genotype matrix must be tranposed to fit into a SnpMatrix object
  int_genotypes <- t(int_genotypes)

  # Return the list with both the SnpMatrix and the metadata
  # The metadata contains information about the markers
  return(list(Genotypes = new("SnpMatrix", .Data = int_genotypes),
              Markers = raw_data[ , 1:11]))
}
