#' Export an snp_data into the hapmap format.
#'
#' Description
#'
#' Details
#'
#' @param snp_data The snp_data to be exported.
#' @param filename Filename for the output file.
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
#'
write_hapmap <- function(snp_data, filename) {

  #  Determine what are the reference/alternate alleles for all positions.
  ref_geno      <- gsub("/.*", "", snp_data$Markers$alleles)
  ref_alleles   <- paste0(ref_geno, ref_geno)
  other_geno    <- gsub(".*/", "", snp_data$Markers$alleles)
  other_alleles <- paste0(other_geno, other_geno)
  het_alleles   <- paste0(ref_geno, other_geno)


  # Convert the numerical matrix into a character matrix with the correct alleles.
  recode_matrix <- apply(t(as.matrix(snp_data$Genotypes)), 2, function(x) {
    return(ifelse(x == 1, ref_alleles,
                  ifelse(x == 2, het_alleles,
                         ifelse(x == 3, other_alleles, "NN"))))
  })
  # Add the Metadata.
  full_df <- cbind(snp_data$Markers, recode_matrix)

  # Export everything to disk.
  write.table(full_df, file = filename, row.names = FALSE,
              col.names = TRUE, sep = "\t", quote = FALSE)
}
