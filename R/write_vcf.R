#' Write VCF output
#'
#' This will generate a VCF file containing the genotypes. This will only work if the
#' initial file was in VCF format.
#'
#' Details
#'
#' @param snp_data The object to be written to the file.
#' @param filename The name of the output file.
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
#'
write_vcf <- function(snp_data, filename) {
  if(!is.null(snp_data$VCF)) {
    VariantAnnotation::writeVcf(snp_data$VCF, filename)
  } else {
    stop("Not implemented")
  }
}
