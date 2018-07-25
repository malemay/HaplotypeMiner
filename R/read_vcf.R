#' Parse a vcf file.
#'
#' Description
#'
#' Details
#'
#' @param filename The file to be parsed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
#'
read_vcf <- function(filename) {

  collapsed_vcf <- VariantAnnotation::readVcf(filename, genome = "Gmax")
  snp_matrix <- VariantAnnotation::genotypeToSnpMatrix(collapsed_vcf)

  row_meta <- SummarizedExperiment::rowRanges(collapsed_vcf)

  metadata <- data.frame(rs = names(row_meta),
                         alleles = paste(as.character(row_meta$REF), as.character(unlist(row_meta$ALT)), sep = "/"),
                         chrom = GenomeInfoDb::seqnames(row_meta),
                         pos = BiocGenerics::start(row_meta),
                         strand = BiocGenerics::strand(row_meta),
                         assembly = NA,
                         center = NA,
                         protLSID = NA,
                         assayLSID = NA,
                         panelLSID = NA,
                         QCcode = NA)

  metadata$rs <- as.character(metadata$rs)
  metadata$alleles <- as.character(metadata$alleles)
  metadata$chrom <- as.character(metadata$chrom)
  metadata$strand <- as.character(metadata$strand)

  return(list(Genotypes = snp_matrix$genotypes,
              Markers = metadata,
              VCF = collapsed_vcf))
}
