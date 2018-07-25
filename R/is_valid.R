#' Title
#'
#' Description
#'
#' Details
#'
#' @param snp_data To be completed
#' @param center_pos To be completed
#'
#' @return To be completed
#'
#' @examples
#' NULL
#'
is_valid <- function(snp_data, center_pos) {

  # First checking that the snp_data is not NULL
  if(is.null(snp_data)) {
    warning("Empty object: all markers have been filtered out. Aborting.")
    return(FALSE)
  }

  # Checking that there are markers on both sides of the gene
  cond2 <- any(snp_data$Markers$pos < center_pos, na.rm = TRUE) && any(snp_data$Markers$pos >= center_pos, na.rm = TRUE)

  if(!cond2) {
    warning("Markers are all on the same side of the gene. Aborting.")
    return(FALSE)
  }

  return(TRUE)
}
