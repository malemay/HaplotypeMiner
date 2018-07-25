#' Title
#'
#' Description
#'
#' Details
#'
#' @param kinship_file To be completed
#'
#' @return To be completed
#' @export
#'
#' @examples
#' NULL
#'
read_kinship <- function(kinship_file) {
    kinship <- read.table(kinship_file)
    return(kinship)
}
