#' Title
#'
#' Description
#'
#' Details
#'
#' @param structure_file To be completed
#'
#' @return To be completed
#' @export
#'
#' @examples
#' NULL
#'
read_structure <- function(structure_file) {
  # There should be some input checking here
    structure_data <- read.table(structure_file)
    return(structure_data)
}
