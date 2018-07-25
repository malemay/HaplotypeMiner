#' Title
#'
#' Description
#'
#' Details
#'
#' @param x To be completed
#'
#' @return To be completed
#'
#' @examples
#' NULL
#'
mac <- function(x) {
  x <- matrix(as.character(x), nrow = nrow(x))
  x[x == "00"] <- NA
  min_count <- apply(x, 2, function(x) sum(ifelse(x == "01", 2, ifelse(x == "02", 1, 0)), na.rm = TRUE))
  maj_count <- apply(x, 2, function(x) sum(ifelse(x == "03", 2, ifelse(x == "02", 1, 0)), na.rm = TRUE))
  pmin(maj_count, min_count)
}
