#' Convert r^2 measure into a printable expression
#'
#' This function is used internally by package \code{autohaplo} to generate
#' pretty labels for graphs.
#'
#' @param r2_measure A character of length one. The R-squared measure used
#'   to represent linkage disequilibrium. Must be either \code{"r2"},
#'   \code{"r2v"}, \code{"r2s"} or \code{"r2vs"}, otherwise an error will
#'   be thrown.
#'
#' @return An object of class \code{expression} which can be used to generate
#'   pretty labels for graphs.
#'
#' @examples
#' NULL
#'
r2expr <- function(r2_measure) {
  if(r2_measure == "r2") {
    return(expression(r^2))
  } else if (r2_measure == "r2v") {
    return(expression(r[V]^2))
  } else if (r2_measure == "r2s") {
    return(expression(r[S]^2))
  } else if (r2_measure == "r2vs") {
    return(expression(r[VS]^2))
  } else {
    stop(r2_measure, " is not a correct r^2 measure.")
  }
}
