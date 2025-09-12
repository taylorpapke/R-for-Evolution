#' Detect GLM Family from Fitness Column
#'
#' If fitness takes values {0,1}, returns \code{binomial(link="logit")};
#' otherwise returns \code{gaussian()}.
#'
#' @param fitness_vector Numeric vector of fitness values.
#' @return An object of class \code{family}.
#' @export

detect_family <- function(fitness_vector) {
  u <- unique(stats::na.omit(fitness_vector))
  if (length(u) <= 2 && all(u %in% c(0, 1))) {
    stats::binomial(link = "logit")
  } else {
    stats::gaussian()
  }
}
