#' Standardized Selection Differential (S)
#'
#' Computes \eqn{S = \mathrm{Cov}(w_{rel}, z)}.
#' If \code{assume_standardized = TRUE}, \eqn{z} is assumed mean 0, sd 1;
#' otherwise \eqn{z} is internally z-scored.
#' If \code{use_relative = TRUE}, fitness is divided by its mean.
#'
#' @param data A data.frame.
#' @param fitness_col Character. Fitness column.
#' @param trait_col Character. Single trait column.
#' @param assume_standardized Logical; if FALSE, z-score inside.
#' @param use_relative Logical; if TRUE, compute relative fitness internally.
#'
#' @return Numeric selection differential S.
#' @export
selection_differential <- function(data, fitness_col, trait_col,
                                   assume_standardized = TRUE,
                                   use_relative = TRUE) {
  # Extract columns
  z <- data[[trait_col]]
  w <- data[[fitness_col]]

  # Drop rows with NA in either
  keep <- stats::complete.cases(z, w)
  z <- z[keep]; w <- w[keep]

  # Standardize trait if needed
  if (!assume_standardized) {
    sd_z <- stats::sd(z)
    if (!is.finite(sd_z) || sd_z == 0) {
      stop("Trait column has zero variance; cannot standardize.")
    }
    z <- as.numeric(scale(z))
  }

  # Relative fitness if requested
  if (use_relative) {
    mu <- mean(w)
    if (is.finite(mu) && mu != 0) {
      w <- w / mu
    } else {
      stop("Cannot compute relative fitness: mean fitness is zero or non-finite.")
    }
  }

  # Cov(w_rel, z) = mean(w_rel * z) because mean(z) = 0 after standardization
  mean(z * w)
}
