#' Prepare Data for Selection Analysis
#'
#' Standardize trait columns (z-score) if requested, and (optionally) add
#' relative fitness. Flexible NA handling with warnings/drop/none.
#'
#'
#' @param data A \code{data.frame}.
#' @param fitness_col Character. Name of the fitness column (0/1 or continuous).
#' @param trait_cols Character vector. Names of trait columns to standardize.
#' @param standardize Logical; whether to z-score traits. Default \code{TRUE}.
#' @param add_relative Logical; whether to add a relative-fitness column.
#' @param na_action One of \code{"warn"}, \code{"drop"}, \code{"none"}.
#'                  \itemize{
#'                    \item \code{"warn"}: keep rows, emit warnings (default)
#'                    \item \code{"drop"}: drop rows with NA in fitness or traits
#'                    \item \code{"none"}: keep rows silently
#'                  }
#' @param name_relative Character. Name of the relative-fitness column to add.
#'
#' @return A \code{data.frame} with standardized traits (if requested) and,
#'         optionally, a relative-fitness column.
#' @export
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(fitness = rbinom(10, 1, 0.4),
#'                  size = rnorm(10, 5, 1),
#'                  speed = rnorm(10, 10, 2))
#' out <- prepare_selection_data(df, "fitness", c("size", "speed"),
#'                               standardize = TRUE,
#'                               add_relative = TRUE,
#'                               na_action = "warn")
#' head(out)

prepare_selection_data <- function(data,
                                   fitness_col,
                                   trait_cols,
                                   standardize = TRUE,
                                   add_relative = TRUE,
                                   na_action = c("warn", "drop", "none"),
                                   name_relative = "relative_fitness") {
  na_action <- match.arg(na_action)

  # Copy and basic checks
  df <- data.frame(data, check.names = FALSE)
  if (!is.character(fitness_col) || length(fitness_col) != 1L)
    stop("`fitness_col` must be a single column name (character).")
  if (!all(c(fitness_col, trait_cols) %in% names(df)))
    stop("`fitness_col` or some `trait_cols` not found in data.")

  # Coerce fitness to numeric if it's logical/integer
  if (is.logical(df[[fitness_col]]) || is.integer(df[[fitness_col]])) {
    df[[fitness_col]] <- as.numeric(df[[fitness_col]])
  }

  # NA handling
  cols_check <- c(fitness_col, trait_cols)
  na_rows <- !stats::complete.cases(df[, cols_check, drop = FALSE])

  if (any(na_rows)) {
    n_bad <- sum(na_rows)
    msg <- sprintf("Detected %d row(s) with NA in fitness/traits: %s",
                   n_bad, paste(cols_check, collapse = ", "))
    if (na_action == "warn") {
      warning(msg)
    } else if (na_action == "drop") {
      df <- df[!na_rows, , drop = FALSE]
    } # if "none": do nothing
  }

  # Standardize traits
  if (standardize) {
    for (t in trait_cols) {
      # Keep vectorized scale, strip attributes
      df[[t]] <- as.numeric(scale(df[[t]]))
    }
  }

  # Add relative fitness
  if (add_relative) {
    mean_fit <- mean(df[[fitness_col]], na.rm = TRUE)
    if (!is.finite(mean_fit) || mean_fit == 0) {
      warning("Mean fitness is zero or non-finite; cannot compute relative fitness. Skipping.")
    } else {
      df[[name_relative]] <- df[[fitness_col]] / mean_fit
    }
  }

  df
}
