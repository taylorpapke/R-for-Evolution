#' Lande and Arnold (1983) Selection Analysis Framework
#'
#' Implementation of the Lande & Arnold (1983) regression approach
#' to quantify natural selection on quantitative traits.
#'
#' @param data A data.frame containing fitness and traits.
#' @param fitness_col Character. Name of the fitness column.
#' @param trait_cols Character vector. Names of trait columns.
#' @param standardize Logical; whether to z-score traits (default TRUE).
#' @param family Optional GLM family. If NULL, auto-detected from fitness.
#'
#' @return A list with:
#' \item{coefficients_table}{Tidy table of linear/quadratic/correlational coefficients}
#' \item{prepared_data}{Preprocessed data (traits standardized if requested)}
#' \item{fitness_type}{"binary" or "continuous"}
#' \item{analysis_date}{Date of analysis}
#' \item{models}{List with fitted linear/nonlinear models}
#' @export
selection_analysis <- function(data,
                               fitness_col,
                               trait_cols,
                               standardize = TRUE,
                               family = NULL) {
  if (length(trait_cols) < 1) {
    stop("At least 1 trait is required for selection analysis")
  }

  df <- prepare_selection_data(data, fitness_col, trait_cols, standardize)
  fam <- if (is.null(family)) detect_family(df[[fitness_col]]) else family
  fitness_type <- ifelse(fam$family == "binomial", "binary", "continuous")

  lin <- analyze_linear_selection(df, fitness_col, trait_cols, fitness_type)
  nonlin <- analyze_nonlinear_selection(df, fitness_col, trait_cols, fitness_type)

  coeff_tab <- extract_linear_coefficients(trait_cols, lin)
  if (length(trait_cols) >= 2) {
    coeff_tab <- rbind(
      coeff_tab,
      extract_quadratic_coefficients(trait_cols, nonlin),
      extract_interaction_coefficients(trait_cols, nonlin)
    )
  }

  list(
    coefficients_table = coeff_tab,
    prepared_data = df,
    fitness_type = fitness_type,
    analysis_date = Sys.Date(),
    models = list(linear = lin$model, nonlinear = nonlin$model)
  )
}
