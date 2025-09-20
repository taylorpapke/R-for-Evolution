#' Calculate Selection Coefficients (Lande & Arnold, 1983)
#'
#' Computes directional (linear, beta), quadratic (diagonal gamma),
#' and correlational (off-diagonal gamma) selection coefficients.
#' Uses standardized traits and (optionally) relative fitness prepared upstream.
#'
#' @param data A data.frame with fitness and traits.
#' @param fitness_col Character. Name of the fitness column.
#' @param trait_cols Character vector. Names of trait columns.
#' @param fitness_type Either "binary" or "continuous".
#' @param standardize Logical; whether to z-score traits (default TRUE).
#' @param use_relative_for_fit Logical; if TRUE and fitness_type="continuous",
#'   models use relative_fitness (recommended). Ignored for binary GLM p-values.
#'
#' @return A tidy data.frame with columns:
#'   Term, Type (Linear/Quadratic/Correlational), Beta_Coefficient, Standard_Error, P_Value
#' @export

selection_coefficients <- function(data, fitness_col, trait_cols,
                                   fitness_type = c("binary","continuous"),
                                   standardize = TRUE,
                                   use_relative_for_fit = TRUE) {
  fitness_type <- match.arg(fitness_type)
  
  # --- 1) Prepare data  ---
  df <- prepare_selection_data(
    data = data,
    fitness_col = fitness_col,
    trait_cols = trait_cols,
    standardize = standardize,
    add_relative = TRUE,  
    na_action = "warn"
  )
  
  # --- 2) call the functions ---
  # linear
  linear_result <- analyze_linear_selection(
    data = df,
    fitness_col = fitness_col,
    trait_cols = trait_cols,
    fitness_type = fitness_type
  )
  
  # non-linear
  nonlinear_result <- analyze_nonlinear_selection(
    data = df,
    fitness_col = fitness_col,
    trait_cols = trait_cols,
    fitness_type = fitness_type
  )
  
  # --- 3) extract all coefficients ---
  linear_coefs <- extract_linear_coefficients(trait_cols, linear_result)
  quadratic_coefs <- extract_quadratic_coefficients(trait_cols, nonlinear_result)
  interaction_coefs <- extract_interaction_coefficients(trait_cols, nonlinear_result)
  
  # --- 4) combine all results ---
  all_coefs <- rbind(linear_coefs, quadratic_coefs, interaction_coefs)
  
  # add variance 
  all_coefs$Variance <- all_coefs$Standard_Error^2
  
  return(all_coefs)
}