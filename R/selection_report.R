#' Full Selection Report: Coefficients (+Variance) and Plots
#'
#' Computes linear, quadratic (gamma_ii), and correlational (gamma_ij) coefficients
#' with SE, p-values, and Variance; generates standardized univariate GAM splines
#' and thin-plate spline correlational surfaces.
#'
#' @param data data.frame
#' @param fitness_col character, name of fitness column
#' @param trait_cols character vector, names of trait columns (>=1)
#' @param fitness_type "binary" or "continuous"
#' @param standardize logical; z-score traits via prepare_selection_data (default TRUE)
#' @param spline_k integer; basis dimension for univariate GAM splines
#' @param bins integer; number of filled contour bins for correlational surface
#' @param save_dir optional character; if provided, ggsave all plots there (PNG)
#'
#' @return list with:
#'   - table: tidy data.frame (Term, Type, Beta_Coefficient, Standard_Error, P_Value, Variance)
#'   - uni_models: named list of GAM models (per trait)
#'   - uni_plots:  named list of ggplot objects (per trait)
#'   - tps_models: named list of TPS models (per trait pair)
#'   - tps_plots:  named list of ggplot objects (per trait pair)
#'   - prepared_data: standardized data (and relative_fitness if applicable)
#' @export
selection_report <- function(data, fitness_col, trait_cols,
                             fitness_type = c("binary","continuous"),
                             standardize = TRUE,
                             spline_k = 10,
                             bins = 12,
                             save_dir = NULL) {
  fitness_type <- match.arg(fitness_type)

  df <- prepare_selection_data(
    data = data,
    fitness_col = fitness_col,
    trait_cols = trait_cols,
    standardize = standardize,
    add_relative = TRUE,
    na_action = "warn"
  )

  coef_tab <- selection_coefficients(
    data = df,
    fitness_col = fitness_col,
    trait_cols = trait_cols,
    fitness_type = fitness_type,
    standardize = FALSE,
    use_relative_for_fit = (fitness_type == "continuous")
  )
  coef_tab$Variance <- coef_tab$Standard_Error^2

  uni_models <- list()
  uni_plots  <- list()
  for (t in trait_cols) {
    um <- univariate_spline(
      data = df,
      fitness_col = fitness_col,
      trait_col = t,
      fitness_type = fitness_type,
      k = spline_k
    )
    uni_models[[t]] <- um$model
    up <- univariate_surface(um, t, title = paste("Univariate spline:", t))
    uni_plots[[t]] <- up
  }

  tps_models <- list()
  tps_plots  <- list()
  if (length(trait_cols) >= 2) {
    pairs <- utils::combn(trait_cols, 2, simplify = FALSE)
    for (pr in pairs) {
      key <- paste(pr, collapse = "_")
      tp <- correlational_tps(
        data = df,
        fitness_col = fitness_col,
        trait_cols = pr,
        use_relative = (fitness_type == "continuous")
      )
      tps_models[[key]] <- tp$model
      tp_plot <- correlation_surface(tp, pr, bins = bins)
      tps_plots[[key]] <- tp_plot
    }
  }

  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

    for (nm in names(uni_plots)) {
      ggplot2::ggsave(
        filename = file.path(save_dir, paste0("univariate_", nm, ".png")),
        plot = uni_plots[[nm]], width = 6, height = 4, dpi = 300
      )
    }

    for (nm in names(tps_plots)) {
      ggplot2::ggsave(
        filename = file.path(save_dir, paste0("correlational_", nm, ".png")),
        plot = tps_plots[[nm]], width = 6, height = 5, dpi = 300
      )
    }
  }

  list(
    table = coef_tab,
    uni_models = uni_models,
    uni_plots  = uni_plots,
    tps_models = tps_models,
    tps_plots  = tps_plots,
    prepared_data = df
  )
}
