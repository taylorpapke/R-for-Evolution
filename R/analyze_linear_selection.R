
# internal: fit linear model for directional selection
analyze_linear_selection <- function(data, fitness_col, trait_cols, fitness_type) {

  rhs <- paste(trait_cols, collapse = " + ")

  if (fitness_type == "continuous") {
    # Prefer to use relative_fitness if prepared by prepare_selection_data()
    resp <- if ("relative_fitness" %in% names(data)) {
      "relative_fitness"
    } else {
      fitness_col
    }
    fml <- stats::as.formula(paste(resp, "~", rhs))
    fit <- stats::lm(fml, data = data, na.action = stats::na.omit)
    sm  <- summary(fit)

    # Type III ANOVA if car package is available
    if (requireNamespace("car", quietly = TRUE)) {
      an <- car::Anova(fit, type = "III")
    } else {
      an <- NULL
    }

  } else {
    # Binary fitness: logistic regression
    fml <- stats::as.formula(paste(fitness_col, "~", rhs))
    fit <- stats::glm(fml, data = data,
                      family = stats::binomial("logit"),
                      na.action = stats::na.omit)
    sm  <- summary(fit)

    # Wald chi-square test for logistic regression
    if (requireNamespace("car", quietly = TRUE)) {
      an <- car::Anova(fit, type = "III", test.statistic = "Wald")
    } else {
      an <- NULL
    }
  }

  list(model = fit, summary = sm, anova = an)
}


