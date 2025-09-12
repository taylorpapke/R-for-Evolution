
# internal: fit quadratic and interaction model for nonlinear selection
analyze_nonlinear_selection <- function(data, fitness_col, trait_cols, fitness_type) {
  # Quadratic terms
  quad <- paste0("I(", trait_cols, "^2)")

  # Interaction terms (pairwise products)
  inter <- if (length(trait_cols) >= 2) {
    pairs <- utils::combn(trait_cols, 2, simplify = FALSE)
    vapply(pairs, function(p) paste0(p[1], ":", p[2]), "")
  } else character(0)

  # Construct full right-hand side
  rhs <- paste(c(trait_cols, quad, inter), collapse = " + ")

  if (fitness_type == "continuous") {
    # Prefer relative_fitness if available
    resp <- if ("relative_fitness" %in% names(data)) {
      "relative_fitness"
    } else {
      fitness_col
    }
    fml <- stats::as.formula(paste(resp, "~", rhs))
    fit <- stats::lm(fml, data = data, na.action = stats::na.omit)
    sm  <- summary(fit)

    # Type III ANOVA if car is installed
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

    # Wald chi-square tests if car is installed
    if (requireNamespace("car", quietly = TRUE)) {
      an <- car::Anova(fit, type = "III", test.statistic = "Wald")
    } else {
      an <- NULL
    }
  }

  list(model = fit, summary = sm, anova = an)
}
