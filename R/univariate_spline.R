
# purpose: estimate one trait fitness function using 
# nonparametric spline-based smoother 
# the funciton fits a GAM with a univariate smooth term to describe how fitness varaies as a function of a single trait
# continuous fitness is using Gaussian GAM 
# Binary fitness is using GAM with binomial family 

univariate_spline <- function(data, fitness_col, trait_col,
                              fitness_type = c("binary","continuous"),
                              k = 10) {
  fitness_type <- match.arg(fitness_type)

  if (length(trait_col) != 1L)
    stop("`trait_col` must be a single column name.")
  if (!trait_col %in% names(data))
    stop("`trait_col` not found in `data`.")
  if (!is.numeric(data[[trait_col]]))
    stop("`trait_col` must be numeric (standardize upstream if needed).")

  # Use relative_fitness if present for continuous fitness, 
  # else compute on the fly
  if (fitness_type == "continuous") {
    if ("relative_fitness" %in% names(data)) {
      y <- data[["relative_fitness"]]
    } else {
      y <- data[[fitness_col]] / mean(data[[fitness_col]], na.rm = TRUE)
    }
    fam <- stats::gaussian()
  } else {
    y <- data[[fitness_col]]
    fam <- stats::binomial("logit")
  }

  df <- data
  df[[".y"]] <- y

  fml <- stats::as.formula(paste0(".y ~ s(", trait_col, ", k = ", k, ")"))

  # fit GAM; omit rows with missing values
  fit <- mgcv::gam(fml, data = df, family = fam, method = "REML",
                   na.action = stats::na.omit)

  # prediction grid across observed trait range
  rng <- range(df[[trait_col]], na.rm = TRUE)
  grid <- data.frame(seq(rng[1], rng[2], length.out = 200))
  names(grid) <- trait_col

  # predict on link scale, then transform to response via linkinv
  pr <- stats::predict(fit, newdata = grid, se.fit = TRUE, type = "link")
  linkinv <- fit$family$linkinv
  grid$fit <- linkinv(pr$fit)
  grid$lwr <- linkinv(pr$fit - 1.96 * pr$se.fit)
  grid$upr <- linkinv(pr$fit + 1.96 * pr$se.fit)

  list(model = fit, grid = grid)
}
