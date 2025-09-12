#' Bootstrap Confidence Intervals for Selection Coefficients
#'
#' Percentile bootstrap over individuals for Lande & Arnold coefficients.
#'
#' @param data A data.frame (raw; traits will be standardized inside).
#' @param fitness_col Character. Fitness column.
#' @param trait_cols Character vector. Trait columns.
#' @param fitness_type "binary" or "continuous".
#' @param B Integer; number of bootstrap replicates (default 1000).
#' @param conf Numeric; confidence level (default 0.95).
#' @param seed Optional integer seed for reproducibility (default NULL).
#'
#' @return A list with:
#' \item{draws}{All bootstrap coefficient draws (long format)}
#' \item{ci}{Percentile CI per term (Term, Type, est, lwr, upr)}
#' @export
bootstrap_selection <- function(data, fitness_col, trait_cols,
                                fitness_type = c("binary","continuous"),
                                B = 1000, conf = 0.95, seed = NULL) {
  fitness_type <- match.arg(fitness_type)

  if (!is.null(seed)) set.seed(seed)

  # 1) Clean + standardize once, up front (avoid silent NA drops inside fits)
  cols_needed <- c(fitness_col, trait_cols)
  ok <- stats::complete.cases(data[, cols_needed, drop = FALSE])
  df0 <- data[ok, , drop = FALSE]
  if (nrow(df0) == 0L) stop("No complete cases for fitness and trait_cols.")
  # z-score traits
  for (t in trait_cols) df0[[t]] <- as.numeric(scale(df0[[t]]))

  n <- nrow(df0)
  alpha <- (1 - conf) / 2

  # 2) One pass to define the output schema
  get_coeffs <- function(df) {
    lin   <- analyze_linear_selection(df, fitness_col, trait_cols, fitness_type)
    nonlin<- analyze_nonlinear_selection(df, fitness_col, trait_cols, fitness_type)

    # Use your SIMPLE extractors (summary p-values only)
    lin_tab  <- extract_linear_coefficients(trait_cols, lin)
    quad_tab <- extract_quadratic_coefficients(trait_cols, nonlin)
    int_tab  <- extract_interaction_coefficients(trait_cols, nonlin)

    # Bind (handle empties)
    out <- rbind(lin_tab, quad_tab, int_tab)
    if (is.null(out)) {
      out <- data.frame(Term=character(), Type=character(),
                        Beta_Coefficient=numeric(), Standard_Error=numeric(),
                        P_Value=numeric(), check.names = FALSE)
    }
    out
  }

  # 3) Bootstrap
  draws_list <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    res <- get_coeffs(df0[idx, , drop = FALSE])
    if (nrow(res)) {
      res$replicate <- b
    } else {
      # Keep structure even if empty
      res$replicate <- integer(0)
    }
    draws_list[[b]] <- res
  }
  draws <- do.call(rbind, draws_list)
  rownames(draws) <- NULL

  # 4) Percentile CI
  if (nrow(draws) == 0L) {
    ci <- data.frame(Term=character(), Type=character(),
                     est=numeric(), lwr=numeric(), upr=numeric(),
                     check.names = FALSE)
  } else {
    key <- unique(draws[, c("Term","Type"), drop = FALSE])
    ci_rows <- vector("list", nrow(key))
    for (i in seq_len(nrow(key))) {
      sel <- draws$Term == key$Term[i] & draws$Type == key$Type[i]
      vals <- draws$Beta_Coefficient[sel]
      vals <- vals[is.finite(vals)]
      if (!length(vals)) {
        ci_rows[[i]] <- data.frame(Term = key$Term[i], Type = key$Type[i],
                                   est = NA_real_, lwr = NA_real_, upr = NA_real_,
                                   check.names = FALSE)
      } else {
        est <- stats::median(vals, na.rm = TRUE)
        qs  <- stats::quantile(vals, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
        ci_rows[[i]] <- data.frame(Term = key$Term[i], Type = key$Type[i],
                                   est = est, lwr = qs[1], upr = qs[2],
                                   check.names = FALSE)
      }
    }
    ci <- do.call(rbind, ci_rows)
    rownames(ci) <- NULL
  }

  list(draws = draws, ci = ci)
}
