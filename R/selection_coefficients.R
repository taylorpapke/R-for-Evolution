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

  # --- 1) Prepare data (uses your fusion prepare_selection_data) ---
  df <- prepare_selection_data(
    data = data,
    fitness_col = fitness_col,
    trait_cols = trait_cols,
    standardize = standardize,
    add_relative = TRUE,      # we want relative fitness available
    na_action = "warn"
  )

  # choose response for modeling
  resp_lin <- if (fitness_type == "continuous" && use_relative_for_fit &&
                  "relative_fitness" %in% names(df)) "relative_fitness" else fitness_col

  # --- 2) Build formulas ---
  # Linear: resp ~ x1 + x2 + ...
  lin_terms <- paste(trait_cols, collapse = " + ")
  f_lin <- stats::as.formula(paste0(resp_lin, " ~ ", lin_terms))

  # Quadratic+Correlational: resp ~ x1 + ... + I(x1^2) + ... + x1:x2 + ...
  quad_terms <- paste(sprintf("I(%s^2)", trait_cols), collapse = " + ")
  int_pairs <- if (length(trait_cols) >= 2) {
    utils::combn(trait_cols, 2, FUN = function(v) paste0(v[0+1],"*",v[1+1]))
  } else character(0)
  rhs_full <- paste(c(lin_terms, quad_terms, int_pairs), collapse = " + ")
  f_full <- stats::as.formula(paste0(resp_lin, " ~ ", rhs_full))

  # --- 3) Fit models and extract stats ---
  # For p-values with binary fitness, GLM(binomial) is appropriate.
  fit_fun_lin  <- if (fitness_type == "binary") stats::glm else stats::lm
  fit_fun_full <- if (fitness_type == "binary") stats::glm else stats::lm
  fam <- if (fitness_type == "binary") stats::binomial("logit") else NULL

  fit_lin  <- if (is.null(fam)) fit_fun_lin(f_lin,  data = df)
  else fit_fun_lin(f_lin,  data = df, family = fam)
  fit_full <- if (is.null(fam)) fit_fun_full(f_full, data = df)
  else fit_fun_full(f_full, data = df, family = fam)

  summ_lin  <- summary(fit_lin)$coefficients
  summ_full <- summary(fit_full)$coefficients

  # helper to tidy one row
  tidy_row <- function(term, type, est, se, p) {
    data.frame(
      Term = term,
      Type = type,
      Beta_Coefficient = unname(est),
      Standard_Error = unname(se),
      P_Value = unname(p),
      row.names = NULL,
      check.names = FALSE
    )
  }

  # --- 3a) Linear betas ---
  out_lin <- list()
  for (t in trait_cols) {
    if (t %in% rownames(summ_lin)) {
      est <- summ_lin[t, "Estimate"]; se <- summ_lin[t, "Std. Error"]; p <- summ_lin[t, "Pr(>|t|)"]
      # in GLM(summary.glm), p-value col name differs:
      if (is.na(p) && "Pr(>|z|)" %in% colnames(summ_lin)) p <- summ_lin[t, "Pr(>|z|)"]
      out_lin[[t]] <- tidy_row(t, "Linear", est, se, p)
    }
  }
  out_lin <- if (length(out_lin)) do.call(rbind, out_lin) else
    data.frame(Term=character(), Type=character(),
               Beta_Coefficient=numeric(), Standard_Error=numeric(), P_Value=numeric())

  # --- 3b) Quadratic (γ_ii = 2 * b_ii on I(x^2)) ---
  out_quad <- list()
  for (t in trait_cols) {
    sq_term <- paste0("I(", t, "^2)")
    if (sq_term %in% rownames(summ_full)) {
      est_b <- summ_full[sq_term, "Estimate"]
      se_b  <- summ_full[sq_term, "Std. Error"]
      p_b   <- summ_full[sq_term, "Pr(>|t|)"]
      if (is.na(p_b) && "Pr(>|z|)" %in% colnames(summ_full)) p_b <- summ_full[sq_term, "Pr(>|z|)"]

      # Lande & Arnold: gamma_ii = 2 * b_ii
      out_quad[[t]] <- tidy_row(t, "Quadratic", 2*est_b, 2*se_b, p_b)
    }
  }
  out_quad <- if (length(out_quad)) do.call(rbind, out_quad) else
    data.frame(Term=character(), Type=character(),
               Beta_Coefficient=numeric(), Standard_Error=numeric(), P_Value=numeric())

  # --- 3c) Correlational (γ_ij = b_ij on x_i:x_j) ---
  out_int <- list()
  if (length(trait_cols) >= 2) {
    pairs <- utils::combn(trait_cols, 2, simplify = FALSE)
    for (pr in pairs) {
      lab <- paste0(pr[[1]], ":", pr[[2]])
      # model.matrix uses ":" for interaction; our formula used "*", which expands to ":" and main effects.
      if (lab %in% rownames(summ_full)) {
        est <- summ_full[lab, "Estimate"]; se <- summ_full[lab, "Std. Error"]; p <- summ_full[lab, "Pr(>|t|)"]
        if (is.na(p) && "Pr(>|z|)" %in% colnames(summ_full)) p <- summ_full[lab, "Pr(>|z|)"]
        out_int[[lab]] <- tidy_row(lab, "Correlational", est, se, p)
      }
    }
  }
  out_int <- if (length(out_int)) do.call(rbind, out_int) else
    data.frame(Term=character(), Type=character(),
               Beta_Coefficient=numeric(), Standard_Error=numeric(), P_Value=numeric())

  # --- 4) Bind tidy tables ---
  res <- rbind(out_lin, out_quad, out_int)
  rownames(res) <- NULL
  res$Variance <- res$Standard_Error^2
  res
}

