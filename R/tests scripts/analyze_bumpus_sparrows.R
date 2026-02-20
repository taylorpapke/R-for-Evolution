
# ---- Packages ----
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
if (!requireNamespace("mgcv", quietly = TRUE)) {
  message("Note: 'mgcv' not installed; 'correlational_tps' with method='gam' will need it.")
}
if (!requireNamespace("fields", quietly = TRUE)) {
  message("Note: 'fields' not installed; only needed for method='tps'.")
}

output_dir <- "complete_test_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("1. Loading core function files...\n")
core_function_files <- c(
  "prepare_selection_data.R", "analyze_linear_selection.R", "analyze_nonlinear_selection.R",
  "extract_results.R", "selection_coefficients.R", "detect_family.R", "selection_differential.R", 
  "univariate_spline.R", "univariate_surface.R", "correlational_tps.R", "bootstrap_selection.R"
)

loaded_files <- character(0)
for (file in core_function_files) {
  if (file.exists(file)) {
    source(file)
    loaded_files <- c(loaded_files, file)
    cat("Loaded:", file, "\n")
  } else if (file.exists(file.path("..", file))) {
    source(file.path("..", file))
    loaded_files <- c(loaded_files, file)
    cat("Loaded:", file.path("..", file), "\n")
  } else {
    cat("File not found:", file, "\n")
  }
}


if (!exists("detect_family")) {
  detect_family <- function(fitness) {
    uniq <- unique(stats::na.omit(fitness))
    if (length(uniq) == 2 && all(sort(uniq) == c(0, 1))) {
      list(family = stats::binomial("logit"), type = "binary")
    } else {
      list(family = stats::gaussian(), type = "continuous")
    }
  }
  cat("Fallback detect_family() defined in script\n")
}

# ---- correlation_surface (for TPS/GAM grid) ----
if (!exists("correlation_surface")) {
  correlation_surface <- function(tps, trait_cols, bins = 12) {
    stopifnot(is.list(tps), "grid" %in% names(tps))
    df <- tps$grid
    if (length(trait_cols) != 2L) stop("trait_cols must be length 2.")
    if (!all(trait_cols %in% names(df))) stop("Trait columns not found in tps$grid.")
    if (!".fit" %in% names(df)) stop("tps$grid must contain '.fit'.")
    if (!is.numeric(df[[trait_cols[1]]]) || !is.numeric(df[[trait_cols[2]]])) stop("Traits must be numeric.")
    if (!is.numeric(df$.fit)) df$.fit <- as.numeric(df$.fit)
    
    ggplot(
      df,
      aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[[".fit"]]
      )
    ) +
      geom_contour_filled(bins = bins) +
      geom_contour(color = "black", alpha = 0.3) +
      labs(x = trait_cols[1], y = trait_cols[2], fill = "Probability") +
      theme_minimal(base_size = 12)
  }
  cat("Created correlation_surface function\n")
}


overlay_glm_on_tps <- function(data_df, tps_result, trait_cols,
                               bins = 15, grid_n = 30,
                               levels = c(0.25, 0.5, 0.75),
                               line_color = "white", line_size = 0.8, line_alpha = 0.9) {
  stopifnot(is.list(tps_result), "grid" %in% names(tps_result))
  x_col <- trait_cols[1]; y_col <- trait_cols[2]
  stopifnot(all(c("Survived", x_col, y_col) %in% names(data_df)))
  tps_grid <- tps_result$grid
  stopifnot(all(c(x_col, y_col, ".fit") %in% names(tps_grid)))
  
  form <- stats::reformulate(c(x_col, y_col), response = "Survived")
  fit_glm <- stats::glm(form, family = stats::binomial(), data = data_df)
  b <- stats::coef(fit_glm)
  if (any(!is.finite(b))) stop("GLM coefficients not finite")
  
  # HL = (logit(p) - b0 - b1*TL)/b2
  logit <- function(p) log(p/(1-p))
  TL_seq <- seq(min(data_df[[x_col]], na.rm = TRUE),
                max(data_df[[x_col]], na.rm = TRUE), length.out = 200)
  
  contour_df_list <- list()
  for (p in levels) {
    eta <- logit(p)
    if (abs(b[[y_col]]) < .Machine$double.eps^0.5) next
    HL_vals <- (eta - b[["(Intercept)"]] - b[[x_col]] * TL_seq) / b[[y_col]]
    dfp <- data.frame(
      x = TL_seq,
      y = HL_vals,
      p = p
    )
    y_min <- min(data_df[[y_col]], na.rm = TRUE)
    y_max <- max(data_df[[y_col]], na.rm = TRUE)
    dfp <- dfp[dfp$y >= y_min & dfp$y <= y_max, , drop = FALSE]
    if (nrow(dfp)) contour_df_list[[length(contour_df_list) + 1]] <- dfp
  }
  contour_df <- if (length(contour_df_list)) do.call(rbind, contour_df_list) else
    data.frame(x = numeric(0), y = numeric(0), p = numeric(0))
  
  base_df <- data.frame(
    x = tps_grid[[x_col]],
    y = tps_grid[[y_col]],
    f = tps_grid$.fit
  )
  
  ggplot() +
    ggplot2::geom_contour_filled(
      data = base_df,
      mapping = ggplot2::aes(x = x, y = y, z = f),
      bins = bins
    ) +
    ggplot2::geom_contour(
      data = base_df,
      mapping = ggplot2::aes(x = x, y = y, z = f),
      colour = "black", alpha = 0.3
    ) +
    ggplot2::geom_path(
      data = contour_df,
      mapping = ggplot2::aes(x = x, y = y, group = p),
      colour = line_color, linewidth = line_size, alpha = line_alpha
    ) +
    ggplot2::labs(
      x = x_col, y = y_col, fill = "Probability",
      subtitle = paste0("Overlay: linear-GLM contours at p = ",
                        paste(levels, collapse = ", "))
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


save_overlays <- function(data_df, tps_results, pairs, out_dir,
                          bins = 15, grid_n = 30, levels = c(0.25, 0.5, 0.75)) {
  for (nm in names(pairs)) {
    pair <- pairs[[nm]]
    if (!is.null(tps_results[[nm]]) && !is.null(tps_results[[nm]]$tps)) {
      p <- overlay_glm_on_tps(
        data_df    = data_df,
        tps_result = tps_results[[nm]]$tps,
        trait_cols = pair,
        bins = bins, grid_n = grid_n, levels = levels
      )
      print(p)
      ggplot2::ggsave(
        file.path(out_dir, paste0("tps_", nm, "_overlay.png")),
        p, width = 10, height = 8, dpi = 300
      )
      cat("Overlay saved for", nm, "\n")
    } else {
      cat("Skipped overlay:", nm, "not successful\n")
    }
  }
}


cat("\n2. Loading Bumpus data...\n")
data_path <- "../test_data/Bumpus_data.csv"
if (!file.exists(data_path)) {
  # Fallback if running from project root
  data_path <- "test_data/Bumpus_data.csv"
}
if (!file.exists(data_path)) stop("Cannot find Bumpus_data.csv")

bumpus <- utils::read.csv(data_path)
bumpus$Survived <- as.numeric(bumpus$surv == "alive")
bumpus$TL <- bumpus$totlen
bumpus$HL <- bumpus$humer
bumpus$WT <- bumpus$wgt
bumpus$KL <- bumpus$stern
cat("Bumpus data loaded:", paste(dim(bumpus), collapse = " x "), "\n")
cat("Survival rate:", mean(bumpus$Survived), "\n")
cat("Columns:", paste(names(bumpus), collapse = ", "), "\n")

# ---- selection_coefficients (original scale to match slides) ----
cat("\n3. Testing selection_coefficients...\n")
model_configs <- list(
  "TL_only"    = list(traits = c("TL"),                   fitness_type = "binary"),
  "TL_HL"      = list(traits = c("TL", "HL"),             fitness_type = "binary"),
  "TL_HL_WT"   = list(traits = c("TL", "HL", "WT"),       fitness_type = "binary"),
  "Full_model" = list(traits = c("TL", "HL", "WT", "KL"), fitness_type = "binary")
)

selection_results <- list()
for (model_name in names(model_configs)) {
  config <- model_configs[[model_name]]
  cat("   Testing:", model_name, "- Traits:", paste(config$traits, collapse = ", "), "\n")
  tryCatch({
    result <- selection_coefficients(
      data = bumpus,
      fitness_col = "Survived",
      trait_cols = config$traits,
      fitness_type = config$fitness_type,
      standardize = FALSE,        
      use_relative_for_fit = FALSE
    )
    selection_results[[model_name]] <- result
    cat("Coefficients rows:", nrow(result), "\n")
    print(result)
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    selection_results[[model_name]] <- NULL
  })
}

# ---- Independent fns ----
cat("\n4. Testing independent functions...\n")
cat("Testing detect_family independently...\n")
family_info <- detect_family(bumpus$Survived)
cat("Detected family:\n"); print(family_info)

cat("Testing selection_differential for all traits...\n")
sel_diffs <- list()
for (trait in c("TL", "HL", "WT", "KL")) {
  tryCatch({
    sel_diff <- selection_differential(
      data = bumpus,
      fitness_col = "Survived",
      trait_col = trait,
      assume_standardized = FALSE,
      use_relative = FALSE
    )
    sel_diffs[[trait]] <- sel_diff
    cat(trait, "selection differential:", round(sel_diff, 4), "\n")
  }, error = function(e) {
    cat(trait, "error:", e$message, "\n")
  })
}

# ---- Univariate plots ----
cat("\n5. Testing visualization functions...\n")

univariate_plots <- list()
for (trait in c("TL", "HL")) {
  tryCatch({
    spline_result <- univariate_spline(
      data = bumpus,
      fitness_col = "Survived",
      trait_col = trait,
      fitness_type = "binary",
      k = 6
    )
    surface_plot <- univariate_surface(
      uni = spline_result,
      trait_col = trait,
      title = paste("Survival vs", trait)
    )
    univariate_plots[[trait]] <- list(spline = spline_result, plot = surface_plot)
    print(surface_plot)
    ggsave(file.path(output_dir, paste0("univariate_", trait, ".png")),
           surface_plot, width = 8, height = 6, dpi = 300)
    cat("Created & saved univariate surface for", trait, "\n")
  }, error = function(e) {
    cat("Univariate surface for", trait, "failed:", e$message, "\n")
  })
}

# ---- Bootstrap ----
cat("\n6. Testing bootstrap_selection...\n")
bootstrap_result <- NULL
tryCatch({
  set.seed(42)
  bootstrap_result <- bootstrap_selection(
    data = bumpus,
    fitness_col = "Survived",
    trait_cols = c("TL", "HL"),
    fitness_type = "binary",
    B = 100,
    seed = 42
  )
  cat("Bootstrap completed:", nrow(bootstrap_result$draws), "draws\n")
  cat("Confidence intervals:\n"); print(bootstrap_result$ci)
  utils::write.csv(bootstrap_result$ci, file.path(output_dir, "bootstrap_ci.csv"), row.names = FALSE)
  cat("Saved bootstrap results\n")
}, error = function(e) {
  cat("Bootstrap failed:", e$message, "\n")
})

# ---- Bivariate surfaces (GAM auto for binary) ----
cat("\n=== TESTING CORRELATIONAL TPS ===\n\n")
all_trait_pairs <- list(
  c("TL", "HL"),
  c("TL", "WT"),
  c("HL", "KL"),
  c("TL", "KL"),
  c("HL", "WT")
)

tps_results <- list()
for (pair in all_trait_pairs) {
  pair_name <- paste(pair, collapse = "_")
  cat("Testing:", pair_name, "\n")
  if (all(pair %in% names(bumpus))) {
    tryCatch({
      tps_result <- correlational_tps(
        data = bumpus,
        fitness_col = "Survived",
        trait_cols = pair,
        grid_n = 30
        # , method = "gam", k = 3  
      )
      cat("TPS SUCCESS - Method:", tps_result$method, "\n")
      cat("Range:", round(range(tps_result$grid$.fit), 4), "\n")

      corr_plot <- correlation_surface(tps_result, pair, bins = 15)
      print(corr_plot)
      ggsave(file.path(output_dir, paste0("tps_", pair_name, ".png")),
             corr_plot, width = 10, height = 8, dpi = 300)
      cat("Surface plot created & saved\n")

      tps_results[[pair_name]] <- list(tps = tps_result, plot = corr_plot, status = "success")
    }, error = function(e) {
      cat("Failed:", e$message, "\n")
      tps_results[[pair_name]] <- list(tps = NULL, plot = NULL, status = "failed", error = e$message)
    })
  } else {
    cat("Traits not found in data:", paste(pair[!pair %in% names(bumpus)], collapse = ", "), "\n")
    tps_results[[pair_name]] <- list(tps = NULL, plot = NULL, status = "traits_missing")
  }
  cat("\n")
}


pairs_to_overlay <- list(
  "TL_HL" = c("TL","HL"),
  "TL_WT" = c("TL","WT"),
  "TL_KL" = c("TL","KL")
)
save_overlays(bumpus, tps_results, pairs_to_overlay, output_dir,
              bins = 15, grid_n = 30, levels = c(0.25, 0.5, 0.75))

# ---- TPS summary ----
cat("=== TPS ANALYSIS SUMMARY ===\n")
success_count <- sum(vapply(tps_results, function(x) is.list(x) && identical(x$status, "success"), logical(1)))
total_count <- length(tps_results)
cat("Successful TPS analyses:", success_count, "/", total_count, "\n")
if (success_count > 0) {
  successful_pairs <- names(tps_results)[vapply(tps_results, function(x) identical(x$status, "success"), logical(1))]
  cat("Successful pairs:", paste(successful_pairs, collapse = ", "), "\n\n")
}
cat("Key Results:\n")
for (pair_name in names(tps_results)) {
  result <- tps_results[[pair_name]]
  if (is.list(result) && identical(result$status, "success")) {
    tps_data <- result$tps
    fit_range <- tryCatch(range(tps_data$grid$.fit), error = function(e) c(NA, NA))
    cat("  ", pair_name, ":\n")
    if (!any(is.na(fit_range))) {
      cat("    Survival probability range:", round(fit_range[1], 4), "-", round(fit_range[2], 4), "\n")
    }
    max_fit_idx <- tryCatch(which.max(tps_data$grid$.fit), error = function(e) 1)
    optimal_point <- tryCatch(tps_data$grid[max_fit_idx, ], error = function(e) NULL)
    if (!is.null(optimal_point)) {
      traits <- strsplit(pair_name, "_")[[1]]
      cat("    Optimal survival at:",
          traits[1], "=", round(optimal_point[[1]], 2), ",",
          traits[2], "=", round(optimal_point[[2]], 3), "\n")
      cat("    Max survival probability:", round(optimal_point$.fit, 4), "\n\n")
    }
  }
}

# ---- Slide-aligned GLMs + Robust coefficient comparison ----
cat("\n9. COMPARISON WITH STAT 226 RESULTS (robust)…\n")
slide_models <- list()
tryCatch({
  slide_models[["TL_only"]]    <- stats::glm(Survived ~ TL, family = stats::binomial, data = bumpus)
  slide_models[["TL_HL"]]      <- stats::glm(Survived ~ TL + HL, family = stats::binomial, data = bumpus)
  slide_models[["TL_HL_WT"]]   <- stats::glm(Survived ~ TL + HL + WT, family = stats::binomial, data = bumpus)
  slide_models[["Full_model"]] <- stats::glm(Survived ~ TL + HL + WT + KL, family = stats::binomial, data = bumpus)

  guess_col <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit)) hit[1] else NA_character_
  }
  normalize_term <- function(x) {
    x <- as.character(x)
    x <- gsub("`", "", x)
    x <- gsub("\\s+", "", x)
    x <- gsub("^z_", "", x)
    x <- gsub("_std$", "", x)
    x <- gsub("^scale\\(([^)]+)\\)$", "\\1", x)
    x <- gsub("^s\\(([^)]+)\\)$", "\\1", x)
    x <- gsub("^te\\(([^)]+)\\)$", "\\1", x)
    tolower(x)
  }
  lookup_estimate <- function(your_df, term_col, est_col, var) {
    var_norm <- normalize_term(var)
    terms_norm <- normalize_term(your_df[[term_col]])
    idx <- which(terms_norm == var_norm)
    if (length(idx) >= 1) return(as.numeric(your_df[[est_col]][idx[1]]))
    if (var_norm %in% c("(intercept)", "intercept")) {
      idx2 <- grep("(intercept)", terms_norm, fixed = TRUE)
      if (length(idx2) >= 1) return(as.numeric(your_df[[est_col]][idx2[1]]))
    }
    idx3 <- which(grepl(paste0("(^|[^a-z])", var_norm, "($|[^a-z])"), terms_norm))
    if (length(idx3) >= 1) return(as.numeric(your_df[[est_col]][idx3[1]]))
    NA_real_
  }

  for (model_name in names(slide_models)) {
    if (model_name %in% names(selection_results) && !is.null(selection_results[[model_name]])) {
      slide_coef <- stats::coef(slide_models[[model_name]])
      your_df    <- selection_results[[model_name]]

      term_col <- guess_col(your_df, c("trait","term","variable","name","param","feature"))
      est_col  <- guess_col(your_df, c("estimate","coef","beta","value","estimate_beta","coefficient"))
      if (is.na(term_col) || is.na(est_col)) {
        cat("   ", model_name, ": could not guess columns (term/estimate). Available cols:",
            paste(names(your_df), collapse = ", "), "\n")
        next
      }

      cat("   ", model_name, "model (term_col =", term_col, ", est_col =", est_col, ")\n")
      vars <- names(slide_coef)
      comp <- data.frame(
        Variable      = vars,
        STAT226       = round(as.numeric(slide_coef), 6),
        Your_Function = NA_real_,
        Difference    = NA_real_,
        stringsAsFactors = FALSE
      )
      for (i in seq_len(nrow(comp))) {
        v <- comp$Variable[i]
        est <- lookup_estimate(your_df, term_col, est_col, v)
        comp$Your_Function[i] <- if (is.finite(est)) round(est, 6) else NA_real_
        if (!is.na(comp$Your_Function[i])) {
          comp$Difference[i] <- round(abs(comp$STAT226[i] - comp$Your_Function[i]), 6)
        }
      }
      if (all(is.na(comp$Your_Function[comp$Variable == "(Intercept)"]))) {
        comp <- comp[comp$Variable != "(Intercept)", , drop = FALSE]
      }
      print(comp)
      matched <- sum(!is.na(comp$Your_Function))
      total   <- nrow(comp)
      cat("     Matched", matched, "/", total, "coefficients\n\n")
    }
  }
}, error = function(e) {
  cat(message, "\n")
})

# ---- Final saving summary ----
cat("\n10. Saving summary complete. Files written to:", normalizePath(output_dir, winslash = "/"), "\n")

# ---- Final summary ----
cat("\n11. FINAL SUMMARY...\n")
cat("Functions tested:", length(core_function_files), " (loaded:", length(loaded_files), ")\n")
cat("Models analyzed:", length(selection_results), "\n") 
cat("Univariate visualizations:", length(Filter(Negate(is.null), lapply(univariate_plots, `[[`, "plot"))), "\n")
cat("TPS analyses:", success_count, "/", total_count, "successful\n")
cat("Bootstrap analysis:", ifelse(!is.null(bootstrap_result), "COMPLETED", "FAILED"), "\n")
cat("STAT 226 comparison: COMPLETED (see console)\n")
cat("Output saved to:", output_dir, "\n")
cat("\n", strrep("=", 60), "\n")
if (success_count > 0) {
  cat("TESTING COMPLETED WITH SUCCESS\n")
} else {
  cat("TESTING COMPLETED WITH ISSUES\n")
}
cat(strrep("=", 60), "\n")
