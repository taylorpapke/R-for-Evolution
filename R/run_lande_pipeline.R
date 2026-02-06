# ==============================================================================
# LANDE PACKAGE: MASTER PIPELINE (The "Pure Matrix" Fix)
# ==============================================================================

# 1. SETUP & LOADING
r_dir <- "~/Github/R-for-Evolution/R/"
core_modules <- c(
  "prepare_selection_data.R", "detect_family.R", "analyze_linear_selection.R",
  "analyze_nonlinear_selection.R", "extract_results.R", "selection_analysis.R",
  "univariate_spline.R", "univariate_surface.R", 
  "correlational_tps.R", "correlation_surface.R"
)

# Ensure dependencies are active
library(car); library(mgcv); library(fields); library(ggplot2)

for (script in core_modules) {
  path <- paste0(r_dir, script)
  if (file.exists(path)) source(path)
}

# 2. THE FINAL HOT-FIX 
# This version explicitly destroys all "list" attributes to prevent the coercion error
correlational_tps <- function(data, fitness_col, trait_cols, use_relative = TRUE, grid_n = 60) {
  
  # Filter complete cases
  df <- data[stats::complete.cases(data[, c(fitness_col, trait_cols)]), ]
  
  # Extract numeric vectors
  x1 <- as.numeric(df[[trait_cols[1]]])
  x2 <- as.numeric(df[[trait_cols[2]]])
  y  <- as.numeric(df[[fitness_col]])
  
  # Standardize fitness if not binary
  if (use_relative && length(unique(y)) > 2) y <- y / mean(y, na.rm = TRUE)
  
  # Fit the Tps model using a clean matrix
  X_mat <- cbind(x1, x2)
  storage.mode(X_mat) <- "double"
  fit <- fields::Tps(X_mat, y)
  
  # CREATE THE PREDICTION GRID
  x1_seq <- seq(min(x1), max(x1), length.out = grid_n)
  x2_seq <- seq(min(x2), max(x2), length.out = grid_n)
  grid_df <- expand.grid(x1 = x1_seq, x2 = x2_seq)
  names(grid_df) <- trait_cols
  
  # --- THE FIX: Convert to a naked numeric matrix ---
  # We use unlist and matrix to ensure NO list attributes remain
  Z_mat <- matrix(as.numeric(as.matrix(grid_df)), ncol = 2)
  storage.mode(Z_mat) <- "double" 
  
  # Predict and attach to grid
  grid_df$.fit <- as.numeric(stats::predict(fit, Z_mat))
  
  return(list(model = fit, grid = grid_df))
}

# 3. GENERATE TEST DATA
set.seed(42)
n <- 500
trait1 <- rnorm(n, 15, 3); trait2 <- rnorm(n, 5, 2)
# Selection primarily on Trait 1
prob <- 1 / (1 + exp(-(0.6 * scale(trait1) - 0.2 * scale(trait2))))
fitness <- rbinom(n, 1, prob)
df <- data.frame(fitness = fitness, trait1 = trait1, trait2 = trait2)

# 4. EXECUTE PIPELINE
df_ready <- prepare_selection_data(df, fitness_col = "fitness", trait_cols = c("trait1", "trait2"))
results  <- selection_analysis(df_ready, fitness_col = "fitness", trait_cols = c("trait1", "trait2"))

# 5. PRINT RESULTS
cat("\n--- LANDE-ARNOLD SELECTION GRADIENTS ---\n")
print(results$coefficients_table)

# 6. VISUAL VALIDATION
# Univariate
uni_obj <- univariate_spline(df_ready, fitness_col = "fitness", trait_col = "trait1")
print(univariate_surface(uni_obj, trait_col = "trait1"))

# Multivariate (Using the redefined function)
tps_obj <- correlational_tps(df_ready, fitness_col = "fitness", trait_cols = c("trait1", "trait2"))
print(correlation_surface(tps_obj, trait_cols = c("trait1", "trait2")))