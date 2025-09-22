cat("=== COMPLETE SELECTION ANALYSIS WORKFLOW ===\n\n")

# === Dependency Check & Auto-Install ===
required_packages <- c("ggplot2", "fields", "car")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# --- Ensure script runs from its own folder ---
args <- commandArgs(trailingOnly = FALSE)
scriptPath <- NULL
fileArg <- grep("^--file=", args, value = TRUE)
if (length(fileArg) > 0) {
  scriptPath <- normalizePath(sub("^--file=", "", fileArg))
} else {
  scriptPath <- normalizePath(".")
}
this.dir <- dirname(scriptPath)
setwd(this.dir)

# 1. set up environment
cat("1. Setting up environment...\n")
start_time <- Sys.time()

# create results folders
create_directories <- function() {
  dirs <- c("results", "results/data", "results/tables", "results/figures", "results/summary")
  for (dir in dirs) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
}
create_directories()

# 2. load all the functions
cat("2. Loading all functions...\n")
function_files <- c(
  "prepare_selection_data.R", "analyze_linear_selection.R", "analyze_nonlinear_selection.R",
  "extract_results.R", "selection_coefficients.R", "detect_family.R", "selection_differential.R",
  "univariate_spline.R", "univariate_surface.R", "correlational_tps.R", "correlation_surface.R",
  "bootstrap_selection.R"
)

for (file in function_files) {
  if (file.exists(file)) {
    source(file)
    cat(file, "\n")
  } else {
    cat(file, "not found\n")
  }
}

# 3. build dataset 
cat("3. Creating test data...\n")
set.seed(42)

test_data <- list(
  binary = data.frame(
    fitness = rbinom(50, 1, 0.6),
    size = rnorm(50, 5, 1),
    speed = rnorm(50, 10, 2),
    color = rnorm(50, 3, 0.5)
  ),
  
  continuous = data.frame(
    fitness = rpois(50, 3) + 1,
    size = rnorm(50, 5, 1),
    speed = rnorm(50, 10, 2),
    color = rnorm(50, 3, 0.5)
  )
)

cat("Created binary data:", dim(test_data$binary), "\n")
cat("Created continuous data:", dim(test_data$continuous), "\n")

# 4. run core tests 
cat("4. Running core tests...\n")
run_core_tests <- function(data) {
  results <- list()
  
  # Test selection_coefficients
  results$binary_coef <- selection_coefficients(
    data$binary, "fitness", c("size", "speed"), "binary", TRUE, TRUE
  )
  
  results$continuous_coef <- selection_coefficients(
    data$continuous, "fitness", c("size", "speed"), "continuous", TRUE, TRUE
  )
  
  # Test selection_differential
  results$selection_diff <- selection_differential(
    data$continuous, "fitness", "size", FALSE, TRUE
  )
  
  # Test detect_family
  results$family_binary <- detect_family(data$binary$fitness)
  results$family_continuous <- detect_family(data$continuous$fitness)
  
  return(results)
}

core_results <- run_core_tests(test_data)
cat("Core tests completed\n")

# 5. run visualization 
cat("5. Running visualization tests...\n")
run_extended_viz_tests <- function(data) {
  results <- list()
  
  # 5.1 testing uni variate plots for continuous data 
  cat("Testing univariate surfaces for continuous data...\n")
  results$univariate_continuous <- list()
  for (trait in c("size", "speed", "color")) {
    spline <- univariate_spline(data$continuous, "fitness", trait, "continuous", k = 8)
    plot <- univariate_surface(spline, trait)
    results$univariate_continuous[[trait]] <- list(spline = spline, plot = plot)
    cat("Continuous", trait, "surface\n")
  }
  
  # 5.2 testing uni variate plots for binary data 
  cat("Testing univariate surfaces for binary data...\n")
  results$univariate_binary <- list()
  for (trait in c("size", "speed")) {
    spline <- univariate_spline(data$binary, "fitness", trait, "binary", k = 8)
    plot <- univariate_surface(spline, trait)
    results$univariate_binary[[trait]] <- list(spline = spline, plot = plot)
    cat("Binary", trait, "surface\n")
  }
  
  # 5.3 testing correlation surfaces for continuous data
  cat("Testing correlation surfaces for continuous data...\n")
  results$correlation_continuous <- list()
  trait_pairs <- list(c("size", "speed"), c("size", "color"), c("speed", "color"))
  
  for (pair in trait_pairs) {
    tps <- correlational_tps(data$continuous, "fitness", pair, TRUE, 30)
    plot <- correlation_surface(tps, pair)
    results$correlation_continuous[[paste(pair, collapse = "_")]] <- list(tps = tps, plot = plot)
    cat("     ✓ Correlation", paste(pair, collapse = " vs "), "\n")
  }
  
  # 5.4 testing correlation surfaces for binary data
  cat("Testing correlation surfaces for binary data...\n")
  results$correlation_binary <- list()
  for (pair in trait_pairs[1:2]) {  
    tps <- correlational_tps(data$binary, "fitness", pair, FALSE, 25)
    plot <- correlation_surface(tps, pair)
    results$correlation_binary[[paste(pair, collapse = "_")]] <- list(tps = tps, plot = plot)
    cat("     ✓ Binary correlation", paste(pair, collapse = " vs "), "\n")
  }
  
  # 5.5 testing different spline complexities 
  cat("Testing different spline complexities...\n")
  results$spline_complexity <- list()
  k_values <- c(5, 8, 12)
  for (k in k_values) {
    spline <- univariate_spline(data$continuous, "fitness", "size", "continuous", k = k)
    plot <- univariate_surface(spline, "size")
    results$spline_complexity[[paste0("k", k)]] <- list(spline = spline, plot = plot)
    cat("Spline k =", k, "\n")
  }
  
  return(results)
}

viz_results <- run_extended_viz_tests(test_data)
cat("visualization tests completed\n")

# 6. run bootstrap tests 
cat("6. Running bootstrap tests...\n")
run_bootstrap_tests <- function(data) {
  bootstrap_selection(
    data$continuous, "fitness", c("size", "speed"), "continuous", B = 50, seed = 42
  )
}

bootstrap_results <- run_bootstrap_tests(test_data)
cat("Bootstrap tests completed\n")

# 7. save all the results 
cat("7. Saving all results...\n")
save_all_results <- function(core_results, viz_results, bootstrap_results, test_data) {
  
  # save core results 
  saveRDS(core_results, "results/data/core_test_results.rds")
  cat("Core results saved\n")
  
  # save bootstrap results 
  saveRDS(bootstrap_results, "results/data/bootstrap_results.rds")
  cat("Bootstrap results saved\n")
  
  # save test data 
  saveRDS(test_data, "results/data/test_data.rds")
  cat("Test data saved\n")
  
  # save csv 
  write.csv(core_results$binary_coef, "results/tables/binary_coefficients.csv", row.names = FALSE)
  write.csv(core_results$continuous_coef, "results/tables/continuous_coefficients.csv", row.names = FALSE)
  write.csv(bootstrap_results$ci, "results/tables/bootstrap_ci.csv", row.names = FALSE)
  cat("CSV tables saved\n")
  
  # 7.1 save plots 
  cat("Saving univariate plots...\n")
  for (trait in names(viz_results$univariate_continuous)) {
    filename <- paste0("results/figures/univariate_continuous_", trait, ".png")
    ggsave(filename, viz_results$univariate_continuous[[trait]]$plot, width = 8, height = 6)
    cat(filename, "\n")
  }
  
  for (trait in names(viz_results$univariate_binary)) {
    filename <- paste0("results/figures/univariate_binary_", trait, ".png")
    ggsave(filename, viz_results$univariate_binary[[trait]]$plot, width = 8, height = 6)
    cat(filename, "\n")
  }
  
  # 7.2 save correlation plots 
  cat("Saving correlation plots...\n")
  for (pair_name in names(viz_results$correlation_continuous)) {
    filename <- paste0("results/figures/correlation_continuous_", pair_name, ".png")
    ggsave(filename, viz_results$correlation_continuous[[pair_name]]$plot, width = 8, height = 6)
    cat(filename, "\n")
  }
  
  for (pair_name in names(viz_results$correlation_binary)) {
    filename <- paste0("results/figures/correlation_binary_", pair_name, ".png")
    ggsave(filename, viz_results$correlation_binary[[pair_name]]$plot, width = 8, height = 6)
    cat(filename, "\n")
  }
  
  # 7.3 save spline complexity plots 
  cat("Saving spline complexity plots...\n")
  for (k_name in names(viz_results$spline_complexity)) {
    filename <- paste0("results/figures/spline_", k_name, ".png")
    ggsave(filename, viz_results$spline_complexity[[k_name]]$plot, width = 8, height = 6)
    cat(filename, "\n")
  }
  
  cat("All plots saved\n")
}

save_all_results(core_results, viz_results, bootstrap_results, test_data)

# 8. generate test reports 
cat("8. Generating test report...\n")
generate_test_report <- function(core_results, viz_results, bootstrap_results) {
  
  # calc number of plots 
  total_plots <- length(viz_results$univariate_continuous) + 
    length(viz_results$univariate_binary) +
    length(viz_results$correlation_continuous) +
    length(viz_results$correlation_binary) +
    length(viz_results$spline_complexity)
  
  report <- paste(
    "SELECTION ANALYSIS TEST REPORT",
    strrep("=", 50),
    paste("Test completed:", Sys.time()),
    paste("Duration:", round(difftime(Sys.time(), start_time, units = "secs"), 1), "seconds"),
    "",
    "TEST RESULTS SUMMARY:",
    paste("- Binary coefficients:", nrow(core_results$binary_coef), "rows"),
    paste("- Continuous coefficients:", nrow(core_results$continuous_coef), "rows"),
    paste("- Selection differential:", round(core_results$selection_diff, 6)),
    paste("- Bootstrap CIs:", nrow(bootstrap_results$ci), "intervals"),
    "",
    "VISUALIZATION SUMMARY:",
    paste("- Total plots generated:", total_plots),
    paste("- Continuous univariate plots:", length(viz_results$univariate_continuous)),
    paste("- Binary univariate plots:", length(viz_results$univariate_binary)),
    paste("- Continuous correlation plots:", length(viz_results$correlation_continuous)),
    paste("- Binary correlation plots:", length(viz_results$correlation_binary)),
    paste("- Spline complexity plots:", length(viz_results$spline_complexity)),
    "",
    "FAMILY DETECTION:",
    paste("- Binary data:", core_results$family_binary$type),
    paste("- Continuous data:", core_results$family_continuous$type),
    "",
    "SUCCESS: All tests and visualizations completed successfully!",
    sep = "\n"
  )
  
  writeLines(report, "results/summary/test_report.txt")
  return(report)
}

test_report <- generate_test_report(core_results, viz_results, bootstrap_results)
cat("   ✓ Test report generated\n")

# 9. final output 
cat("\n9. Final output...\n")
cat(strrep("=", 60), "\n")
cat("EXTENDED ANALYSIS COMPLETE - ALL TESTS PASSED!\n")
cat(strrep("=", 60), "\n\n")

cat("Generated files:\n")
cat("----------------\n")

# list all the files 
result_files <- list.files("results", recursive = TRUE)
for (file in result_files) {
  cat("•", file, "\n")
}

cat("\nTotal plots generated:", 
    length(viz_results$univariate_continuous) + 
      length(viz_results$univariate_binary) +
      length(viz_results$correlation_continuous) + 
      length(viz_results$correlation_binary) +
      length(viz_results$spline_complexity), "\n\n")

cat("Test summary:\n")
cat("-------------\n")
cat(test_report, "\n")

cat("\n", strrep("=", 60), "\n", sep = "")
cat("COMPLETE ANALYSIS WORKFLOW FINISHED SUCCESSFULLY!\n")
cat(strrep("=", 60), "\n", sep = "")
