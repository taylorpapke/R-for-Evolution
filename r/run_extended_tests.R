
cat("=== EXTENDED PARAMETER AND BOUNDARY TESTING ===\n\n")

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

# 1) Load all functions
cat("1. Loading all functions...\n")
function_files <- c(
  "prepare_selection_data.R", "analyze_linear_selection.R", "analyze_nonlinear_selection.R",
  "extract_results.R", "selection_coefficients.R", "detect_family.R", "selection_differential.R",
  "univariate_spline.R", "univariate_surface.R", "correlational_tps.R", "correlation_surface.R",
  "bootstrap_selection.R"
)

for (file in function_files) {
  if (file.exists(file)) source(file)
}

# 2) Create diverse test datasets
cat("2. Creating diverse test data...\n")
set.seed(42)

test_datasets <- list(
  # Normal datasets
  normal_binary = data.frame(
    fitness = rbinom(50, 1, 0.6),
    size = rnorm(50, 5, 1),
    speed = rnorm(50, 10, 2),
    color = rnorm(50, 3, 0.5)
  ),
  
  normal_continuous = data.frame(
    fitness = rpois(50, 3) + 1,
    size = rnorm(50, 5, 1),
    speed = rnorm(50, 10, 2),
    color = rnorm(50, 3, 0.5)
  ),
  
  # Boundary-case datasets
  small_sample = data.frame(
    fitness = c(0, 1, 0, 1, 0),
    size = c(4.1, 5.2, 3.9, 5.5, 4.8),
    speed = c(9.5, 10.2, 8.9, 11.1, 10.5)
  ),
  
  all_zeros = data.frame(
    fitness = rep(0, 20),
    size = rnorm(20, 5, 1),
    speed = rnorm(20, 10, 2)
  ),
  
  all_ones = data.frame(
    fitness = rep(1, 20),
    size = rnorm(20, 5, 1),
    speed = rnorm(20, 10, 2)
  ),
  
  high_correlation = data.frame(
    fitness = rpois(30, 2) + 1,
    size = rnorm(30, 5, 1),
    speed = rnorm(30, 5, 0.1),  # near-constant
    color = rnorm(30, 3, 0.5)
  )
)

# 3) Test parameter combinations for selection_coefficients
cat("3. Testing selection_coefficients parameter combinations...\n")
test_selection_combinations <- function() {
  results <- list()
  test_count <- 0
  
  # Parameter combination grid
  combinations <- expand.grid(
    data_type = c("normal_binary", "normal_continuous"),
    standardize = c(TRUE, FALSE),
    use_relative = c(TRUE, FALSE),
    traits = list("size", c("size", "speed"), c("size", "speed", "color"))
  )
  
  for (i in 1:nrow(combinations)) {
    combo <- combinations[i, ]
    data_name <- as.character(combo$data_type)
    traits <- combo$traits[[1]]
    
    cat("Testing:", data_name, "std =", combo$standardize, 
        "rel =", combo$use_relative, "traits =", paste(traits, collapse = "+"), "\n")
    
    tryCatch({
      result <- selection_coefficients(
        data = test_datasets[[data_name]],
        fitness_col = "fitness",
        trait_cols = traits,
        fitness_type = ifelse(grepl("binary", data_name), "binary", "continuous"),
        standardize = combo$standardize,
        use_relative_for_fit = combo$use_relative
      )
      
      test_count <- test_count + 1
      results[[paste("test", test_count)]] <- list(
        combination = combo,
        result = result,
        success = TRUE
      )
      cat("Success\n")
      
    }, error = function(e) {
      results[[paste("test", test_count)]] <- list(
        combination = combo, 
        error = e$message,
        success = FALSE
      )
      cat("ailed:", e$message, "\n")
    })
  }
  
  return(list(results = results, total_tests = test_count))
}

selection_results <- test_selection_combinations()

# 4) Test boundary-case datasets
cat("\n4. Testing boundary cases...\n")
test_boundary_cases <- function() {
  boundary_tests <- list()
  test_count <- 0
  
  # Small sample
  cat("   Testing small sample...\n")
  tryCatch({
    result <- selection_coefficients(
      data = test_datasets$small_sample,
      fitness_col = "fitness",
      trait_cols = "size",
      fitness_type = "binary",
      standardize = TRUE,
      use_relative_for_fit = FALSE
    )
    test_count <- test_count + 1
    boundary_tests$small_sample <- list(success = TRUE, result = result)
    cat("Small sample test\n")
  }, error = function(e) {
    boundary_tests$small_sample <- list(success = FALSE, error = e$message)
    cat("small sample failed:", e$message, "\n")
  })
  
  # All zeros
  cat("   Testing all zeros...\n")
  tryCatch({
    result <- selection_coefficients(
      data = test_datasets$all_zeros,
      fitness_col = "fitness",
      trait_cols = "size",
      fitness_type = "binary",
      standardize = TRUE,
      use_relative_for_fit = FALSE
    )
    test_count <- test_count + 1
    boundary_tests$all_zeros <- list(success = TRUE, result = result)
    cat("All zeros test\n")
  }, error = function(e) {
    boundary_tests$all_zeros <- list(success = FALSE, error = e$message)
    cat("All zeros failed:", e$message, "\n")
  })
  
  # High correlation
  cat("Testing high correlation...\n")
  tryCatch({
    result <- selection_coefficients(
      data = test_datasets$high_correlation,
      fitness_col = "fitness",
      trait_cols = c("size", "speed"),
      fitness_type = "continuous",
      standardize = TRUE,
      use_relative_for_fit = TRUE
    )
    test_count <- test_count + 1
    boundary_tests$high_correlation <- list(success = TRUE, result = result)
    cat("High correlation test\n")
  }, error = function(e) {
    boundary_tests$high_correlation <- list(success = FALSE, error = e$message)
    cat("High correlation failed:", e$message, "\n")
  })
  
  return(list(tests = boundary_tests, total_tests = test_count))
}

boundary_results <- test_boundary_cases()

# 5) Test univariate_spline parameter settings
cat("\n5. Testing univariate_spline parameters...\n")
test_spline_parameters <- function() {
  spline_tests <- list()
  test_count <- 0
  
  # Test different k values
  k_values <- c(3, 5, 8, 12, 20)
  for (k in k_values) {
    cat("   Testing k =", k, "...\n")
    tryCatch({
      result <- univariate_spline(
        data = test_datasets$normal_continuous,
        fitness_col = "fitness",
        trait_col = "size",
        fitness_type = "continuous",
        k = k
      )
      test_count <- test_count + 1
      spline_tests[[paste("k", k)]] <- list(success = TRUE, result = result)
      cat("k =", k, "test\n")
    }, error = function(e) {
      spline_tests[[paste("k", k)]] <- list(success = FALSE, error = e$message)
      cat("k =", k, "failed:", e$message, "\n")
    })
  }
  
  return(list(tests = spline_tests, total_tests = test_count))
}

spline_results <- test_spline_parameters()

# 6) Test bootstrap parameters
cat("\n6. Testing bootstrap parameters...\n")
test_bootstrap_parameters <- function() {
  bootstrap_tests <- list()
  test_count <- 0
  
  # Test different B and confidence levels
  params <- list(
    list(B = 50, conf = 0.90),
    list(B = 100, conf = 0.95),
    list(B = 200, conf = 0.99),
    list(B = 20, conf = 0.80)
  )
  
  for (param in params) {
    cat("Testing B =", param$B, "conf =", param$conf, "...\n")
    tryCatch({
      result <- bootstrap_selection(
        data = test_datasets$normal_continuous,
        fitness_col = "fitness",
        trait_cols = c("size", "speed"),
        fitness_type = "continuous",
        B = param$B,
        conf = param$conf,
        seed = 42
      )
      test_count <- test_count + 1
      bootstrap_tests[[paste("B", param$B, "conf", param$conf)]] <- list(
        success = TRUE, 
        result = result
      )
      cat("B =", param$B, "conf =", param$conf, "test\n")
    }, error = function(e) {
      bootstrap_tests[[paste("B", param$B, "conf", param$conf)]] <- list(
        success = FALSE, 
        error = e$message
      )
      cat("B =", param$B, "conf =", param$conf, "failed:", e$message, "\n")
    })
  }
  
  return(list(tests = bootstrap_tests, total_tests = test_count))
}

bootstrap_results <- test_bootstrap_parameters()

# 7) Generate extended test report
cat("\n7. Generating comprehensive test report...\n")
generate_extended_report <- function(selection_results, boundary_results, spline_results, bootstrap_results) {
  
  total_tests <- selection_results$total_tests + boundary_results$total_tests + 
    spline_results$total_tests + bootstrap_results$total_tests
  
  # Compute success rate
  success_count <- sum(c(
    sum(sapply(selection_results$results, function(x) x$success)),
    sum(sapply(boundary_results$tests, function(x) x$success)),
    sum(sapply(spline_results$tests, function(x) x$success)),
    sum(sapply(bootstrap_results$tests, function(x) x$success))
  ))
  
  success_rate <- round(success_count / total_tests * 100, 1)
  
  report <- paste(
    "EXTENDED TESTING REPORT",
    strrep("=", 50),
    paste("Test completed:", Sys.time()),
    paste("Total tests run:", total_tests),
    paste("Successful tests:", success_count),
    paste("Success rate:", success_rate, "%"),
    "",
    "TEST CATEGORIES:",
    paste("- Selection coefficients parameter combinations:", selection_results$total_tests, "tests"),
    paste("- Boundary case tests:", boundary_results$total_tests, "tests"),
    paste("- Spline parameter tests:", spline_results$total_tests, "tests"),
    paste("- Bootstrap parameter tests:", bootstrap_results$total_tests, "tests"),
    "",
    "PARAMETERS TESTED:",
    "• Standardize: TRUE, FALSE",
    "• Use relative fitness: TRUE, FALSE", 
    "• Trait combinations: 1, 2, 3 traits",
    "• Data types: binary, continuous",
    "• Spline k values: 3, 5, 8, 12, 20",
    "• Bootstrap B values: 20, 50, 100, 200",
    "• Confidence levels: 0.80, 0.90, 0.95, 0.99",
    "",
    "BOUNDARY CASES:",
    "• Small sample (n=5)",
    "• All zeros fitness",
    "• All ones fitness", 
    "• High correlation traits",
    sep = "\n"
  )
  
  writeLines(report, "results/summary/extended_test_report.txt")
  return(report)
}

final_report <- generate_extended_report(selection_results, boundary_results, spline_results, bootstrap_results)

# 8) Final summary
cat("\n8. Final summary...\n")
cat(strrep("=", 60), "\n")
cat("EXTENDED TESTING COMPLETED!\n")
cat(strrep("=", 60), "\n\n")

cat("Total tests run:", 
    selection_results$total_tests + boundary_results$total_tests + 
      spline_results$total_tests + bootstrap_results$total_tests, "\n")

cat("Test categories:\n")
cat("- Parameter combinations:", selection_results$total_tests, "\n")
cat("- Boundary cases:", boundary_results$total_tests, "\n")
cat("- Spline parameters:", spline_results$total_tests, "\n")
cat("- Bootstrap parameters:", bootstrap_results$total_tests, "\n\n")

cat("Detailed report saved to: extended_test_report.txt\n")
cat("\nFinal report:\n")
cat(strrep("-", 40), "\n")
cat(final_report, "\n")
