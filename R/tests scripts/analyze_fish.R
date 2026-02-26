
## 1. Set up output directory -------------------------------------------------

results_dir <- "function_test_results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
  cat("   Created output directory:", results_dir, "\n")
} else {
  cat("   Using existing output directory:", results_dir, "\n")
}

## 2. Load data and function definitions -------------------------------------
cat("\n2. LOADING DATA AND FUNCTIONS...\n")

# Load Crescent Pond puppyfish dataset
data_dir <- file.path("..", "test_data")
crescent_path <- NULL

# Try exact filenames (URL-encoded and standard)
possible_names <- c(
  "Crescent+Pond+-+size-corrected+trait+data+++survival+++growth+++d13C+++d15N.csv",
  "Crescent Pond - size-corrected trait data + survival + growth + d13C + d15N.csv"
)

for (f in possible_names) {
  p <- file.path(data_dir, f)
  if (file.exists(p)) { crescent_path <- p; break }
}

# Fallback: fuzzy search in the data directory if not found
if (is.null(crescent_path) && dir.exists(data_dir)) {
  candidates <- list.files(data_dir, pattern = "Crescent.*Pond.*\\.csv$", full.names = TRUE)
  if (length(candidates) > 0) crescent_path <- candidates[1]
}

if (is.null(crescent_path) || !file.exists(crescent_path)) stop("Data file 'Crescent Pond...' not found in 'R/test_data/' directory.")
cat("   Using data file:", crescent_path, "\n")
crescent_data <- read.csv(crescent_path)

# Define fitness variables and trait set
fitness_binary     <- "survival"
fitness_continuous <- "ln.growth"
morphological_traits <- c("jaw", "eye", "body", "nasal", "mouth", "SL")

# Assemble analysis dataset
test_data <- crescent_data[, c(fitness_binary, fitness_continuous, morphological_traits)]
test_data <- test_data[complete.cases(test_data), ]

cat("   Test dataset contains", nrow(test_data), "complete observations\n")

# Load core analysis functions
function_files <- c(
  "prepare_selection_data.R", "selection_coefficients.R",
  "univariate_spline.R", "univariate_surface.R",
  "correlational_tps.R", "correlation_surface.R"
)

for (func_file in function_files) {
  file_path <- file.path("..", func_file)
  if (file.exists(file_path)) {
    source(file_path)
    cat("   Loaded:", func_file, "\n")
  }
}

## 2.1 Load required packages ------------------------------------------------
cat("   2.1 LOADING PACKAGES...\n")

required_packages <- c("ggplot2", "mgcv") # ggplot2 for plotting, mgcv for splines
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
  cat("   Loaded:", pkg, "\n")
}

## 3. Run analyses and record results ----------------------------------------
cat("\n3. RUNNING ANALYSES AND RECORDING RESULTS...\n")

# Initialize master results object
all_results <- list(
  test_info = list(
    date = Sys.time(),
    data_source = "Crescent Pond puppyfish",
    n_observations = nrow(test_data),
    fitness_variables = c(binary = fitness_binary,
                          continuous = fitness_continuous),
    traits_analyzed = morphological_traits
  ),
  data_preparation = list(),
  selection_analysis = list(),
  fitness_surfaces = list(),
  validation_summary = list()
)

## 3.1 Data preparation -------------------------------------------------------
cat("   3.1 DATA PREPARATION...\n")

if (exists("prepare_selection_data")) {
  
  # Binary fitness preparation
  prepared_binary <- tryCatch({
    result <- prepare_selection_data(
      data = test_data,
      fitness_col = fitness_binary,
      trait_cols = morphological_traits,
      standardize = TRUE
    )
    all_results$data_preparation$binary <- list(
      status = "success",
      dimensions = dim(result)
    )
    result
  }, error = function(e) {
    all_results$data_preparation$binary <- list(
      status = "failed",
      error = e$message
    )
    NULL
  })
  
  # Continuous fitness preparation
  prepared_continuous <- tryCatch({
    result <- prepare_selection_data(
      data = test_data,
      fitness_col = fitness_continuous,
      trait_cols = morphological_traits,
      standardize = TRUE
    )
    all_results$data_preparation$continuous <- list(
      status = "success",
      dimensions = dim(result)
    )
    result
  }, error = function(e) {
    all_results$data_preparation$continuous <- list(
      status = "failed",
      error = e$message
    )
    NULL
  })
  
  cat("   Data preparation completed\n")
}

## 3.2 Selection coefficient analysis ----------------------------------------
cat("   3.2 SELECTION COEFFICIENT ANALYSIS...\n")

if (exists("selection_coefficients")) {
  
  # Binary fitness model
  if (!is.null(prepared_binary)) {
    selection_binary <- tryCatch({
      result <- selection_coefficients(
        data = prepared_binary,
        fitness_col = fitness_binary,
        trait_cols = morphological_traits,
        fitness_type = "binary",
        standardize = FALSE
      )
      all_results$selection_analysis$binary <- list(
        status = "success",
        results = result
      )
      result
    }, error = function(e) {
      all_results$selection_analysis$binary <- list(
        status = "failed",
        error = e$message
      )
      NULL
    })
  }
  
  # Continuous fitness model
  if (!is.null(prepared_continuous)) {
    selection_continuous <- tryCatch({
      result <- selection_coefficients(
        data = prepared_continuous,
        fitness_col = fitness_continuous,
        trait_cols = morphological_traits,
        fitness_type = "continuous",
        standardize = FALSE
      )
      all_results$selection_analysis$continuous <- list(
        status = "success",
        results = result
      )
      result
    }, error = function(e) {
      all_results$selection_analysis$continuous <- list(
        status = "failed",
        error = e$message
      )
      NULL
    })
  }
  
  cat("   Selection coefficient analysis completed\n")
}

## 3.3 Fitness surface estimation --------------------------------------------
cat("   3.3 FITNESS SURFACE ESTIMATION...\n")

if (exists("univariate_spline") && !is.null(prepared_binary)) {
  
  # Univariate fitness surfaces
  univariate_results <- list()
  for (trait in morphological_traits[1:3]) {
    spline_result <- tryCatch({
      univariate_spline(
        data = prepared_binary,
        fitness_col = fitness_binary,
        trait_col = trait,
        fitness_type = "binary",
        k = 6
      )
    }, error = function(e) NULL)
    
    if (!is.null(spline_result)) {
      univariate_results[[trait]] <- spline_result
      
      if (exists("univariate_surface")) {
        plot_file <- file.path(
          results_dir,
          paste0("univariate_surface_", trait, ".png")
        )
        png(plot_file, width = 800, height = 600)
        print(
          univariate_surface(spline_result, trait_col = trait) +
            labs(title = paste("Fitness Surface:", trait))
        )
        dev.off()
      }
    }
  }
  
  all_results$fitness_surfaces$univariate <- list(
    traits_tested = names(univariate_results),
    n_surfaces = length(univariate_results)
  )
  
  cat("   Univariate fitness surfaces generated (",
      length(univariate_results), " traits)\n", sep = "")
}

# Bivariate (correlational) fitness surface
if (exists("correlational_tps") && !is.null(prepared_binary)) {
  trait_pair <- morphological_traits[1:2]
  
  tps_result <- tryCatch({
    result <- correlational_tps(
      data = prepared_binary,
      fitness_col = fitness_binary,
      trait_cols = trait_pair,
      method = "auto",
      grid_n = 20
    )
    
    if (exists("correlation_surface")) {
      plot_file <- file.path(results_dir, "correlation_surface.png")
      png(plot_file, width = 1000, height = 800)
      print(
        correlation_surface(result, trait_cols = trait_pair) +
          labs(title = "Bivariate Fitness Surface")
      )
      dev.off()
    }
    
    all_results$fitness_surfaces$correlational <- list(
      status = "success",
      traits = trait_pair,
      grid_dimensions = dim(result$grid)
    )
    result
  }, error = function(e) {
    all_results$fitness_surfaces$correlational <- list(
      status = "failed",
      error = e$message
    )
    NULL
  })
  
  cat("   Correlational fitness surface completed\n")
}
