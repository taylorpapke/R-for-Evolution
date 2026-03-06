
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
crescent_path <- NULL

possible_dirs <- c(
  file.path("..", "data"),
  file.path("..", "test_data"),
  "data",
  "test_data"
)

for (d in possible_dirs) {
  if (dir.exists(d)) {
    candidates <- list.files(d, pattern = "Crescent.*Pond.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
    if (length(candidates) > 0) {
      crescent_path <- candidates[1]
      break
    }
  }
}

if (is.null(crescent_path) || !file.exists(crescent_path)) stop("Data file 'Crescent Pond...' not found. Searched: ", paste(possible_dirs, collapse = ", "))
cat("   Using data file:", crescent_path, "\n")
crescent_data <- read.csv(crescent_path)

# Define fitness variables and trait set
fitness_binary     <- "survival"
fitness_continuous <- "ln.growth"
morphological_traits <- c("jaw", "eye", "body", "nasal", "mouth", "SL")

# Filter for high density to match Martin (2016)
if ("density" %in% names(crescent_data)) {
  cat("   Filtering for High Density (H) enclosure...\n")
  crescent_data <- crescent_data[crescent_data$density == "H", ]
}

# Assemble analysis dataset
test_data <- crescent_data[, c(fitness_binary, fitness_continuous, morphological_traits)]
test_data <- test_data[complete.cases(test_data), ]

cat("   Test dataset contains", nrow(test_data), "complete observations\n")

# Load core analysis functions
function_files <- c(
  "1_prepare_selection_data.R", "2_linear_selection_analysis.R",
  "3_nonlinear_selection_analysis.R", "5_adaptive_landscape_analysis.R",
  "univariate_spline.R", "plot_univariate_fitness.R",
  "plot_adaptive_landscape.R", "plot_correlated_fitness.R"
)

func_search_dirs <- c(
  file.path("..", "functions"), file.path("..", "plotting"), file.path("..", "scripts"),
  file.path("R", "functions"), file.path("R", "plotting"), file.path("R", "scripts"),
  "functions", "plotting", "scripts",
  ".."
)

functions_loaded <- 0

for (func_file in function_files) {
  found <- FALSE
  for (d in func_search_dirs) {
    file_path <- file.path(d, func_file)
    if (file.exists(file_path)) {
      source(file_path)
      cat("   Loaded:", func_file, "\n")
      found <- TRUE
      functions_loaded <- functions_loaded + 1
      break
    }
  }
  if (!found) cat("   Warning: Could not find", func_file, "\n")
}

if (functions_loaded == 0) stop("No analysis functions were loaded. Searched: ", paste(func_search_dirs, collapse = ", "))

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
  
  if (!is.null(prepared_binary)) {
    cat("\n   --- Data Summary (Binary Fitness) ---\n")
    print(table(Survival = prepared_binary[[fitness_binary]]))
    cat("   -------------------------------------\n")
  }

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

if (exists("analyze_linear_selection")) {
  
  # Binary fitness model
  if (!is.null(prepared_binary)) {
    selection_binary <- tryCatch({
      result <- analyze_linear_selection(
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

    if (!is.null(selection_binary)) {
      cat("\n   --- Linear Selection Gradients (Binary) ---\n")
      print(selection_binary)
      cat("   -------------------------------------------\n")
    }
  }
  
  # Continuous fitness model
  if (!is.null(prepared_continuous)) {
    selection_continuous <- tryCatch({
      result <- analyze_linear_selection(
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

    if (!is.null(selection_continuous)) {
      cat("\n   --- Linear Selection Gradients (Continuous) ---\n")
      print(selection_continuous)
      cat("   -----------------------------------------------\n")
    }
  }
  
  cat("   Selection coefficient analysis completed\n")
} else {
  cat("   Skipping selection analysis: 'analyze_linear_selection' function not found.\n")
}

## 3.3 Fitness surface estimation --------------------------------------------
cat("   3.3 FITNESS SURFACE ESTIMATION...\n")

if (exists("univariate_spline") && !is.null(prepared_binary)) {
  
  # Univariate fitness surfaces
  univariate_results <- list()
  for (trait in morphological_traits) {
    spline_result <- tryCatch({
      # Use mgcv::gam directly with Thin Plate Splines (bs='tp') and k=20
      # This matches the methodology used for the bivariate analysis to capture complexity
      fml <- as.formula(paste(fitness_binary, "~ s(", trait, ", bs='tp', k=20)"))
      model <- mgcv::gam(fml, data = prepared_binary, family = binomial)
      
      # 1. Create prediction data
      pred_data <- data.frame(seq(min(prepared_binary[[trait]]), 
                                  max(prepared_binary[[trait]]), 
                                  length.out = 100))
      names(pred_data) <- trait
      
      # 2. Get predictions with standard errors
      preds <- predict(model, newdata = pred_data, se.fit = TRUE)
      
      # 3. Build the specific grid structure the plot function wants
      grid_to_plot <- data.frame(
        trait = as.vector(pred_data[[trait]]),
        fit   = as.vector(plogis(preds$fit)),
        lwr   = as.vector(plogis(preds$fit - (1.96 * preds$se.fit))),
        upr   = as.vector(plogis(preds$fit + (1.96 * preds$se.fit)))
      )
      colnames(grid_to_plot) <- c("trait", "fit", "lwr", "upr")
      
      # Add the specific trait column as a duplicate of 'trait' to satisfy potential aes() mappings
      grid_to_plot[[trait]] <- grid_to_plot$trait
      
      # Keep a copy of the actual trait name if the function needs it elsewhere
      attr(grid_to_plot, "trait_name") <- trait 
      
      res <- list(model = model, grid = grid_to_plot)
      class(res) <- "univariate_fitness"
      res
    }, error = function(e) NULL)
    
    if (!is.null(spline_result)) {
      univariate_results[[trait]] <- spline_result

      cat("\n   --- Univariate Model Summary:", trait, "---\n")
      print(summary(spline_result$model))
      cat("   --------------------------------------\n")
      
      if (exists("plot_univariate_fitness")) {
        tryCatch({
          plot_file <- file.path(
            results_dir,
            paste0("univariate_surface_", trait, ".png")
          )
          png(plot_file, width = 800, height = 600)
          print(
            plot_univariate_fitness(spline_result, trait_col = trait) +
              labs(title = paste("Fitness Surface:", trait))
          )
          dev.off()
        }, error = function(e) {
          cat("   Warning: Plotting failed for", trait, ":", e$message, "\n")
          if (length(dev.list()) > 0) dev.off()
        })
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
tps_result <- NULL

if (!is.null(prepared_binary)) {
  # Use Jaw and Body for the most complex interaction (per Martin 2016)
  trait_pair <- c("jaw", "body")
  
  # 1. Run Analysis
  tps_result <- tryCatch({
    # Use mgcv::gam with pure Thin Plate Regression Splines (bs='tp')
    # This more closely matches the 'fields' package methodology used in Martin (2016)
    fml <- as.formula(paste(fitness_binary, "~ s(", trait_pair[1], ",", trait_pair[2], ", bs='tp', k=20)"))
    
    model <- mgcv::gam(fml, data = prepared_binary, family = binomial)
    
    list(model = model)
  }, error = function(e) {
    cat("   Error in bivariate analysis:", e$message, "\n")
    NULL
  })

  # 2. Print Summary
  if (!is.null(tps_result)) {
    cat("\n   --- Bivariate TPS Model Summary ---\n")
    print(summary(tps_result$model))
    cat("   -----------------------------------\n")
    
    all_results$fitness_surfaces$correlational <- list(
      status = "success",
      traits = trait_pair,
      grid_dimensions = NULL
    )
  }

  # 3. Plotting
  # Note: plot_correlated_fitness requires a pre-calculated grid which we haven't generated here.
  # We will rely on adaptive_landscape (Section 3.4) for visualization.
  
  cat("   Correlational fitness surface completed\n")
}

# 1. Load the landscape analysis function
# (Make sure "5_adaptive_landscape_analysis.R" is in your function_files list)

if (exists("adaptive_landscape") && !is.null(tps_result)) {
  cat("   3.4 GENERATING ADAPTIVE LANDSCAPE (Population Level)...\n")
  
  # This creates the mean fitness surface shown in Figure 1
  landscape_result <- adaptive_landscape(
    data = prepared_binary,
    fitness_model = tps_result$model,
    trait_cols = trait_pair
  )
  
  if (exists("plot_adaptive_landscape")) {
    landscape_file <- file.path(results_dir, "adaptive_landscape_fig1.png")
    png(landscape_file, width = 1000, height = 800)
    print(
      plot_adaptive_landscape(landscape_result, trait_cols = trait_pair) +
        labs(title = "Adaptive Landscape (Population Mean Fitness)")
    )
    dev.off()
    cat("   Adaptive landscape plot saved to:", landscape_file, "\n")
  }
}
