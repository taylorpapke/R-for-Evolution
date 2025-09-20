# test_all_functions_fixed.R
# 修复了字符串重复错误的版本

cat("=== COMPREHENSIVE FUNCTION TESTING ===\n\n")

# 1. 加载所有函数文件
cat("1. Loading all function files...\n")
function_files <- c(
  "prepare_selection_data.R", "analyze_linear_selection.R", "analyze_nonlinear_selection.R",
  "extract_results.R", "selection_coefficients.R", "detect_family.R", "selection_differential.R", 
  "univariate_spline.R", "univariate_surface.R", "correlational_tps.R", "correlation_surface.R", 
  "bootstrap_selection.R"
)

for (file in function_files) {
  if (file.exists(file)) {
    source(file)
    cat("   Loaded:", file, "\n")
  } else {
    cat("   ❗ File not found:", file, "\n")
  }
}

# 2. 创建测试数据
cat("\n2. Creating test data...\n")
set.seed(42)

df_binary <- data.frame(
  fitness = rbinom(50, 1, 0.6),
  size = rnorm(50, 5, 1),
  speed = rnorm(50, 10, 2),
  color = rnorm(50, 3, 0.5)
)

df_continuous <- data.frame(
  fitness = rpois(50, 3) + 1,
  size = rnorm(50, 5, 1),
  speed = rnorm(50, 10, 2),
  color = rnorm(50, 3, 0.5)
)

cat("   Created df_binary:", dim(df_binary), "\n")
cat("   Created df_continuous:", dim(df_continuous), "\n")

# 3. 测试之前的核心函数
cat("\n3. Testing core functions...\n")

# 3.1 测试 selection_coefficients
cat("   Testing selection_coefficients...\n")
test1 <- selection_coefficients(
  data = df_binary,
  fitness_col = "fitness",
  trait_cols = c("size", "speed"),
  fitness_type = "binary",
  standardize = TRUE,
  use_relative_for_fit = TRUE
)

test2 <- selection_coefficients(
  data = df_continuous,
  fitness_col = "fitness",
  trait_cols = c("size", "speed"),
  fitness_type = "continuous",
  standardize = TRUE,
  use_relative_for_fit = TRUE
)

cat("   ✓ Binary result:", nrow(test1), "rows\n")
cat("   ✓ Continuous result:", nrow(test2), "rows\n")

# 3.2 测试 analyze_linear_selection
cat("   Testing analyze_linear_selection...\n")
linear_result <- analyze_linear_selection(
  data = df_continuous,
  fitness_col = "fitness",
  trait_cols = c("size", "speed"),
  fitness_type = "continuous"
)
cat("   ✓ Linear analysis complete\n")

# 3.3 测试 extract_results
cat("   Testing extract_results...\n")
linear_coefs <- extract_linear_coefficients(c("size", "speed"), linear_result)
cat("   ✓ Extracted", nrow(linear_coefs), "linear coefficients\n")

# 4. 测试新函数
cat("\n4. Testing new functions...\n")

# 4.1 测试 detect_family
cat("   Testing detect_family...\n")
family_binary <- detect_family(df_binary$fitness)
family_continuous <- detect_family(df_continuous$fitness)
cat("   ✓ Binary family:", family_binary$type, "\n")
cat("   ✓ Continuous family:", family_continuous$type, "\n")

# 4.2 测试 selection_differential
cat("   Testing selection_differential...\n")
sel_diff <- selection_differential(
  data = df_continuous,
  fitness_col = "fitness",
  trait_col = "size",
  assume_standardized = FALSE,
  use_relative = TRUE
)
cat("   ✓ Selection differential:", round(sel_diff, 4), "\n")

# 4.3 测试 univariate_spline 和 univariate_surface
cat("   Testing univariate_spline and univariate_surface...\n")
spline_result <- univariate_spline(
  data = df_continuous,
  fitness_col = "fitness",
  trait_col = "size",
  fitness_type = "continuous",
  k = 8
)

surface_plot <- univariate_surface(
  uni = spline_result,
  trait_col = "size"
)
cat("   ✓ Spline model:", class(spline_result$model), "\n")
cat("   ✓ Surface plot:", class(surface_plot), "\n")

# 4.4 测试 correlational_tps 和 correlation_surface
cat("   Testing correlational_tps and correlation_surface...\n")
tps_result <- correlational_tps(
  data = df_continuous,
  fitness_col = "fitness",
  trait_cols = c("size", "speed"),
  use_relative = TRUE,
  grid_n = 30
)

correlation_plot <- correlation_surface(
  tps = tps_result,
  trait_cols = c("size", "speed")
)
cat("   ✓ TPS model:", class(tps_result$model), "\n")
cat("   ✓ Correlation plot:", class(correlation_plot), "\n")

# 4.5 测试 bootstrap_selection
cat("   Testing bootstrap_selection...\n")
bootstrap_results <- bootstrap_selection(
  data = df_continuous,
  fitness_col = "fitness",
  trait_cols = c("size", "speed"),
  fitness_type = "continuous",
  B = 25,
  seed = 42
)
cat("   ✓ Bootstrap draws:", nrow(bootstrap_results$draws), "\n")
cat("   ✓ Bootstrap CI:", nrow(bootstrap_results$ci), "confidence intervals\n")

# 5. 验证所有函数输出
cat("\n5. Validating all function outputs...\n")

validate_output <- function(name, object, expected_type) {
  if (is.null(object)) {
    return(paste("❌", name, "is NULL"))
  }
  
  if (expected_type == "data.frame" && !is.data.frame(object)) {
    return(paste("❌", name, "is not a data.frame"))
  }
  
  if (expected_type == "list" && !is.list(object)) {
    return(paste("❌", name, "is not a list")) 
  }
  
  if (expected_type == "numeric" && !is.numeric(object)) {
    return(paste("❌", name, "is not numeric"))
  }
  
  if (expected_type == "plot" && !inherits(object, "gg")) {
    return(paste("❌", name, "is not a ggplot"))
  }
  
  return(paste("✅", name, "PASS"))
}

# 验证所有输出
results <- list(
  validate_output("selection_coefficients binary", test1, "data.frame"),
  validate_output("selection_coefficients continuous", test2, "data.frame"),
  validate_output("detect_family binary", family_binary, "list"),
  validate_output("selection_differential", sel_diff, "numeric"),
  validate_output("univariate_spline", spline_result, "list"),
  validate_output("univariate_surface", surface_plot, "plot"),
  validate_output("correlational_tps", tps_result, "list"),
  validate_output("correlation_surface", correlation_plot, "plot"),
  validate_output("bootstrap_selection", bootstrap_results, "list")
)

cat("Validation results:\n")
for (result in results) {
  cat("   ", result, "\n")
}

# 6. 保存测试结果
cat("\n6. Saving test results...\n")

# 保存图形
if (exists("surface_plot")) {
  tryCatch({
    ggsave("test_univariate_surface.png", surface_plot, width = 8, height = 6)
    cat("   ✓ Saved: test_univariate_surface.png\n")
  }, error = function(e) cat("   ❗ Failed to save surface plot:", e$message, "\n"))
}

if (exists("correlation_plot")) {
  tryCatch({
    ggsave("test_correlation_surface.png", correlation_plot, width = 8, height = 6)
    cat("   ✓ Saved: test_correlation_surface.png\n")
  }, error = function(e) cat("   ❗ Failed to save correlation plot:", e$message, "\n"))
}

# 保存数据结果
test_output <- list(
  selection_coefficients_binary = test1,
  selection_coefficients_continuous = test2,
  selection_differential = sel_diff,
  bootstrap_results = bootstrap_results$ci,
  univariate_spline_grid = spline_result$grid,
  tps_grid = tps_result$grid
)

tryCatch({
  saveRDS(test_output, "test_output_results.rds")
  cat("   ✓ Saved: test_output_results.rds\n")
}, error = function(e) cat("   ❗ Failed to save RDS:", e$message, "\n"))

# 7. 完整的可视化测试
cat("\n7. Comprehensive visualization testing...\n")

# 7.1 测试所有性状的单变量图
cat("   Testing univariate surfaces for all traits...\n")

univariate_plots <- list()

# 测试连续数据的所有性状
for (trait in c("size", "speed", "color")) {
  cat("   Testing", trait, "univariate surface...\n")
  
  spline <- univariate_spline(
    data = df_continuous,
    fitness_col = "fitness",
    trait_col = trait,
    fitness_type = "continuous",
    k = 8
  )
  
  plot <- univariate_surface(
    uni = spline,
    trait_col = trait
  )
  
  univariate_plots[[paste("continuous", trait)]] <- plot
  cat("     ✓ Created", trait, "surface plot\n")
}

# 测试二分类数据的单变量图
for (trait in c("size", "speed")) {
  cat("   Testing binary", trait, "univariate surface...\n")
  
  spline <- univariate_spline(
    data = df_binary,
    fitness_col = "fitness",
    trait_col = trait,
    fitness_type = "binary",
    k = 8
  )
  
  plot <- univariate_surface(
    uni = spline,
    trait_col = trait
  )
  
  univariate_plots[[paste("binary", trait)]] <- plot
  cat("     ✓ Created binary", trait, "surface plot\n")
}

# 7.2 测试所有性状组合的相关表面图
cat("   Testing correlation surfaces for all trait pairs...\n")

correlation_plots <- list()

# 测试所有可能的性状组合
trait_combinations <- list(
  c("size", "speed"),
  c("size", "color"), 
  c("speed", "color")
)

for (traits in trait_combinations) {
  cat("   Testing correlation surface for", paste(traits, collapse = " vs "), "...\n")
  
  tps <- correlational_tps(
    data = df_continuous,
    fitness_col = "fitness",
    trait_cols = traits,
    use_relative = TRUE,
    grid_n = 30
  )
  
  plot <- correlation_surface(
    tps = tps,
    trait_cols = traits
  )
  
  plot_name <- paste(traits, collapse = "_")
  correlation_plots[[plot_name]] <- plot
  cat("     ✓ Created", plot_name, "correlation surface\n")
}

# 7.3 测试二分类数据的相关表面图
cat("   Testing binary data correlation surfaces...\n")

for (traits in trait_combinations[1:2]) {  # 只测试前两个组合，避免过多计算
  cat("   Testing binary correlation surface for", paste(traits, collapse = " vs "), "...\n")
  
  tps <- correlational_tps(
    data = df_binary,
    fitness_col = "fitness",
    trait_cols = traits,
    use_relative = FALSE,
    grid_n = 25
  )
  
  plot <- correlation_surface(
    tps = tps,
    trait_cols = traits
  )
  
  plot_name <- paste("binary", paste(traits, collapse = "_"))
  correlation_plots[[plot_name]] <- plot
  cat("     ✓ Created binary", plot_name, "correlation surface\n")
}

# 7.4 保存所有图形
cat("\n8. Saving all visualization outputs...\n")

# 保存单变量图
for (plot_name in names(univariate_plots)) {
  filename <- paste0("univariate_", gsub(" ", "_", plot_name), ".png")
  tryCatch({
    ggsave(filename, univariate_plots[[plot_name]], width = 8, height = 6, dpi = 300)
    cat("   ✓ Saved:", filename, "\n")
  }, error = function(e) {
    cat("   ❗ Failed to save", filename, ":", e$message, "\n")
  })
}

# 保存相关表面图
for (plot_name in names(correlation_plots)) {
  filename <- paste0("correlation_", plot_name, ".png")
  tryCatch({
    ggsave(filename, correlation_plots[[plot_name]], width = 8, height = 6, dpi = 300)
    cat("   ✓ Saved:", filename, "\n")
  }, error = function(e) {
    cat("   ❗ Failed to save", filename, ":", e$message, "\n")
  })
}

# 7.5 验证图形属性
cat("\n9. Validating plot properties...\n")

validate_plot <- function(plot, plot_name) {
  cat("   Validating", plot_name, "...\n")
  
  checks <- c(
    "Is ggplot" = inherits(plot, "gg"),
    "Has layers" = length(plot$layers) > 0,
    "Has x label" = !is.null(plot$labels$x),
    "Has y label" = !is.null(plot$labels$y),
    "Has data" = !is.null(plot$data)
  )
  
  if (all(checks)) {
    cat("     ✅", plot_name, "PASS\n")
    return(TRUE)
  } else {
    cat("     ❌", plot_name, "FAIL - Issues:", paste(names(checks)[!checks], collapse = ", "), "\n")
    return(FALSE)
  }
}

# 验证所有图形
plot_validation <- list()
for (name in names(univariate_plots)) {
  plot_validation[[name]] <- validate_plot(univariate_plots[[name]], name)
}

for (name in names(correlation_plots)) {
  plot_validation[[name]] <- validate_plot(correlation_plots[[name]], name)
}

# 7.6 图形测试总结
cat("\n10. Visualization test summary...\n")

total_plots <- length(univariate_plots) + length(correlation_plots)
successful_plots <- sum(unlist(plot_validation))

cat("   Total plots created:", total_plots, "\n")
cat("   Successful plots:", successful_plots, "\n")
cat("   Plot success rate:", round(successful_plots/total_plots * 100, 1), "%\n")

cat("   Univariate plots:", length(univariate_plots), "\n")
cat("   Correlation plots:", length(correlation_plots), "\n")

# 更新最终总结
cat("\n11. Updated final summary...\n")

cat("   Functions tested:", length(function_files), "\n")
cat("   Tests passed:", passed, "/", total, "\n")
cat("   Plots created:", total_plots, "\n")
cat("   Successful plots:", successful_plots, "\n")
cat("   Output files:", 3 + total_plots, "\n")  # 3个数据文件 + 所有图形文件

if (passed == total && successful_plots == total_plots) {
  cat("\n", strrep("=", 50), "\n")
  cat("🎉 ALL TESTS AND VISUALIZATIONS PASSED SUCCESSFULLY!\n")
  cat(strrep("=", 50), "\n")
} else {
  cat("\n", strrep("=", 50), "\n")
  cat("⚠️  SOME ISSUES DETECTED:\n")
  if (passed < total) cat("   -", total - passed, "function tests failed\n")
  if (successful_plots < total_plots) cat("   -", total_plots - successful_plots, "plots failed validation\n")
  cat(strrep("=", 50), "\n")
}