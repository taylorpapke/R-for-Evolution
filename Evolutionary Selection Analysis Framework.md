# Evolutionary Selection Analysis Framework

## 1. Overview

This project implements a comprehensive quantitative genetic framework for analyzing natural selection on phenotypic traits. The framework follows the Lande & Arnold (1983) approach while extending it with modern nonparametric methods for visualizing adaptive landscapes.

**Key capabilities:**

| Analysis Type | Description |
|---------------|-------------|
| Selection differentials ($\mathbf{S}$) | Total selection acting on traits |
| Linear selection gradients ($\boldsymbol{\beta}$) | Direct directional selection |
| Nonlinear selection gradients ($\boldsymbol{\gamma}$) | Stabilizing/disruptive and correlational selection |
| Univariate fitness surfaces | Flexible spline-based fitness functions |
| Multivariate adaptive landscapes | Thin-plate spline (TPS) fitness surfaces |
| Binary fitness handling | Logistic regression for survival/reproduction data |

---

## 2. Mathematical Framework

### 2.1 Relative Fitness

Individual fitness is converted to relative fitness to enable comparisons across populations and generations:

$$w_i = \frac{W_i}{\bar{W}}$$

where:
- $W_i$ = absolute fitness of individual $i$ (e.g., offspring count, survival status)
- $\bar{W} = \frac{1}{n}\sum_{i=1}^n W_i$ = mean population fitness
- $w_i$ = relative fitness of individual $i$

**Why relative fitness?** Selection operates on relative, not absolute, differences in fitness. This standardization ensures that selection gradients are comparable across studies.

### 2.2 Trait Standardization

Traits are standardized to mean 0 and variance 1 before analysis:

$$z_i = \frac{x_i - \bar{x}}{s_x}$$

where:
- $x_i$ = original trait value
- $\bar{x}$ = sample mean
- $s_x$ = sample standard deviation
- $z_i$ = standardized trait value

**Purpose:** Standardization puts all traits on the same scale, allowing direct comparison of selection gradients. A $\beta = 0.5$ means that a one-standard-deviation increase in the trait increases relative fitness by 0.5 units.

### 2.3 Selection Differential

The selection differential ($S$) measures the total selection acting on a trait, including both direct and indirect effects through correlated traits:

$$S = \text{Cov}(z, w) = \mathbb{E}[(z - \bar{z})(w - \bar{w})]$$

For multiple traits, the vector of selection differentials is:

$$\mathbf{S} = \text{Cov}(\mathbf{z}, w) = \begin{bmatrix} \text{Cov}(z_1, w) \\ \text{Cov}(z_2, w) \\ \vdots \\ \text{Cov}(z_m, w) \end{bmatrix}$$

**Interpretation:** $S$ represents the shift in trait mean after one generation of selection. A positive $S$ indicates the trait mean increases.

### 2.4 Linear Selection Gradients

Linear selection gradients ($\boldsymbol{\beta}$) measure the direct directional selection on each trait, controlling for correlations with other traits:

$$w = \alpha + \boldsymbol{\beta}^T \mathbf{z} + \varepsilon$$

In matrix form:

$$\boldsymbol{\beta} = \mathbf{P}^{-1} \mathbf{S}$$

where $\mathbf{P}$ is the phenotypic variance-covariance matrix.

For a model with two traits:

$$w = \alpha + \beta_1 z_1 + \beta_2 z_2 + \varepsilon$$

**Interpretation:** $\beta_i$ is the partial regression coefficient—the change in relative fitness for a one-standard-deviation increase in trait $i$, holding all other traits constant.

### 2.5 Nonlinear Selection Gradients

Quadratic selection gradients ($\boldsymbol{\gamma}$) measure nonlinear selection, including stabilizing/disruptive selection (diagonal elements) and correlational selection (off-diagonal elements):

$$w = \alpha + \boldsymbol{\beta}^T \mathbf{z} + \frac{1}{2} \mathbf{z}^T \boldsymbol{\gamma} \mathbf{z} + \varepsilon$$

For two traits, this expands to:

$$w = \alpha + \beta_1 z_1 + \beta_2 z_2 + \frac{1}{2}\gamma_{11} z_1^2 + \frac{1}{2}\gamma_{22} z_2^2 + \gamma_{12} z_1 z_2 + \varepsilon$$

**Important:** The $\frac{1}{2}$ factor in the quadratic terms follows Lande & Arnold's convention, making $\gamma_{ii}$ the second derivative of the fitness surface. In practice, most implementations (including this one) report $2 \times$ the quadratic coefficient from regression, which equals $\gamma_{ii}$.

**Interpretation:**

| $\gamma$ value | Type | Interpretation |
|----------------|------|----------------|
| $\gamma_{ii} < 0$ | Stabilizing selection | Intermediate trait values have highest fitness |
| $\gamma_{ii} > 0$ | Disruptive selection | Extreme trait values have highest fitness |
| $\gamma_{ij} > 0$ | Positive correlational | Selection favors positive correlation between traits |
| $\gamma_{ij} < 0$ | Negative correlational | Selection favors negative correlation between traits |

### 2.6 Nonparametric Adaptive Landscapes

Traditional quadratic models impose a fixed shape on the fitness function. Nonparametric methods allow the data to reveal unexpected patterns.

**Univariate spline model:**

$$w = \alpha + s(z) + \varepsilon$$

where $s(z)$ is a smooth function estimated from the data using penalized regression splines.

**Multivariate thin-plate spline:**

$$w = \alpha + f(z_1, z_2) + \varepsilon$$

where $f(z_1, z_2)$ is a two-dimensional smooth surface. This creates a visual "adaptive landscape" showing how fitness varies across trait combinations.

---

## 3. Framework Structure

### 3.1 Core Scripts and Their Roles

| Script | Purpose | Key Functions |
|--------|---------|---------------|
| `prepare_selection_data.R` | Data preprocessing | Standardization, relative fitness calculation, NA handling |
| `detect_family.R` | Automatic fitness type detection | Determines binary vs continuous fitness |
| `selection_differential.R` | Selection differential calculation | $S = \text{Cov}(z, w)$ |
| `selection_coefficients.R` | Main wrapper function | Runs complete selection analysis |
| `analyze_linear_selection.R` | Linear gradient estimation | Fits $w \sim \beta z$ models with diagnostics |
| `analyze_nonlinear_selection.R` | Quadratic gradient estimation | Fits models with $z^2$ and interactions |
| `analyze_disruptive_selection.R` | Single-trait quadratic analysis | Tests stabilizing/disruptive selection |
| `correlational_tps.R` | 2D fitness surface estimation | TPS for continuous, GAM for binary |
| `univariate_spline.R` | 1D fitness function | GAM with univariate smooth |
| `correlation_surface.R` | Visualization | Contour plots of fitness surfaces |
| `univariate_surface.R` | Visualization | Line plots with confidence bands |
| `extract_results.R` | Results compilation | Coefficient extraction from model objects |

### 3.2 Function Dependencies

```
                  selection_coefficients()
                           |
                    prepare_selection_data()
                    (standardizes traits,
                     computes relative fitness)
                           +
                      detect_family()
                    (auto-detects binary vs
                     continuous fitness)
                    /         |         \
                   ↓          ↓          ↓
    fitness_type =    fitness_type =    fitness_type = 
       "auto"           "binary"         "continuous"
        (uses               |                 |
      detect_family         ↓                 ↓
       to choose)    analyze_linear_    analyze_linear_
                          selection()        selection()
                           |                 |
                           ↓                 ↓
                    analyze_nonlinear_  analyze_nonlinear_
                         selection()         selection()
                           |                 |
                           ↓                 ↓
                    extract_linear_     extract_linear_
                      coefficients()       coefficients()
                           ↓                 ↓
                    extract_quadratic_  extract_quadratic_
                      coefficients()       coefficients()
                           ↓                 ↓
                    extract_interaction_ extract_interaction_
                      coefficients()        coefficients()
                           \                 /
                            \               /
                             ↓             ↓
                          selection_coefficients()
                                    ↓
                         Returns combined results
                         data frame with all
                         selection gradients
```

---

## 4. Installation and Requirements

### 4.1 Required Packages

```r
# Install required packages
packages <- c("mgcv", "fields", "car", "ggplot2", "tidyr", "dplyr")

install.packages(packages)
```

### 4.2 Package Descriptions

| Package | Purpose |
|---------|---------|
| `mgcv` | Generalized additive models for spline fitting |
| `fields` | Thin-plate splines for bivariate surfaces |
| `car` | Type III ANOVA and variance inflation factors |
| `ggplot2` | Visualization of fitness surfaces |
| `tidyr`/`dplyr` | Data manipulation (recommended) |

---

## 5. Detailed Function Reference

### 5.1 `prepare_selection_data()`

**Purpose:** Preprocesses raw data for selection analysis.

```r
prepare_selection_data(
  data,                    # Data frame
  fitness_col,             # Name of fitness column (character)
  trait_cols,              # Vector of trait column names
  standardize = TRUE,       # Standardize traits to mean 0, SD 1
  add_relative = TRUE,      # Add relative fitness column
  na_action = c("warn", "drop", "none"),  # How to handle NAs
  name_relative = "relative_fitness"       # Name for relative fitness column
)
```

**Returns:** Modified data frame with:
- Standardized traits (if `standardize = TRUE`)
- Relative fitness column (if `add_relative = TRUE`)
- Potentially dropped NA rows (if `na_action = "drop"`)

**Example:**

```r
# Original data
raw_data <- data.frame(
  survival = c(1, 0, 1, 1, NA, 1),
  size = c(10.2, 8.5, NA, 11.3, 9.8, 10.7),
  color = c(2.1, 3.2, 1.8, 2.5, 3.5, 2.0)
)

# Prepare data
prepared <- prepare_selection_data(
  data = raw_data,
  fitness_col = "survival",
  trait_cols = c("size", "color"),
  na_action = "warn"
)
# Warning: Detected 2 row(s) with NA in fitness/traits: survival, size, color
```

### 5.2 `detect_family()`

**Purpose:** Automatically determines whether fitness is binary or continuous.

```r
detect_family(y)  # y = fitness vector
```

**Logic:**
- **Binary** if: ≤ 2 unique values, all in {0,1}, and n ≥ 10
- **Continuous** otherwise

**Returns:** List with:
- `$type`: "binary" or "continuous"
- `$family`: Appropriate GLM family object

**Example:**
```r
# Binary fitness
fitness_binary <- c(1, 0, 1, 1, 0, 1)
detect_family(fitness_binary)
# $type: "binary"
# $family: binomial(logit)

# Continuous fitness
fitness_cont <- c(2.3, 4.1, 3.2, 5.6, 1.8)
detect_family(fitness_cont)
# $type: "continuous"
# $family: gaussian
```

### 5.3 `selection_differential()`

**Purpose:** Calculates the selection differential ($S$) for a single trait.

```r
selection_differential(
  data,
  fitness_col,
  trait_col,
  assume_standardized = TRUE,   # Whether trait is already standardized
  use_relative = TRUE            # Use relative fitness
)
```

**Formula:** $
S = \text{Cov}(z, w)
$

When traits are standardized ($\bar{z} = 0$), this simplifies to:
$$S = \mathbb{E}(z \times w) = \frac{1}{n}\sum_{i=1}^n z_i w_i$$

**Example:**
```r
# Calculate selection differential for size
S_size <- selection_differential(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_col = "size",
  assume_standardized = TRUE
)
# Returns: 0.234 (positive directional selection)
```

### 5.4 `analyze_linear_selection()`

**Purpose:** Estimates linear selection gradients ($\boldsymbol{\beta}$) with comprehensive diagnostics.

```r
analyze_linear_selection(
  data,
  fitness_col,
  trait_cols,
  fitness_type        # "binary" or "continuous"
)
```

**Features:**
- Sample size checks (warning if n < 10)
- Missing data handling
- Convergence checks for GLM
- Type III ANOVA (if `car` package available)
- Separation detection for binary models

**Returns:** List containing:
- `$model`: Fitted model object (lm or glm)
- `$summary`: Model summary
- `$anova`: Type III ANOVA table (if available)

**Example:**
```r
linear_results <- analyze_linear_selection(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_cols = c("size", "color"),
  fitness_type = "continuous"
)

# View coefficients
summary(linear_results$model)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  1.02345    0.08765  11.674  < 2e-16 ***
# size         0.15678    0.04532   3.459  0.00123 **
# color       -0.08934    0.05123  -1.744  0.08945 .
```

### 5.5 `analyze_nonlinear_selection()`

**Purpose:** Estimates quadratic ($\gamma_{ii}$) and correlational ($\gamma_{ij}$) selection gradients.

```r
analyze_nonlinear_selection(
  data,
  fitness_col,
  trait_cols,
  fitness_type
)
```

**Model specification:** Automatically constructs:
- Linear terms: $z_i$
- Quadratic terms: $z_i^2$
- Interaction terms: $z_i \times z_j$ (for i < j)

**Diagnostics:**
- Sample size checks (n < 20 warning)
- Variance inflation factors (VIF > 10 warning)
- Convergence monitoring for GLM

**Example:**
```r
nonlinear_results <- analyze_nonlinear_selection(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_cols = c("size", "color"),
  fitness_type = "continuous"
)

# Extract quadratic term for size (multiply by 2 for gamma)
coef(summary(nonlinear_results$model))["I(size^2)", ]
#              Estimate Std. Error   t value     Pr(>|t|)
# I(size^2)   -0.078234   0.034521 -2.266327  0.0289345
# gamma_size = 2 * (-0.078234) = -0.156 (stabilizing selection)
```

### 5.6 `analyze_disruptive_selection()`

**Purpose:** Specialized function for detecting disruptive or stabilizing selection on a single trait.

```r
analyze_disruptive_selection(
  data,
  fitness_col,
  trait_col,
  fitness_type = c("binary", "continuous"),
  standardize = TRUE
)
```

**Model:** $w \sim z + z^2$

**Returns:** Data frame with:
- Linear term ($\beta$): Directional selection
- Quadratic term ($\gamma = 2 \times b_2$): Disruptive (positive) or stabilizing (negative)

**Example:**
```r
disruptive_test <- analyze_disruptive_selection(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_col = "size",
  fitness_type = "continuous"
)

print(disruptive_test)
#     Term Type       Beta_Coefficient Standard_Error P_Value Variance
# 1   size Linear              0.1567         0.0453  0.0012  0.00205
# 2 size²  Quadratic           -0.1565         0.0690  0.0289  0.00476
# Interpretation: Negative quadratic = stabilizing selection
```

### 5.7 `univariate_spline()`

**Purpose:** Fits a flexible nonparametric fitness function using GAMs.

```r
univariate_spline(
  data,
  fitness_col,
  trait_col,
  fitness_type = c("binary", "continuous"),
  k = 10          # Basis dimension (smoothness parameter)
)
```

**Technical details:**
- Uses penalized regression splines with REML smoothing parameter selection
- For binary fitness: binomial GAM with logit link
- For continuous: Gaussian GAM
- Returns predictions with 95% confidence intervals

**Example:**
```r
spline_fit <- univariate_spline(
  data = prepared,
  fitness_col = "survival",
  trait_col = "size",
  fitness_type = "binary",
  k = 8
)

# Access prediction grid
head(spline_fit$grid)
#       size       fit       lwr       upr
# 1 7.890123 0.1234567 0.0891234 0.1678901
# 2 8.123456 0.1456789 0.1123456 0.1890123
# ...
```

### 5.8 `correlational_tps()`

**Purpose:** Estimates 2D fitness surfaces using thin-plate splines (continuous) or GAMs (binary).

```r
correlational_tps(
  data,
  fitness_col,
  trait_cols,        # Exactly 2 traits
  grid_n = 60,       # Grid resolution (grid_n × grid_n points)
  method = "auto",   # "auto", "gam", or "tps"
  scale_traits = TRUE,
  k = 30             # Basis dimension for GAM method
)
```

**Method selection:**
- `"tps"`: Uses `fields::Tps()` for continuous fitness (spline-based)
- `"gam"`: Uses `mgcv::gam()` with tensor product smooth (works for both)
- `"auto"`: Chooses based on fitness type (GAM for binary, TPS for continuous)

**Returns:** List containing:
- `$model`: Fitted model object
- `$grid`: Data frame with trait values and predicted fitness
- `$method`: Method used
- `$data_type`: "binary" or "continuous"

**Example:**
```r
surface <- correlational_tps(
  data = prepared,
  fitness_col = "relative_fitness",
  trait_cols = c("size", "color"),
  method = "gam",
  grid_n = 50
)

# Find optimum (maximum fitness)
optimum <- surface$grid[which.max(surface$grid$.fit), ]
print(optimum)
#      size   color      .fit
# 345 12.34   2.56     1.2345
```

### 5.9 `correlation_surface()` and `correlation_surface_enhanced()`

**Purpose:** Visualize 2D fitness surfaces with contour plots.

```r
correlation_surface(
  tps,              # Output from correlational_tps()
  trait_cols,       # Two trait names
  bins = 12         # Number of contour bins
)

correlation_surface_enhanced(
  tps,
  trait_cols,
  original_data = NULL,   # Original data points
  fitness_col = NULL,      # For coloring points by outcome
  bins = 12
)
```

**Example:**
```r
# Basic contour plot
p1 <- correlation_surface(surface, c("size", "color"))

# Enhanced plot with data points and optimum
p2 <- correlation_surface_enhanced(
  surface,
  c("size", "color"),
  original_data = prepared,
  fitness_col = "survival"
)

# Display the plot
print(p2)
```

### 5.10 `univariate_surface()`

**Purpose:** Visualize univariate fitness functions with confidence bands.

```r
univariate_surface(
  uni,              # Output from univariate_spline()
  trait_col,
  title = NULL
)
```

**Example:**
```r
p <- univariate_surface(
  spline_fit,
  "size",
  title = "Fitness Function for Body Size"
)
print(p)
```

### 5.11 `selection_coefficients()`

**Purpose:** Main wrapper function that runs a complete selection analysis.

```r
selection_coefficients(
  data,
  fitness_col,
  trait_cols,
  fitness_type = c("auto", "binary", "continuous"),
  standardize = TRUE,
  use_relative_for_fit = TRUE
)
```

**Workflow:**
1. Prepares data (standardization, relative fitness)
2. Detects fitness type (if `auto`)
3. Runs linear selection analysis
4. Runs nonlinear selection analysis
5. Extracts all coefficients into a single table

**Returns:** Data frame with columns:
- `Term`: Coefficient name (e.g., "size", "size²", "size×color")
- `Type`: "Linear", "Quadratic", or "Correlational"
- `Beta_Coefficient`: Estimated selection gradient
- `Standard_Error`: Standard error of estimate
- `P_Value`: Statistical significance
- `Variance`: Square of standard error

**Example:**
```r
# Complete analysis
results <- selection_coefficients(
  data = my_data,
  fitness_col = "survival",
  trait_cols = c("size", "color", "mass")
)

# View results
print(results)
#        Term          Type Beta_Coefficient Standard_Error P_Value Variance
# 1      size        Linear           0.1567         0.0453  0.0012  0.00205
# 2     color        Linear          -0.0893         0.0512  0.0895  0.00262
# 3      mass        Linear           0.2134         0.0489  0.0001  0.00239
# 4    size²     Quadratic          -0.1565         0.0690  0.0289  0.00476
# 5   color²     Quadratic           0.0876         0.0589  0.1456  0.00347
# 6    mass²     Quadratic          -0.2345         0.0721  0.0023  0.00520
# 7 size×color Correlational           0.0456         0.0389  0.2456  0.00151
# 8  size×mass Correlational           0.1234         0.0412  0.0045  0.00170
# 9 color×mass Correlational          -0.0678         0.0398  0.0923  0.00158
```

### 5.12 Extraction Functions

Internal functions used by `selection_coefficients()`:

- `extract_linear_coefficients(trait_cols, results)`
- `extract_quadratic_coefficients(trait_cols, results)`
- `extract_interaction_coefficients(trait_cols, results)`

These functions safely extract coefficients from model objects, handling different column names and missing terms.

---

## 6. Complete Analysis Workflow

### 6.1 Step-by-Step Pipeline

```r
# ==================================================
# COMPLETE SELECTION ANALYSIS WORKFLOW
# ==================================================

# 1. LOAD REQUIRED PACKAGES
# --------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(mgcv)
library(fields)
library(car)

# 2. SOURCE ALL FUNCTIONS
# --------------------------------------------------
source("prepare_selection_data.R")
source("detect_family.R")
source("selection_differential.R")
source("selection_coefficients.R")
source("analyze_linear_selection.R")
source("analyze_nonlinear_selection.R")
source("analyze_disruptive_selection.R")
source("correlational_tps.R")
source("univariate_spline.R")
source("univariate_surface.R")
source("correlation_surface.R")
source("extract_results.R")

# 3. LOAD AND EXPLORE DATA
# --------------------------------------------------
# Example: Darwin's finch data (body size, beak depth, survival)
data <- read.csv("finch_data.csv")

# Examine structure
str(data)
# 'data.frame': 200 obs. of 4 variables:
#  $ survival: int  1 0 1 1 0 1 1 1 0 1 ...
#  $ size    : num  10.2 8.5 11.3 9.8 7.9 10.7 11.2 9.3 8.1 10.5 ...
#  $ beak    : num  2.1 3.2 1.8 2.5 3.5 2.0 2.3 2.8 3.1 1.9 ...
#  $ mass    : num  5.2 4.1 5.8 4.9 3.8 5.5 5.9 4.7 4.2 5.3 ...

# Check for missing values
colSums(is.na(data))
# survival    size    beak    mass 
#        0       3       1       2 

# 4. PREPARE DATA (STANDARDIZE, RELATIVE FITNESS)
# --------------------------------------------------
prepared <- prepare_selection_data(
  data = data,
  fitness_col = "survival",
  trait_cols = c("size", "beak", "mass"),
  standardize = TRUE,
  add_relative = TRUE,
  na_action = "drop"  # Remove rows with NAs
)

# Check prepared data
head(prepared)
#   survival  size   beak   mass relative_fitness
# 1        1  0.54   0.23   0.67        1.333333
# 2        0 -1.22   1.32  -0.89        0.000000
# 3        1  1.08  -0.56   0.95        1.333333
# ...

# 5. RUN COMPLETE SELECTION ANALYSIS
# --------------------------------------------------
# Get all selection coefficients
results <- selection_coefficients(
  data = prepared,
  fitness_col = "survival",
  trait_cols = c("size", "beak", "mass"),
  fitness_type = "binary"  # We know it's binary
)

# View results
print(results)

# 6. CHECK INDIVIDUAL TRAITS FOR DISRUPTIVE SELECTION
# --------------------------------------------------
# Test each trait separately
disruptive_size <- analyze_disruptive_selection(
  data = prepared,
  fitness_col = "survival",
  trait_col = "size",
  fitness_type = "binary"
)

disruptive_beak <- analyze_disruptive_selection(
  data = prepared,
  fitness_col = "survival",
  trait_col = "beak",
  fitness_type = "binary"
)

# 7. VISUALIZE UNIVARIATE FITNESS FUNCTIONS
# --------------------------------------------------
# Size fitness function
size_spline <- univariate_spline(
  data = prepared,
  fitness_col = "survival",
  trait_col = "size",
  fitness_type = "binary",
  k = 8
)

p_size <- univariate_surface(
  size_spline,
  "size",
  title = "Survival Probability by Body Size"
)

# Beak fitness function
beak_spline <- univariate_spline(
  data = prepared,
  fitness_col = "survival",
  trait_col = "beak",
  fitness_type = "binary",
  k = 8
)

p_beak <- univariate_surface(
  beak_spline,
  "beak",
  title = "Survival Probability by Beak Depth"
)

# Arrange plots side by side
library(patchwork)
p_size + p_beak

# 8. CREATE BIVARIATE FITNESS SURFACES
# --------------------------------------------------
# Size × Beak surface
surface_sb <- correlational_tps(
  data = prepared,
  fitness_col = "survival",
  trait_cols = c("size", "beak"),
  method = "gam",  # Use GAM for binary data
  grid_n = 50,
  k = 25
)

# Visualize with enhanced plot
p_surface <- correlation_surface_enhanced(
  surface_sb,
  c("size", "beak"),
  original_data = prepared,
  fitness_col = "survival",
  bins = 15
)

print(p_surface)

# 9. QUANTIFY OPTIMAL PHENOTYPE
# --------------------------------------------------
# Find the trait combination with highest survival
optimum <- surface_sb$grid[which.max(surface_sb$grid$.fit), ]
cat("Optimal phenotype:\n")
cat("  Size:", optimum$size, "\n")
cat("  Beak:", optimum$beak, "\n")
cat("  Predicted survival:", optimum$.fit, "\n")

# 10. EXPORT RESULTS
# --------------------------------------------------
# Save coefficients
write.csv(results, "selection_coefficients.csv", row.names = FALSE)

# Save model objects
saveRDS(list(
  linear = linear_model,
  nonlinear = nonlinear_model,
  surfaces = list(size_beak = surface_sb)
), "selection_models.rds")

# Generate report
sink("selection_report.txt")
cat("SELECTION ANALYSIS REPORT\n")
cat("========================\n\n")
cat("Data:", nrow(prepared), "individuals\n")
cat("Traits:", paste(c("size", "beak", "mass"), collapse=", "), "\n\n")

cat("LINEAR SELECTION GRADIENTS\n")
print(results[results$Type == "Linear", ])

cat("\n\nNONLINEAR SELECTION\n")
print(results[results$Type != "Linear", ])

cat("\n\nINTERPRETATION\n")
for(i in 1:nrow(results)) {
  term <- results$Term[i]
  type <- results$Type[i]
  beta <- results$Beta_Coefficient[i]
  p <- results$P_Value[i]
  
  if(type == "Quadratic") {
    if(beta < 0 && p < 0.05) cat("  - ", term, ": Significant stabilizing selection\n")
    if(beta > 0 && p < 0.05) cat("  - ", term, ": Significant disruptive selection\n")
  }
  
  if(type == "Correlational") {
    if(p < 0.05) cat("  - ", term, ": Significant correlational selection\n")
  }
}
sink()
```

### 6.2 Binary Fitness Example with Interpretation

```r
# Example: Survival analysis in a plant species
plant_data <- data.frame(
  survival = c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1),
  height = c(25.3, 18.7, 32.1, 15.4, 28.9, 31.2, 16.8, 29.5, 14.2, 27.8,
             30.1, 26.4, 17.9, 33.2, 19.5, 28.7, 31.8, 15.9, 29.3, 27.1),
  width = c(12.1, 8.3, 15.7, 7.2, 13.5, 14.8, 7.9, 13.9, 6.8, 13.2,
            14.3, 12.8, 8.7, 16.1, 9.2, 13.8, 15.2, 7.5, 14.1, 12.9),
  flowers = c(8, 3, 12, 2, 9, 11, 4, 10, 2, 8,
              10, 7, 4, 13, 5, 9, 11, 3, 10, 8)
)

# Run analysis
plant_results <- selection_coefficients(
  data = plant_data,
  fitness_col = "survival",
  trait_cols = c("height", "width", "flowers"),
  fitness_type = "binary"
)

print(plant_results)
#         Term          Type Beta_Coefficient Standard_Error   P_Value   Variance
# 1     height        Linear           0.8234         0.2345   0.0004   0.05499
# 2      width        Linear           0.4567         0.1987   0.0215   0.03948
# 3    flowers        Linear           1.2345         0.3123   0.0001   0.09753
# 4   height²     Quadratic          -0.3345         0.1456   0.0215   0.02120
# 5    width²     Quadratic           0.1234         0.1234   0.3178   0.01523
# 6  flowers²     Quadratic          -0.5678         0.2345   0.0156   0.05499
# 7 height×width Correlational           0.2345         0.1123   0.0367   0.01261
# 8 height×flowers Correlational           0.4567         0.1456   0.0018   0.02120
# 9 width×flowers Correlational           0.1234         0.0987   0.2115   0.00974

# Interpretation:
# - All traits show positive directional selection (higher values increase survival)
# - Height shows stabilizing selection (γ = -0.33, p = 0.02)
# - Flowers shows stabilizing selection (γ = -0.57, p = 0.02)
# - Positive correlational selection between height and flowers
#   (taller plants with more flowers have especially high survival)
```

### 6.3 Continuous Fitness Example

```r
# Example: Reproductive success in birds
bird_data <- data.frame(
  fitness = c(2, 4, 0, 3, 5, 2, 3, 1, 4, 3, 2, 4, 3, 2, 5, 1, 3, 4, 2, 3),
  wing = c(45.2, 52.3, 38.7, 48.9, 55.1, 44.8, 49.2, 40.3, 51.7, 47.5,
           43.9, 50.2, 46.8, 42.1, 53.4, 39.8, 48.3, 52.9, 44.5, 47.1),
  tail = c(22.1, 25.3, 18.9, 23.4, 26.7, 21.8, 24.1, 19.5, 25.8, 23.1,
           21.2, 24.7, 22.8, 20.1, 26.2, 18.7, 23.9, 25.4, 21.5, 23.2),
  beak = c(3.2, 3.8, 2.9, 3.5, 4.1, 3.1, 3.6, 2.8, 3.9, 3.4,
           3.0, 3.7, 3.3, 2.9, 4.0, 2.7, 3.5, 3.9, 3.1, 3.4)
)

# Run analysis (auto-detect continuous)
bird_results <- selection_coefficients(
  data = bird_data,
  fitness_col = "fitness",
  trait_cols = c("wing", "tail", "beak")
)

# Note: For continuous fitness, models use relative fitness
print(bird_results)
```

---

## 7. Interpretation Guide

### 7.1 Statistical Significance

| P-value | Evidence |
|---------|----------|
| < 0.001 | Strong evidence of selection |
| 0.001-0.01 | Moderate evidence |
| 0.01-0.05 | Weak evidence |
| > 0.05 | No significant selection detected |

### 7.2 Selection Types

**Linear Selection ($\beta$):**

| $\beta$ sign | Interpretation |
|--------------|----------------|
| Positive (+) | Selection favors larger trait values |
| Negative (-) | Selection favors smaller trait values |
| Zero (0) | No directional selection |

**Quadratic Selection ($\gamma$):**

| $\gamma$ sign | Interpretation | Shape |
|---------------|----------------|-------|
| Negative (-) | Stabilizing selection | ∩-shaped |
| Positive (+) | Disruptive selection | ∪-shaped |
| Zero (0) | No nonlinear selection | Linear |

**Correlational Selection ($\gamma_{ij}$):**

| $\gamma_{ij}$ sign | Interpretation |
|--------------------|----------------|
| Positive (+) | Selection favors positive correlation (both traits high or both low) |
| Negative (-) | Selection favors negative correlation (one high, one low) |
| Zero (0) | Traits evolve independently |

### 7.3 Effect Sizes

**For standardized traits ($\mu=0, \sigma=1$):**

| $|\beta|$ | Effect size |
|-----------|-------------|
| < 0.1 | Weak selection |
| 0.1 - 0.3 | Moderate selection |
| > 0.3 | Strong selection |

### 7.4 Example Interpretations

```
Scenario 1: β_size = 0.25, p = 0.001
→ "Body size experiences significant positive directional selection.
   A one-standard-deviation increase in size increases relative fitness by 0.25 units."

Scenario 2: γ_size = -0.15, p = 0.02
→ "Body size shows significant stabilizing selection.
   Intermediate sizes have highest fitness; extremes are disadvantageous."

Scenario 3: γ_size×beak = 0.32, p = 0.004
→ "Significant positive correlational selection between size and beak depth.
   Selection favors individuals that are either both large or both small."
```

---

## 8. Troubleshooting and Common Issues

### 8.1 Sample Size Recommendations

| Analysis | Minimum n | Recommended n |
|----------|-----------|---------------|
| Linear selection (continuous) | 10 | > 30 |
| Linear selection (binary) | 20 per group | > 50 per group |
| Nonlinear selection | 20 | > 100 |
| Bivariate TPS | 30 | > 200 |

### 8.2 Convergence Issues in Binary Models

**Symptoms:**
- Very large coefficients (> 10)
- Extremely large standard errors
- Warning about "complete separation"

**Solutions:**
```r
# 1. Check for complete separation
table(data$fitness, cut(data$trait, 5))

# 2. Use Firth's penalized likelihood
if (!require("logistf")) install.packages("logistf")
library(logistf)
fit <- logistf(fitness ~ trait, data = data)

# 3. Collect more data
# 4. Simplify model (fewer traits)
```

### 8.3 Missing Data Handling

```r
# Check missing data pattern
library(naniar)
gg_miss_upset(data)

# Option 1: Drop rows with missing (if < 5%)
prepared <- prepare_selection_data(..., na_action = "drop")

# Option 2: Impute missing values
library(mice)
imp <- mice(data, m = 5)
data_complete <- complete(imp, 1)

# Option 3: Use na_action = "warn" to investigate
prepared <- prepare_selection_data(..., na_action = "warn")
```

### 8.4 Multicollinearity

High correlations between traits can inflate standard errors:

```r
# Check correlations
cor(data[, trait_cols])

# Calculate VIF
if (requireNamespace("car", quietly = TRUE)) {
  vif_vals <- car::vif(model)
  # VIF > 10 indicates problematic collinearity
}
```

**Solutions:**
- Combine correlated traits (PCA)
- Use ridge regression
- Collect more data

---

## 9. Output Formats and Reporting

### 9.1 Table Format for Publications

```r
# Create publication-ready table
library(knitr)

results_table <- results %>%
  mutate(
    Beta = sprintf("%.3f", Beta_Coefficient),
    SE = sprintf("%.3f", Standard_Error),
    P = case_when(
      P_Value < 0.001 ~ "<0.001",
      P_Value < 0.01 ~ sprintf("%.3f", P_Value),
      TRUE ~ sprintf("%.3f", P_Value)
    ),
    Sig = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(Term, Type, Beta, SE, P, Sig)

kable(results_table, caption = "Table 1. Selection gradients (±SE)")
```

### 9.2 Visualization for Reports

```r
# Create summary plot
library(ggplot2)

ggplot(results, aes(x = Term, y = Beta_Coefficient, fill = Type)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = Beta_Coefficient - 1.96*Standard_Error,
                    ymax = Beta_Coefficient + 1.96*Standard_Error),
                width = 0.2, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Trait", y = "Selection Gradient (β or γ)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

---

## 10. References

1. Lande, R., & Arnold, S. J. (1983). The measurement of selection on correlated characters. *Evolution*, 37(6), 1210-1226.

2. Arnold, S. J., & Wade, M. J. (1984). On the measurement of natural and sexual selection: Theory. *Evolution*, 38(4), 709-719.

3. Brodie, E. D., Moore, A. J., & Janzen, F. J. (1995). Visualizing and quantifying natural selection. *Trends in Ecology & Evolution*, 10(8), 313-318.

4. Schluter, D. (1988). Estimating the form of natural selection on a quantitative trait. *Evolution*, 42(5), 849-861.

5. Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R* (2nd ed.). Chapman and Hall/CRC.

6. Kingsolver, J. G., et al. (2001). The strength of phenotypic selection in natural populations. *The American Naturalist*, 157(3), 245-261.

7. Hereford, J., Hansen, T. F., & Houle, D. (2004). Comparing strengths of directional selection: How strong is strong? *Evolution*, 58(10), 2133-2143.

## 
