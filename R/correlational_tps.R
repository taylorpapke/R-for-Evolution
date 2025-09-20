#' Thin-Plate Spline Surface for Correlational Selection (2 Traits)
#'
#' Fits \code{fields::Tps} on two traits and fitness.
#' For continuous fitness, you may use relative fitness (divide by mean).
#' For binary fitness (0/1), \code{use_relative} is ignored.
#'
#' @param data A data.frame.
#' @param fitness_col Character. Fitness column.
#' @param trait_cols Character vector of length 2. Two numeric trait columns.
#' @param use_relative Logical; if TRUE and fitness is continuous, divide by mean (default TRUE).
#' @param grid_n Integer; resolution per trait for the prediction grid (default 60).
#'
#' @return A list with \code{model} and \code{grid} (data.frame of gridded predictions).
#' @export

correlational_tps <- function(data, fitness_col, trait_cols, use_relative = TRUE, grid_n = 60) {
  if (length(trait_cols) != 2L) stop("`trait_cols` must be length 2.")
  if (!all(c(fitness_col, trait_cols) %in% names(data))) stop("Columns not found.")
  
  # keep complete rows
  need <- c(fitness_col, trait_cols)
  keep <- stats::complete.cases(data[, need, drop = FALSE])
  df <- data[keep, , drop = FALSE]
  if (!nrow(df)) stop("No complete cases for fitness and traits.")
  
  # strict numeric coercion with clear errors
  to_num <- function(x, nm) {
    if (is.numeric(x)) return(as.numeric(x))
    x2 <- suppressWarnings(as.numeric(x))
    if (any(is.na(x2) & !is.na(x))) stop(sprintf("Column '%s' cannot be coerced to numeric.", nm))
    x2
  }
  x1 <- to_num(df[[trait_cols[1]]], trait_cols[1])
  x2 <- to_num(df[[trait_cols[2]]], trait_cols[2])
  y  <- to_num(df[[fitness_col]],  fitness_col)
  
  # relative fitness only if not binary
  u <- sort(unique(stats::na.omit(y)))
  is_bin <- length(u) <= 2 && all(u %in% c(0,1))
  if (use_relative && !is_bin) {
    mu <- mean(y, na.rm = TRUE)
    if (!is.finite(mu) || mu == 0) stop("Mean fitness not finite/zero; cannot compute relative fitness.")
    y <- y / mu
  }
  
  # enforce storage mode
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  y <- as.numeric(y)
  storage.mode(x1) <- "double"
  storage.mode(x2) <- "double"
  storage.mode(y)  <- "double"
  
  ok2 <- is.finite(x1) & is.finite(x2) & is.finite(y)
  x1 <- x1[ok2]; x2 <- x2[ok2]; y <- y[ok2]
  if (!length(y)) stop("No finite observations after filtering.")
  
  X <- cbind(x1, x2)
  colnames(X) <- trait_cols
  storage.mode(X) <- "double"
  
  fit <- tryCatch(
    fields::Tps(X, y),
    error = function(e) {
      stop("fields::Tps failed. Underlying error: ", conditionMessage(e))
    }
  )
  
  # Fix: Create numerical sequences directly to avoid expand.grid issues
  x1_seq <- as.numeric(seq(min(x1), max(x1), length.out = grid_n))
  x2_seq <- as.numeric(seq(min(x2), max(x2), length.out = grid_n))
  storage.mode(x1_seq) <- "double"
  storage.mode(x2_seq) <- "double"
  
  grid_df <- expand.grid(
    x1_seq, 
    x2_seq
  )
  names(grid_df) <- trait_cols
  
  # Fix: Make sure to convert to a numeric matrix
  Z_matrix <- as.matrix(grid_df)
  Z_matrix <- matrix(as.numeric(Z_matrix), nrow = nrow(Z_matrix), ncol = ncol(Z_matrix))
  storage.mode(Z_matrix) <- "double"
  colnames(Z_matrix) <- trait_cols
  
  zhat <- tryCatch(
    predict(fit, Z_matrix),
    error = function(e) {
      stop("predict(Tps, ...) failed. Underlying error: ", conditionMessage(e))
    }
  )
  grid_df$.fit <- as.numeric(zhat)
  
  list(model = fit, grid = grid_df)
}
