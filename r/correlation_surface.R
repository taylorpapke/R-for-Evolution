#' Standardized Correlational Selection Surface (Contour Plot)
#'
#' Plot thin-plate spline surface as filled contours.
#'
#' @param tps A list returned by \code{correlational_tps()}.
#' @param trait_cols Character vector (length 2). Trait names.
#' @param bins Integer; number of filled contour bins (default 12).
#'
#' @return A \code{ggplot} object.
#' @export

correlation_surface <- function(tps, trait_cols, bins = 12) {
  stopifnot(is.list(tps), "grid" %in% names(tps))
  df <- tps$grid
  if (length(trait_cols) != 2L) stop("trait_cols must be length 2.")
  if (!all(trait_cols %in% names(df))) stop("Trait columns not found in tps$grid.")
  if (!".fit" %in% names(df)) stop("tps$grid must contain '.fit'.")
  if (!is.numeric(df[[trait_cols[1]]]) || !is.numeric(df[[trait_cols[2]]])) stop("Traits must be numeric.")
  if (!is.numeric(df$.fit)) df$.fit <- as.numeric(df$.fit)

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[trait_cols[1]]],
                                   y = .data[[trait_cols[2]]],
                                   z = .data[[".fit"]])) +
    ggplot2::geom_contour_filled(bins = bins) +
    ggplot2::geom_contour(color = "black", alpha = 0.3) +
    ggplot2::labs(x = trait_cols[1], y = trait_cols[2], fill = "Fitness") +
    ggplot2::theme_minimal(base_size = 12)
}
