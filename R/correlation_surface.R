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


correlation_surface_enhanced <- function(tps, trait_cols, original_data = NULL, fitness_col = NULL, bins = 12) {
  p <- correlation_surface(tps, trait_cols, bins)
  
  if (!is.null(original_data) && !is.null(fitness_col)) {
    if (all(trait_cols %in% names(original_data)) && fitness_col %in% names(original_data)) {
      plot_data <- original_data[, c(trait_cols, fitness_col)]
      plot_data <- plot_data[complete.cases(plot_data), ]
      
      p <- p + 
        ggplot2::geom_point(
          data = plot_data,
          ggplot2::aes(x = .data[[trait_cols[1]]], 
                       y = .data[[trait_cols[2]]],
                       color = as.factor(.data[[fitness_col]])),
          alpha = 0.6, size = 1.5
        ) +
        ggplot2::scale_color_manual(
          values = c("0" = "red", "1" = "green"),
          labels = c("0" = "Perished", "1" = "Survived"),
          name = "Outcome"
        )
    }
  }
  
  df <- tps$grid
  optimal_point <- df[which.max(df$.fit), ]
  p <- p + 
    ggplot2::geom_point(
      data = optimal_point,
      ggplot2::aes(x = .data[[trait_cols[1]]], y = .data[[trait_cols[2]]]),
      color = "yellow", size = 3, shape = 18
    )
  
  return(p)
}