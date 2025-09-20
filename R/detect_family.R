#' Detect Family from Fitness Vector
#'
#' Heuristically picks the appropriate GLM family based on fitness values.
#' Robust to edge cases like all zeros, all ones, and small samples.
#'
#' @param y Numeric vector of fitness values
#' @return List with type ("binary" or "continuous") and family object
#' @export

detect_family <- function(y) {
  # remove NA 
  y_clean <- y[!is.na(y)]
  if (length(y_clean) == 0) {
    stop("No non-NA values in fitness vector")
  }
  
  # get unique value 
  unique_vals <- sort(unique(y_clean))
  
  # check if it is binary 
  is_binary <- length(unique_vals) <= 2 && 
    all(unique_vals %in% c(0, 1)) &&
    length(y_clean) >= 10  
  
  if (is_binary) {
    # check separation 
    if (all(y_clean == 0) || all(y_clean == 1)) {
      warning("Complete separation detected: all fitness values are ", 
              ifelse(all(y_clean == 0), "0", "1"),
              ". Consider using Firth regression or more data.")
    }
    
    return(list(
      type = "binary",
      family = binomial("logit")
    ))
  } else {
    return(list(
      type = "continuous", 
      family = gaussian()
    ))
  }
}