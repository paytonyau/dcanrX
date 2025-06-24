#' Calculate Biweight Midvariance
#'
#' A robust measure of variance based on the median absolute deviation.
#' This is a helper function for ROS-DET's weighting step.
#'
#' @param x A numeric vector.
#' @return The biweight midvariance of the vector.
#' @keywords internal
biweight_midvariance <- function(x) {
  x <- x[!is.na(x)]
  median_x <- stats::median(x)
  # Use median absolute deviation as a robust measure of spread
  mad_x <- stats::mad(x)

  # Avoid division by zero if spread is zero
  if (mad_x == 0) return(0)

  # Calculate weights based on distance from the median
  u <- (x - median_x) / (9 * mad_x)
  w <- (1 - u^2)^2
  w[abs(u) > 1] <- 0

  # Weighted sum of squares
  numerator <- sum(w * (x - median_x)^2)
  denominator <- sum(w)

  # Handle cases where all weights are zero
  if (denominator == 0) return(0)

  # Return the final variance value
  return(length(x) * numerator / (denominator^2))
}
