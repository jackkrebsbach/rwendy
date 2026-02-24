#' @importFrom stats smooth.spline approx spline loess lm.fit ksmooth predict
#' @importFrom utils head tail
NULL

# Apply one interpolation method to map rows of U from tt_vec onto tt_target.
# Returns an interpolated matrix with nrow = length(tt_target), ncol = D.
# @keywords internal
interpolate_to_grid <- function(U, tt_vec, tt_target, method) {
  result <- switch(method,
    spline = {
      fits <- apply(U, 2, function(col) smooth.spline(tt_vec, col))
      sapply(fits, function(fit) predict(fit, tt_target)$y)
    },
    linear   = apply(U, 2, function(col) approx(tt_vec, col, xout = tt_target)$y),
    cubic    = apply(U, 2, function(col) spline(tt_vec, col, xout = tt_target, method = "natural")$y),
    cubic_ls = {
      X     <- cbind(1, tt_vec, tt_vec^2, tt_vec^3)
      X_new <- cbind(1, tt_target, tt_target^2, tt_target^3)
      apply(U, 2, function(col) drop(X_new %*% lm.fit(X, col)$coefficients))
    },
    loess  = apply(U, 2, function(col) predict(loess(col ~ tt_vec), newdata = tt_target)),
    kernel = {
      bw <- diff(range(tt_vec)) / sqrt(length(tt_vec))
      apply(U, 2, function(col) ksmooth(tt_vec, col, kernel = "normal", bandwidth = bw, x.points = tt_target)$y)
    },
    stop("Unknown interpolation_method: ", method)
  )

  # Substitute original observed values at their exact time points
  tol <- sqrt(.Machine$double.eps)
  for (i in seq_along(tt_vec)) {
    j <- which(abs(tt_target - tt_vec[i]) < tol)
    if (length(j) == 1L) result[j, ] <- U[i, ]
  }

  result
}

# Apply one named interpolation method to (U, tt), handling non-uniform spacing
# and minimum-number-of-points densification.  Returns list(U, tt).
# @keywords internal
interpolate_data <- function(U, tt, method, control) {
  diff_dt <- diff(as.vector(tt))
  dt      <- mean(diff_dt, na.rm = TRUE)

  if (max(abs(diff(tt) - dt)) > sqrt(.Machine$double.eps)) {
    cat("Non uniform spacing detected, interpolating data using", method, "...\n")
    n      <- max(floor((max(tt) - min(tt)) / dt), control$min_number_points)
    tt_new <- seq(min(tt), max(tt), length.out = n)
    U      <- interpolate_to_grid(U, as.vector(tt), tt_new, method)
    tt     <- matrix(tt_new, ncol = 1)
  }

  if (nrow(U) < control$min_number_points) {
    tt_vec   <- as.vector(tt)
    tt_dense <- tt_vec
    while (2 * length(tt_dense) - 1 <= control$min_number_points * 2) {
      mids     <- (head(tt_dense, -1) + tail(tt_dense, -1)) / 2
      tt_dense <- sort(c(tt_dense, mids))
    }
    U  <- interpolate_to_grid(U, tt_vec, tt_dense, method)
    tt <- matrix(tt_dense, ncol = 1)
  }

  list(U = U, tt = tt)
}
