# https://faculty.washington.edu/rjl/fdmbook/matlab/fdcoeffF.m
fdcoeffF <- function(k, xbar, x) {
  n <- length(x)

  if (k >= n) {
    stop("length of x must be larger than k")
  }

  m <- k  # change to m=n-1 if you want coefficients for all derivatives

  c1 <- 1
  c4 <- x[1] - xbar
  C <- matrix(0, nrow = n, ncol = m + 1)
  C[1, 1] <- 1
  for (i in 1:(n-1)) {
    i1 <- i+1
    mn <- min(i, m)
    c2 <- 1
    c5 <- c4
    c4 <- x[i1] - xbar
    for (j in 0:(i-1)) {
      j1 <- j+1
      c3 <- x[i1] - x[j1]
      c2 <- c2 * c3
      if (j == (i - 1)) {
          for (s in seq(mn, 1, by = -1)) {
            s1 <- s + 1
            C[i1, s1] <- c1 * (s * C[i1 - 1, s1 - 1] - c5 * C[i1 - 1, s1]) / c2
          }
        C[i1, 1] <- -c1 * c5 * C[i1 - 1, 1] / c2
      }
      for (s in seq(mn, 1, by = -1)) {
          s1 <- s + 1
          C[j1, s1] <- (c4 * C[j1, s1] - s * C[j1, s1 - 1]) / c3
      }
      C[j1, 1] <- c4 * C[j1, 1] / c3
    }
    c1 <- c2
  }

  return(C)
}

#' Estimate noise standard deviation from a high-order finite-difference
#' (high-pass) residual.
#'
#' The order-\code{k} FD filter annihilates smooth signal and normalizes white
#' noise to unit variance; the per-state noise SD is the scale of the filtered
#' residual. \code{robust = FALSE} (default) uses the RMS
#' (\eqn{\sqrt{\mathrm{mean}(\cdot^2)}}): the efficient choice, and on coarse
#' grids the resulting signal-curvature leakage usefully inflates
#' \eqn{\hat\sigma} so the GLS weights cover the discretization error (better
#' \eqn{\hat p}). \code{robust = TRUE} uses the MAD: its ~50\% breakdown ignores
#' the leakage-dominated stencil positions and recovers the TRUE measurement
#' noise on coarse grids -- preferred for UQ (covariance calibration), where the
#' inflated RMS over-covers. Identical on resolved grids and at high noise; the
#' two diverge only when curvature leaks (coarse + low noise). See
#' examples/validation/noise_mad_phat_regression.R for the \eqn{\hat p} trade-off
#' that motivates keeping RMS for estimation and MAD only for UQ.
#'
#' @param U Numeric matrix of observations (rows = time points, cols = states).
#' @param k Integer; order of the finite-difference filter (default 6).
#' @param robust Logical; \code{TRUE} uses the MAD (robust) scale instead of the
#'   RMS. Default \code{FALSE}.
#' @return Numeric vector of length \code{ncol(U)}: per-state noise SD estimates.
#' @keywords internal
estimate_std <- function(U, k = 6, robust = FALSE) {
  D <- ncol(U)
  std_vec <- numeric(D)

  for (d in 1:D) {
    f <- U[, d]
    x <- seq(-k-2, k+2)
    C <- fdcoeffF(k, 0, x)

    filter <- C[, ncol(C)]

    filter <- filter / sqrt(sum(filter^2))

    convolved <- stats::convolve(f, filter, type = "filter")

    std_vec[d] <- if (robust) stats::mad(convolved) else sqrt(mean(convolved^2))
  }

  return(std_vec)
}

# For multiplicative log normal noise we can't have negative data or data at 0
preprocess_data <- function(U, tt){
  mp1 <- nrow(U)
  indices <- logical(mp1)
  for(i in seq(mp1)){
    indices[i] <- all(U[i,] > 0)
  }
  return(list(U = log(U[indices,  ,drop = FALSE]), tt = tt[indices, drop = FALSE]))
}

