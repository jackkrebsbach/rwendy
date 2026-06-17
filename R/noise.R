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

#' Estimate noise standard deviation from the robust scale of a high-order
#' finite-difference (high-pass) residual.
#'
#' The order-\code{k} FD filter annihilates smooth signal and normalizes white
#' noise to unit variance; the per-state noise SD is then the MAD (median
#' absolute deviation, a robust scale) of the filtered residual rather than its
#' RMS. On coarse grids residual signal curvature (\eqn{\propto \Delta t^k
#' f^{(k)}}) leaks through the filter: the RMS adds it in quadrature at every
#' point and floors \eqn{\hat\sigma} as the true noise shrinks (mistaking
#' curvature for noise on under-resolved grids), whereas the MAD's ~50\%
#' breakdown ignores the minority of leakage-dominated stencil positions and
#' tracks the true noise down. Identical to the RMS on resolved grids and at
#' high noise; 3-14x more accurate at low-noise/coarse (validated, see
#' examples/validation). Past 50\% contamination (very coarse grids) the data is
#' under-sampled and no aggregator recovers \eqn{\sigma}.
#'
#' @param U Numeric matrix of observations (rows = time points, cols = states).
#' @param k Integer; order of the finite-difference filter (default 6).
#' @return Numeric vector of length \code{ncol(U)}: per-state noise SD estimates.
#' @keywords internal
estimate_std <- function(U, k = 6) {
  D <- ncol(U)
  std_vec <- numeric(D)

  for (d in 1:D) {
    f <- U[, d]
    x <- seq(-k-2, k+2)
    C <- fdcoeffF(k, 0, x)

    filter <- C[, ncol(C)]

    filter <- filter / sqrt(sum(filter^2))

    convolved <- stats::convolve(f, filter, type = "filter")

    # MAD (robust scale), not RMS: on coarse grids signal curvature leaks through
    # the FD filter and mean(convolved^2) adds it in quadrature, flooring sigma_hat
    # as the true noise shrinks. MAD's ~50% breakdown ignores the leakage-dominated
    # minority of positions; 1.4826 (mad default) makes it consistent for Gaussians.
    std_vec[d] <- stats::mad(convolved)
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

