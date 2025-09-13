
convolve1d_inner <- function(signal, kernel) {
  signal_len <- length(signal)
  kernel_len <- length(kernel)

  output_len <- signal_len - kernel_len + 1

  if (output_len <= 0) {
    stop("Signal too short for convolution with this kernel")
  }

  result <- numeric(output_len)

  for (i in 1:output_len) {
    result[i] <- sum(signal[i:(i + kernel_len - 1)] * kernel)
  }

  return(result)
}

fdcoeffF <- function(k, xbar, x) {
  n <- length(x)

  if (k >= n) {
    stop("length of x must be larger than k")
  }

  m <- k  # change to m=n-1 if you want coefficients for all derivatives

  c1 <- 1.0
  c4 <- x[1] - xbar
  C <- matrix(0, nrow = n, ncol = m + 1)

  C[1, 1] <- 1

  for (i in 2:n) {
    mn <- min(i - 1, m)
    c2 <- 1.0
    c5 <- c4
    c4 <- x[i] - xbar

    for (j in 1:(i-1)) {
      c3 <- x[i] - x[j]
      c2 <- c2 * c3

      if (j == (i - 1)) {
        if (mn >= 1) {
          for (s in (mn + 1):2) {
            C[i, s] <- c1 * ((s-1) * C[i - 1, s - 1] - c5 * C[i - 1, s]) / c2
          }
        }
        C[i, 1] <- -c1 * c5 * C[i - 1, 1] / c2
      }

      if (mn >= 1) {
        for (s in (mn + 1):2) {
          C[j, s] <- (c4 * C[j, s] - (s-1) * C[j, s - 1]) / c3
        }
      }
      C[j, 1] <- c4 * C[j, 1] / c3
    }
    c1 <- c2
  }

  return(C)
}

estimate_std <- function(U, k = 6) {
  D <- ncol(U)
  std_vec <- numeric(D)

  for (d in 1:D) {
    f <- U[, d]
    x <- seq(-k-2, k+2)
    C <- fdcoeffF(k, 0, x)

    filter <- C[, ncol(C)]

    filter <- filter / sqrt(sum(filter^2))

    convolved <- convolve1d_inner(f, filter)
    squared <- convolved^2
    mean_squared <- mean(squared)

    std_vec[d] <- sqrt(mean_squared)
  }

  return(std_vec)
}

