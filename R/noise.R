
convolve1d_inner <- function(signal, kernel) {
  signal_len <- length(signal)
  kernel_len <- length(kernel)

  output_len <- signal_len - kernel_len + 1

  if (output_len <= 0) {
    stop("Signal too short for convolution with this kernel")
  }

  result <- numeric(output_len)

  for (i in 1:output_len) {
     sum <- 0
     for(j in 1:kernel_len){
       sum <- sum +  signal[i + j - 1] * kernel[j];
     }
    result[i] <- sum
  }

  return(result)
}


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

