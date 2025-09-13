
find_last <- function(arr, pred) {
  n <- length(arr)
  for (i in n:1) {  # Search from end to beginning
    if (pred(arr[i])) {
      return(i)
    }
  }
  return(-1)
}

get_corner_index <- function(y, xx_in = NULL) {
  N <- length(y)
  M <- N - 1

  if (M <= 0) {
    return(1)  # R is 1-indexed
  }

  E <- numeric(N)

  # Create x values if not provided
  if (is.null(xx_in)) {
    x <- 1:N
  } else {
    x <- xx_in
  }

  # Calculate error for each potential corner point
  for (k in 2:(M)) {  # R indexing: k from 2 to M (equivalent to 1 to M-1 in C++)
    x0 <- x[1]
    xk <- x[k]
    xM <- x[N]  # x[M+1] in C++ becomes x[N] in R

    y0 <- y[1]
    yk <- y[k]
    yM <- y[N]

    # Calculate slopes for two line segments
    slope1 <- (yk - y0) / (xk - x0)
    slope2 <- (yM - yk) / (xM - xk)

    # Define line functions
    L1 <- function(x_val) {
      slope1 * (x_val - x0) + y0
    }

    L2 <- function(x_val) {
      slope2 * (x_val - xk) + yk
    }

    # Calculate sum of squared relative errors for first segment
    sum1 <- 0.0
    for (m in 1:k) {  # m from 1 to k (0 to k in C++)
      err <- (L1(x[m]) - y[m]) / y[m]
      sum1 <- sum1 + err^2
    }

    # Calculate sum of squared relative errors for second segment
    sum2 <- 0.0
    for (m in k:N) {  # m from k to N (k to M in C++)
      err <- (L2(x[m]) - y[m]) / y[m]
      sum2 <- sum2 + err^2
    }

    E[k] <- sqrt(sum1 + sum2)
  }

  # Set boundary values to high number to avoid selection
  INF_APPROX <- 1e300
  E[1] <- INF_APPROX
  E[N] <- INF_APPROX

  # Return index of minimum error
  return(which.min(E))
}

phi <- function(t, eta) {
  exp(-eta / (1 - t^2))
}

dphi_dt <- function(t, eta) {
  denom <- 1 - t^2
  phi_val <- exp(-eta / denom)
  phi_val * (-eta * 2 * t) / (denom^2)
}

phi_vec <- function(t_vec, eta) {
  sapply(t_vec, function(t) exp(-eta / (1 - t^2)))
}

test_function_derivative <- function(radius, dt, order) {
  scale_factor <- (radius * dt)^(-order)

  function(t) {
    if (order == 0) {
      return(phi(t, 9))
    } else {
      return(dphi_dt(t, 9) * (dt * radius)^(-1))
    }
  }
}

get_test_function_support_indices <- function(radius, len_tt) {
  diameter <- 2 * radius + 1
  len_interior <- len_tt - 2

  n <- len_interior - diameter + 1  # number of test functions

  if (diameter > len_interior) {
    stop("diameter must be less than len_interior")
  }

  indices_list <- vector("list", n)

  for (i in 1:n) {
    start <- i
    end <- start + diameter - 1
    indices_list[[i]] <- start:end
  }

  return(indices_list)
}

build_test_function_matrix <- function(tt, radius, order = 0) {
  len_tt <- length(tt)
  dt <- mean(diff(tt))

  # Diameter cannot be larger than the interior of the domain
  diameter <- 2 * radius + 1

  if (len_tt < diameter) {
    warning("diameter outside of domain")
    radius <- floor((len_tt - 2) / 2)
    diameter <- 2 * radius + 1
  }

  indices <- get_test_function_support_indices(radius, len_tt)

  lin <- seq(-1, 1, length.out = diameter)
  xx <- lin[2:(diameter-1)]

  if (order == 0) {
    f_vals <- phi_vec(xx, 9)
  } else {
    f_vals <- sapply(xx, function(t) dphi_dt(t, 9) * (dt * radius)^(-1))
  }

  norm_vals <- phi_vec(xx, 9)
  norm_factor <- sqrt(sum(norm_vals^2))
  v_row <- f_vals / norm_factor

  v_row_padded <- c(0, v_row, 0)

  V <- matrix(0, nrow = length(indices), ncol = len_tt)

  for (i in seq_along(indices)) {
    support_indices <- indices[[i]]
    n_support <- length(support_indices)
    V[i, support_indices] <- v_row_padded[1:n_support]
  }

  return(V)
}

build_full_test_function_matrix <- function(tt, radii, order = 0) {
  test_matrices <- vector("list", length(radii))

  for (i in seq_along(radii)) {
    test_matrices[[i]] <- build_test_function_matrix(tt, radii[i], order)
  }

  V_full <- do.call(rbind, test_matrices)

  return(V_full)
}

find_min_radius_int_error <- function(U, tt, radius_min, radius_max, num_radii, sub_sample_rate) {
  Mp1 <- nrow(U)  # Number of data points
  D <- ncol(U)    # Dimension of the system

  step <- max(1, ceiling((radius_max - radius_min) / num_radii))
  radii <- seq(radius_min, radius_max - 1, by = step)

  errors <- numeric(length(radii))

  IX <- floor((Mp1 - 1) / sub_sample_rate)

  for (i in seq_along(radii)) {
    radius <- radii[i]
    V_r <- build_test_function_matrix(tt, radius)
    K <- nrow(V_r)  # Number of test functions for a given radius

    # Create G tensor (K, Mp1, D)
    # Element wise for each dimension phi(t_i) and u(t_i) for all phi_k
    G <- array(0, dim = c(K, Mp1, D))
    for (k in 1:K) {
      for (d in 1:D) {
        G[k, , d] <- V_r[k, ] * U[, d]
      }
    }

    GT <- aperm(G, c(1, 3, 2))

    GT_reshaped <- matrix(GT, nrow = K * D, ncol = Mp1)

    f_hat_G <- apply(GT_reshaped, 1, fft)

    if (is.matrix(f_hat_G)) {
      f_hat_G_imag <- Im(f_hat_G[IX + 1, ])  # R is 1-indexed
    } else {
      f_hat_G_imag <- Im(f_hat_G[[IX + 1]])
    }

    errors[i] <- sqrt(sum(f_hat_G_imag^2))  # L2 norm
  }

  log_errors <- log(errors)
  ix <- get_corner_index(log_errors)

  return(list(index = ix, errors = errors, radii = radii))
}

build_full_test_function_matrices <- function(U, tt, test_function_params, compute_svd = TRUE) {
  # cat("<< Building test matrices >>\n")

  dt <- mean(diff(tt))
  mp1 <- nrow(U)

  radii <- test_function_params$radius_params
  radius_min_time <- test_function_params$radius_min_time
  radius_max_time <- test_function_params$radius_max_time

  min_radius <- max(ceiling(radius_min_time / dt), 2)
  max_radius <- floor(radius_max_time / dt)

  radius_min_max <- floor(max_radius / max(radii))

  if (radius_min_max < min_radius) {
    radius_min_max <- min_radius * 10
  }

  max_radius_for_interior <- floor((mp1 - 2) / 2)
  if (max_radius > max_radius_for_interior) {
    max_radius <- max_radius_for_interior
  }

  if (radius_min_max <= min_radius) {
    radius_min_max <- min_radius * 10
  }

  #cat(sprintf("  Min radius: %d\n", min_radius))
  #cat(sprintf("  Max radius: %d\n", max_radius))
  #cat(sprintf("  Minmax radius: %d\n", radius_min_max))

  result <- find_min_radius_int_error(U, tt, min_radius, radius_min_max,
                                      num_radii = 10, sub_sample_rate = 1)
  ix <- result$index
  errors <- result$errors
  radii_sweep <- result$radii

  min_radius_int_error <- radii_sweep[ix]

  min_radius_errors <- errors
  min_radius_radii <- radii_sweep
  min_radius_ix <- ix
  min_radius <- min_radius_int_error

  #cat(sprintf("  Integral Error min radius: %d\n", min_radius_int_error))

  radii <- test_function_params$radius_params * min_radius_int_error

  radii_filtered <- radii[radii < max_radius]

  if (length(radii_filtered) == 0) {
    radii_filtered <- max_radius
  }

  #cat(sprintf("  Radii [%s]\n", paste(radii_filtered, collapse = ", ")))

  V_ <- build_full_test_function_matrix(tt, radii_filtered, order = 0)
  V_prime_ <- build_full_test_function_matrix(tt, radii_filtered, order = 1)

  if (!compute_svd) {
    return(list(
      V = V_,
      V_prime = V_prime_,
      radii = radii_filtered,
      min_radius_errors = min_radius_errors,
      min_radius_radii = min_radius_radii,
      min_radius_ix = min_radius_ix,
      min_radius = min_radius
    ))
  }

  k_full <- nrow(V_)
  #cat(sprintf("  K Full: %d\n", k_full))

  SVD <- svd(V_)
  U_svd <- SVD$u
  singular_values <- SVD$d
  Vt <- t(SVD$v)

  sum_singular_values <- sum(singular_values)
  k_max <- min(k_full, mp1, test_function_params$k_max)

  info_numbers <- numeric(k_max)
  for (i in 2:k_max) {
    info_numbers[i] <- sum(singular_values[1:(i-1)]) / sum_singular_values
  }

  condition_numbers <- singular_values[1] / singular_values

  k1 <- find_last(condition_numbers, function(x) {
    x < test_function_params$max_test_fun_condition_number
  })
  k2 <- find_last(info_numbers, function(x) {
    x < test_function_params$min_test_fun_info_number
  })

  if (k1 == -1) k1 <- .Machine$integer.max
  if (k2 == -1) k2 <- .Machine$integer.max

  K <- min(k1, k2, k_max)

  #cat(sprintf("  Condition Number is now: %.4f\n", condition_numbers[K]))
  #cat(sprintf("  Info Number is now: %.4f\n", info_numbers[K]))
  #cat(sprintf("  K is: %d\n", K))

  V_final <- Vt[1:K, , drop = FALSE]

  #cat("  Calculating Vprime\n")

  U_T <- t(U_svd)[1:K, , drop = FALSE]
  UV <- U_T %*% V_prime_
  inv_s <- 1.0 / singular_values[1:K]
  V_prime_final <- UV * inv_s

  return(list(
    V = V_final,
    V_prime = V_prime_final,
    radii = radii_filtered,
    min_radius_errors = min_radius_errors,
    min_radius_radii = min_radius_radii,
    min_radius_ix = min_radius_ix,
    min_radius = min_radius,
    singular_values = singular_values[1:K],
    condition_numbers = condition_numbers[1:K],
    info_numbers = info_numbers[1:K],
    K = K
  ))
}


