
phi <- function(t, eta = 9) {
  exp(-eta / (1 - t^2))
}

test_function_derivative <- function(radius, dt, order) {
  # Chain rule to account for (t/a)^2 we get factors of (1/a), a=dt*radius
  scale_factor <- (radius * dt)^(-order)
  t <- symengine::S("t")
  phi_sym <- phi(t)

  phi_deriv_sym <- phi_sym
  for (i in seq_len(order)) {
    phi_deriv_sym <- symengine::D(phi_deriv_sym, "t")
  }

  symengine::lambdify(scale_factor * phi_deriv_sym, t)
}

get_test_function_support_indices <- function(radius, len_tt) {
  diameter <- 2 * radius + 1
  len_interior <- len_tt - 2

  n <- len_interior - diameter + 2 # Number of test functions to build dense sweep

  if (diameter > len_interior) {
    stop("diameter must be less than len_interior")
  }

  indices <- lapply(2:n, function(i) {
    end <- i + diameter - 1
    i:end
  })
  return(indices)
}

build_test_function_matrix <- function(tt, radius, order = 0) {
  len_tt <- length(tt)
  dt <- mean(diff(tt))

  diameter <- 2 * radius + 1

  if (len_tt < diameter) {
    warning("diameter outside of domain")
    radius <- floor((len_tt - 2) / 2)
    diameter <- 2 * radius + 1
  }

  # For one radius get the support indices for all phi_k (different centers)
  indices <- get_test_function_support_indices(radius, len_tt)

  # Don't include the endpoints (outside of support)
  lin <- seq(-1, 1, length.out = diameter)
  xx <- lin[2:(diameter-1)]

  ph <- test_function_derivative(radius, dt, order)
  f_vals <- ph(xx)
  norm_factor <- sqrt(sum(phi(xx)^2))
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
  Mp1 <- nrow(U)
  D <- ncol(U)

  step <- max(1, ceiling((radius_max - radius_min) / num_radii))
  radii <- seq(radius_min, radius_max - 1, by = step)

  errors <- numeric(length(radii))

  IX <- floor((Mp1 - 1) / sub_sample_rate)

  for (i in seq_along(radii)) {
    radius <- radii[i]
    V_r <- build_test_function_matrix(tt, radius)
    K <- nrow(V_r)

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
      f_hat_G_imag <- Im(f_hat_G[IX + 1,])
    } else {
      f_hat_G_imag <- Im(f_hat_G[[IX + 1]])
    }

    errors[i] <- sqrt(sum(f_hat_G_imag^2))
  }

  log_errors <- log(errors)
  ix <- get_corner_index(log_errors)

  return(list(index = ix, errors = errors, radii = radii))
}

build_full_test_function_matrices <- function(U, tt, test_function_params, compute_svd = TRUE) {
  #cat("<< Building test matrices >>\n")
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

  # cat(sprintf("  Min radius: %d\n", min_radius))
  # cat(sprintf("  Max radius: %d\n", max_radius))
  # cat(sprintf("  Minmax radius: %d\n", radius_min_max))

  result <- find_min_radius_int_error(U, tt, min_radius, radius_min_max,
                                      num_radii = 100, sub_sample_rate = 1)
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

  # cat(sprintf("  Condition Number is now: %.4f\n", condition_numbers[K]))
  # cat(sprintf("  Info Number is now: %.4f\n", info_numbers[K]))
  # cat(sprintf("  K is: %d\n", K))

  V_final <- Vt[1:K, , drop = FALSE]

  # cat("  Calculating Vprime\n")

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

find_last <- function(arr, pred) {
  n <- length(arr)
  for (i in n:1) {
    if (pred(arr[i])) {
      return(i)
    }
  }
  return(-1)
}

get_corner_index <- function(y, xx_in = NULL) {
  N <- length(y)
  M <- N - 1

  if (M <= 0) { return(1) }

  E <- numeric(N)

  if (is.null(xx_in)) {
    x <- 1:N
  } else {
    x <- xx_in
  }

  for (k in 2:(M)) {
    x0 <- x[1]
    xk <- x[k]
    xM <- x[N]

    y0 <- y[1]
    yk <- y[k]
    yM <- y[N]

    slope1 <- (yk - y0) / (xk - x0)
    slope2 <- (yM - yk) / (xM - xk)

    L1 <- function(x_val) {
      slope1 * (x_val - x0) + y0
    }

    L2 <- function(x_val) {
      slope2 * (x_val - xk) + yk
    }

    sum1 <- 0.0
    for (m in 1:k) {
      err <- (L1(x[m]) - y[m]) / y[m]
      sum1 <- sum1 + err^2
    }

    sum2 <- 0.0
    for (m in k:N) {
      err <- (L2(x[m]) - y[m]) / y[m]
      sum2 <- sum2 + err^2
    }

    E[k] <- sqrt(sum1 + sum2)
  }

  # Don't want to choose the first or last index
  INF_APPROX <- 1e300
  E[1] <- INF_APPROX
  E[N] <- INF_APPROX

  return(which.min(E))
}
