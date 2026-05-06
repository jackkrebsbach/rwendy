# Piecewise polynomial test function 
# Advantageous because we have access to analytic form of Fourier coefficients
# 𝛹(t; r, p) = [(1- (t/r)²)]ᵖ₊
psi <- function(t, r, p = 16){
  (1 - ( t / r)^2)^p
} 

# 𝐂^{∞} Bump test function
# 𝜱(t; r, n) = exp(η/ [(1 - (t/r)²)]₊)
phi <- function(t, r, eta = 9) {
  exp(-eta / (1 - (t / r)^2))
}

# Fourier coefficient of 𝛹(t; r, p)  r = dt * radius
psi_hat <- function(freqs, radius, dt, T, C, p = 16){
   psih <- complex(length(freqs))

   n_max <- floor(length(freqs) /  2)  
   r <- dt * radius

   n0 <- 2 * r^(2 * p + 1) * sum(sapply(0:p, function(j){choose(p, j) * (-1)^j / (2 * j + 1)}))
   ng1 <- sapply(1:n_max, function(n){
             sqrt(pi) * (r * T /(n * pi))^(p + 1/2) * gamma(p+1) * besselJ(2 * pi * n * r / T, p + 1/2)
          })
   coeffs <- (C / sqrt(T)) *  c(n0, ng1)
    
   psih[1:(n_max + 1)] <- coeffs
   n_neg <- length(psih) - n_max - 1
   psih[(n_max + 2):length(psih)] <- rev(coeffs[2:(n_neg + 1)])
   
   return(psih)
}

test_function_derivative <- function(test_function, radius, dt, order) {
  t <- symengine::S("t")
  r <- radius * dt
  phi_sym <- test_function(t, r)

  phi_deriv_sym <- phi_sym
  for (i in seq_len(order)) {
    phi_deriv_sym <- symengine::D(phi_deriv_sym, "t")
  }

  symengine::lambdify(phi_deriv_sym, t)
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

build_test_function_matrix <- function(test_function, tt, radius, order = 0) {
  len_tt <- length(tt)
  dt <- mean(diff(tt))

  diameter <- 2 * radius + 1

  if (diameter > len_tt - 2) {
    radius <- floor((len_tt - 3) / 2)
    diameter <- 2 * radius + 1
  }

  # For one radius get the support indices for all phi_k (different centers)
  indices <- get_test_function_support_indices(radius, len_tt)
  r <- dt * radius
  lin <- seq(-r, r, length.out = diameter)
  # Don't include the endpoints (outside of support for 𝚽 or zero for 𝚿) 
  xx <- lin[2:(diameter-1)]
  ph <- test_function_derivative(test_function, radius, dt, order)
  f_vals <- ph(xx)
  norm_factor <- sqrt(sum(test_function(xx, r)^2))
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

build_full_test_function_matrix <- function(test_function, tt, radii, order = 0) {
  do.call(rbind, lapply(radii, function(r) build_test_function_matrix(test_function, tt, r, order)))
}

find_min_radius_int_error <- function(U, tt, radius_min, radius_max, num_radii, sub_sample_rate = 2) {

  D   <- ncol(U)
  Mp1 <- nrow(U)

  step  <- max(1, ceiling((radius_max - radius_min) / num_radii))
  radii <- seq(radius_min, radius_max - 1, by = step)
  
  errors  <- numeric(length(radii))

  IX <- floor(Mp1 / sub_sample_rate) 

  for (i in seq_along(radii)) {
    radius <- radii[i]
    V_r <- build_test_function_matrix(phi, tt, radius)
    K   <- nrow(V_r)

    GT <- do.call(rbind, lapply(seq_len(D), function(d){
      sweep(V_r, 2, U[, d], `*`)
    }))

    f_hat_G_imag <- Im(mvfft(t(GT))[IX, ])
    
    errors[i] <- (4 * pi * IX / Mp1) * sqrt(sum(f_hat_G_imag^2) / K)
    
  }

  log_errors <- log(errors)
  ix <- get_corner_index(log_errors)

  return(list(index = ix, errors = errors, radii = radii))
}

build_full_test_function_matrices_ssl <- function(U, tt, control) {

  dt <- mean(diff(tt))
  mp1 <- nrow(U)

  if (!is.null(control$fixed_radius)) {
    radius_c <- control$fixed_radius
    rc_errors <- NULL
    rc_radii <- NULL
  } else {
    data <- compute_r_c_hat(U, tt, control$S, control$p)
    radius_c <- data$rc
    rc_errors <- data$ehat
    rc_radii <- data$radii
  }

  V  <- build_full_test_function_matrix(psi, tt, c(radius_c), order = 0)
  Vp <- build_full_test_function_matrix(psi, tt, c(radius_c), order = 1)

  return(list(V = V, V_prime = Vp, radius_c = radius_c, rc_errors = rc_errors, rc_radii = rc_radii))

}

build_full_test_function_matrices_msg <- function(U, tt, control, compute_svd = TRUE) {
  # cat("<< Building test matrices >>\n")
  dt <- mean(diff(tt))
  mp1 <- nrow(U)

  radii <- control$radius_params
  radius_min_time <- control$radius_min_time
  radius_max_time <- control$radius_max_time

  test_fun <- get(control$test_fun)

  min_radius <- max(ceiling(radius_min_time / dt), 2)
  max_radius <- floor(radius_max_time / dt)

  max_radius_for_interior <- floor((mp1 - 2) / 2)
  if (max_radius > max_radius_for_interior) {
    max_radius <- max_radius_for_interior
  }

  # cat(sprintf("  Min radius: %d\n", min_radius))
  # cat(sprintf("  Max radius: %d\n", max_radius))

  if (!is.null(control$fixed_radius)) {
    min_radius_int_error <- control$fixed_radius
    min_radius_errors <- NA
    min_radius_radii <- NA
    min_radius_ix <- NA
  } else {
    if(control$test_fun == "phi"){
      # 𝜱
      result <- find_min_radius_int_error(U, tt, min_radius, max_radius, num_radii = 100, sub_sample_rate = 2)
      min_radius_int_error <- result$radii[result$index]
      min_radius_errors <- result$errors
      min_radius_radii <- result$radii
      min_radius_ix <- result$index
    } else {
      # 𝚿
      result <- compute_r_c_hat(U, tt, control$S, control$p)
      min_radius_int_error <- result$rc
      min_radius_errors    <- result$errors
      min_radius_radii     <- result$radii
      min_radius_ix        <- result$ix
    }
  }

  min_radius <- min_radius_int_error

  # cat(sprintf("  Integral Error min radius: %d\n", min_radius_int_error))

  radii <- control$radius_params * min_radius_int_error

  radii_filtered <- radii[radii < max_radius]

  if (length(radii_filtered) == 0) {
    radii_filtered <- max_radius
  }

  # cat(sprintf("  Radii [%s]\n", paste(radii_filtered, collapse = ", ")))

  V_ <- build_full_test_function_matrix(test_fun, tt, radii_filtered, order = 0)
  V_prime_ <- build_full_test_function_matrix(test_fun, tt, radii_filtered, order = 1)

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

  # cat(sprintf("  K Full: %d\n", k_full))
  # cat(sprintf(" Computing SVD"))

  SVD <- svd(V_)
  U_svd <- SVD$u
  singular_values <- SVD$d
  Vt <- t(SVD$v)

  sum_singular_values <- sum(singular_values)
  k_max <- min(k_full, mp1, control$k_max)

  info_numbers <- numeric(k_max)
  for (i in 2:k_max) {
    info_numbers[i] <- sum(singular_values[1:(i-1)]) / sum_singular_values
  }
  condition_numbers <- singular_values[1] / singular_values

  k1 <- find_last(condition_numbers, function(x) {
    x < control$max_test_fun_condition_number
  })

  k2 <- find_last(info_numbers, function(x) {
    x < control$min_test_fun_info_number
  })

  if (k1 == -1) k1 <- .Machine$integer.max
  if (k2 == -1) k2 <- .Machine$integer.max

  K <- min(k1, k2, k_max)

  # cat(sprintf("  Condition Number is now: %.4f\n", condition_numbers[K]))
  # cat(sprintf("  Info Number is now: %.4f\n", info_numbers[K]))
  # cat(sprintf("  K is: %d\n", K))

  V_final <- Vt[1:K, , drop = FALSE] * dt

  # cat("  Calculating Vprime\n")

  U_T <- t(U_svd)[1:K, , drop = FALSE]
  UV <- U_T %*% V_prime_
  inv_s <- 1.0 / singular_values[1:K]
  V_prime_final <- UV * inv_s * dt

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

  if (N <= 2) { return(1) }

  if (is.null(xx_in)) {
    x <- 1:N
  } else {
    x <- xx_in
  }

  x0 <- x[1]; y0 <- y[1]
  xM <- x[N]; yM <- y[N]

  E <- numeric(N - 2)

  for (k in 2:(N-1)) {
    xk <- x[k]
    yk <- y[k]

    slope1 <- (yk - y0) / (xk - x0)
    L1 <- function(x_val) slope1 * (x_val - x0) + y0

    sum1 <- 0.0
    for (m in 1:k) {
      err <- (L1(x[m]) - y[m]) / y[m]
      sum1 <- sum1 + err^2
    }

    slope2 <- (yM - yk) / (xM - xk)
    L2 <- function(x_val) slope2 * (x_val - xk) + yk

    sum2 <- 0.0
    for (m in k:N) {
      err <- (L2(x[m]) - y[m]) / y[m]
      sum2 <- sum2 + err^2
    }

    E[k - 1] <- sqrt(sum1 + sum2)
  }

  return(which.min(E) + 1)
}

endpoint_derivative_finite_difference_approx <- function(f, dt, l, mu){
  # Compute endpoints finite difference approximation of l-th derivative
  # with order mu accurate
  # f data: 1D vector 
  # dt: timestep
  # l: derivative order, l = 0, 1, ...
  # mu: order of accuracy mu = 1, 2, 3... 

  if(l == 0){
    D0 <- f[1]
    DT <- tail(f, n=1) 

  } else {
    # Forward difference weights
    A <- matrix(0, nrow = mu + l, ncol = mu + l)
    b <- rep(0, mu + l)
    
    for(i in seq(0, mu + l - 1)){
      for(m in seq(0, mu + l - 1)){
        A[(i+1),(m+1)] <- m^i
        if(i == l){
          b[i+1] = base::factorial(l)
        }
      }
    }
    
    fwd_diff_weights <- solve(A,b)
    bck_diff_weights <-  if(l %% 2 == 0) fwd_diff_weights else -rev(fwd_diff_weights)

    D0 <- as.numeric(fwd_diff_weights %*% f[1:(mu + l)] / dt^l)
    DT <- as.numeric(bck_diff_weights %*% tail(f, mu + l) / dt^l)
  }
    return(DT - D0)
}

compute_endpoint_derivatives <- function(U, dt, mus = c(1,2,1), S = 1){
  D <- ncol(U)
  
  endpoints <- lapply(1:D, function(d) {
    Ud <- U[,d]
    
    if(S == 0){
      c(endpoint_derivative_finite_difference_approx(Ud, dt, 0, mus[1]),
        endpoint_derivative_finite_difference_approx(Ud, dt, 1, mus[2]))
    } else {
      sapply(0:(2*S), function(l) {
        endpoint_derivative_finite_difference_approx(Ud, dt, l, mus[l+1])
      })
    }
  })

  return(endpoints)
}

buildI <- function(freqs, S, dt, T, endpoint_derivatives_d) {
  I <- complex(length(freqs))

  for(idx in seq_along(freqs)) {
    n <- freqs[idx] 
    In <- endpoint_derivatives_d[1] +
      sum(
        sapply(1:S, function(s) {
          inner_sum <- sum(
            sapply(0:(2*s), function(l) {
              choose(2*s, l) * ( (2 * pi * n * 1i) / T )^(2 * s - l) * endpoint_derivatives_d[l+1]
            })
          )
          bn <- bernoulli_numbers(2*s)
          last_bn <- bn[2 * s + 1, 1] / bn[2 * s + 1, 2]
          bernoulli_term <- last_bn / factorial(2*s) * dt^(2*s)
          
          return(inner_sum * bernoulli_term)
        })
      )
    
    I[idx] <- In
  }
  return(I)
}

fftfreq <- function(n, d = 1.0) {
  val <- 1.0 / (n * d)
  if (n %% 2 == 0) {
    freq <- c(0:(n/2 - 1), -(n/2):-1)
  } else {
    freq <- c(0:((n-1)/2), -((n-1)/2):-1)
  }
  return(freq * val)
}

# Change point in the integration error as a function of test function suppport
compute_r_c_hat <- function(U, tt, S, p){
  D <- ncol(U)
  mp1 <- nrow(U)
  dt <- mean(diff(tt))
  T <- as.numeric(tail(tt, n = 1) - tt[1])

  endpoint_derivatives <- compute_endpoint_derivatives(U, dt)

  M_tilde <- nrow(U) - 1
  freqs <- fftfreq(M_tilde, d = 1 / M_tilde) # Fourier frequencies
  Is <- lapply(seq(D), function(d){ buildI(freqs, S, dt, T, endpoint_derivatives[[d]]) })
  
  radii <- seq(2, floor(mp1 / 2))
  
  ehat <- numeric(length(radii))

  for(r_ix in seq(length(radii))){ 
    radius <- radii[r_ix]
    centers <- seq(radius + 1, mp1 - radius)
    K <- length(centers)

    C  <- 1 / (sqrt(sum(sapply(0:(2*p), function(k){
          choose(2*p, k) * (-1)^(k) / (2*k +1)
    }))) * (dt * radius)^(2 * p) * sqrt(2 * dt * radius))

    psi_hat_eval <- psi_hat(freqs, radius, dt, T, C, p = 16) # 𝚿̂̂
    
    e_int_hat <- t(sapply(seq(D), function(d){
       fft(psi_hat_eval * Is[[d]])[centers] / sqrt(T)
    }))
    
    e_int_hat_norm <- norm(e_int_hat, "F") / sqrt(K)
    ehat[r_ix] <- e_int_hat_norm
  }

  ix <- get_corner_index(log(ehat))
  
  rc <- radii[ix]

  return(list(rc = rc, ehat = ehat, radii = radii, ix = ix))
}
