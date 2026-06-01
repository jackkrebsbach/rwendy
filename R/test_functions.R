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
  t <- sym_symbol("t")
  r <- radius * dt
  phi_sym <- test_function(t, r)

  phi_deriv_sym <- phi_sym
  for (i in seq_len(order)) {
    phi_deriv_sym <- sym_diff(phi_deriv_sym, t)
  }

  sym_lambdify(phi_deriv_sym, t)
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

# Build the boundary-layer rows for ONE side of ONE radius at a given derivative order.
# n_bl and step follow the estimate_u0 heuristic. Rows are placed so the peak of phi
# lands at column 1+(i-1)*step (left) or mp1-(i-1)*step (right); the support is
# clipped to [1, mp1]. Returned matrix is RAW (no trapezoid endpoint weight, no dt).
build_boundary_layer_block <- function(test_function, tt, radius, order = 0,
                                       side = c("left", "right"), n_bl = NULL) {
  side <- match.arg(side)
  dt   <- mean(diff(tt))
  mp1  <- length(tt)
  diameter <- 2 * radius + 1
  r <- dt * radius

  lin <- seq(-r, r, length.out = diameter)
  xx  <- lin[2:(diameter - 1)]
  norm_factor <- sqrt(sum(test_function(xx, r)^2))
  ph_d <- test_function_derivative(test_function, radius, dt, order)
  row_full <- c(0, ph_d(xx) / norm_factor, 0)   # length = diameter, zero at ends

  n_bl <- if (!is.null(n_bl)) max(1L, as.integer(n_bl)) else max(3L, as.integer(ceiling(radius / 4)))
  step <- max(1L, floor(radius / n_bl) - 2L)

  V <- matrix(0, nrow = n_bl, ncol = mp1)
  for (i in seq_len(n_bl)) {
    if (side == "left") {
      shift_start <- max(1L, radius + 1L - (i - 1L) * step)
      slice       <- row_full[shift_start:diameter]
      s           <- min(length(slice), mp1)
      V[i, 1:s]   <- slice[1:s]
    } else {
      shift_end <- min(diameter, radius + 1L + (i - 1L) * step)
      slice     <- row_full[1:shift_end]
      s         <- min(length(slice), mp1)
      V[i, (mp1 - s + 1L):mp1] <- slice[(length(slice) - s + 1L):length(slice)]
    }
  }
  V
}

find_min_radius_int_error <- function(U, tt, radius_min, radius_max, num_radii, sub_sample_rate = 2) {

  D   <- ncol(U)
  Mp1 <- nrow(U)
  T <- tt[Mp1] - tt[1]
  dt <- mean(diff(tt))
  
  step  <- max(1, ceiling((radius_max - radius_min) / num_radii))
  hi    <- max(radius_max - 1L, radius_min)  # guard against collapsed range (very small mp1)
  radii <- seq(radius_min, hi, by = step)
  # radii <- radii[1:min(50, length(radii))]
  
  errors  <- numeric(length(radii))

  IX <- floor(Mp1 / sub_sample_rate) 

  for (i in seq_along(radii)) {
    radius <- radii[i]
    V_r <- build_test_function_matrix(phi, tt, radius)
    K   <- nrow(V_r)

    GT <- do.call(rbind, lapply(seq_len(D), function(d){
      sweep(V_r, 2, U[, d], `*`)
    }))

    f_hat_G_imag <- dt / sqrt(T) * Im(mvfft(t(GT))[IX, ])
    
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

  # psi shape parameter p (= control$p): higher p -> sharper, more concentrated
  # test functions; lower p -> wider, flatter. Drives the SVD modes of the basis.
  p_shape <- control$p %||% 16
  psi_p   <- function(t, r) psi(t, r, p = p_shape)

  V   <- build_full_test_function_matrix(psi_p, tt, c(radius_c), order = 0)
  Vp  <- build_full_test_function_matrix(psi_p, tt, c(radius_c), order = 1)
  Vpp <- build_full_test_function_matrix(psi_p, tt, c(radius_c), order = 2)

  K_interior <- nrow(V)
  K_bl       <- 0L
  bl_phi_t1  <- NULL
  bl_phi_tM  <- NULL
  # SSL uses the "sum" convention (no dt in V), so the IBP boundary terms and
  # the EM correction (both absolute units) need to be scaled by 1/dt to match.
  bdry_scale <- 1.0 / dt

  if (isTRUE(control$include_boundary_layer)) {
    max_em_order <- 4L
    bl_per_order <- lapply(0:max_em_order, function(ord) {
      rbind(
        build_boundary_layer_block(psi_p, tt, radius_c, order = ord, side = "left",  n_bl = control$n_bl),
        build_boundary_layer_block(psi_p, tt, radius_c, order = ord, side = "right", n_bl = control$n_bl)
      )
    })
    names(bl_per_order) <- paste0("phi", 0:max_em_order)
    K_bl <- nrow(bl_per_order$phi0)

    bl_phi_t1 <- vapply(bl_per_order, function(B) B[, 1],   numeric(K_bl))
    bl_phi_tM <- vapply(bl_per_order, function(B) B[, mp1], numeric(K_bl))
    colnames(bl_phi_t1) <- colnames(bl_phi_tM) <- paste0("phi", 0:max_em_order)

    # Apply trapezoid endpoint weights (x0.5 at columns 1 and mp1) to BL rows
    # only. This is needed so that V_BL %*% F matches the trapezoid sum that
    # the Euler-Maclaurin correction is derived for; without it the residual
    # picks up an O(h) artifact from the missing endpoint half-weight. Interior
    # SSL phi vanishes at the endpoints, so the same operation on interior rows
    # would be a no-op (we skip it).
    apply_trap_weights <- function(B) {
      B[, 1]   <- B[, 1]   * 0.5
      B[, mp1] <- B[, mp1] * 0.5
      B
    }

    # Append BL rows in the SSL convention (no dt scaling on the integration
    # block). The boundary terms and EM correction are still scaled by 1/dt via
    # bdry_scale below to compensate for SSL's missing h.
    V   <- rbind(V,   apply_trap_weights(bl_per_order$phi0))
    Vp  <- rbind(Vp,  apply_trap_weights(bl_per_order$phi1))
    Vpp <- rbind(Vpp, apply_trap_weights(bl_per_order$phi2))
  }

  return(list(
    V = V, V_prime = Vp, V_pp = Vpp, radius_c = radius_c, rc_errors = rc_errors, rc_radii = rc_radii,
    K_interior = K_interior, K_bl = K_bl,
    bl_phi_t1 = bl_phi_t1, bl_phi_tM = bl_phi_tM,
    bdry_scale = bdry_scale
  ))
}

# build_full_test_function_matrices_ssl <- function(U, tt, control) {

#   dt <- mean(diff(tt))
#   mp1 <- nrow(U)

#   if (!is.null(control$fixed_radius)) {
#     radius_c <- control$fixed_radius
#     rc_errors <- NULL
#     rc_radii <- NULL
#   } else {
#     data <- compute_r_c_hat(U, tt, control$S, control$p)
#     radius_c <- data$rc
#     rc_errors <- data$ehat
#     rc_radii <- data$radii
#   }

#   V  <- build_full_test_function_matrix(psi, tt, c(radius_c), order = 0)
#   Vp <- build_full_test_function_matrix(psi, tt, c(radius_c), order = 1)

#   # Integral convention: V %*% F approximates ∫ phi(t) F(t) dt directly.
#   V  <- V  * dt
#   Vp <- Vp * dt

#   list(V = V, V_prime = Vp, radius_c = radius_c,
#        rc_errors = rc_errors, rc_radii = rc_radii)
# }

build_full_test_function_matrices_msg <- function(U, tt, control, compute_svd = TRUE) {
  # cat("<< Building test matrices >>\n")
  dt <- mean(diff(tt))
  mp1 <- nrow(U)

  # Radius selection is data-driven. The test-function matrices themselves are
  # grid-based.
  U_rad <- U

  radii <- control$radius_params
  radius_min_time <- control$radius_min_time
  radius_max_time <- control$radius_max_time

  test_fun <- get(control$test_fun)

  min_radius <- max(ceiling(radius_min_time / dt), 2)
  max_radius <- floor(radius_max_time / dt / 2)

  max_radius_for_interior <- floor((mp1 - 2) / 2)
  if (max_radius > max_radius_for_interior) {
    max_radius <- floor(max_radius_for_interior / 2)
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
      result <- find_min_radius_int_error(U_rad, tt, min_radius, max_radius, num_radii = 100, sub_sample_rate = 2)
      min_radius_int_error <- result$radii[result$index]
      min_radius_errors <- result$errors
      min_radius_radii <- result$radii
      min_radius_ix <- result$index
    } else {
      # 𝚿
      result <- compute_r_c_hat(U_rad, tt, control$S, control$p)
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

  # cat("Computing SVD")

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

compute_last_bernoulli_num <- function(n) {
  B <- numeric(n + 1)
  B[1] <- 1
  for (m in 1:n) {
    B[m + 1] <- -1 /
      (m + 1) *
      sum(sapply(0:(m - 1), function(k) choose(m + 1, k) * B[k + 1]))
  }
  return(B[n + 1])
}

buildI <- function(freqs, S, dt, T, endpoint_derivatives_d) {
  I <- complex(length(freqs))

  bernoulli_terms <- vapply(1:S, function(s) {
    compute_last_bernoulli_num(2*s) / factorial(2*s) * dt^(2*s)
  }, numeric(1))

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

          return(inner_sum * bernoulli_terms[s])
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
