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

# Boundary-layer augmentation shared by the SSL and MSG builders.
#
# Builds the left+right BL blocks at orders 0..4 for a single radius, applies the
# trapezoid endpoint weight (x0.5 at columns 1 and mp1) to the orders 0..2 blocks
# that get appended as rows of V / V_prime / V_pp, and extracts the RAW (un-trap-
# weighted) phi^(0..4) at the two boundaries that the integration-by-parts
# boundary term (phi0) and the Euler-Maclaurin correction (phi0..4) consume.
#
# Both builders append these rows AFTER their interior basis. The interior rows
# vanish at the endpoint columns -- raw SSL rows by construction, and MSG SVD
# modes because they are linear combinations of such rows -- so the appended BL
# block is the only carrier of boundary content, and it is never subject to the
# MSG SVD truncation. The downstream consumer (build_wendy_problem) folds the
# boundary term into V_prime (-> b and L0) and adds the EM defect to g / Jp_r.
#
# n_bl: per-side BL count. Defaults to control$n_bl (NULL -> build_boundary_layer_block's
# own radius-based heuristic). The MSG builder passes an explicit value capped against
# K_interior so the boundary layer cannot dominate the interior basis.
#
# Returns: V0/V1/V2 (trap-weighted BL rows for orders 0/1/2), K_bl, and the
# (K_bl x 5) raw endpoint matrices bl_phi_t1 / bl_phi_tM with columns phi0..phi4.
build_boundary_layer_augmentation <- function(test_function, tt, radius, control, n_bl = control$n_bl) {
  mp1 <- length(tt)
  max_em_order <- 4L

  bl_per_order <- lapply(0:max_em_order, function(ord) {
    rbind(
      build_boundary_layer_block(test_function, tt, radius, order = ord, side = "left",  n_bl = n_bl),
      build_boundary_layer_block(test_function, tt, radius, order = ord, side = "right", n_bl = n_bl)
    )
  })
  names(bl_per_order) <- paste0("phi", 0:max_em_order)
  K_bl <- nrow(bl_per_order$phi0)

  bl_phi_t1 <- vapply(bl_per_order, function(B) B[, 1],   numeric(K_bl))
  bl_phi_tM <- vapply(bl_per_order, function(B) B[, mp1], numeric(K_bl))
  colnames(bl_phi_t1) <- colnames(bl_phi_tM) <- paste0("phi", 0:max_em_order)

  apply_trap_weights <- function(B) {
    B[, 1]   <- B[, 1]   * 0.5
    B[, mp1] <- B[, mp1] * 0.5
    B
  }

  list(
    V0 = apply_trap_weights(bl_per_order$phi0),
    V1 = apply_trap_weights(bl_per_order$phi1),
    V2 = apply_trap_weights(bl_per_order$phi2),
    K_bl = K_bl,
    bl_phi_t1 = bl_phi_t1,
    bl_phi_tM = bl_phi_tM
  )
}

find_min_radius_int_error <- function(U, tt, radius_min, radius_max, num_radii, sub_sample_rate = 2) {

  D   <- ncol(U)
  Mp1 <- nrow(U)
  T <- tt[Mp1] - tt[1]
  dt <- mean(diff(tt))
  
  step  <- max(1, ceiling((radius_max - radius_min) / num_radii))
  hi    <- max(radius_max - 1L, radius_min)  # guard against collapsed range (very small mp1)
  radii <- seq(radius_min, hi, by = step)
  radii <- radii[1:min(50, length(radii))]
  
  errors  <- numeric(length(radii))

  IX <- floor(Mp1 / sub_sample_rate)

  # Only the single Fourier coefficient at frequency (IX-1), and only its
  # imaginary part, is used. For real GT that coefficient is
  # GT %*% Im(exp(-2i*pi*(IX-1)*m/Mp1)), so the full per-radius mvfft collapses to
  # one real matmul. Folding in the per-state column scaling of GT, the imaginary
  # coefficient for state d is V_r %*% (U[,d] * sinvec), i.e. V_r %*% M below.
  sinvec <- -sin(2 * pi * (IX - 1) * (0:(Mp1 - 1)) / Mp1)
  M      <- U * sinvec                                   # Mp1 x D (sinvec recycled per column)
  scale  <- (4 * pi * IX / Mp1) * (dt / sqrt(T))

  for (i in seq_along(radii)) {
    V_r <- build_test_function_matrix(phi, tt, radii[i])
    VM  <- V_r %*% M                                     # K x D, block d = Im Fourier coef for state d
    errors[i] <- scale * sqrt(sum(VM^2) / nrow(V_r))
  }

  log_errors <- log(errors)
  ix <- get_corner_index(log_errors)

  return(list(index = ix, errors = errors, radii = radii))
}

build_full_test_function_matrices_ssl <- function(U, tt, control) {

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
  bl_phi_t1  <- NULL   # (K_bl x 5) raw phi^(0..4) at t_1, per BL test function
  bl_phi_tM  <- NULL   # (K_bl x 5) raw phi^(0..4) at t_M

  if (isTRUE(control$include_boundary_layer)) {
    # Build + trap-weight the BL rows and extract the raw endpoint derivatives
    # the boundary term / EM correction need (see build_boundary_layer_augmentation).
    bl        <- build_boundary_layer_augmentation(psi_p, tt, radius_c, control)
    K_bl      <- bl$K_bl
    bl_phi_t1 <- bl$bl_phi_t1
    bl_phi_tM <- bl$bl_phi_tM

    # Append the boundary-layer rows. Like the interior rows these are in the
    # raw (sum) convention; the caller applies the dt quadrature weight uniformly.
    V   <- rbind(V,   bl$V0)
    Vp  <- rbind(Vp,  bl$V1)
    Vpp <- rbind(Vpp, bl$V2)
  }

  return(list(
    V = V, V_prime = Vp, V_pp = Vpp, radius_c = radius_c,
    rc_errors = rc_errors, rc_radii = rc_radii,
    K_interior = K_interior, K_bl = K_bl,
    bl_phi_t1 = bl_phi_t1, bl_phi_tM = bl_phi_tM
  ))
}


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

  bl_radius <- max(radii_filtered)

  make_bl <- function() {
    if (!isTRUE(control$include_boundary_layer)) return(NULL)
    n_side <- control$n_bl
    if (is.null(n_side)) n_side <- max(3L, as.integer(ceiling(bl_radius / 4)))
    build_boundary_layer_augmentation(test_fun, tt, bl_radius, control, n_bl = n_side)
  }

  if (!compute_svd) {
    bl     <- make_bl()
    V_out  <- V_
    Vp_out <- V_prime_
    if (!is.null(bl)) {
      V_out  <- rbind(V_out,  bl$V0)
      Vp_out <- rbind(Vp_out, bl$V1)
    }
    return(list(
      V = V_out,
      V_prime = Vp_out,
      radii = radii_filtered,
      min_radius_errors = min_radius_errors,
      min_radius_radii = min_radius_radii,
      min_radius_ix = min_radius_ix,
      min_radius = min_radius,
      K_interior = nrow(V_),
      K_bl       = if (!is.null(bl)) bl$K_bl      else 0L,
      bl_phi_t1  = if (!is.null(bl)) bl$bl_phi_t1 else NULL,
      bl_phi_tM  = if (!is.null(bl)) bl$bl_phi_tM else NULL
    ))
  }

  bl           <- make_bl()
  V_pool       <- if (!is.null(bl)) rbind(V_,       bl$V0) else V_
  V_prime_pool <- if (!is.null(bl)) rbind(V_prime_, bl$V1) else V_prime_
  n_int        <- nrow(V_)
  n_bl_pool    <- if (!is.null(bl)) bl$K_bl else 0L
  bl_pool_idx  <- if (n_bl_pool > 0L) (n_int + 1L):(n_int + n_bl_pool) else integer(0)

  k_full <- nrow(V_pool)

  # Truncated SVD via the smaller Gram matrix is about 1.5x faster than a full svd()
  left_gram <- (k_full <= mp1)
  eg <- eigen(if (left_gram) tcrossprod(V_pool) else crossprod(V_pool), symmetric = TRUE)
  singular_values <- sqrt(pmax(eg$values, 0))

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

  K0 <- min(k1, k2, k_max)

  # Truncation is purely spectral: K is set by the energy / condition-number
  # cutoffs above. Whatever boundary content survives in the top-K modes is what
  # the boundary correction acts on; there is no special-casing to retain a
  # "boundary mode", consistent with boundary content not partitioning cleanly.
  K <- K0
  inv_s <- 1.0 / singular_values[1:K]

  if (left_gram) {
    U_K     <- eg$vectors[, 1:K, drop = FALSE]        # top-K left singular vectors
    V_final <- t(crossprod(V_pool, U_K)) * inv_s      # right vectors v_k = V^T u_k / s_k (rows)
    U_T     <- t(U_K)
  } else {
    V_K     <- eg$vectors[, 1:K, drop = FALSE]        # top-K right singular vectors
    V_final <- t(V_K)
    U_T     <- t(V_pool %*% V_K) * inv_s              # left vectors u_k = V v_k / s_k (rows)
  }
  UV <- U_T %*% V_prime_pool
  V_prime_final <- UV * inv_s
  K_bl            <- 0L
  bl_phi_t1       <- NULL
  bl_phi_tM       <- NULL
  if (n_bl_pool > 0L) {
    coef_bl   <- U_T[, bl_pool_idx, drop = FALSE] * inv_s     # K x n_bl_pool
    bl_phi_t1 <- coef_bl %*% bl$bl_phi_t1                     # K x 5 (phi0..phi4)
    bl_phi_tM <- coef_bl %*% bl$bl_phi_tM
    colnames(bl_phi_t1) <- colnames(bl_phi_tM) <- colnames(bl$bl_phi_t1)
    K_bl <- K
  }

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
    K = K,
    # BL rows are mixed into the SVD pool, so there is no separate interior block:
    # all K modes carry the exact boundary correction (K_interior = 0, K_bl = K),
    # and build_wendy_problem's bl_rows = K_interior + 1:K_bl resolves to 1:K.
    K_interior      = 0L,
    K_bl            = K_bl,
    bl_phi_t1       = bl_phi_t1,
    bl_phi_tM       = bl_phi_tM
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
