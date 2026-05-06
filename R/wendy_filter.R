
# estimate_u0 <- function(U, f_, J_u, J_t, tt, p, test_function_params){
#   # Dimension of system
#   mp1 <- nrow(U)
#   D <- ncol(U)
#   J <- length(p)

#   # Build boundary layer test functions 
#   data <- compute_r_c_hat(U, tt, test_function_params$S, test_function_params$p)
#   radius_c <- data$rc
#   diameter <- 2 * radius_c + 1
#   dt <- mean(diff(tt))
#   r <- dt * radius_c
#   lin <- seq(-r, r, length.out = diameter)
#   xx <- lin[2:(diameter-1)] # Don't include the endpoints (outside of support for 𝚽 or zero for 𝚿) 

#   ph <- test_function_derivative(psi, radius_c, dt, order = 0)
#   php <- test_function_derivative(psi, radius_c, dt, order = 1)
#   phpp <- test_function_derivative(psi, radius_c, dt, order = 2)
#   phppp <- test_function_derivative(psi, radius_c, dt, order = 3)
#   phpppp <- test_function_derivative(psi, radius_c, dt, order = 4)

#   phxx <- ph(xx)
#   phpxx <- php(xx)
#   phppxx <- phpp(xx)
#   phpppxx <- phppp(xx)
#   phppppxx <- phpppp(xx)

#   scale <- sqrt(sum(phxx^2))

#   p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
#   input <- rbind(p_mat, t(U), matrix(tt, nrow = 1L))
#   Fp <- f_(input)
  
#   n_bl <- max(1L, floor(radius_c / 5))

#   v_row <- c(0, phxx / scale, 0)
#   vp_row <- c(0, phpxx / scale, 0)
#   vpp_row <- c(0, phppxx / scale, 0)
#   vppp_row <- c(0, phpppxx / scale, 0)
#   vpppp_row <- c(0, phppppxx / scale, 0)

#   step <- max(1L, floor(radius_c / n_bl) - 2)

#   V <- matrix(0, nrow = n_bl, ncol = mp1) 
#   Vp <- matrix(0, nrow = n_bl, ncol = mp1)
#   Vpp <- matrix(0, nrow = n_bl, ncol = mp1)
#   Vppp <- matrix(0, nrow = n_bl, ncol = mp1)
#   Vpppp <- matrix(0, nrow = n_bl, ncol = mp1)

#   for (i in seq(nrow(V))) {
#     shift <- radius_c + 1 - (i - 1) * step
    
#     l <- length(v_row)

#     v_trimmed <- v_row[shift:l]
#     vp_trimmed <- vp_row[shift:l]
#     vpp_trimmed <- vpp_row[shift:l]
#     vppp_trimmed <- vppp_row[shift:l]
#     vpppp_trimmed <- vpppp_row[shift:l]

#     s <- length(v_trimmed)

#     if (s > 0) {
#       V[i, 1:s] <- v_trimmed
#       Vp[i, 1:s] <- vp_trimmed
#       Vpp[i, 1:s] <- vpp_trimmed
#       Vppp[i, 1:s] <- vppp_trimmed
#       Vpppp[i, 1:s] <- vpppp_trimmed
#     }
#   }

#   # Solve uncorrected system -ϕ(0)u(0) = ∫g(t)dt for u0 and k boundary layer functions
#   B <- V[, 1] # Boundary layer test functions
#   Q <- diag(c(0.5, rep(1, mp1 - 2), 0.5) * dt) # Quadrature matrix for Trapezoidal Rule
#   residual <- -(V %*% Q %*% Fp  + Vp %*% Q %*% U) # -(Weak residual)
#   u00 <- solve(crossprod(B), t(B) %*% residual) # Least squares solution to over determined system

#   # Iterate with the Euler-Maclaurin formula
#   un0 <- u00
#   t0 <- tt[1]

#   for (iter in seq_len(3)) {
#     input  <- c(p, un0, t0)
#     f_val  <- matrix(f_(input), ncol = 1)
#     Ju_val <- t(matrix(as.vector(J_u(input)), nrow = D, ncol = D))
#     Jt_val <- as.vector(J_t(input))
    
#     # g(t) := φ(t)f(p,u,t) + φ′(t)u(t)
#     # g'(t) = 2φ′(t)f(p,u,t) + φ(t)( ∇ᵤf(p,u,t) ⋅ f(p,u,t) + ∇ₜf(p,u,t)) + φ′′(t)u(t) 
#     gprime0 <- matrix(t(vapply(seq(length(B)), function(row){
#           return(2 * Vp[row, 1] * f_val + V[row, 1] * (as.vector(Ju_val %*% f_val) + Jt_val) + Vpp[row,1] * as.vector(un0))
#     }, numeric(D))), nrow = length(B), ncol = D)

#     # EM Correction is h^2 /12 (g'(t0) - g'(tf)) but g'(tf) = 0 to estimate initial condition
#     un0 <- solve(crossprod(B), t(B) %*% (residual - dt^2 / 12  * gprime0)) # Least squares solution to over determined system
#   }

#   return(as.numeric(un0))
# }

#' Estimate the initial condition using Euler-MacLaurin Trapezoid correction on the function g(t) := φ(t)f(p,u,t) + φ′(t)u(t)
#' @export
estimate_u0 <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, test_function_params){
  mp1 <- nrow(U)
  D <- ncol(U)
  J <- length(p)

  data <- compute_r_c_hat(U, tt, test_function_params$S, test_function_params$p)
  radius_c <- data$rc
  diameter <- 2 * radius_c + 1
  dt <- mean(diff(tt))
  r <- dt * radius_c
  lin <- seq(-r, r, length.out = diameter)
  xx <- lin[2:(diameter-1)]

  ph <- test_function_derivative(psi, radius_c, dt, order = 0)
  php <- test_function_derivative(psi, radius_c, dt, order = 1)
  phpp <- test_function_derivative(psi, radius_c, dt, order = 2)
  phppp <- test_function_derivative(psi, radius_c, dt, order = 3)
  phpppp <- test_function_derivative(psi, radius_c, dt, order = 4)

  phxx     <- ph(xx)
  scale    <- sqrt(sum(phxx^2))

  # Store scaled phi derivatives as a named list — extensible to any order
  phi_rows <- list(
    phi0 = c(0, phxx               / scale, 0),
    phi1 = c(0, php(xx)            / scale, 0),
    phi2 = c(0, phpp(xx)           / scale, 0),
    phi3 = c(0, phppp(xx)          / scale, 0),
    phi4 = c(0, phpppp(xx)         / scale, 0)
  )

  n_bl <- max(1L, floor(radius_c / 5))
  step <- max(1L, floor(radius_c / n_bl) - 2)
  l <- length(phi_rows$phi0)

  # Build test function derivative matrices
  make_V <- function(row_vals) {
    V <- matrix(0, nrow = n_bl, ncol = mp1)
    for (i in seq_len(n_bl)) {
      shift     <- radius_c + 1 - (i - 1) * step
      trimmed   <- row_vals[shift:l]
      s         <- length(trimmed)
      if (s > 0) V[i, 1:s] <- trimmed
    }
    V
  }

  Vs <- lapply(phi_rows, make_V)

  B        <- Vs$phi0[, 1]
  Q        <- diag(c(0.5, rep(1, mp1 - 2), 0.5) * dt)
  p_mat    <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
  input    <- rbind(p_mat, t(U), matrix(tt, nrow = 1L))
  Fp       <- f_(input)
  residual <- -(Vs$phi0 %*% Q %*% Fp + Vs$phi1 %*% Q %*% U)
  u00      <- solve(crossprod(B), t(B) %*% residual)

  # g derivative at t0 — takes named list of phi scalars and f-derivative quantities
  # g^(n) = Σ_{k=0}^{n} C(n,k) φ^(k+1) f^(n-k)  +  φ^(n+1) u (from product rule on φf, plus derivative of φ′u term)
  g_deriv_at_t0 <- function(phi_scalars, f_derivs, u0_vec, order) {
    # phi_scalars[[k+1]] = φ^(k);  f_derivs[[m+1]] = d^m f/dt^m
    # g = φF + φ'u, to calcualte g^(n) use Leibniz on each term:
    #   d^n[φF]/dt^n = Σ_{k=0}^n C(n,k) φ^(k) F^(n-k)
    #   d^n[φ'u]/dt^n = φ^(n+1)u + Σ_{k=0}^{n-1} C(n,k) φ^(k+1) F^(n-k-1) (using u^(0)=u, u^(m)=F^(m-1) for m≥1)
    n  <- order
    result <- phi_scalars[[n + 2]] * u0_vec          # φ^(n+1)·u
    for (k in seq(0, n)) {
      result <- result + choose(n, k) * phi_scalars[[k + 1]] * f_derivs[[n - k + 1]] # d^n[φF]/dt^n
    }
    if (n >= 1) { 
      for (k in seq(0, n - 1)){
        result <- result + choose(n, k) * phi_scalars[[k + 2]] * f_derivs[[n - k]] # remaining d^n[φ'u]/dt^n terms
      }
    }
    result
  }

  un0 <- u00
  t0  <- tt[1]

  for (iter in seq_len(3)) {
    input <- c(p, un0, t0)

    f_derivs <- list(
      as.vector(f_(input)),
      as.vector(dF_dt_(input)),
      as.vector(d2F_dt2_(input)),
      as.vector(d3F_dt3_(input))
    )

    make_g_deriv <- function(ord)
      matrix(t(vapply(seq_along(B), function(row) {
        phi_scalars <- lapply(Vs, function(V) V[row, 1])
        g_deriv_at_t0(phi_scalars, f_derivs, as.vector(un0), order = ord)
      }, numeric(D))), nrow = length(B), ncol = D)

    gprime0  <- make_g_deriv(1)
    g3prime0 <- make_g_deriv(3)

    un0 <- solve(crossprod(B), t(B) %*% (residual - dt^2/12 * gprime0 + dt^4/720 * g3prime0))
  }

  return(as.numeric(un0))
}

# Evaluate ψ and its first four scaled derivatives at interior points of [-r, r].
build_phi_rows <- function(radius_c, dt) {
  diameter <- 2L * radius_c + 1L
  r <- dt * radius_c
  xx <- seq(-r, r, length.out = diameter)[2:(diameter - 1L)]

  ph  <- test_function_derivative(psi, radius_c, dt, order = 0L)
  phxx  <- ph(xx)
  scale <- sqrt(sum(phxx^2))

  make_row <- function(ord){
    c(0, test_function_derivative(psi, radius_c, dt, order = ord)(xx) / scale, 0)
  }

  list(
    phi0 = c(0, phxx / scale, 0),
    phi1 = make_row(1L),
    phi2 = make_row(2L),
    phi3 = make_row(3L),
    phi4 = make_row(4L)
  )
}

# Build boundary-layer matrices anchored at column k.
#   "right": support [k, k+r], φ(k) ≠ 0, φ(k+r) = 0  — right half of ψ
#   "left" : support [k-r, k], φ(k) ≠ 0, φ(k-r) = 0  — left  half of ψ
build_bl_matrices <- function(phi_rows, k, mp1, n_bl, step, orientation) {
  radius_c <- (length(phi_rows$phi0) - 1L) %/% 2L

  lapply(phi_rows, function(row_vals) {
    V <- matrix(0, nrow = n_bl, ncol = mp1)
    for (i in seq_len(n_bl)) {
      if (orientation == "right") {
        shift   <- radius_c + 1L - (i - 1L) * step
        trimmed <- row_vals[shift:length(row_vals)]
        cols    <- k:(k + length(trimmed) - 1L)
        ok      <- cols <= mp1
        if (any(ok)) V[i, cols[ok]] <- trimmed[ok]
      } else {
        end     <- radius_c + 1L - (i - 1L) * step
        trimmed <- row_vals[1L:end]
        cols    <- (k - length(trimmed) + 1L):k
        ok      <- cols >= 1L
        if (any(ok)) V[i, cols[ok]] <- trimmed[ok]
      }
    }
    V
  })
}

# g^(n)(t_k): Leibniz rule on g = φF + φ'u along the trajectory.
#   phi_scalars[[j+1]] = φ^(j)(t_k),  f_derivs[[m+1]] = d^m f/dt^m at (p, u_k, t_k)
g_deriv_at_tk <- function(phi_scalars, f_derivs, uk_vec, order) {
  n      <- order
  result <- phi_scalars[[n + 2L]] * uk_vec
  for (j in seq(0L, n))
    result <- result + choose(n, j) * phi_scalars[[j + 1L]] * f_derivs[[n - j + 1L]]
  if (n >= 1L)
    for (j in seq(0L, n - 1L))
      result <- result + choose(n, j) * phi_scalars[[j + 2L]] * f_derivs[[n - j]]
  result
}

# EM-corrected LS solve for u_k, iterating 3 times.
em_corrected_solve <- function(B, residual, Vs, k, f_, dF_dt_, d2F_dt2_, d3F_dt3_, p, tk, u_init, dt, D) {
  un <- u_init
  for (iter in seq_len(3L)) {
    input    <- c(p, as.vector(un), tk)
    f_derivs <- list(
      as.vector(f_(input)),
      as.vector(dF_dt_(input)),
      as.vector(d2F_dt2_(input)),
      as.vector(d3F_dt3_(input))
    )
    make_gn <- function(ord)
      matrix(t(vapply(seq_along(B), function(row) {
        phi_scalars <- lapply(Vs, function(V) V[row, k])
        g_deriv_at_tk(phi_scalars, f_derivs, as.vector(un), order = ord)
      }, numeric(D))), nrow = length(B), ncol = D)

    un <- solve(crossprod(B), t(B) %*% (residual - dt^2/12 * make_gn(1L) + dt^4/720 * make_gn(3L)))
  }
  as.numeric(un)
}

#' Estimate u(t_k) using boundary-layer test functions + Euler-MacLaurin correction.
#'
#' Orientation is chosen automatically:
#'   right (k + radius_c <= mp1): support to the right of k, φ(k+r) = 0
#'   left  (k + radius_c  > mp1): support to the left  of k, φ(k-r) = 0
#' The EM correction formula is identical for both orientations; only the
#' residual sign flips.
#' @export
estimate_u_k <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p,
                          test_function_params, k) {
  tt   <- as.vector(tt)
  mp1  <- nrow(U)
  D    <- ncol(U)
  J    <- length(p)
  dt   <- mean(diff(tt))

  data     <- compute_r_c_hat(U, tt, test_function_params$S, test_function_params$p)
  radius_c <- data$rc
  # Determine orientation with the unclamped radius, then clamp only for the
  # relevant side.  Clamping to both (mp1-k) AND (k-1) is wrong: for k near
  # the start (right-oriented), k-1 could be tiny even though there is ample
  # room to the right, producing a useless 1-point support.
  orientation <- if (k + radius_c <= mp1) "right" else "left"
  radius_c    <- if (orientation == "right") min(radius_c, mp1 - k) else min(radius_c, k - 1L)
  radius_c    <- max(radius_c, 1L)

  phi_rows <- build_phi_rows(radius_c, dt)
  n_bl     <- max(1L, floor(radius_c / 5L))
  step     <- max(1L, floor(radius_c / n_bl) - 2L)

  Vs <- build_bl_matrices(phi_rows, k, mp1, n_bl, step, orientation)
  B  <- Vs$phi0[, k]

  # Support-local Q: half-weight at the anchor k (the only endpoint where φ ≠ 0),
  # full weight elsewhere.  The global Q is wrong for interior k because it gives
  # dt instead of 0.5*dt at the anchor, producing an O(h) error that swamps the
  # EM corrections.
  q_weights    <- rep(dt, mp1)
  q_weights[k] <- 0.5 * dt
  Q            <- diag(q_weights)

  p_mat <- matrix(rep(p, mp1), nrow = J, ncol = mp1)
  Fp    <- f_(rbind(p_mat, t(U), matrix(tt, nrow = 1L)))

  sign     <- if (orientation == "right") -1 else +1
  residual <- sign * (Vs$phi0 %*% Q %*% Fp + Vs$phi1 %*% Q %*% U)
  u_init   <- solve(crossprod(B), t(B) %*% residual)

  em_corrected_solve(B, residual, Vs, k, f_, dF_dt_, d2F_dt2_, d3F_dt3_,
                     p, tt[k], u_init, dt, D)
}

#' Estimate u(t_k) for a sequence of indices and return a corrected trajectory.
#'
#' @param indices Integer vector of row indices to estimate (default: all rows).
#' @export
step_through <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p,
                          test_function_params, indices = NULL) {
  if (is.null(indices)) indices <- seq_len(nrow(U))
  D   <- ncol(U)
  raw <- vapply(indices, function(k)
    estimate_u_k(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, test_function_params, k),
    numeric(D)
  )
  if (D == 1L) matrix(raw, ncol = 1L) else t(raw)
}

estimate_U_star <- function(U, f_, J_u, J_t, tt, p, test_function_params, sigma = NULL) {
  tt   <- as.vector(tt)
  mp1  <- nrow(U)
  D    <- ncol(U)
  I_D  <- diag(D)

  # Noise covariance 
  noise_sd <- as.numeric(if (!is.null(sigma)) sigma else estimate_std(U, k = 6))
  noise_sd <- mean(noise_sd)   # collapse to scalar if estimate_std returns per-column vector
  R_mat    <- noise_sd^2 * I_D

  # Process noise: small so the filter trusts the ODE model (maybe based on the covariance estiamtes of p^\star)
  Q_mat <- (noise_sd * 1e-3)^2 * I_D

  u0 <- U[1,]
  P0 <- noise_sd^2 * I_D

  # Fuse u0 prior with first observation
  K0 <- P0 %*% solve(P0 + R_mat)
  u_filt <- matrix(0, mp1, D)
  u_filt[1, ] <- u0 + K0 %*% (U[1, ] - u0)
  P_filt <- array(0, c(mp1, D, D))
  P_filt[1,,] <- (I_D - K0) %*% P0

  # Storage for smoother
  u_pred <- matrix(0, mp1, D)
  P_pred <- array(0,  c(mp1, D, D))
  F_store <- array(0,  c(mp1 - 1L, D, D))  # linearised transition Jacobians

  # Forward Pass Extend Kalman Filter 
  for (k in seq_len(mp1 - 1L)) {
    dt_k <- tt[k + 1L] - tt[k]
    uk   <- u_filt[k, ]

    # RK4 state prediction
    k1 <- f_(c(p, uk,                           tt[k]               ))
    k2 <- f_(c(p, uk + 0.5 * dt_k * k1,        tt[k] + 0.5 * dt_k ))
    k3 <- f_(c(p, uk + 0.5 * dt_k * k2,        tt[k] + 0.5 * dt_k ))
    k4 <- f_(c(p, uk +       dt_k * k3,        tt[k] +       dt_k ))
    u_pred[k + 1L, ] <- uk + (dt_k / 6) * (k1 + 2*k2 + 2*k3 + k4)

    # Linearised Covariance Propogation: F ≈ I + dt * J_u  (first-order)
    Ju_k <- t(matrix(as.vector(J_u(c(p, uk, tt[k]))), D, D))
    Fk  <- I_D + dt_k * Ju_k
    F_store[k,,] <- Fk

    Pk_pred <- Fk %*% P_filt[k,,] %*% t(Fk) + Q_mat
    P_pred[k + 1L,,] <- Pk_pred

    # Update with observation y_{k+1}
    S_k              <- Pk_pred + R_mat          # innovation covariance (H = I)
    Kk               <- Pk_pred %*% solve(S_k)
    u_filt[k + 1L, ] <- u_pred[k + 1L, ] + Kk %*% (U[k + 1L, ] - u_pred[k + 1L, ])
    P_filt[k + 1L,,] <- (I_D - Kk) %*% Pk_pred
  }

  # Backward pass RTS smoother
  u_smooth      <- matrix(0, mp1, D)
  P_smooth      <- array(0,  c(mp1, D, D))
  u_smooth[mp1, ]  <- u_filt[mp1, ]
  P_smooth[mp1,,]  <- P_filt[mp1,,]

  for (k in seq(mp1 - 1L, 1L)) {
    Pk   <- P_filt[k,,]
    Pp   <- P_pred[k + 1L,,]
    Fk   <- F_store[k,,]

    # Smoother gain
    Gk <- Pk %*% t(Fk) %*% solve(Pp)
    u_smooth[k, ] <- u_filt[k, ] + Gk %*% (u_smooth[k + 1L, ] - u_pred[k + 1L, ])
    P_smooth[k,,] <- Pk           + Gk %*% (P_smooth[k + 1L,,] - Pp) %*% t(Gk)
  }

  list(
    U_star   = u_smooth,
    P_smooth = P_smooth,
    u_filt   = u_filt,
    P_filt   = P_filt,
    u_pred   = u_pred,
    P_pred   = P_pred
  )
}

#' EKF + RTS smoother using the wendy step_through estimator as the predictor.
#'
#' Two-pass approach to avoid double-counting the noisy observations:
#'   Pass 1 — RK4-EKF + RTS smoother on the raw data; the smoothed trajectory
#'             (bidirectional) is used as the wendy input in Pass 2.
#'   Pass 2 — EKF whose prediction at each step comes from \code{estimate_u_k}
#'             applied to the Pass-1 smoothed trajectory.  The Kalman update
#'             still fuses against the original observations.
#'
#' Q_mat for Pass 2 is set to \code{noise_sd^2 * I_D} (same order as R_mat) so
#' the filter does not blindly trust wendy predictions that are degraded near
#' the boundaries or by imperfect parameters.
#'
#' @export
estimate_U_star_wendy <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_,
                                  J_u, tt, p, test_function_params,
                                  sigma = NULL) {
  tt   <- as.vector(tt)
  mp1  <- nrow(U)
  D    <- ncol(U)
  I_D  <- diag(D)

  noise_sd <- as.numeric(if (!is.null(sigma)) sigma else estimate_std(U, k = 6))
  noise_sd <- mean(noise_sd)
  R_mat    <- noise_sd^2 * I_D
  Q_mat1   <- (noise_sd * 1e-3)^2 * I_D   # tight process noise for RK4 pass
  Q_mat2   <- noise_sd^2 * I_D             # wendy prediction error ~ measurement noise

  .init_filter <- function() {
    u0      <- U[1, ]
    P0      <- noise_sd^2 * I_D
    K0      <- P0 %*% solve(P0 + R_mat)
    uf      <- matrix(0, mp1, D)
    uf[1, ] <- u0 + K0 %*% (U[1, ] - u0)
    Pf      <- array(0, c(mp1, D, D))
    Pf[1,,] <- (I_D - K0) %*% P0
    list(u_filt = uf, P_filt = Pf)
  }

  # ── Pass 1: RK4-EKF forward + RTS backward to get a smoothed trajectory ──────
  s1          <- .init_filter()
  u_filt1     <- s1$u_filt
  P_filt1     <- s1$P_filt
  u_pred1     <- matrix(0, mp1, D)
  P_pred1     <- array(0, c(mp1, D, D))
  F_store1    <- array(0, c(mp1 - 1L, D, D))

  for (k in seq_len(mp1 - 1L)) {
    dt_k <- tt[k + 1L] - tt[k]
    uk   <- u_filt1[k, ]

    k1 <- f_(c(p, uk,                    tt[k]             ))
    k2 <- f_(c(p, uk + 0.5*dt_k*k1,     tt[k] + 0.5*dt_k ))
    k3 <- f_(c(p, uk + 0.5*dt_k*k2,     tt[k] + 0.5*dt_k ))
    k4 <- f_(c(p, uk +     dt_k*k3,     tt[k] +     dt_k  ))
    up_k <- uk + (dt_k / 6) * (k1 + 2*k2 + 2*k3 + k4)

    Ju_k               <- t(matrix(as.vector(J_u(c(p, uk, tt[k]))), D, D))
    Fk                 <- I_D + dt_k * Ju_k
    F_store1[k,,]      <- Fk
    Pk_pred            <- Fk %*% P_filt1[k,,] %*% t(Fk) + Q_mat1
    P_pred1[k + 1L,,]  <- Pk_pred
    u_pred1[k + 1L, ]  <- up_k

    S_k             <- Pk_pred + R_mat
    Kk              <- Pk_pred %*% solve(S_k)
    u_filt1[k+1L, ] <- up_k + Kk %*% (U[k+1L, ] - up_k)
    P_filt1[k+1L,,] <- (I_D - Kk) %*% Pk_pred
  }

  # RTS backward smoother for Pass 1
  u_smooth1        <- matrix(0, mp1, D)
  P_smooth1        <- array(0, c(mp1, D, D))
  u_smooth1[mp1, ] <- u_filt1[mp1, ]
  P_smooth1[mp1,,] <- P_filt1[mp1,,]

  for (k in seq(mp1 - 1L, 1L)) {
    Pk              <- P_filt1[k,,]
    Pp              <- P_pred1[k + 1L,,]
    Fk              <- F_store1[k,,]
    Gk              <- Pk %*% t(Fk) %*% solve(Pp)
    u_smooth1[k, ]  <- u_filt1[k, ] + Gk %*% (u_smooth1[k+1L, ] - u_pred1[k+1L, ])
    P_smooth1[k,,]  <- Pk            + Gk %*% (P_smooth1[k+1L,,] - Pp) %*% t(Gk)
  }

  # ── Pass 2: Wendy-EKF using the Pass-1 smoothed trajectory ──────────────────
  s2      <- .init_filter()
  u_filt  <- s2$u_filt
  P_filt  <- s2$P_filt
  u_pred  <- matrix(0, mp1, D)
  P_pred  <- array(0, c(mp1, D, D))
  F_store <- array(0, c(mp1 - 1L, D, D))

  for (k in seq_len(mp1 - 1L)) {
    dt_k <- tt[k + 1L] - tt[k]
    uk   <- u_filt[k, ]

    # Wendy predictor on RTS-smoothed (bidirectional) data
    u_pred[k + 1L, ] <- estimate_u_k(
      u_smooth1, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, test_function_params, k + 1L
    )

    Ju_k         <- t(matrix(as.vector(J_u(c(p, uk, tt[k]))), D, D))
    Fk           <- I_D + dt_k * Ju_k
    F_store[k,,] <- Fk

    Pk_pred          <- Fk %*% P_filt[k,,] %*% t(Fk) + Q_mat2
    P_pred[k + 1L,,] <- Pk_pred

    S_k              <- Pk_pred + R_mat
    Kk               <- Pk_pred %*% solve(S_k)
    u_filt[k + 1L, ] <- u_pred[k + 1L, ] + Kk %*% (U[k + 1L, ] - u_pred[k + 1L, ])
    P_filt[k + 1L,,] <- (I_D - Kk) %*% Pk_pred
  }

  # ── Backward pass: RTS smoother ──────────────────────────────────────────────
  u_smooth        <- matrix(0, mp1, D)
  P_smooth        <- array(0, c(mp1, D, D))
  u_smooth[mp1, ] <- u_filt[mp1, ]
  P_smooth[mp1,,] <- P_filt[mp1,,]

  for (k in seq(mp1 - 1L, 1L)) {
    Pk   <- P_filt[k,,]
    Pp   <- P_pred[k + 1L,,]
    Fk   <- F_store[k,,]

    Gk             <- Pk %*% t(Fk) %*% solve(Pp)
    u_smooth[k, ]  <- u_filt[k, ] + Gk %*% (u_smooth[k + 1L, ] - u_pred[k + 1L, ])
    P_smooth[k,,]  <- Pk           + Gk %*% (P_smooth[k + 1L,,] - Pp) %*% t(Gk)
  }

  list(
    U_star   = u_smooth,
    P_smooth = P_smooth,
    u_filt   = u_filt,
    P_filt   = P_filt,
    u_pred   = u_pred,
    P_pred   = P_pred
  )
}