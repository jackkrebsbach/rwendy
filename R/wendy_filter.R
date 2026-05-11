
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

  un0 <- U[1,]
  t0  <- tt[1]

  for (iter in seq_len(5)) {
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

#' Estimate the state using the RTS smoother
#' @export
wendy_erts <- function(U, f_, J_u, tt, p, test_function_params,
                       sigma = NULL) {
  tt   <- as.vector(tt)
  mp1  <- nrow(U)
  D    <- ncol(U)
  I_D  <- diag(D)

  noise_sd <- as.numeric(if (!is.null(sigma)) sigma else estimate_std(U, k = 6))
  noise_sd <- mean(noise_sd)
  R_obs    <- noise_sd^2 * I_D
  Q_mat    <- (noise_sd * 0.1)^2 * I_D

  u0 <- as.vector(U[1, ])

  P0           <- noise_sd^2 * I_D
  S0           <- P0 + R_obs
  K0           <- P0 %*% solve(S0)
  u_filt       <- matrix(0, mp1, D)
  u_filt[1, ]  <- u0 + K0 %*% (U[1, ] - u0)
  P_filt       <- array(0, c(mp1, D, D))
  P_filt[1,,]  <- (I_D - K0) %*% P0 %*% t(I_D - K0) + K0 %*% R_obs %*% t(K0)

  u_pred  <- matrix(0, mp1, D)
  P_pred  <- array(0, c(mp1, D, D))
  F_store <- array(0, c(mp1 - 1L, D, D))

  for (k in seq_len(mp1 - 1L)) {
    dt_k <- tt[k + 1L] - tt[k]
    uk   <- u_filt[k, ]

    k1 <- as.vector(f_(matrix(c(p, uk,                    tt[k]          ), ncol = 1)))
    k2 <- as.vector(f_(matrix(c(p, uk + 0.5*dt_k*k1,     tt[k]+0.5*dt_k ), ncol = 1)))
    k3 <- as.vector(f_(matrix(c(p, uk + 0.5*dt_k*k2,     tt[k]+0.5*dt_k ), ncol = 1)))
    k4 <- as.vector(f_(matrix(c(p, uk +     dt_k*k3,     tt[k]+    dt_k  ), ncol = 1)))
    u_pred[k + 1L, ] <- uk + (dt_k / 6) * (k1 + 2*k2 + 2*k3 + k4)

    Ju_k <- t(matrix(as.vector(J_u(c(p, uk, tt[k]))), D, D))
    Fk   <- I_D + dt_k * Ju_k
    F_store[k,,] <- Fk

    Pk_pred <- Fk %*% P_filt[k,,] %*% t(Fk) + Q_mat
    P_pred[k + 1L,,] <- Pk_pred

    Sk               <- Pk_pred + R_obs
    Kk               <- Pk_pred %*% solve(Sk)
    innov            <- U[k + 1L, ] - u_pred[k + 1L, ]
    u_filt[k + 1L, ] <- u_pred[k + 1L, ] + Kk %*% innov
    P_filt[k + 1L,,] <- (I_D - Kk) %*% Pk_pred %*% t(I_D - Kk) + Kk %*% R_obs %*% t(Kk)
  }

  u_smooth         <- matrix(0, mp1, D)
  P_smooth         <- array(0,  c(mp1, D, D))
  u_smooth[mp1, ]  <- u_filt[mp1, ]
  P_smooth[mp1,,]  <- P_filt[mp1,,]

  for (k in seq(mp1 - 1L, 1L)) {
    Pk <- P_filt[k,,]
    Pp <- P_pred[k + 1L,,]
    Fk <- F_store[k,,]

    Gk <- Pk %*% t(Fk) %*% solve(Pp + 1e-10 * I_D)
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

# Joint state-parameter Ensemble Kalman Filter (stochastic EnKF).
# Augmented state: [u(t); p]. Parameters are treated as constant between steps.
# Returns list(U_star, p, p_ens) — smoothed state mean, parameter mean, full parameter ensemble.
wendy_enkf <- function(p, U, tt, f_, sigma = NULL, N_e = 100, u0 = NULL) {
  mp1 <- nrow(U)
  D   <- ncol(U)
  J   <- length(p)
  tt  <- as.vector(tt)

  noise_sd <- mean(as.numeric(if (!is.null(sigma)) sigma else estimate_std(U, k = 6)))
  R_obs    <- noise_sd^2 * diag(D)

  u0_init <- if (!is.null(u0)) as.numeric(u0) else as.vector(U[1, ])
  if (length(u0_init) != D)
    stop(sprintf("u0 must have length D = %d", D))

  p_sd  <- pmax(abs(p) * 0.1, 1e-6)
  u0_sd <- rep(noise_sd, D)

  ens <- rbind(
    u0_init + matrix(rnorm(D * N_e) * u0_sd, D, N_e),
    p       + matrix(rnorm(J * N_e, sd = p_sd), J, N_e)
  )

  rk4_step <- function(u, pi, t, dt) {
    k1 <- as.vector(f_(matrix(c(pi, u,                  t         ), ncol = 1)))
    k2 <- as.vector(f_(matrix(c(pi, u + 0.5*dt*k1,     t + 0.5*dt), ncol = 1)))
    k3 <- as.vector(f_(matrix(c(pi, u + 0.5*dt*k2,     t + 0.5*dt), ncol = 1)))
    k4 <- as.vector(f_(matrix(c(pi, u +     dt*k3,     t +     dt), ncol = 1)))
    u + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
  }

  solve_safe <- function(A, b) {
    tryCatch(solve(A, b), error = function(e) MASS::ginv(A) %*% b)
  }

  U_mean       <- matrix(0, mp1, D)
  U_mean[1, ]  <- rowMeans(ens[seq_len(D), , drop = FALSE])

  for (k in seq_len(mp1 - 1L)) {
    dt_k <- tt[k + 1L] - tt[k]

    ens_f <- ens
    for (i in seq_len(N_e)) {
      u_i <- ens[seq_len(D), i]
      p_i <- ens[D + seq_len(J), i]
      ens_f[seq_len(D), i] <- tryCatch(
        rk4_step(u_i, p_i, tt[k], dt_k),
        error = function(e) u_i
      )
    }

    y_k    <- as.vector(U[k + 1L, ])
    mean_f <- rowMeans(ens_f)
    A_f    <- ens_f - mean_f
    H_ens  <- ens_f[seq_len(D), , drop = FALSE]
    A_H    <- H_ens - rowMeans(H_ens)

    C_uy <- (A_f %*% t(A_H)) / (N_e - 1)
    C_yy <- (A_H %*% t(A_H)) / (N_e - 1) + R_obs

    K      <- C_uy %*% solve_safe(C_yy, diag(D))
    y_pert <- y_k + matrix(rnorm(D * N_e, sd = noise_sd), D, N_e)
    ens    <- ens_f + K %*% (y_pert - H_ens)

    U_mean[k + 1L, ] <- rowMeans(ens[seq_len(D), , drop = FALSE])
  }

  list(
    U_star = U_mean,
    p      = rowMeans(ens[D + seq_len(J), , drop = FALSE]),
    p_ens  = t(ens[D + seq_len(J), , drop = FALSE])
  )
}
