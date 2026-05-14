# Leibniz expansion of g^(n)(t*) where g(t) = phi(t) F(p,u(t),t) + phi'(t) u(t)
# and trajectory derivatives F^(m) are passed in (precomputed via dF_dt_ etc.).
#
# phi_scalars: list of length (order+2) with phi^(0..order+1) at the endpoint.
# f_derivs:    list of length (order+1) with F^(0..order) at the endpoint.
# u_vec:       state at the endpoint (numeric vector of length D).
# order:       derivative order n.
g_deriv_at_endpoint <- function(phi_scalars, f_derivs, u_vec, order) {
  n      <- order
  result <- phi_scalars[[n + 2]] * u_vec                        # φ^(n+1)·u
  for (k in seq(0, n)) {
    result <- result + choose(n, k) * phi_scalars[[k + 1]] * f_derivs[[n - k + 1]]
  }
  if (n >= 1) {
    for (k in seq(0, n - 1)) {
      result <- result + choose(n, k) * phi_scalars[[k + 2]] * f_derivs[[n - k]]
    }
  }
  result
}

# Build the EM correction closure for the augmented BL residual.
#
# bl_phi_t1, bl_phi_tM: (K_bl × 5) matrices with raw phi^(0..4) at the
#   left/right boundary for each BL test function (columns = derivative orders).
# f_, dF_dt_, d2F_dt2_, d3F_dt3_: trajectory-derivative callables of the RHS.
#
# Returns function(U, p, tt) -> (K_bl × D) torch tensor, or NULL when K_bl == 0.
# EM_k = -dt²/12·(g_k'(t_M) - g_k'(t_1)) + dt⁴/720·(g_k'''(t_M) - g_k'''(t_1))
# with g(t) = phi_k(t) F(p,u(t),t) + phi_k'(t) u(t).
build_em_correction <- function(bl_phi_t1, bl_phi_tM,
                                f_, dF_dt_, d2F_dt2_, d3F_dt3_,
                                dt, scale = 1.0, device = torch::torch_device("cpu")) {
  if (is.null(bl_phi_t1) || nrow(bl_phi_t1) == 0L) return(NULL)
  dtype <- torch::torch_float64()
  K_bl  <- nrow(bl_phi_t1)
  c2    <- scale * dt^2 / 12
  c4    <- scale * dt^4 / 720

  function(U, p, tt) {
    D    <- ncol(U)
    M    <- nrow(U)
    u_t1 <- as.numeric(U[1L, ])
    u_tM <- as.numeric(U[M, ])
    t1   <- tt[1L]
    tM   <- tt[M]

    eval_derivs <- function(u_pt, t_pt) {
      input <- matrix(c(p, u_pt, t_pt), ncol = 1L)
      list(
        as.vector(f_(input)),
        as.vector(dF_dt_(input)),
        as.vector(d2F_dt2_(input)),
        as.vector(d3F_dt3_(input))
      )
    }

    f_derivs_t1 <- eval_derivs(u_t1, t1)
    f_derivs_tM <- eval_derivs(u_tM, tM)

    em <- matrix(0, nrow = K_bl, ncol = D)
    for (k in seq_len(K_bl)) {
      phi_t1 <- as.list(bl_phi_t1[k, ])   # phi^(0..4) at t_1
      phi_tM <- as.list(bl_phi_tM[k, ])   # phi^(0..4) at t_M

      g1_t1 <- g_deriv_at_endpoint(phi_t1, f_derivs_t1, u_t1, order = 1L)
      g1_tM <- g_deriv_at_endpoint(phi_tM, f_derivs_tM, u_tM, order = 1L)
      g3_t1 <- g_deriv_at_endpoint(phi_t1, f_derivs_t1, u_t1, order = 3L)
      g3_tM <- g_deriv_at_endpoint(phi_tM, f_derivs_tM, u_tM, order = 3L)

      em[k, ] <- -c2 * (g1_tM - g1_t1) + c4 * (g3_tM - g3_t1)
    }

    torch::torch_tensor(em, dtype = dtype, device = device)
  }
}

#' Estimate u0 and uM by minimising the augmented BL weak residual.
#'
#' Builds SSL boundary-layer test functions at the SSL change-point radius r_c,
#' assembles the augmented residual (trapezoid + 2-term Euler-Maclaurin) and
#' optimises only the two corner rows U[1, ] and U[M, ] via Levenberg-Marquardt.
#' The interior of U (rows 2..M-1) is held at its observed values throughout.
#' With K_bl*D ~ 20+ residual equations and 2*D unknowns, the LS is comfortably
#' overdetermined and needs no regularisation; warm-started from the observed
#' corner values.
#'
#' @return Named list with U_hat (full M x D state, interior verbatim, corners
#'   replaced), u0hat (= U_hat[1, ]), uMhat (= U_hat[M, ]), and r_c.
#' @export
estimate_u0 <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, test_function_params) {
  M  <- nrow(U)
  D  <- ncol(U)
  J  <- length(p)
  tt_vec <- as.vector(tt)
  dt <- mean(diff(tt_vec))

  # SSL change-point radius (sets the BL test-function support).
  data_rc <- compute_r_c_hat(U, tt, test_function_params$S, test_function_params$p)
  r_c     <- data_rc$rc

  # Build BL test-function blocks at orders 0..4 for left and right sides.
  bl_left  <- lapply(0:4, function(ord) build_boundary_layer_block(psi, tt_vec, r_c, order = ord, side = "left"))
  bl_right <- lapply(0:4, function(ord) build_boundary_layer_block(psi, tt_vec, r_c, order = ord, side = "right"))
  n_bl <- nrow(bl_left[[1]])
  K_bl <- 2L * n_bl

  # Trapezoid endpoint weights on BL rows (so V_BL %*% F is the trapezoid sum
  # the EM correction is derived for).
  apply_trap <- function(B) { B[, 1] <- B[, 1] * 0.5; B[, M] <- B[, M] * 0.5; B }
  V_BL  <- rbind(apply_trap(bl_left[[1]]), apply_trap(bl_right[[1]]))   # K_bl x M
  Vp_BL <- rbind(apply_trap(bl_left[[2]]), apply_trap(bl_right[[2]]))   # K_bl x M

  # Endpoint phi^(0..4) per BL row (raw, used by boundary terms and EM).
  bl_phi_t1 <- matrix(0, nrow = K_bl, ncol = 5)
  bl_phi_tM <- matrix(0, nrow = K_bl, ncol = 5)
  for (ord in 0:4) {
    bl_phi_t1[1:n_bl,           ord + 1] <- bl_left [[ord + 1]][, 1]
    bl_phi_tM[1:n_bl,           ord + 1] <- bl_left [[ord + 1]][, M]
    bl_phi_t1[(n_bl + 1L):K_bl, ord + 1] <- bl_right[[ord + 1]][, 1]
    bl_phi_tM[(n_bl + 1L):K_bl, ord + 1] <- bl_right[[ord + 1]][, M]
  }

  reconstruct_U <- function(theta) {
    U_hat <- U
    U_hat[1, ] <- theta[seq_len(D)]
    U_hat[M, ] <- theta[(D + 1L):(2L * D)]
    U_hat
  }

  residual_fn <- function(theta) {
    U_hat <- reconstruct_U(theta)

    # F(p, U_hat, t) on the t-grid: build_fn returns M x D (do_transpose = TRUE).
    input  <- rbind(matrix(rep(p, M), nrow = J), t(U_hat), matrix(tt_vec, nrow = 1L))
    F_eval <- f_(input)

    trap_phiF  <- dt * (V_BL  %*% F_eval)        # K_bl x D
    trap_phipU <- dt * (Vp_BL %*% U_hat)         # K_bl x D
    bdry       <- bl_phi_t1[, 1] %o% U_hat[1, ] - bl_phi_tM[, 1] %o% U_hat[M, ]

    u_t1 <- as.vector(U_hat[1, ])
    u_tM <- as.vector(U_hat[M, ])
    inp_t1 <- matrix(c(p, u_t1, tt_vec[1]), ncol = 1L)
    inp_tM <- matrix(c(p, u_tM, tt_vec[M]), ncol = 1L)
    fd_t1 <- list(as.vector(f_(inp_t1)), as.vector(dF_dt_(inp_t1)),
                  as.vector(d2F_dt2_(inp_t1)), as.vector(d3F_dt3_(inp_t1)))
    fd_tM <- list(as.vector(f_(inp_tM)), as.vector(dF_dt_(inp_tM)),
                  as.vector(d2F_dt2_(inp_tM)), as.vector(d3F_dt3_(inp_tM)))

    EM <- matrix(0, nrow = K_bl, ncol = D)
    for (k in seq_len(K_bl)) {
      phi_t1_k <- as.list(bl_phi_t1[k, ])
      phi_tM_k <- as.list(bl_phi_tM[k, ])
      g1_t1 <- g_deriv_at_endpoint(phi_t1_k, fd_t1, u_t1, order = 1L)
      g1_tM <- g_deriv_at_endpoint(phi_tM_k, fd_tM, u_tM, order = 1L)
      g3_t1 <- g_deriv_at_endpoint(phi_t1_k, fd_t1, u_t1, order = 3L)
      g3_tM <- g_deriv_at_endpoint(phi_tM_k, fd_tM, u_tM, order = 3L)
      EM[k, ] <- -(dt^2 / 12) * (g1_tM - g1_t1) + (dt^4 / 720) * (g3_tM - g3_t1)
    }

    r_aug <- trap_phiF + trap_phipU + bdry + EM   # K_bl x D
    as.vector(r_aug)
  }

  # Warm start from the observed corner values.
  theta0 <- c(as.vector(U[1, ]), as.vector(U[M, ]))
  fit    <- nls.lm(par = theta0, fn = residual_fn)

  u0hat <- as.numeric(fit$par[seq_len(D)])
  uMhat <- as.numeric(fit$par[(D + 1L):(2L * D)])
  U_hat <- U
  U_hat[1, ] <- u0hat
  U_hat[M, ] <- uMhat

  list(
    U_hat = U_hat,
    u0hat = u0hat,
    uMhat = uMhat,
    r_c   = r_c
  )
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
