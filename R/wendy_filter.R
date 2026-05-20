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
estimate_u0 <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, r_c) {
  M  <- nrow(U)
  D  <- ncol(U)
  J  <- length(p)
  tt_vec <- as.vector(tt)
  dt <- mean(diff(tt_vec))

  # r_c sets the BL test-function support (passed in by caller: SSL change-point
  # radius if test_fun_type == "SSL", else the MSG min radius).  Clamp so the
  # left/right windows don't overlap.
  r_c <- min(r_c, floor((M - 1L) / 2L))

  # Variable layout: theta = [vec(U_left); vec(U_right)], each r_c x D.
  l_idx   <- seq_len(r_c)
  r_idx   <- seq(M - r_c + 1L, M)
  int_idx <- seq(r_c + 1L, M - r_c)
  n_l     <- length(l_idx) * D
  n_r     <- length(r_idx) * D

  # Use n_bl = r_c BL test functions per side so K_bl * D = 2 * r_c * D matches
  # the number of variables (the full BL window).  With the default n_bl ≈ r_c/5
  # the system is severely underdetermined once we expand from 2 corner rows to
  # the full r_c rows on each side.
  n_bl_per_side <- max(1L, as.integer(r_c))

  # Build BL test-function blocks at orders 0..4 for left and right sides.
  bl_left  <- lapply(0:4, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_c, order = ord, side = "left",  n_bl = n_bl_per_side))
  bl_right <- lapply(0:4, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_c, order = ord, side = "right", n_bl = n_bl_per_side))
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

  U_interior <- U[int_idx, , drop = FALSE]

  reconstruct_U <- function(theta) {
    U_l <- matrix(theta[seq_len(n_l)],       nrow = length(l_idx), ncol = D)
    U_r <- matrix(theta[n_l + seq_len(n_r)], nrow = length(r_idx), ncol = D)
    U_hat <- U
    U_hat[l_idx, ] <- U_l
    U_hat[r_idx, ] <- U_r
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

  # Warm start from the input rows in the BL window.
  theta0 <- c(as.vector(U[l_idx, ]), as.vector(U[r_idx, ]))
  fit    <- nls.lm(par = theta0, fn = residual_fn)

  U_hat       <- reconstruct_U(fit$par)
  U_left_hat  <- U_hat[l_idx, , drop = FALSE]
  U_right_hat <- U_hat[r_idx, , drop = FALSE]

  list(
    U_hat       = U_hat,
    U_left_hat  = U_left_hat,           # r_c x D — refined left BL window
    U_right_hat = U_right_hat,          # r_c x D — refined right BL window
    u0hat       = U_left_hat[1L, ],     # back-compat: corner of left window
    uMhat       = U_right_hat[r_c, ],   # back-compat: corner of right window
    r_c         = r_c
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