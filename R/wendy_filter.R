# Leibniz expansion of g^(n)(t*) where g(t) = phi(t) F(p,u(t),t) + phi'(t) u(t)
# and trajectory derivatives F^(m) are passed in (precomputed via dF_dt_ etc.).
#
# phi_scalars: list of length (order+2) with phi^(0..order+1) at the endpoint.
# f_derivs:    list of length (order+1) with F^(0..order) at the endpoint.
# u_vec:       state at the endpoint (numeric vector of length D).
# order:       derivative order n.
g_deriv_at_endpoint <- function(phi_scalars, f_derivs, u_vec, order) {
  n      <- order
  result <- phi_scalars[[n + 2]] * u_vec # φ^(n+1)·u
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
# Returns function(U, p, tt) -> (K_bl × D) matrix, or NULL when K_bl == 0.
# EM_k = -dt²/12·(g_k'(t_M) - g_k'(t_1)) + dt⁴/720·(g_k'''(t_M) - g_k'''(t_1))
# with g(t) = phi_k(t) F(p,u(t),t) + phi_k'(t) u(t).
build_em_correction <- function(bl_phi_t1, bl_phi_tM,
                                f_, dF_dt_, d2F_dt2_, d3F_dt3_,
                                dt, scale = 1.0) {
  if (is.null(bl_phi_t1) || nrow(bl_phi_t1) == 0L) return(NULL)
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

    em
  }
}

#' Estimate u(0) via iterative defect-correction on left BL test functions
#'
#' Iterates the linear least-squares system
#' \deqn{B \, u_0^{(n+1)} = r_{\text{trap}} - \Delta_{EM}(u_0^{(n)})}
#' where \eqn{B} is the column vector of \eqn{\psi_k(0)} for K_bl left BL
#' test functions, \eqn{r_{\text{trap}}} is the fixed trapezoidal residual
#' \eqn{-T_h[f(u,\hat\theta)\psi] - T_h[u\,\psi']} (evaluated once on observed
#' U), and \eqn{\Delta_{EM}} is the analytic Euler-Maclaurin defect using
#' total time derivatives of f along the ODE. Contraction rate
#' \eqn{\kappa = O(h^2)}; typically 3-5 iterations.
#'
#' EM(2) keeps only the \eqn{h^2/12} correction (\eqn{O(h^4)} accuracy);
#' EM(4) adds \eqn{h^4/720} (\eqn{O(h^6)}). The LEFT-BL boundary term is
#' \eqn{B u_0}; right-side EM contributions vanish because left BL test
#' functions and all their derivatives are zero at \eqn{t=T}.
#'
#' @param U Numeric matrix (M x D) of observed states.
#' @param f_,dF_dt_,d2F_dt2_,d3F_dt3_ Callable RHS and total time-derivative
#'   evaluators built from the symbolic engine.
#' @param tt Numeric vector (length M) of time points.
#' @param p Numeric parameter vector \eqn{\hat\theta} (held fixed).
#' @param r_c Integer; left BL window half-width.
#' @param n_bl Optional integer; number of left BL test functions
#'   (default \code{max(3, ceiling(r_c/4))}).
#' @param max_iter,tol Fixed-point iteration controls.
#' @param em_order Either 2 or 4.
#' @param update_trap_u0 Logical. If \code{FALSE} (slide-literal default), the
#'   trapezoidal residual \eqn{r_{\text{trap}}} is computed once on observed
#'   \code{U} and held fixed across iterations. If \code{TRUE}, row 1 of
#'   \code{U} is replaced by the current \eqn{u_0^{(n)}} estimate at each
#'   iteration and \eqn{r_{\text{trap}}} is recomputed; the integrand no
#'   longer contains the noisy boundary observation but the contraction map
#'   changes.
#' @param J_u Callable state Jacobian \eqn{\partial f/\partial u}
#'   (\code{matrix(as.vector(J_u(c(p,u,t))), D, D)} with entry
#'   \eqn{[a,b] = \partial f_a/\partial u_b}). Used together with \code{sigma}
#'   (when \code{update_trap_u0 = FALSE}) so that \code{cov_u0} is the
#'   noise-channel propagation \eqn{\sum_m DU_m\,\mathrm{diag}(\sigma^2)\,DU_m^T}
#'   rather than the over-conservative LS-residual variance.
#' @param sigma Data-noise standard deviation, scalar or length-\eqn{D} (per
#'   state), used to scale the noise-channel \code{cov_u0}.
#' @return Named list with \code{U_hat} (U with row 1 replaced by
#'   \code{u0hat}), \code{u0hat}, \code{cov_u0} (D x D covariance of
#'   \eqn{\hat u_0}; the noise-channel propagation, falling back to the
#'   LS-residual variance \eqn{s_d^2/B^TB} only when \code{update_trap_u0 = TRUE}
#'   or \code{sigma} is degenerate),
#'   \code{cov_u0_resid} (the LS-residual variance, always returned for
#'   reference), \code{cov_method} (\code{"noise_propagation"} or
#'   \code{"ls_residual"}), plus \code{iters}, \code{converged},
#'   \code{diverged}, \code{u0_history}, \code{r_c}, \code{n_bl}, \code{K_bl},
#'   \code{em_order}, \code{update_trap_u0}.
#' @export
estimate_IC <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, r_c, J_u, sigma,
                        n_bl           = NULL,
                        max_iter       = 20L,
                        tol            = 1e-12,
                        em_order       = c(4L, 2L),
                        update_trap_u0 = FALSE) {
  em_order <- as.integer(em_order[1])
  if (!em_order %in% c(2L, 4L)) {
    stop("em_order must be 2 or 4", call. = FALSE)
  }

  M      <- nrow(U)
  D      <- ncol(U)
  J      <- length(p)
  tt_vec <- as.vector(tt)
  dt     <- mean(diff(tt_vec))

  r_c  <- min(r_c, floor((M - 1L) / 2L))
  n_bl <- if (is.null(n_bl)) max(3L, as.integer(ceiling(r_c / 4)))
          else                max(1L, as.integer(n_bl))

  bl_left <- lapply(0:4, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_c, order = ord,
                               side = "left", n_bl = n_bl))
  K_bl <- nrow(bl_left[[1]])

  apply_trap <- function(B) { B[, 1] <- B[, 1] * 0.5; B[, M] <- B[, M] * 0.5; B }
  V_BL  <- apply_trap(bl_left[[1]])
  Vp_BL <- apply_trap(bl_left[[2]])

  bl_phi_t1 <- matrix(0, nrow = K_bl, ncol = 5)
  for (ord in 0:4) bl_phi_t1[, ord + 1] <- bl_left[[ord + 1]][, 1]

  B   <- bl_phi_t1[, 1]
  BtB <- sum(B * B)
  if (!is.finite(BtB) || BtB < .Machine$double.eps) {
    u0_obs <- as.numeric(U[1, ])
    U_hat <- U; U_hat[1, ] <- u0_obs
    return(list(U_hat = U_hat, u0hat = u0_obs, cov_u0 = NULL,
                iters = 0L, converged = FALSE, diverged = FALSE,
                u0_history = matrix(u0_obs, nrow = 1),
                r_c = r_c, n_bl = n_bl, K_bl = K_bl, em_order = em_order,
                update_trap_u0 = update_trap_u0))
  }

  compute_r_trap <- function(U_in) {
    input  <- rbind(matrix(rep(p, M), nrow = J), t(U_in), matrix(tt_vec, nrow = 1L))
    F_eval <- f_(input)
    -dt * (V_BL %*% F_eval) - dt * (Vp_BL %*% U_in)
  }
  r_trap <- compute_r_trap(U)

  c2 <- dt^2 / 12
  c4 <- if (em_order >= 4L) dt^4 / 720 else 0

  em_correction <- function(u0_curr) {
    u_t1   <- as.vector(u0_curr)
    inp_t1 <- matrix(c(p, u_t1, tt_vec[1]), ncol = 1L)
    fd_t1  <- list(
      as.vector(f_(inp_t1)),
      as.vector(dF_dt_(inp_t1)),
      as.vector(d2F_dt2_(inp_t1)),
      as.vector(d3F_dt3_(inp_t1))
    )
    EM <- matrix(0, nrow = K_bl, ncol = D)
    for (k in seq_len(K_bl)) {
      phi_t1_k <- as.list(bl_phi_t1[k, ])
      g1 <- g_deriv_at_endpoint(phi_t1_k, fd_t1, u_t1, order = 1L)
      g3 <- if (em_order >= 4L)
              g_deriv_at_endpoint(phi_t1_k, fd_t1, u_t1, order = 3L)
            else rep(0, D)
      EM[k, ] <- c2 * g1 - c4 * g3
    }
    EM
  }

  # Track the best iterate by residual norm; the iteration contracts only when
  # the EM-correction map has spectral radius < 1 (kappa = O(h^2)). For
  # under-sampled stiff systems the contraction can fail and iterates blow up.
  # Returning the best-seen iterate keeps us no worse than the uncorrected LS.
  # Residual: ||B u0^T - (r_trap - EM(u0))||_F over the K_bl x D system.
  residual_norm <- function(u0_curr) {
    EM   <- em_correction(u0_curr)
    pred <- outer(B, as.numeric(u0_curr))       # K_bl x D
    sqrt(sum((pred - (r_trap - EM))^2))
  }

  u0       <- as.numeric(crossprod(B, r_trap) / BtB)
  u0_hist  <- list(u0)
  best_u0  <- u0
  best_res <- if (all(is.finite(u0))) residual_norm(u0) else Inf
  iters    <- 0L
  converged <- FALSE
  diverged  <- FALSE

  for (it in seq_len(max_iter)) {
    iters  <- it
    if (update_trap_u0) {
      U_curr <- U; U_curr[1, ] <- u0
      r_trap <- compute_r_trap(U_curr)
    }
    EM     <- em_correction(u0)
    rhs    <- r_trap - EM
    u0_new <- as.numeric(crossprod(B, rhs) / BtB)
    u0_hist[[it + 1L]] <- u0_new

    if (!all(is.finite(u0_new))) { diverged <- TRUE; break }

    res_new <- residual_norm(u0_new)
    if (is.finite(res_new) && res_new < best_res) {
      best_res <- res_new
      best_u0  <- u0_new
    }

    delta <- sqrt(sum((u0_new - u0)^2))
    u0    <- u0_new
    if (is.finite(delta) && delta < tol) {
      converged <- TRUE
      break
    }
  }

  u0 <- best_u0

  # --- Covariance of u0hat ------------------------------------------------
  #
  # Residual-based closed-form LS variance (fallback): for the K_bl x 1
  # per-state regression B u0_d = rhs_d, rhs_d = (r_trap - EM(u0))[, d], the
  # per-state residual variance is s_d^2 = ||B u0_d - rhs_d||^2 / (K_bl - 1)
  # and Var(u0_d) = s_d^2 / B^T B. This OVER-states Var(u0hat): the residual
  # norm absorbs the DETERMINISTIC Euler-Maclaurin / trapezoidal truncation
  # mismatch on top of the propagated noise (empirically ~1.8x too wide in
  # SE, ~100% coverage of a nominal-95% interval on the logistic problem).
  cov_u0_resid <- tryCatch({
    EM_final <- em_correction(u0)
    rhs      <- r_trap - EM_final            # K_bl x D
    e        <- outer(B, u0) - rhs            # K_bl x D residuals
    df       <- max(K_bl - 1L, 1L)
    s2       <- colSums(e * e) / df           # length-D
    diag(s2 / BtB, nrow = D, ncol = D)
  }, error = function(err) NULL)

  # Noise-channel propagation (preferred when J_u and sigma are supplied).
  # u0 = B^T(r_trap - EM(u0)) / B^T B is a deterministic function of the data
  # U. Differentiating the converged fixed-point relation
  #   (B^T B) I u0 + B^T EM(u0) = B^T r_trap(U)
  # w.r.t. U gives, for data entry (m, c),  (B^T B I + A) du0 = c^{(m,c)} with
  #   A[d,e]      = sum_k B[k] dEM[k,d]/du0_e        (EM hits U only via u0),
  #   c^{(m,c)}_d = -dt (B^T V_BL)[m] J_u(t_m)[d,c] - dt (B^T Vp_BL)[m] d_{dc},
  # because r_trap = -dt V_BL f(p,U,t) - dt Vp_BL U and f(p, U[m,], t_m) is
  # local in m with state Jacobian J_u. For per-state data noise with column
  # variances diag(sigma^2) (sigma scalar or length-D),
  #   Cov_noise(u0|p) = sum_m DU_m diag(sigma^2) DU_m^T,
  #   DU_m = (B^TB I + A)^-1 C_m,
  #   C_m  = -dt (B^T V_BL)[m] J_u(t_m) - dt (B^T Vp_BL)[m] I.
  # Only the m within the (compact) BL window contribute. This reflects only
  # sampling variability (calibrated ~95% coverage; validated against a
  # numeric Jacobian to <0.5%). The parameter contribution J_p Cov(phat) J_p^T
  # is left to the caller (negligible when p is well identified; the library
  # reports the Fisher param cov separately). Only valid when r_trap depends on
  # the raw data, i.e. !update_trap_u0.
  # J_u and sigma are required; the noise-channel form is used whenever it is
  # valid (raw-data residual, well-formed per-state sigma) and falls back to
  # the LS-residual variance only for update_trap_u0 = TRUE or degenerate sigma.
  I_D <- diag(D)
  use_noiseprop <- length(sigma) %in% c(1L, D) && all(is.finite(sigma)) &&
                   !update_trap_u0

  cov_u0_noise <- if (use_noiseprop) tryCatch({
    sig_vec <- if (length(sigma) == 1L) rep(sigma, D) else as.numeric(sigma)
    Sig2    <- diag(sig_vec^2, nrow = D, ncol = D)
    h <- 1e-6 * max(1, sqrt(sum(u0^2)))
    A <- matrix(0, D, D)
    for (e_i in seq_len(D)) {
      up <- u0; up[e_i] <- up[e_i] + h
      dn <- u0; dn[e_i] <- dn[e_i] - h
      dEM_e    <- (em_correction(up) - em_correction(dn)) / (2 * h)  # K_bl x D
      A[, e_i] <- as.numeric(crossprod(B, dEM_e))
    }
    Minv <- solve(BtB * I_D + A)
    bV   <- as.numeric(crossprod(B, V_BL))   # length M  (= B^T V_BL)
    bVp  <- as.numeric(crossprod(B, Vp_BL))  # length M  (= B^T Vp_BL)
    acc  <- matrix(0, D, D)
    for (m in seq_len(M)) {
      if (bV[m] == 0 && bVp[m] == 0) next
      Ju_m <- matrix(as.vector(J_u(c(p, U[m, ], tt_vec[m]))), D, D)
      C_m  <- -dt * bV[m] * Ju_m - dt * bVp[m] * I_D
      DU_m <- Minv %*% C_m
      acc  <- acc + DU_m %*% Sig2 %*% t(DU_m)
    }
    acc
  }, error = function(err) NULL) else NULL

  cov_u0     <- if (!is.null(cov_u0_noise)) cov_u0_noise else cov_u0_resid
  cov_method <- if (!is.null(cov_u0_noise)) "noise_propagation" else "ls_residual"

  U_hat <- U; U_hat[1, ] <- u0

  list(
    U_hat          = U_hat,
    u0hat          = u0,
    cov_u0         = cov_u0,
    cov_u0_resid   = cov_u0_resid,
    cov_method     = cov_method,
    iters          = iters,
    converged      = converged,
    diverged       = diverged,
    u0_history     = do.call(rbind, u0_hist),
    r_c            = r_c,
    n_bl           = n_bl,
    K_bl           = K_bl,
    em_order       = em_order,
    update_trap_u0 = update_trap_u0
  )
}

#' Estimate the state using the RTS smoother
#'
#' Extended Kalman filter forward pass followed by a Rauch-Tung-Striebel
#' backward smoother on the parameter-conditioned dynamics.
#'
#' When \code{param_cov} and \code{J_p} are both supplied, the process
#' covariance at step k is propagated from the WENDy parameter covariance:
#' \deqn{Q_k = (\Delta t_k)^2 \, \nabla_p f(\hat{p}, u_k, t_k) \, \hat{C} \,
#'        \nabla_p f(\hat{p}, u_k, t_k)^T,}
#' which absorbs parameter uncertainty into the predict step. Without them,
#' \code{Q_k} falls back to the naive \code{(0.1 sigma)^2 I_D}.
#'
#' @param U Numeric matrix (mp1 x D) of noisy observations.
#' @param f_ Callable f(p, u, t) RHS evaluator built from the symbolic engine.
#' @param J_u Callable Jacobian d f / d u evaluator from the symbolic engine.
#' @param tt Numeric vector of time points (length mp1).
#' @param p Numeric vector of parameter estimates.
#' @param test_function_params Unused; retained for API compatibility.
#' @param sigma Optional scalar or per-state noise SD; if NULL it is estimated
#'   from U via \code{estimate_std}.
#' @param u0_init Optional length-D vector to seed the filter at \code{tt[1]}
#'   (e.g. \code{estimate_IC()$u0hat}). Defaults to \code{U[1, ]}.
#' @param P0_init Optional D x D prior covariance for \code{u0_init}. When
#'   \code{u0_init} is supplied, defaults to \code{0.1 * sigma^2 I_D} (more
#'   confident than the raw observation); otherwise \code{sigma^2 I_D}.
#' @param param_cov Optional J x J parameter covariance \eqn{\hat{C}}.
#' @param J_p Optional callable returning the D x J Jacobian
#'   \eqn{\nabla_p f(p, u, t)} (flattened d-fast, like the rest of the package).
#' @return Named list with U_star (smoothed state), P_smooth (posterior
#'   covariance), and the intermediate filter/predictor states.
#' @export
wendy_erts <- function(U, f_, J_u, tt, p, test_function_params,
                       sigma = NULL,
                       u0_init = NULL,
                       P0_init = NULL,
                       param_cov = NULL,
                       J_p = NULL) {
  tt   <- as.vector(tt)
  mp1  <- nrow(U)
  D    <- ncol(U)
  J    <- length(p)
  I_D  <- diag(D)

  noise_sd <- as.numeric(if (!is.null(sigma)) sigma else estimate_std(U, k = 6))
  noise_sd <- mean(noise_sd)
  R_obs    <- noise_sd^2 * I_D

  use_param_Q <- !is.null(param_cov) && !is.null(J_p)
  Q_default   <- (noise_sd * 0.1)^2 * I_D

  u0 <- if (!is.null(u0_init)) as.vector(u0_init) else as.vector(U[1, ])
  P0 <- if (!is.null(P0_init)) P0_init
        else if (!is.null(u0_init)) 0.1 * noise_sd^2 * I_D
        else noise_sd^2 * I_D

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

    Ju_k <- matrix(as.vector(J_u(c(p, uk, tt[k]))), D, D)  # J[a, b] = df_a/du_b
    Fk   <- I_D + dt_k * Ju_k
    F_store[k,,] <- Fk

    Q_k <- if (use_param_Q) {
      Jp_k <- matrix(as.vector(J_p(c(p, uk, tt[k]))), D, J)  # J[a, j] = df_a/dp_j
      (dt_k^2) * (Jp_k %*% param_cov %*% t(Jp_k))
    } else {
      Q_default
    }

    Pk_pred <- Fk %*% P_filt[k,,] %*% t(Fk) + Q_k
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