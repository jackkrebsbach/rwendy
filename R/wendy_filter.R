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

#' Estimate u0 and/or uM by minimising the augmented BL weak residual.
#'
#' Builds SSL boundary-layer test functions at the SSL change-point radius r_c,
#' assembles the augmented residual (trapezoid + 2-term Euler-Maclaurin) and
#' optimises the corner rows of U via Levenberg-Marquardt. Which corners are
#' treated as unknowns is controlled by `corners`:
#'   - `"both"` (default): both U[1, ] and U[M, ] are optimised (2*D unknowns).
#'   - `"u0"`: only U[1, ] is optimised; U[M, ] is held at its observed value.
#'   - `"uM"`: only U[M, ] is optimised; U[1, ] is held at its observed value.
#' The interior of U (rows 2..M-1) is held at its observed values throughout.
#' With K_bl*D ~ 20+ residual equations and up to 2*D unknowns, the LS is
#' overdetermined and needs no regularisation; warm-started from the observed
#' corner values.
#'
#' @param U Numeric matrix (M x D) of observed states.
#' @param f_ Callable f(p, u, t) RHS evaluator built from the symbolic engine.
#' @param dF_dt_,d2F_dt2_,d3F_dt3_ Callable evaluators for the first, second,
#'   and third total time derivatives of f along the trajectory.
#' @param tt Numeric vector or column matrix of time points (length M).
#' @param p Numeric vector of parameter estimates.
#' @param r_c Integer; SSL change-point radius (BL window size).
#' @param corners One of `"both"`, `"u0"`, `"uM"`; selects which corner rows
#'   are treated as unknowns. Defaults to `"both"`.
#' @param n_bl_per_side Optional integer; number of BL test functions per side.
#'   Defaults to `max(3, ceiling(r_c / 4))` — empirically near-optimal across
#'   logistic and Lotka-Volterra benchmarks under uniform placement.
#' @return Named list with `U_hat` (full M x D state, interior verbatim, corner
#'   rows replaced by the LS estimate where applicable), `u0hat` (= `U_hat[1, ]`),
#'   `uMhat` (= `U_hat[M, ]`), `r_c`, and `corners`.
#' @export
estimate_u0 <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, r_c,
                        corners = c("both", "u0", "uM"),
                        n_bl_per_side = NULL) {
  corners <- match.arg(corners)
  M  <- nrow(U)
  D  <- ncol(U)
  J  <- length(p)
  tt_vec <- as.vector(tt)
  dt <- mean(diff(tt_vec))

  # Largest allowed radius so left/right windows never overlap.
  r_c <- min(r_c, floor((M - 1L) / 2L))

  # Default count empirically optimal at r_c/4 across logistic and LV; lower
  # bound at 3 because anything fewer underdetermines the placement spread.
  n_bl_per_side <- if (is.null(n_bl_per_side)) {
    max(3L, as.integer(ceiling(r_c / 4)))
  } else {
    max(1L, as.integer(n_bl_per_side))
  }

  # Uniform placement: peaks evenly spaced across the BL window so the basis
  # spans the whole window instead of collapsing onto the corner.
  step <- max(1L, as.integer(floor((r_c - 1L) / max(1L, n_bl_per_side - 1L))))

  bl_left  <- lapply(0:4, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_c, order = ord, side = "left",
                               n_bl = n_bl_per_side, step = step))
  bl_right <- lapply(0:4, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_c, order = ord, side = "right",
                               n_bl = n_bl_per_side, step = step))
  n_bl <- nrow(bl_left[[1]])
  K_bl <- 2L * n_bl

  # Trapezoid endpoint weights so V_BL %*% F matches the trapezoid sum
  # the Euler-Maclaurin correction is derived for.
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

  # theta layout depends on `corners`: "both" => c(U[1, ], U[M, ]) (2*D),
  # "u0" => U[1, ] only (D), "uM" => U[M, ] only (D). Whichever corner is not
  # being optimised stays at its observed value; the interior of U is always
  # held at its observed values.
  residual_fn <- function(theta) {
    U_hat <- U
    if (corners == "both") {
      U_hat[1, ] <- theta[1:D]
      U_hat[M, ] <- theta[D + 1:D]
    } else if (corners == "u0") {
      U_hat[1, ] <- theta[1:D]
    } else {
      U_hat[M, ] <- theta[1:D]
    }

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

    as.vector(trap_phiF + trap_phipU + bdry + EM)
  }

  theta0 <- switch(corners,
                   both = c(U[1, ], U[M, ]),
                   u0   = as.numeric(U[1, ]),
                   uM   = as.numeric(U[M, ]))

  fit <- nls.lm(par = theta0, fn = residual_fn,
                control = nls.lm.control(maxiter = 100))

  U_hat <- U
  if (corners == "both") {
    U_hat[1, ] <- fit$par[1:D]
    U_hat[M, ] <- fit$par[D + 1:D]
  } else if (corners == "u0") {
    U_hat[1, ] <- fit$par[1:D]
  } else {
    U_hat[M, ] <- fit$par[1:D]
  }

  # Asymptotic NLS covariance: Cov(theta_hat) = s^2 * (J^T J)^{-1}, where
  # fit$hessian = J^T J (verified empirically for minpack.lm) and s^2 is the
  # residual variance at the optimum. Slice into u0 / uM blocks per `corners`.
  n_res <- length(fit$fvec)
  n_par <- length(fit$par)
  df    <- max(1L, n_res - n_par)
  s2    <- sum(fit$fvec^2) / df
  cov_full <- tryCatch(
    s2 * solve(fit$hessian + 1e-12 * diag(n_par)),
    error = function(e) NULL
  )
  cov_u0 <- NULL
  cov_uM <- NULL
  if (!is.null(cov_full)) {
    if (corners == "both") {
      cov_u0 <- cov_full[1:D, 1:D, drop = FALSE]
      cov_uM <- cov_full[D + 1:D, D + 1:D, drop = FALSE]
    } else if (corners == "u0") {
      cov_u0 <- cov_full
    } else {
      cov_uM <- cov_full
    }
  }

  list(
    U_hat   = U_hat,
    u0hat   = U_hat[1, ],
    uMhat   = U_hat[M, ],
    cov_u0  = cov_u0,
    cov_uM  = cov_uM,
    r_c     = r_c,
    corners = corners
  )
}

#' Iterative defect-correction estimator for u(0) via left BL test functions
#'
#' Alternative to \code{estimate_u0}: instead of NLS over corner unknowns,
#' iterate the linear least-squares system
#' \deqn{B \, u_0^{(n+1)} = r_{\text{trap}} - \Delta_{EM}(u_0^{(n)})}
#' where \eqn{B} is the column vector of \eqn{\psi_k(0)} for K_bl left BL
#' test functions, \eqn{r_{\text{trap}}} is the fixed trapezoidal residual
#' \eqn{-T_h[f(u,\hat\theta)\psi] - T_h[u\,\psi']} (evaluated once on observed
#' U), and \eqn{\Delta_{EM}} is the analytic Euler-Maclaurin defect using
#' total time derivatives of f along the ODE. Contraction rate
#' \eqn{\kappa = O(h^2)}; typically 3-5 iterations.
#'
#' EM(2) keeps only the \eqn{h^2/12} correction (\eqn{O(h^4)} accuracy);
#' EM(4) adds \eqn{h^4/720} (\eqn{O(h^6)}). Sign convention matches
#' \code{estimate_u0}: the LEFT-BL boundary term is \eqn{B u_0}; right-side
#' EM contributions vanish because left BL test functions and all their
#' derivatives are zero at \eqn{t=T}.
#'
#' @param U Numeric matrix (M x D) of observed states.
#' @param f_,dF_dt_,d2F_dt2_,d3F_dt3_ Callable RHS and total time-derivative
#'   evaluators built from the symbolic engine.
#' @param tt Numeric vector (length M) of time points.
#' @param p Numeric parameter vector \eqn{\hat\theta} (held fixed).
#' @param r_c Integer; left BL window half-width.
#' @param n_bl Optional integer; number of left BL test functions
#'   (default \code{max(3, ceiling(r_c/4))}, same heuristic as
#'   \code{estimate_u0}).
#' @param max_iter,tol Fixed-point iteration controls.
#' @param em_order Either 2 or 4.
#' @param update_trap_u0 Logical. If \code{FALSE} (slide-literal default), the
#'   trapezoidal residual \eqn{r_{\text{trap}}} is computed once on observed
#'   \code{U} and held fixed across iterations. If \code{TRUE}, row 1 of
#'   \code{U} is replaced by the current \eqn{u_0^{(n)}} estimate at each
#'   iteration and \eqn{r_{\text{trap}}} is recomputed; the integrand no
#'   longer contains the noisy boundary observation but the contraction map
#'   changes.
#' @return Named list matching \code{estimate_u0}'s interface where applicable:
#'   \code{U_hat} (U with row 1 replaced by \code{u0hat}), \code{u0hat},
#'   \code{uMhat} (pass-through of observed \code{U[M, ]}; bldc does not
#'   estimate the right corner), \code{cov_u0} (D x D diagonal from the
#'   closed-form LS variance \eqn{s_d^2 / B^T B}), \code{cov_uM} (NULL),
#'   plus method-specific fields \code{iters}, \code{converged},
#'   \code{diverged}, \code{u0_history}, \code{r_c}, \code{n_bl}, \code{K_bl},
#'   \code{em_order}, \code{update_trap_u0}.
#' @export
estimate_u0_bldc <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, r_c,
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
  step <- max(1L, as.integer(floor((r_c - 1L) / max(1L, n_bl - 1L))))

  bl_left <- lapply(0:4, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_c, order = ord, side = "left",
                               n_bl = n_bl, step = step))
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
    return(list(U_hat = U_hat, u0hat = u0_obs, uMhat = as.numeric(U[M, ]),
                cov_u0 = NULL, cov_uM = NULL,
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

  # Closed-form LS covariance for the K_bl x 1 per-state regression
  # B * u0_d = rhs_d, where rhs_d = (r_trap - EM(u0))[, d]. With residuals
  # e_d = B u0_d - rhs_d, the per-state residual variance is s_d^2 = ||e_d||^2
  # / max(K_bl - 1, 1) and Var(u0_d) = s_d^2 / B^T B. Cross-state correlation
  # is left as zero (data noise iid across states); off-diagonals would
  # require the residual cross-covariance which we don't track here.
  cov_u0 <- tryCatch({
    EM_final <- em_correction(u0)
    rhs      <- r_trap - EM_final            # K_bl x D
    e        <- outer(B, u0) - rhs            # K_bl x D residuals
    df       <- max(K_bl - 1L, 1L)
    s2       <- colSums(e * e) / df           # length-D
    diag(s2 / BtB, nrow = D, ncol = D)
  }, error = function(err) NULL)

  U_hat <- U; U_hat[1, ] <- u0

  list(
    U_hat          = U_hat,
    u0hat          = u0,
    uMhat          = as.numeric(U[M, ]),
    cov_u0         = cov_u0,
    cov_uM         = NULL,
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

#' Propagate parameter sensitivity along a trajectory
#'
#' Solves the linear sensitivity ODE
#' \deqn{\frac{d S}{d t} = J_u(u(t), \hat\theta, t)\, S + J_p(u(t), \hat\theta, t),
#'        \quad S(0) = 0}
#' along the supplied state trajectory using Heun's method (RK2). Returns the
#' \eqn{D \times J} sensitivity matrix \eqn{S(t_k) = \partial u(t_k) / \partial
#' \theta} at each grid time. No deSolve call — the right-hand side is
#' evaluated via the analytic Jacobians \code{J_u} and \code{J_p} the
#' symbolic engine already produced.
#'
#' When \code{param_cov} is supplied, also returns the trajectory variance
#' contribution \eqn{P_{\text{param}}(t_k) = S(t_k)\, \hat C\, S(t_k)^T} as
#' a \eqn{D \times D \times M} array; this stacks additively with the ERTS
#' posterior \code{P_smooth} to give the total state band.
#'
#' @param U Numeric matrix (M x D) of state values along the trajectory; the
#'   smoothed state \code{wendy_erts()$U_star} is the natural input.
#' @param tt Numeric vector (length M) of time points.
#' @param p Numeric parameter vector \eqn{\hat\theta}.
#' @param J_u Callable Jacobian \eqn{\partial f / \partial u} (D*D flat,
#'   d-fast within each row).
#' @param J_p Callable Jacobian \eqn{\partial f / \partial p} (D*J flat,
#'   d-fast within each row).
#' @param param_cov Optional J x J parameter covariance \eqn{\hat C}; if
#'   supplied, \code{P_param} is returned.
#' @return List with \code{S} (length-M list of D x J matrices) and, when
#'   \code{param_cov} is provided, \code{P_param} (M x D x D array).
#' @export
propagate_param_sensitivity <- function(U, tt, p, J_u, J_p, param_cov = NULL) {
  tt <- as.vector(tt)
  M  <- nrow(U)
  D  <- ncol(U)
  J  <- length(p)

  S <- vector("list", M)
  S[[1L]] <- matrix(0, D, J)

  eval_AB <- function(u_pt, t_pt) {
    inp <- c(p, as.numeric(u_pt), t_pt)
    list(A = matrix(as.vector(J_u(inp)), D, D),
         B = matrix(as.vector(J_p(inp)), D, J))
  }

  ab_prev <- eval_AB(U[1L, ], tt[1L])
  for (k in seq_len(M - 1L)) {
    dt_k    <- tt[k + 1L] - tt[k]
    ab_next <- eval_AB(U[k + 1L, ], tt[k + 1L])
    k1      <- ab_prev$A %*% S[[k]] + ab_prev$B
    S_pred  <- S[[k]] + dt_k * k1
    k2      <- ab_next$A %*% S_pred + ab_next$B
    S[[k + 1L]] <- S[[k]] + 0.5 * dt_k * (k1 + k2)
    ab_prev <- ab_next
  }

  out <- list(S = S)
  if (!is.null(param_cov)) {
    P_param <- array(0, c(M, D, D))
    for (k in seq_len(M)) {
      Sk <- S[[k]]
      P_param[k,,] <- Sk %*% param_cov %*% t(Sk)
    }
    out$P_param <- P_param
  }
  out
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
#'   (e.g. \code{estimate_u0()$u0hat}). Defaults to \code{U[1, ]}.
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