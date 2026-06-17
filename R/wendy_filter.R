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

# Euler-Maclaurin defect for the boundary-layer weak residual.
#
# A trapezoidal weak integral int (phi F + phi' u) over a BL test function that
# does NOT vanish at the endpoints carries an O(h^2) Euler-Maclaurin error. This
# returns that defect so the caller can add it to the BL residual, making it an
# unbiased weak form (O(h^6) with both terms below). With g(t) = phi(t) F + phi'(t) u,
#   EM_k = -dt^2/12 (g_k'(t_M) - g_k'(t_1)) + dt^4/720 (g_k'''(t_M) - g_k'''(t_1)).
# Endpoint derivatives g^(n) come from the Leibniz expansion (g_deriv_at_endpoint)
# using the total time derivatives F^(0..3) of the RHS along the trajectory.
#
# bl_phi_t1, bl_phi_tM: (K_bl x 5) raw phi^(0..4) at t_1 / t_M per BL test function.
# f_, dF_dt_, d2F_dt2_, d3F_dt3_: RHS and total time-derivative callables.
# Returns function(U, p, tt) -> (K_bl x D) matrix, or NULL when there are no BL rows.
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

# Analytic p-Jacobian of the boundary-layer EM defect (build_em_correction).
#
# EM_k(p) is a linear combination (via g_deriv_at_endpoint) of the trajectory
# derivatives F^(0..3) = f, dF/dt, d2F/dt2, d3F/dt3 at the two boundaries, with
# constant phi coefficients, plus a phi*u term with NO p-dependence. Hence
# dEM_k/dp_j is the SAME Leibniz combination applied to the p-derivatives
# dF^(m)/dp_j (and the phi*u term drops, i.e. u_vec = 0). No finite differences:
# the dF^(m)/dp callables come straight from the symbolic engine.
#
# dfdp_, dF1dp_, dF2dp_, dF3dp_: callables returning dF^(0..3)/dp at a single
#   point, each reshaping to D x J (d-fast, exactly like J_p).
# Returns function(U, p, tt) -> (K_bl x D x J), or NULL when there are no BL rows.
build_em_jacobian <- function(bl_phi_t1, bl_phi_tM,
                              dfdp_, dF1dp_, dF2dp_, dF3dp_,
                              dt, D, J, scale = 1.0) {
  if (is.null(bl_phi_t1) || nrow(bl_phi_t1) == 0L) return(NULL)
  K_bl  <- nrow(bl_phi_t1)
  c2    <- scale * dt^2 / 12
  c4    <- scale * dt^4 / 720
  zeroD <- numeric(D)

  function(U, p, tt) {
    M  <- nrow(U)
    t1 <- tt[1L]; tM <- tt[M]

    # dF^(0..3)/dp at a boundary point -> list of four D x J matrices.
    dp_at <- function(u_pt, t_pt) {
      inp <- matrix(c(p, u_pt, t_pt), ncol = 1L)
      list(matrix(as.vector(dfdp_(inp)),  D, J),
           matrix(as.vector(dF1dp_(inp)), D, J),
           matrix(as.vector(dF2dp_(inp)), D, J),
           matrix(as.vector(dF3dp_(inp)), D, J))
    }
    A1 <- dp_at(as.numeric(U[1L, ]), t1)
    AM <- dp_at(as.numeric(U[M,  ]), tM)

    jac <- array(0, c(K_bl, D, J))
    for (k in seq_len(K_bl)) {
      phi_t1 <- as.list(bl_phi_t1[k, ])
      phi_tM <- as.list(bl_phi_tM[k, ])
      for (j in seq_len(J)) {
        fd1 <- list(A1[[1]][, j], A1[[2]][, j], A1[[3]][, j], A1[[4]][, j])  # dF^(0..3)/dp_j @ t1
        fdM <- list(AM[[1]][, j], AM[[2]][, j], AM[[3]][, j], AM[[4]][, j])  # @ tM
        dg1_t1 <- g_deriv_at_endpoint(phi_t1, fd1, zeroD, order = 1L)
        dg1_tM <- g_deriv_at_endpoint(phi_tM, fdM, zeroD, order = 1L)
        dg3_t1 <- g_deriv_at_endpoint(phi_t1, fd1, zeroD, order = 3L)
        dg3_tM <- g_deriv_at_endpoint(phi_tM, fdM, zeroD, order = 3L)
        jac[k, , j] <- -c2 * (dg1_tM - dg1_t1) + c4 * (dg3_tM - dg3_t1)
      }
    }
    jac
  }
}

# ---- Boundary-layer IC system helpers ---------------------------------------

# Build the left boundary-layer system for one (r_c, n_bl) design: trap-weighted
# order-0/1 rows, the boundary vector B = psi_k(t_1), the raw endpoint
# derivatives phi^(orders) at t_1, and the data window the rows touch.
# orders = 0:1 suffices for the design-stage (no-EM) covariance; the estimator
# itself needs 0:4 for the Euler-Maclaurin correction.
build_ic_bl_system <- function(tt_vec, r_c, n_bl, orders = 0:4) {
  M <- length(tt_vec)
  bl_left <- lapply(orders, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_c, order = ord,
                               side = "left", n_bl = n_bl))
  K_bl <- nrow(bl_left[[1]])

  apply_trap <- function(V) {
    V[, 1] <- V[, 1] * 0.5
    V[, M] <- V[, M] * 0.5
    V
  }

  bl_phi_t1 <- matrix(0, nrow = K_bl, ncol = length(orders))
  for (i in seq_along(orders)) bl_phi_t1[, i] <- bl_left[[i]][, 1]

  win_cols <- which(colSums(abs(bl_left[[1]]) + abs(bl_left[[2]])) > 0)

  B <- bl_phi_t1[, 1]
  list(V_BL = apply_trap(bl_left[[1]]), Vp_BL = apply_trap(bl_left[[2]]),
       B = B, BtB = sum(B * B), K_bl = K_bl, bl_phi_t1 = bl_phi_t1,
       win_cols = win_cols)
}

# Per-data-point sensitivity of the K_bl boundary-layer equations,
#   X[(k,d),(m,c)] = d r_trap[k,d] / d U[m,c]
#                  = -dt (V_BL[k,m] J_u(t_m)[d,c] + Vp_BL[k,m] delta_dc),
# restricted to the window columns the BL rows actually touch. Layout is
# column-major vec: rows (d-1)*K_bl + k, columns (c-1)*n_win + i with i
# indexing win_cols; s2 carries the matching per-column noise variances.
# Bbold = I_D (x) B is the design matrix of the stacked system
# Bbold u0 = vec(r). X is shared by the GLS weights (Omega = X diag(s2) X^T)
# and, collapsed through the fixed-point Jacobian, by the noise-channel
# covariance of u0hat.
build_ic_noise_sensitivity <- function(bl, U, tt_vec, p, J_u, sig_vec, dt) {
  D    <- ncol(U)
  K_bl <- bl$K_bl
  KD   <- K_bl * D
  win  <- bl$win_cols
  nw   <- length(win)
  rowblk <- function(d) ((d - 1L) * K_bl + 1L):(d * K_bl)

  Bbold <- matrix(0, KD, D)
  for (d in seq_len(D)) Bbold[rowblk(d), d] <- bl$B

  X  <- matrix(0, KD, nw * D)
  s2 <- numeric(nw * D)
  for (i in seq_len(nw)) {
    m    <- win[i]
    Ju_m <- matrix(as.vector(J_u(c(p, U[m, ], tt_vec[m]))), D, D)
    for (cc in seq_len(D)) {
      col     <- (cc - 1L) * nw + i
      s2[col] <- sig_vec[cc]^2
      for (d in seq_len(D)) {
        X[rowblk(d), col] <- -dt * (bl$V_BL[, m] * Ju_m[d, cc] +
                                    (if (d == cc) bl$Vp_BL[, m] else 0))
      }
    }
  }
  list(X = X, s2 = s2, Bbold = Bbold)
}

# GLS (BLUE) weights for the stacked BL system. Omega = X diag(s2) X^T is the
# delta-method covariance of vec(r_trap): the errors of the K_bl equations are
# strongly correlated because they integrate the SAME noisy samples near the
# boundary, which the unweighted (OLS) combine ignores. W = Omega^{-1} (tiny
# ridge for near-duplicate rows). cov_design = (Bbold^T W Bbold)^{-1} is the
# no-EM a-priori GLS covariance of u0hat: it depends only on the design
# (r_c, n_bl), sigma, and J_u along the observed trajectory -- not on the
# solved u0 -- so it doubles as the design-selection criterion.
build_ic_gls_weights <- function(sens) {
  KD    <- nrow(sens$X)
  Omega <- sens$X %*% (sens$s2 * t(sens$X))
  Wm    <- solve(Omega + 1e-10 * mean(diag(Omega)) * diag(KD))
  BtW   <- crossprod(sens$Bbold, Wm)        # D x KD
  BtWB  <- BtW %*% sens$Bbold               # D x D
  list(Wm = Wm, BtW = BtW, BtWB = BtWB, cov_design = solve(BtWB))
}

# Analytic O(sigma^2) bias of the feasible-GLS fixed point (both channels).
#
# A second-order M-estimator expansion of the estimating equation
#   Bb' W(U) (vec r(U) - Bb u0 - vec EM(u0)) = 0
# gives E[u0hat] - u0* = b1 + b2 + O(sigma^4):
#   b1 = P E[dr]                      "f'' mean" channel: the residual is
#        quadratic in the noise through f, E[dr]_(k,d) =
#        -dt sum_m V_BL[k,m] * 0.5 sum_c J_uu[d,c,c](u_m) sigma_c^2;
#   b2 = -P c                         "weight feedback" channel: Omega is built
#        from the SAME noisy data as the residuals, so the weights correlate
#        with the errors they weight,
#        c = sum_n sigma_n^2 (dOmega/dU_n) (W R X)[, n],
#        R = I - (Bb + EMp) P', dOmega/dU_n = T_n S X' + X S T_n',
#        T_n = dX/dU_n through J_uu (central differences on J_u).
# The two channels PARTIALLY CANCEL (opposite signs on the validated systems);
# correcting either alone makes coverage worse — always subtract the sum.
# Validated (examples/validation/tmp_bias_debias_race.R, tmp_bias_o2_fd_check.R):
# each channel matches its MC counterpart to 1-2% on logistic at 5% noise
# (b1 +1.14e-2 vs +1.15e-2, b2 -8.7e-3 vs -8.6e-3 at (23,8)), and the D=3
# tensor algebra matches a dense FD-of-Omega implementation to ~4e-3.
# Plug-in evaluation at the noisy data costs only O(sigma^3).
#
# P is the IFT projection (Mbar^{-1} Bb' W with Mbar = Bb' W (Bb + EMp)),
# shared with the covariance channels. Returns list(b1, b2, b = b1 + b2).
build_ic_bias_o2 <- function(bl, sens, gls, P, EMp, U, tt_vec, p, J_u,
                             sig_vec, dt) {
  D    <- length(sig_vec)
  K    <- bl$K_bl
  KD   <- K * D
  win  <- bl$win_cols
  nw   <- length(win)
  X    <- sens$X
  s2   <- sens$s2
  Yvec <- sens$Bbold + EMp
  A    <- gls$Wm %*% (X - Yvec %*% (P %*% X))   # W (I - (Bb+EMp) P) X
  Eq   <- matrix(0, K, D)
  cvec <- numeric(KD)
  for (i in seq_len(nw)) {
    m <- win[i]
    if (all(bl$V_BL[, m] == 0)) next            # Vp-only column: linear in U
    u_m <- as.numeric(U[m, ])
    h0  <- 1e-5 * max(1, sqrt(sum(u_m^2)))
    Juu <- array(0, c(D, D, D))                 # [d,c,c'] = d2 f_d / du_c du_c'
    for (cp in seq_len(D)) {
      up <- u_m; up[cp] <- up[cp] + h0
      dn <- u_m; dn[cp] <- dn[cp] - h0
      Juu[, , cp] <- (matrix(as.vector(J_u(c(p, up, tt_vec[m]))), D, D) -
                      matrix(as.vector(J_u(c(p, dn, tt_vec[m]))), D, D)) / (2 * h0)
    }
    for (d in seq_len(D)) {
      tr_d <- 0
      for (cc in seq_len(D)) tr_d <- tr_d + Juu[d, cc, cc] * sig_vec[cc]^2
      Eq[, d] <- Eq[, d] - dt * bl$V_BL[, m] * (0.5 * tr_d)
    }
    for (cp in seq_len(D)) {
      n   <- (cp - 1L) * nw + i
      a_n <- A[, n]
      w   <- s2 * as.vector(crossprod(X, a_n))  # S X' a_n
      Tw  <- numeric(KD)
      XSy <- numeric(KD)
      for (d in seq_len(D)) {
        coef <- 0
        for (cc in seq_len(D)) coef <- coef + w[(cc - 1L) * nw + i] * Juu[d, cc, cp]
        Tw[((d - 1L) * K + 1L):(d * K)] <- -dt * bl$V_BL[, m] * coef
      }
      for (cc in seq_len(D)) {
        sdc <- 0
        for (d in seq_len(D))
          sdc <- sdc + Juu[d, cc, cp] *
                 sum(bl$V_BL[, m] * a_n[((d - 1L) * K + 1L):(d * K)])
        col <- (cc - 1L) * nw + i
        XSy <- XSy + X[, col] * (s2[col] * (-dt * sdc))
      }
      cvec <- cvec + sig_vec[cp]^2 * (Tw + XSy)
    }
  }
  b1 <- as.numeric(P %*% as.vector(Eq))
  b2 <- -as.numeric(P %*% cvec)
  list(b1 = b1, b2 = b2, b = b1 + b2)
}

# A-priori (r_c, n_bl) selection: sweep the candidate grid and keep the design
# minimizing the no-EM GLS variance, summed over states relative to the data
# noise (sum_d Var_d / sigma_d^2, scale-free across states). Computable before
# estimating anything -- no Monte-Carlo, no refitting; each candidate costs one
# small (K_bl D)^2 solve. Under GLS the variance is monotone (weakly) improving
# in n_bl and window width (unlike the OLS combine, where clustered test
# functions actively inflate it), so a small grid suffices.
select_ic_design <- function(U, tt_vec, p, J_u, sig_vec, dt,
                             r_c_grid, n_bl_grid, rc_cap) {
  r_cs <- sort(unique(pmin(pmax(as.integer(r_c_grid), 2L), rc_cap)))
  nbls <- sort(unique(pmax(as.integer(n_bl_grid), 1L)))
  rows <- vector("list", length(r_cs) * length(nbls))
  i <- 0L
  for (r_c in r_cs) for (n_bl in nbls) {
    res <- tryCatch({
      bl   <- build_ic_bl_system(tt_vec, r_c, n_bl, orders = 0:1)
      sens <- build_ic_noise_sensitivity(bl, U, tt_vec, p, J_u, sig_vec, dt)
      cov  <- build_ic_gls_weights(sens)$cov_design
      c(sum(diag(cov) / sig_vec^2), mean(sqrt(pmax(diag(cov), 0))))
    }, error = function(err) c(NA_real_, NA_real_))
    i <- i + 1L
    rows[[i]] <- data.frame(r_c = r_c, n_bl = n_bl,
                            obj = res[1], se_mean = res[2])
  }
  tab <- do.call(rbind, rows)
  if (!any(is.finite(tab$obj))) return(NULL)
  best <- which.min(tab$obj)
  list(table = tab, r_c = tab$r_c[best], n_bl = tab$n_bl[best])
}

#' Estimate u(0) via iterative defect-correction on left BL test functions
#'
#' Iterates the linear system
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
#' The K_bl equations share the same noisy samples, so their errors are
#' correlated with covariance \eqn{\Omega = X \mathrm{diag}(\sigma^2) X^T}
#' (\eqn{X} the equation/data sensitivity). With the default
#' \code{combine = "gls"} the equations are combined by GLS (the BLUE),
#' \deqn{u_0 = (\mathbf{B}^T \Omega^{-1} \mathbf{B})^{-1}
#'       \mathbf{B}^T \Omega^{-1} \mathrm{vec}(r),}
#' which is minimum-variance among all linear combinations of the equations
#' (validated at ~1.05-1.2x the window Cramer-Rao bound, vs up to ~25x for the
#' unweighted combine; examples/validation/). Additionally, when \code{n_bl}
#' is \code{NULL}, the design \code{(r_c, n_bl)} is selected a priori by
#' sweeping a small grid and minimizing the no-EM GLS variance
#' \eqn{(\mathbf{B}^T \Omega^{-1} \mathbf{B})^{-1}}, which depends only on the
#' design, \code{sigma}, and \code{J_u} along the observed trajectory.
#' \code{combine = "ols"} restores the legacy unweighted combine and its
#' \code{max(3, ceiling(r_c/8))} count heuristic.
#'
#' The feasible GLS estimator carries an \eqn{O(\sigma^2)} bias with two
#' partially cancelling channels (noise nonlinearity through \eqn{f''}, and
#' the feedback of the noise into the weights through \eqn{\Omega(U)}); at
#' 5\% noise on aggressive designs the net bias is ~0.3-0.5 of the (much
#' smaller) SE. With \code{debias = TRUE} (default) both channels are computed
#' analytically (see \code{build_ic_bias_o2}) and subtracted, restoring
#' centering at the cost of an \eqn{O(\sigma^3)} plug-in error (validated:
#' bias/SD 0.3-0.4 -> ~0.03 on logistic at 5\% noise with unchanged SD;
#' examples/validation/tmp_bias_debias_race.R). The correction is skipped
#' (with \code{debias_applied = FALSE}) when any component exceeds twice its
#' SE — a correction that large signals a regime where the expansion itself
#' is suspect.
#'
#' @param U Numeric matrix (M x D) of observed states.
#' @param f_,dF_dt_,d2F_dt2_,d3F_dt3_ Callable RHS and total time-derivative
#'   evaluators built from the symbolic engine.
#' @param tt Numeric vector (length M) of time points.
#' @param p Numeric parameter vector \eqn{\hat\theta} (held fixed).
#' @param r_c Integer; left BL window half-width. Under \code{combine = "gls"}
#'   with \code{n_bl = NULL} it seeds the design grid (see \code{r_c_grid});
#'   otherwise it is used as given.
#' @param n_bl Optional integer; number of left BL test functions. When
#'   \code{NULL} (default) and \code{combine = "gls"}, both \code{r_c} and
#'   \code{n_bl} are selected by the a-priori minimum-variance design sweep;
#'   when \code{NULL} under \code{combine = "ols"} the legacy heuristic
#'   \code{max(3, ceiling(r_c/8))} applies (small count, wide placement --
#'   under the unweighted combine, clustered test functions inflate
#'   Var(u0hat) by ~7-19\%; see examples/validation/). An explicit value
#'   always skips the sweep.
#' @param max_iter,tol Fixed-point iteration controls.
#' @param em_order Either 2 or 4.
#' @param combine \code{"gls"} (default) for the minimum-variance GLS combine
#'   of the BL equations, or \code{"ols"} for the legacy unweighted combine.
#'   GLS requires a valid \code{sigma} (and \code{J_u}) to build
#'   \eqn{\Omega}; it degrades to \code{"ols"} otherwise.
#' @param debias Logical (default \code{TRUE}). On the GLS path, subtract the
#'   analytic \eqn{O(\sigma^2)} bias (both channels; see Details). Ignored on
#'   the OLS path and on diverged solves.
#' @param r_c_grid,n_bl_grid Optional integer vectors of design candidates for
#'   the a-priori sweep (defaults \code{c(1, 2, 4) * r_c} capped at
#'   \code{floor((M-1)/2)} and deduped, and \code{c(3, 8, 16)}). Only used
#'   when \code{combine = "gls"} and \code{n_bl} is \code{NULL}.
#' @param J_u Callable state Jacobian \eqn{\partial f/\partial u}
#'   (\code{matrix(as.vector(J_u(c(p,u,t))), D, D)} with entry
#'   \eqn{[a,b] = \partial f_a/\partial u_b}). Used together with \code{sigma}
#'   to build the GLS weights and the noise-channel propagation
#'   \code{cov_u0} rather than the over-conservative LS-residual variance.
#' @param sigma Data-noise standard deviation, scalar or length-\eqn{D} (per
#'   state), used to scale the noise-channel \code{cov_u0}.
#' @param param_cov Optional J x J covariance \eqn{\hat C} of \eqn{\hat p}
#'   (e.g. the Fisher form \eqn{(G^T S^{-1} G)^{-1}}). When supplied,
#'   \code{cov_u0} additionally includes the parameter-uncertainty channel
#'   \eqn{S_p \hat C S_p^T} with \eqn{S_p = \partial \hat u_0 / \partial p}
#'   (the explained part of the law of total variance); see Details in the
#'   covariance comments below.
#' @return Named list with \code{U_hat} (U with row 1 replaced by
#'   \code{u0hat}), \code{u0hat}, \code{cov_u0} (D x D covariance of
#'   \eqn{\hat u_0}; the noise-channel propagation plus, when
#'   \code{param_cov} is supplied, the parameter channel
#'   \eqn{S_p \hat C S_p^T}, falling back to the
#'   LS-residual variance \eqn{s_d^2/B^TB} only when \code{sigma} is
#'   degenerate, and set to \code{NULL} when the iteration \code{diverged}
#'   (the covariance is then meaningless, so callers should fall back to the
#'   raw observation)),
#'   \code{cov_u0_resid} (the LS-residual variance, always returned for
#'   reference), \code{cov_u0_param} (the parameter channel alone, or
#'   \code{NULL} when not computed), \code{cov_method}
#'   (\code{"noise_propagation"}, \code{"ls_residual"}, either with a
#'   \code{"+param"} suffix when the parameter channel is included, or
#'   \code{"diverged"}), plus \code{combine} (the combine actually used,
#'   after any degradation), \code{design} (the design-sweep table with
#'   columns \code{r_c}, \code{n_bl}, \code{obj}, \code{se_mean}, or
#'   \code{NULL} when no sweep ran), \code{bias_o2} (the analytic
#'   \eqn{O(\sigma^2)} bias estimate of the UNcorrected solve, length-D, or
#'   \code{NULL} when not computed), \code{debias_applied} (whether
#'   \code{bias_o2} was subtracted from \code{u0hat}), \code{iters},
#'   \code{converged}, \code{diverged}, \code{u0_history}, \code{r_c},
#'   \code{n_bl}, \code{K_bl}, \code{em_order}.
#' @export
estimate_IC <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, r_c, J_u, sigma,
                        param_cov      = NULL,
                        n_bl           = NULL,
                        max_iter       = 20L,
                        tol            = 1e-12,
                        em_order       = c(4L, 2L),
                        combine        = c("gls", "ols"),
                        r_c_grid       = NULL,
                        n_bl_grid      = NULL,
                        debias         = TRUE) {
  em_order <- as.integer(em_order[1])
  if (!em_order %in% c(2L, 4L)) {
    stop("em_order must be 2 or 4", call. = FALSE)
  }
  combine <- match.arg(combine)

  M      <- nrow(U)
  D      <- ncol(U)
  J      <- length(p)
  tt_vec <- as.vector(tt)
  dt     <- mean(diff(tt_vec))

  rc_cap <- floor((M - 1L) / 2L)
  r_c    <- min(r_c, rc_cap)
  r_c_in <- r_c

  use_noiseprop <- length(sigma) %in% c(1L, D) && all(is.finite(sigma))
  sig_vec <- if (use_noiseprop) {
    if (length(sigma) == 1L) rep(sigma, D) else as.numeric(sigma)
  } else NULL
  # GLS needs a valid sigma (with J_u) to build Omega; degrade to the legacy
  # unweighted combine otherwise.
  if (combine == "gls" && !use_noiseprop) combine <- "ols"

  # A-priori design selection (GLS only): when no explicit n_bl is given,
  # sweep (r_c, n_bl) candidates and keep the design with the minimum no-EM
  # GLS variance (select_ic_design). The passed r_c (pipeline radius) seeds
  # the grid. Under the OLS combine the legacy heuristic below applies.
  design_table <- NULL
  if (combine == "gls" && is.null(n_bl)) {
    if (is.null(r_c_grid))  r_c_grid  <- c(1L, 2L, 3L, 4L) * r_c
    if (is.null(n_bl_grid)) n_bl_grid <- c(3L, 8L, 16L)
    sel <- tryCatch(
      select_ic_design(U, tt_vec, p, J_u, sig_vec, dt,
                       r_c_grid, n_bl_grid, rc_cap),
      error = function(err) NULL)
    if (!is.null(sel)) {
      design_table <- sel$table
      r_c          <- sel$r_c
      n_bl         <- sel$n_bl
    }
  }

  # Legacy count heuristic (OLS combine, or sweep unavailable): small count,
  # wide placement. Under the UNWEIGHTED combine, clustered test functions add
  # redundant info with full noise and inflate Var(u0hat) by ~7-19%
  # (examples/validation/); under GLS the sweep above governs instead.
  n_bl <- if (is.null(n_bl)) max(3L, as.integer(ceiling(r_c / 8)))
          else                max(1L, as.integer(n_bl))

  build_system <- function(r_c_use, n_bl_use) {
    bl <- build_ic_bl_system(tt_vec, r_c_use, n_bl_use, orders = 0:4)
    sens <- if (use_noiseprop) tryCatch(
      build_ic_noise_sensitivity(bl, U, tt_vec, p, J_u, sig_vec, dt),
      error = function(err) NULL) else NULL
    gls <- if (combine == "gls" && !is.null(sens)) tryCatch(
      build_ic_gls_weights(sens),
      error = function(err) NULL) else NULL
    list(bl = bl, sens = sens, gls = gls)
  }

  sys <- build_system(r_c, n_bl)
  if (combine == "gls" && is.null(sys$gls)) {
    # GLS weights unavailable (singular Omega, J_u failure, ...). Revert to
    # the FULL legacy path: a sweep-selected design (large n_bl, wide window)
    # is GLS-specific and would be inefficient under the unweighted combine.
    combine <- "ols"
    if (!is.null(design_table)) {
      design_table <- NULL
      r_c  <- r_c_in
      n_bl <- max(3L, as.integer(ceiling(r_c / 8)))
      sys  <- build_system(r_c, n_bl)
    }
  }

  bl        <- sys$bl
  sens      <- sys$sens
  gls       <- sys$gls
  V_BL      <- bl$V_BL
  Vp_BL     <- bl$Vp_BL
  bl_phi_t1 <- bl$bl_phi_t1
  B         <- bl$B
  BtB       <- bl$BtB
  K_bl      <- bl$K_bl

  if (!is.finite(BtB) || BtB < .Machine$double.eps) {
    u0_obs <- as.numeric(U[1, ])
    U_hat <- U; U_hat[1, ] <- u0_obs
    return(list(U_hat = U_hat, u0hat = u0_obs, cov_u0 = NULL,
                iters = 0L, converged = FALSE, diverged = FALSE,
                u0_history = matrix(u0_obs, nrow = 1),
                r_c = r_c, n_bl = n_bl, K_bl = K_bl, em_order = em_order,
                combine = combine, design = design_table,
                bias_o2 = NULL, debias_applied = FALSE
              ))
  }

  compute_r_trap <- function(U_in, p_use = p) {
    input  <- rbind(matrix(rep(p_use, M), nrow = J), t(U_in), matrix(tt_vec, nrow = 1L))
    F_eval <- f_(input)
    -dt * (V_BL %*% F_eval) - dt * (Vp_BL %*% U_in)
  }
  r_trap <- compute_r_trap(U)

  c2 <- dt^2 / 12
  c4 <- if (em_order >= 4L) dt^4 / 720 else 0

  em_correction <- function(u0_curr, p_use = p) {
    u_t1   <- as.vector(u0_curr)
    inp_t1 <- matrix(c(p_use, u_t1, tt_vec[1]), ncol = 1L)
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

  # Projection of the stacked K_bl x D system onto u0. GLS (BLUE) weights the
  # equations by W = Omega^{-1}:  u0 = (Bb' W Bb)^{-1} Bb' W vec(rhs); the OLS
  # branch is the legacy unweighted combine B^T rhs / B^T B. vec layout is
  # column-major (k fast, d slow), matching Bbold in build_ic_noise_sensitivity.
  proj <- if (!is.null(gls)) {
    function(rhs) as.numeric(solve(gls$BtWB, gls$BtW %*% as.vector(rhs)))
  } else {
    function(rhs) as.numeric(crossprod(B, rhs) / BtB)
  }

  # Track the best iterate by residual norm; the iteration contracts only when
  # the EM-correction map has spectral radius < 1 (kappa = O(h^2)). For
  # under-sampled stiff systems the contraction can fail and iterates blow up.
  # Returning the best-seen iterate keeps us no worse than the uncorrected LS.
  # Residual over the K_bl x D system, in the same metric the projection
  # minimizes: W-weighted for GLS, Frobenius for OLS.
  residual_norm <- function(u0_curr) {
    EM <- em_correction(u0_curr)
    e  <- outer(B, as.numeric(u0_curr)) - (r_trap - EM)   # K_bl x D
    if (!is.null(gls)) {
      ev <- as.vector(e)
      sqrt(sum(ev * (gls$Wm %*% ev)))
    } else {
      sqrt(sum(e * e))
    }
  }

  u0       <- proj(r_trap)
  u0_hist  <- list(u0)
  best_u0  <- u0
  best_res <- if (all(is.finite(u0))) residual_norm(u0) else Inf
  iters    <- 0L
  converged <- FALSE
  diverged  <- FALSE

  for (it in seq_len(max_iter)) {
    iters  <- it
    EM     <- em_correction(u0)
    rhs    <- r_trap - EM
    u0_new <- proj(rhs)
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

  # Implicit-function-theorem machinery shared by the noise and parameter
  # channels. u0hat is a deterministic function of the data: differentiating
  # the converged fixed-point relation
  #   (C Bbold) u0 + C vec(EM(u0)) = C vec(r_trap(U)),
  # with collapse matrix C = Bbold^T (OLS) or C = Bbold^T W (GLS), gives
  #   P        = (C Bbold + C EMp)^{-1} C,        EMp = dvec(EM)/du0,
  #   du0/dU   = P X            (X from build_ic_noise_sensitivity),
  #   du0/dp_j = P dvec(r_trap - EM)/dp_j.
  # EMp is taken by central differences (as before). For OLS this reproduces
  # the previous (B^T B I + A)^{-1} form exactly: C Bbold = BtB I_D and
  # C EMp = A.
  I_D <- diag(D)
  KD  <- K_bl * D
  EMp <- tryCatch({
    EMp <- matrix(0, KD, D)
    h <- 1e-6 * max(1, sqrt(sum(u0^2)))
    for (e_i in seq_len(D)) {
      up <- u0
      up[e_i] <- up[e_i] + h
      dn <- u0
      dn[e_i] <- dn[e_i] - h
      EMp[, e_i] <- as.vector((em_correction(up) - em_correction(dn)) / (2 * h))
    }
    EMp
  }, error = function(err) NULL)
  P <- if (!is.null(EMp)) tryCatch({
    if (!is.null(gls)) {
      solve(gls$BtWB + gls$BtW %*% EMp, gls$BtW)
    } else {
      Bbold <- if (!is.null(sens)) sens$Bbold else {
        Bb <- matrix(0, KD, D)
        for (d in seq_len(D)) Bb[((d - 1L) * K_bl + 1L):(d * K_bl), d] <- B
        Bb
      }
      solve(BtB * I_D + crossprod(Bbold, EMp), t(Bbold))
    }
  }, error = function(err) NULL) else NULL

  # Noise channel (preferred when J_u and sigma are supplied):
  #   Cov_noise(u0|p) = (P X) diag(s2) (P X)^T,
  # the delta-method propagation of the per-point data noise through the
  # estimator -- only the BL-window samples contribute (X is built on
  # win_cols). Calibrated ~95% coverage; for OLS identical to the previous
  # per-m accumulation, for GLS additionally the (near-)minimum-variance
  # combine: with EMp -> 0 it telescopes to (Bbold^T Omega^{-1} Bbold)^{-1}
  # (Gauss-Markov; see build_ic_gls_weights). NOTE: the weights W are treated
  # as fixed (feasible GLS); their dependence on U through J_u is one order
  # down in sigma (MC-calibrated 0.97-1.10).
  # Falls back to the LS-residual variance only for degenerate sigma.
  cov_u0_noise <- if (!is.null(sens) && !is.null(P)) tryCatch({
    G <- P %*% sens$X
    G %*% (sens$s2 * t(G))
  }, error = function(err) NULL) else NULL

  # Parameter channel: the noise propagation above is conditional on p = phat,
  # i.e. only the "unexplained" part of the law of total variance
  #   Var(u0hat) = E_p[Var(u0hat | p)] + Var_p(E[u0hat | p]).
  # The "explained" part comes from the dependence of u0hat on phat.
  # Differentiating the fixed-point relation w.r.t. p_j at fixed data U gives
  #   S_p[, j] = P dvec(r_trap - EM)/dp_j,
  # with the p-derivatives taken by central differences (EM's p-dependence
  # runs through total time derivatives of f up to order 3 — messy
  # analytically, trivial numerically), and
  #   Cov_param(u0) = S_p Cov(phat) S_p^T.
  # The cross term between the two channels (phat is estimated from the same
  # data U) is ignored, as in wendy_erts's parameter fold.
  cov_u0_param <- if (!is.null(param_cov) && !is.null(P)) tryCatch({
    rhs_of_p <- function(p_use)
      as.vector(compute_r_trap(U, p_use) - em_correction(u0, p_use))
    S_p <- matrix(0, D, J)
    for (j in seq_len(J)) {
      hj <- 1e-6 * max(1, abs(p[j]))
      pj_up <- p; pj_up[j] <- pj_up[j] + hj
      pj_dn <- p; pj_dn[j] <- pj_dn[j] - hj
      S_p[, j] <- P %*% ((rhs_of_p(pj_up) - rhs_of_p(pj_dn)) / (2 * hj))
    }
    S_p %*% param_cov %*% t(S_p)
  }, error = function(err) NULL) else NULL

  cov_u0     <- if (!is.null(cov_u0_noise)) cov_u0_noise else cov_u0_resid
  cov_method <- if (!is.null(cov_u0_noise)) "noise_propagation" else "ls_residual"
  if (!is.null(cov_u0) && !is.null(cov_u0_param)) {
    cov_u0     <- cov_u0 + cov_u0_param
    cov_method <- paste0(cov_method, "+param")
  }

  # A diverged fixed point means u0hat is only the best-seen iterate, not a
  # solution of the BL system; its covariance — which assumes the converged
  # fixed-point relation — is then meaningless (empirically coverage collapses
  # to ~5-19% on fast-spiking stiff systems such as Hindmarsh-Rose, vs ~95% for
  # non-diverged solves, including coarse-n cases that merely fail to hit `tol`).
  # Invalidate cov_u0 so the caller (solveWendy) falls back to the raw
  # observation instead of seeding the state filter from a bad prior. We gate on
  # `diverged` only, NOT `!converged`: non-converged-but-bounded solves stay well
  # calibrated, so discarding them would needlessly throw away good estimates.
  if (diverged) {
    cov_u0       <- NULL
    cov_u0_param <- NULL
    cov_method   <- "diverged"
  }

  # --- O(sigma^2) debias (GLS path) ----------------------------------------
  # Subtract the analytic second-order bias b1 + b2 (build_ic_bias_o2; both
  # channels — correcting either alone WORSENS coverage because they partially
  # cancel). The covariance is left unchanged: the correction shifts the mean
  # at O(sigma^2) and perturbs the variance only at higher order. Gated on a
  # sane magnitude (each |b_d| <= 2 SE): in every validated regime the net
  # bias is well under one SE, so a larger value signals a regime (stiff,
  # under-resolved) where the expansion itself is no longer trustworthy and
  # the uncorrected estimate is safer.
  bias_o2        <- NULL
  debias_applied <- FALSE
  if (isTRUE(debias) && !diverged && !is.null(gls) && !is.null(sens) &&
      !is.null(P) && !is.null(EMp) && !is.null(cov_u0_noise)) {
    bias_o2 <- tryCatch(
      build_ic_bias_o2(bl, sens, gls, P, EMp, U, tt_vec, p, J_u,
                       sig_vec, dt)$b,
      error = function(err) NULL)
    if (!is.null(bias_o2) && all(is.finite(bias_o2))) {
      se_gate <- sqrt(pmax(diag(cov_u0_noise), 0))
      if (all(abs(bias_o2) <= 2 * se_gate)) {
        u0             <- u0 - bias_o2
        debias_applied <- TRUE
      }
    }
  }

  U_hat <- U; U_hat[1, ] <- u0

  list(
    U_hat          = U_hat,
    u0hat          = u0,
    cov_u0         = cov_u0,
    cov_u0_resid   = cov_u0_resid,
    cov_u0_param   = cov_u0_param,
    cov_method     = cov_method,
    combine        = combine,
    design         = design_table,
    bias_o2        = bias_o2,
    debias_applied = debias_applied,
    iters          = iters,
    converged      = converged,
    diverged       = diverged,
    u0_history     = do.call(rbind, u0_hist),
    r_c            = r_c,
    n_bl           = n_bl,
    K_bl           = K_bl,
    em_order       = em_order
  )
}

# ---- Pointwise weak-form GLS state estimation --------------------------------

# 5 x D Euler-Maclaurin coefficient matrix at one window endpoint. The EM
# defect is LINEAR in phi^(0..4) (g_deriv_at_endpoint is linear in the phi
# scalars), so a K_bl x D defect block is bl_phi_t1 %*% C: one matmul replaces
# K_bl rowwise Leibniz loops (exact; validated to ~3e-17 against the rowwise
# form in examples/validation/tmp_state_gls_parallel.R). EM(4) coefficients.
em_coef_matrix <- function(fd, u_vec, dt, D) {
  C <- matrix(0, 5L, D)
  for (j in 1:5) {
    ph <- as.list(numeric(5))
    ph[[j]] <- 1
    C[j, ] <- dt^2 / 12 * g_deriv_at_endpoint(ph, fd, u_vec, 1L) -
              dt^4 / 720 * g_deriv_at_endpoint(ph, fd, u_vec, 3L)
  }
  C
}

# Analytic O(sigma^2) bias of one anchor's joint (two-sided) GLS system —
# the same two channels as build_ic_bias_o2 (f'' mean + weight feedback),
# generalized to the stacked side systems of wendy_gls_state: rows live in
# per-side blocks (offsets by side, K rows per state within a side), columns
# in the union-window layout (cc-1)*nU + gi, and each side's Jacobian carries
# its time-reversal sign, which applies to J_uu exactly as to J_u. Only the
# psi-channel (V_BL) differentiates w.r.t. the data; the psi' channel is
# linear. Juu_all is the per-sample d2f/du2 array precomputed once per
# trajectory. Returns the length-D net bias (b1 + b2).
build_state_bias_o2 <- function(sides, gidx_all, V, s2, Wm, Yvec, P,
                                Juu_all, sig_vec, dt, D) {
  nU   <- length(gidx_all)
  KDt  <- nrow(V)
  A    <- Wm %*% (V - Yvec %*% (P %*% V))
  offs <- integer(length(sides)); off <- 0L
  for (si in seq_along(sides)) {
    offs[si] <- off; off <- off + nrow(sides[[si]]$Bb)
  }
  Eq   <- numeric(KDt)
  cvec <- numeric(KDt)
  for (gi in seq_len(nU)) {
    g <- gidx_all[gi]
    touch <- list()
    for (si in seq_along(sides)) {
      i_s <- match(g, sides[[si]]$gidx)
      if (!is.na(i_s)) touch[[length(touch) + 1L]] <- c(si = si, i = i_s)
    }
    if (!length(touch)) next
    for (tc in touch) {                       # channel 1: E[dr]
      s <- sides[[tc[["si"]]]]
      K <- nrow(s$Bb) %/% D; o <- offs[tc[["si"]]]
      vcol <- s$V_BL[, tc[["i"]]]
      for (d in seq_len(D)) {
        tr_d <- 0
        for (cc in seq_len(D))
          tr_d <- tr_d + s$sgn * Juu_all[g, d, cc, cc] * sig_vec[cc]^2
        rows <- o + (d - 1L) * K + seq_len(K)
        Eq[rows] <- Eq[rows] - dt * vcol * (0.5 * tr_d)
      }
    }
    for (cp in seq_len(D)) {                  # channel 2: c vector
      a_n <- A[, (cp - 1L) * nU + gi]
      w_g <- numeric(D)                       # (S V' a_n) at the gi entries
      for (cc in seq_len(D)) {
        col <- (cc - 1L) * nU + gi
        w_g[cc] <- s2[col] * sum(V[, col] * a_n)
      }
      Tw  <- numeric(KDt)
      ycc <- numeric(D)
      for (tc in touch) {
        s <- sides[[tc[["si"]]]]
        K <- nrow(s$Bb) %/% D; o <- offs[tc[["si"]]]
        vcol <- s$V_BL[, tc[["i"]]]
        for (d in seq_len(D)) {
          rows <- o + (d - 1L) * K + seq_len(K)
          coef <- 0
          for (cc in seq_len(D))
            coef <- coef + w_g[cc] * (s$sgn * Juu_all[g, d, cc, cp])
          Tw[rows] <- Tw[rows] - dt * vcol * coef
          va <- sum(vcol * a_n[rows])
          for (cc in seq_len(D))
            ycc[cc] <- ycc[cc] - dt * (s$sgn * Juu_all[g, d, cc, cp]) * va
        }
      }
      XSy <- numeric(KDt)
      for (cc in seq_len(D)) {
        col <- (cc - 1L) * nU + gi
        XSy <- XSy + V[, col] * (s2[col] * ycc[cc])
      }
      cvec <- cvec + sig_vec[cp]^2 * (Tw + XSy)
    }
  }
  as.numeric(P %*% Eq) - as.numeric(P %*% cvec)
}

#' Pointwise weak-form GLS state estimation ("weak-form local smoother")
#'
#' Estimates the full state trajectory by solving, at EVERY time point
#' \eqn{t_k}, the same joint GLS boundary-layer system that
#' \code{\link{estimate_IC}} solves at \eqn{t_1}: a forward window
#' \eqn{[t_k, t_k+W]} of left-BL test functions anchored at \eqn{t_k}, plus a
#' backward window \eqn{[t_k-W, t_k]} mapped onto the same machinery by time
#' reversal (\eqn{s = t_k - t}: \eqn{f \to -f}, \eqn{J_u \to -J_u}, total time
#' derivatives \eqn{F^{(m)} \to (-1)^{m+1} F^{(m)}}; valid for non-autonomous
#' \eqn{f} because every evaluation uses the original sample time). Both sides
#' are stacked into ONE GLS system with the joint error covariance
#' \eqn{\Omega} over the union window, so the correlation through the shared
#' anchor sample is exact (no independence approximation), and each anchor's
#' \eqn{(r_c, n_{bl})} design is selected a priori per side by the minimum
#' no-EM GLS variance.
#'
#' There is NO recursion: each anchor is an independent solve, so errors do
#' not accumulate, failures stay local (a failed anchor falls back to the raw
#' observation), and the reported per-anchor covariance is the calibrated
#' delta-method noise propagation (plus the parameter channel when
#' \code{param_cov} is supplied). Validated (examples/validation/
#' tmp_state_gls_boost.R, tmp_state_gls_lorenz.R): matches the ERTS smoother's
#' RMSE on logistic while ERTS over-covers ~2x, and beats it on Lorenz where
#' the EKF's first-order covariance propagation turns 2-3x overconfident
#' (coverage 30-64%) while this estimator stays ~calibrated. Trade-offs: no
#' process-noise tuning knob, but it assumes the dynamics hold exactly within
#' each window, requires a uniform grid, and costs ~2M small GLS solves
#' (anchors are independent and parallelize trivially; kept serial here).
#'
#' @param U Numeric matrix (M x D) of noisy observations.
#' @param f_,J_u,dF_dt_,d2F_dt2_,d3F_dt3_ Callable RHS, state Jacobian and
#'   total time-derivative evaluators built from the symbolic engine.
#' @param tt Numeric vector (length M) of UNIFORM time points.
#' @param p Numeric parameter vector \eqn{\hat\theta} (held fixed).
#' @param sigma Data-noise SD, scalar or length-D; required (the GLS weights
#'   and the covariance are built from it).
#' @param param_cov Optional J x J covariance of \eqn{\hat p}; adds the
#'   parameter channel \eqn{S_p \hat C S_p^T} per anchor exactly as in
#'   \code{\link{estimate_IC}}.
#' @param r_c_grid,n_bl_grid Design candidates for the per-anchor sweep
#'   (defaults \code{c(8, 12, 16, 20, 24)} and \code{c(8, 16)}). r_c is
#'   floored at 8: below the EM(4)-validity floor the defect correction may
#'   not contract and the ill-conditioned \eqn{\Omega} makes the a-priori
#'   criterion spuriously optimistic.
#' @param design_stride Integer. In the interior (where every candidate fits)
#'   the design is selected at every \code{design_stride}-th anchor and reused
#'   between knots (J_u varies smoothly); near the edges, where the feasible
#'   set changes with every anchor, selection is per-anchor.
#' @param design_objective \code{"total"} (default) or \code{"noise"}. With
#'   \code{"total"} and a supplied \code{param_cov}, the per-anchor design is
#'   chosen by minimizing the noise PLUS parameter-channel variance: where
#'   \eqn{|\partial f/\partial p|} is large (attractors, plateaus) a wide
#'   dynamics-leaning window is noise-optimal but imports parameter error, and
#'   the total objective backs off to a shorter window. \code{"noise"}
#'   reproduces the noise-only criterion (and is the behavior whenever
#'   \code{param_cov} is \code{NULL}).
#' @param max_iter,tol Fixed-point iteration controls (one EM evaluation per
#'   iteration; best-iterate divergence guard as in \code{estimate_IC}).
#' @param debias Logical (default \code{TRUE}). Subtract the analytic
#'   \eqn{O(\sigma^2)} bias (f'' mean + weight-feedback channels, both — they
#'   partially cancel) from each anchor's estimate, exactly as in
#'   \code{\link{estimate_IC}}; gated per anchor at twice the noise-channel
#'   SE. \eqn{\partial^2 f/\partial u^2} is taken by central differences on
#'   \code{J_u}, once per trajectory.
#' @return Named list with \code{U_star} (M x D state estimate),
#'   \code{P_smooth} (M x D x D per-anchor covariance, noise channel plus
#'   parameter channel when folded), \code{se} (M x D), \code{fallback}
#'   (logical M; anchors that failed and returned the raw observation with
#'   \eqn{\mathrm{diag}(\sigma^2)}), \code{debiased} (logical M; anchors where
#'   the \eqn{O(\sigma^2)} correction was applied), \code{cov_method}, and
#'   \code{design} (per-anchor selected \code{(r_c, n_bl)} per side, NA where
#'   a side has no feasible window).
#' @export
wendy_gls_state <- function(U, f_, J_u, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, sigma,
                            param_cov        = NULL,
                            r_c_grid         = NULL,
                            n_bl_grid        = NULL,
                            design_stride    = 8L,
                            design_objective = c("total", "noise"),
                            max_iter         = 12L,
                            tol              = 1e-8,
                            debias           = TRUE) {
  M      <- nrow(U)
  D      <- ncol(U)
  J      <- length(p)
  tt_vec <- as.vector(tt)
  dt     <- mean(diff(tt_vec))

  design_objective <- match.arg(design_objective)
  if (!(length(sigma) %in% c(1L, D)) || !all(is.finite(sigma)))
    stop("wendy_gls_state requires a finite sigma (scalar or per-state).",
         call. = FALSE)
  sig_vec <- if (length(sigma) == 1L) rep(sigma, D) else as.numeric(sigma)

  rc_cap <- floor((M - 1L) / 2L)
  if (rc_cap < 8L)
    stop("wendy_gls_state needs M >= 17 samples (r_c floor of 8).", call. = FALSE)
  if (is.null(r_c_grid))  r_c_grid  <- c(8L, 12L, 16L, 20L, 24L)
  if (is.null(n_bl_grid)) n_bl_grid <- c(8L, 16L)
  r_cs  <- sort(unique(pmin(pmax(as.integer(r_c_grid), 8L), rc_cap)))
  nbls  <- sort(unique(pmax(as.integer(n_bl_grid), 1L)))
  cands <- expand.grid(r_c = r_cs, n_bl = nbls)

  # BL systems per candidate on a synthetic uniform grid (only dt and the
  # window length matter), truncated to the touched window.
  bls <- lapply(seq_len(nrow(cands)), function(ci) {
    r_c <- cands$r_c[ci]; n_bl <- cands$n_bl[ci]
    ttw <- seq(0, by = dt, length.out = min(M, 2L * r_c + 2L * n_bl + 2L))
    bl  <- build_ic_bl_system(ttw, r_c, n_bl, orders = 0:4)
    bl$win   <- max(bl$win_cols)
    bl$V_BL  <- bl$V_BL[,  seq_len(bl$win), drop = FALSE]
    bl$Vp_BL <- bl$Vp_BL[, seq_len(bl$win), drop = FALSE]
    bl
  })
  wins    <- vapply(bls, `[[`, integer(1), "win")
  min_win <- min(wins)
  max_win <- max(wins)

  # Backward windows via time reversal s = t_k - t. The reversed trajectory
  # solves u' = -f, its total time derivatives pick up signs (-1)^(m+1) and
  # J_u flips sign; with all evaluations at ORIGINAL sample times this is
  # exact for non-autonomous f as well.
  cal_fwd <- list(f_ = f_, J_u = J_u, dF_dt_ = dF_dt_,
                  d2F_dt2_ = d2F_dt2_, d3F_dt3_ = d3F_dt3_)
  cal_rev <- list(f_ = function(x) -f_(x), J_u = function(x) -J_u(x),
                  dF_dt_ = dF_dt_,
                  d2F_dt2_ = function(x) -d2F_dt2_(x),
                  d3F_dt3_ = function(x) -d3F_dt3_(x))

  # Per-sample f and J_u are anchor-independent: evaluate once per trajectory.
  input_all <- rbind(matrix(rep(p, M), nrow = J), t(U), matrix(tt_vec, nrow = 1L))
  F_all  <- f_(input_all)                                   # M x D
  Ju_all <- array(0, c(M, D, D))
  for (m in seq_len(M))
    Ju_all[m, , ] <- matrix(as.vector(J_u(c(p, U[m, ], tt_vec[m]))), D, D)

  # Per-sample d2f/du2 (central differences on J_u), anchor-independent like
  # Ju_all; feeds the per-anchor O(sigma^2) debias (build_state_bias_o2).
  Juu_all <- if (isTRUE(debias)) tryCatch({
    arr <- array(0, c(M, D, D, D))
    for (m in seq_len(M)) {
      u_m <- as.numeric(U[m, ])
      h0  <- 1e-5 * max(1, sqrt(sum(u_m^2)))
      for (cp in seq_len(D)) {
        up <- u_m; up[cp] <- up[cp] + h0
        dn <- u_m; dn[cp] <- dn[cp] - h0
        arr[m, , , cp] <-
          (matrix(as.vector(J_u(c(p, up, tt_vec[m]))), D, D) -
           matrix(as.vector(J_u(c(p, dn, tt_vec[m]))), D, D)) / (2 * h0)
      }
    }
    arr
  }, error = function(err) NULL) else NULL

  # One side (forward or time-reversed window) of the anchored system.
  side_system <- function(k, side, ci) {
    bl  <- bls[[ci]]; win <- bl$win; K <- bl$K_bl
    sgn <- if (side == "R") 1 else -1
    idx <- if (side == "R") k:(k + win - 1L) else k:(k - win + 1L)
    Uw  <- U[idx, , drop = FALSE]
    F_w <- sgn * F_all[idx, , drop = FALSE]
    r_mat <- -dt * (bl$V_BL %*% F_w) - dt * (bl$Vp_BL %*% Uw)
    # X[(d-1)K+q, (c-1)win+i] = -dt (V[q,i] Ju_i[d,c] + Vp[q,i] delta_dc)
    X <- matrix(0, K * D, win * D)
    for (cc in seq_len(D)) for (d in seq_len(D)) {
      ju  <- sgn * Ju_all[idx, d, cc]
      blk <- -dt * (bl$V_BL * matrix(ju, K, win, byrow = TRUE))
      if (d == cc) blk <- blk - dt * bl$Vp_BL
      X[((d - 1L) * K + 1L):(d * K), ((cc - 1L) * win + 1L):(cc * win)] <- blk
    }
    Bb <- matrix(0, K * D, D)
    for (d in seq_len(D)) Bb[((d - 1L) * K + 1L):(d * K), d] <- bl$B
    cal <- if (side == "R") cal_fwd else cal_rev
    em <- function(u0c, p_use = p) {
      inp <- matrix(c(p_use, u0c, tt_vec[k]), ncol = 1L)
      fd  <- list(as.vector(cal$f_(inp)), as.vector(cal$dF_dt_(inp)),
                  as.vector(cal$d2F_dt2_(inp)), as.vector(cal$d3F_dt3_(inp)))
      bl$bl_phi_t1 %*% em_coef_matrix(fd, u0c, dt, D)
    }
    r_of_p <- function(p_use) {     # parameter channel: r at perturbed p
      input <- rbind(matrix(rep(p_use, win), nrow = J), t(Uw),
                     matrix(tt_vec[idx], nrow = 1L))
      -dt * (bl$V_BL %*% (sgn * f_(input))) - dt * (bl$Vp_BL %*% Uw)
    }
    list(Bb = Bb, r = as.vector(r_mat), X = X, gidx = idx, em = em,
         r_of_p = r_of_p, V_BL = bl$V_BL, sgn = sgn)
  }

  want_param <- !is.null(param_cov)

  # A-priori one-sided design objective (no EM): noise-relative TOTAL variance
  # of the side's GLS solve. The noise term is estimate_IC's select_ic_design
  # criterion; when param_cov is supplied (and design_objective = "total") the
  # parameter channel S_p C S_p^T is added, because it is the dominant term
  # where |df/dp| is large (e.g. near attractors/plateaus): there a wide,
  # dynamics-leaning window is noise-optimal but imports maximal parameter
  # error, and the total-variance criterion correctly backs off to a shorter
  # window closer to pure local averaging (which has no p-sensitivity).
  design_obj <- function(k, side, ci) {
    s   <- side_system(k, side, ci)
    win <- length(s$gidx)
    s2  <- rep(sig_vec^2, each = win)
    Omega <- s$X %*% (s2 * t(s$X))
    Wm    <- solve(Omega + 1e-10 * mean(diag(Omega)) * diag(nrow(s$X)))
    BtW   <- crossprod(s$Bb, Wm)
    BtWB  <- BtW %*% s$Bb
    v <- diag(solve(BtWB))
    if (want_param && identical(design_objective, "total")) {
      P   <- solve(BtWB, BtW)
      S_p <- matrix(0, D, J)
      for (j in seq_len(J)) {
        hj <- 1e-6 * max(1, abs(p[j]))
        pj_up <- p; pj_up[j] <- pj_up[j] + hj
        pj_dn <- p; pj_dn[j] <- pj_dn[j] - hj
        S_p[, j] <- P %*% (as.vector(s$r_of_p(pj_up) - s$r_of_p(pj_dn)) / (2 * hj))
      }
      v <- v + diag(S_p %*% param_cov %*% t(S_p))
    }
    sum(v / sig_vec^2)
  }
  select_side <- function(k, side) {
    avail <- if (side == "R") M - k + 1L else k
    feas  <- which(wins <= avail)
    if (!length(feas)) return(NA_integer_)
    best_ci <- NA_integer_; best_v <- Inf
    for (ci in feas) {
      v <- tryCatch(design_obj(k, side, ci), error = function(err) Inf)
      if (is.finite(v) && v < best_v) { best_v <- v; best_ci <- ci }
    }
    best_ci
  }
  # Near the edges the feasible set changes with every anchor -> select
  # per-anchor; in the interior every candidate fits and J_u varies smoothly
  # -> select at stride knots and reuse the nearest knot's design.
  sel_table <- function(side) {
    out   <- rep(NA_integer_, M)
    avail <- if (side == "R") M - seq_len(M) + 1L else seq_len(M)
    edge  <- which(avail >= min_win & avail < max_win)
    intr  <- which(avail >= max_win)
    for (k in edge) out[k] <- select_side(k, side)
    if (length(intr)) {
      stride <- max(1L, as.integer(design_stride))
      knots  <- intr[unique(c(seq(1L, length(intr), by = stride), length(intr)))]
      kn_ci  <- vapply(knots, select_side, integer(1), side = side)
      for (k in intr) out[k] <- kn_ci[which.min(abs(knots - k))]
    }
    out
  }
  ciL <- sel_table("L")
  ciR <- sel_table("R")

  joint_fit <- function(k) {
    sides <- list()
    if (!is.na(ciR[k])) sides$R <- side_system(k, "R", ciR[k])
    if (!is.na(ciL[k])) sides$L <- side_system(k, "L", ciL[k])
    if (!length(sides)) return(NULL)
    # joint system over the union window; the shared anchor sample couples the
    # two sides exactly through Omega (no independence approximation)
    gidx_all <- sort(unique(unlist(lapply(sides, `[[`, "gidx"))))
    nU  <- length(gidx_all)
    KDt <- sum(vapply(sides, function(s) nrow(s$Bb), integer(1)))
    V <- matrix(0, KDt, nU * D)
    off <- 0L
    for (s in sides) {
      loc <- match(s$gidx, gidx_all); win <- length(s$gidx)
      for (cc in seq_len(D))
        V[off + seq_len(nrow(s$X)), (cc - 1L) * nU + loc] <-
          s$X[, ((cc - 1L) * win + 1L):(cc * win)]
      off <- off + nrow(s$X)
    }
    s2    <- rep(sig_vec^2, each = nU)
    Omega <- V %*% (s2 * t(V))
    Wm <- tryCatch(solve(Omega + 1e-10 * mean(diag(Omega)) * diag(KDt)),
                   error = function(err) NULL)
    if (is.null(Wm)) return(NULL)
    Bb   <- do.call(rbind, lapply(sides, `[[`, "Bb"))
    BtW  <- crossprod(Bb, Wm)
    BtWB <- BtW %*% Bb

    r  <- unlist(lapply(sides, `[[`, "r"))
    em <- function(u0c, p_use = p)
      unlist(lapply(sides, function(s) as.vector(s$em(u0c, p_use))))
    proj <- function(rhs) as.numeric(solve(BtWB, BtW %*% rhs))

    # fixed point, ONE EM evaluation per iteration (cached for both the
    # projection and the best-iterate divergence guard)
    u0 <- proj(r)
    if (!all(is.finite(u0))) return(NULL)
    best <- u0; best_res <- Inf
    for (it in seq_len(max_iter)) {
      EMc <- em(u0)
      e   <- as.numeric(Bb %*% u0) - (r - EMc)
      rn  <- sqrt(sum(e * (Wm %*% e)))
      if (is.finite(rn) && rn < best_res) { best_res <- rn; best <- u0 }
      un <- proj(r - EMc)
      if (!all(is.finite(un))) break
      delta <- sqrt(sum((un - u0)^2))
      u0 <- un
      if (is.finite(delta) && delta < tol) break
    }
    EMc <- em(u0)
    e   <- as.numeric(Bb %*% u0) - (r - EMc)
    rn  <- sqrt(sum(e * (Wm %*% e)))
    if (!(is.finite(rn) && rn <= best_res)) u0 <- best

    res <- tryCatch({
      h   <- 1e-6 * max(1, sqrt(sum(u0^2)))
      EMp <- matrix(0, KDt, D)
      for (e_i in seq_len(D)) {
        up <- u0; up[e_i] <- up[e_i] + h
        dn <- u0; dn[e_i] <- dn[e_i] - h
        EMp[, e_i] <- (em(up) - em(dn)) / (2 * h)
      }
      P <- solve(BtWB + BtW %*% EMp, BtW)
      G <- P %*% V
      covn <- G %*% (s2 * t(G))
      se_noise <- sqrt(pmax(diag(covn), 0))
      if (want_param) {
        rhs_of_p <- function(p_use)
          unlist(lapply(sides, function(s) as.vector(s$r_of_p(p_use)))) -
            em(u0, p_use)
        S_p <- matrix(0, D, J)
        for (j in seq_len(J)) {
          hj <- 1e-6 * max(1, abs(p[j]))
          pj_up <- p; pj_up[j] <- pj_up[j] + hj
          pj_dn <- p; pj_dn[j] <- pj_dn[j] - hj
          S_p[, j] <- P %*% ((rhs_of_p(pj_up) - rhs_of_p(pj_dn)) / (2 * hj))
        }
        covn <- covn + S_p %*% param_cov %*% t(S_p)
      }
      # O(sigma^2) debias: both channels (f'' mean + weight feedback) on the
      # joint system, gated at 2x the noise-channel SE exactly as in
      # estimate_IC (a larger correction flags an untrustworthy expansion).
      bias_b <- if (!is.null(Juu_all)) tryCatch({
        b <- build_state_bias_o2(sides, gidx_all, V, s2, Wm, Bb + EMp, P,
                                 Juu_all, sig_vec, dt, D)
        if (all(is.finite(b)) && all(abs(b) <= 2 * se_noise)) b else NULL
      }, error = function(err) NULL) else NULL
      list(cov = covn, bias = bias_b)
    }, error = function(err) NULL)
    if (is.null(res)) return(NULL)
    if (!is.null(res$bias)) u0 <- u0 - res$bias
    list(u0 = u0, cov = res$cov, debiased = !is.null(res$bias))
  }

  U_star   <- U
  P_arr    <- array(0, c(M, D, D))
  se       <- matrix(NA_real_, M, D)
  fallback <- logical(M)
  debiased <- logical(M)
  for (k in seq_len(M)) {
    fit <- tryCatch(joint_fit(k), error = function(err) NULL)
    if (is.null(fit) || !all(is.finite(fit$u0))) {
      # raw-observation fallback: a failed anchor must not poison the
      # trajectory (same philosophy as estimate_IC's divergence guard)
      fallback[k]  <- TRUE
      P_arr[k, , ] <- diag(sig_vec^2, nrow = D, ncol = D)
      se[k, ]      <- sig_vec
      next
    }
    U_star[k, ]  <- fit$u0
    P_arr[k, , ] <- fit$cov
    se[k, ]      <- sqrt(pmax(diag(fit$cov), 0))
    debiased[k]  <- isTRUE(fit$debiased)
  }

  list(
    U_star     = U_star,
    P_smooth   = P_arr,
    se         = se,
    fallback   = fallback,
    debiased   = debiased,
    cov_method = if (want_param) "noise_propagation+param" else "noise_propagation",
    design     = data.frame(k = seq_len(M),
                            r_c_L = cands$r_c[ciL], n_bl_L = cands$n_bl[ciL],
                            r_c_R = cands$r_c[ciR], n_bl_R = cands$n_bl[ciR])
  )
}

# ---- Coupled weak-form GLS state estimation ----------------------------------

# Legendre P_0..P_3 as ascending power-series coefficients in xi on [-1, 1],
# plus tiny polynomial helpers (coefficient differentiation / Horner eval).
.coupled_leg <- list(c(1), c(0, 1), c(-0.5, 0, 1.5), c(0, -1.5, 0, 2.5))
.poly_deriv  <- function(cf) if (length(cf) <= 1L) 0 else cf[-1L] * seq_len(length(cf) - 1L)
.poly_eval   <- function(cf, x) { y <- 0; for (c0 in rev(cf)) y <- y * x + c0; y }

#' Coupled weak-form GLS state estimation (global two-point system)
#'
#' Estimates the full state trajectory by ONE global GLS solve in all
#' \eqn{M \times D} unknowns \eqn{z = (u(t_1), \dots, u(t_M))}, coupling the
#' unknowns through the dynamics with TWO-POINT weak equations: for a window
#' \eqn{[t_a, t_b]} and any smooth test function \eqn{\psi}, integration by
#' parts + trapezoid + two-endpoint Euler-Maclaurin gives the exact equation
#' \deqn{\psi(b) z_b - \psi(a) z_a = T - EM(z_a, z_b),}
#' where \eqn{T = \Delta t \sum_m w_m [\psi_m f(U_m) + \psi'_m U_m]} is a data
#' integral and the \eqn{O(h^2), O(h^4)} EM defect involves only the endpoint
#' unknowns (the one-sided \eqn{\psi(b) = 0} case reproduces
#' \code{\link{estimate_IC}}). Two row families: SHORT \eqn{[t_k, t_{k+1}]}
#' with \eqn{\psi = 1} (adjacent dynamics chain; EM-safe at any spacing since
#' \eqn{\psi' = 0}) and LONG \eqn{[t_k, t_k + L - 1]} every \code{s_long}
#' samples with Legendre \eqn{P_1..P_3} mapped to the window (level
#' information). \eqn{P_0} long rows are dropped when shorts are present: a
#' \eqn{P_0} row is EXACTLY the telescoped sum of the shorts it spans (same
#' data integral, same EM), and the duplicate makes \eqn{\Omega} singular.
#'
#' All rows share one global \eqn{\Omega = X \mathrm{diag}(\sigma^2) X^T}
#' plus a SECOND-ORDER variance floor: rows whose first-order sensitivity
#' vanishes (\eqn{\psi}-channel with \eqn{J_u \approx 0}) are not noiseless —
#' their error is the second-order \eqn{\frac12 \epsilon^T f'' \epsilon} term
#' the delta method drops, and without the floor GLS enforces them as exact
#' constraints and accumulates unmodeled \eqn{O(\sigma^2)} error along the
#' chain. The \eqn{O(\sigma^2)} MEAN bias (same two channels as
#' \code{\link{estimate_IC}}: \eqn{f''} noise nonlinearity + feasible-GLS
#' weight feedback, partially cancelling) is subtracted analytically when
#' \code{debias = TRUE}.
#'
#' Compared to the pointwise \code{\link{wendy_gls_state}}: dominates on
#' chaotic coarse grids (Lorenz n=128: RMSE ~2x lower than both pointwise and
#' ERTS, smoothest, best-calibrated) and approaches the exact-dynamics
#' \eqn{q=0} CRLB at the trajectory EDGES; RMSE parity on mild systems.
#' Dense implementation: the \eqn{\Omega} factorization costs
#' \eqn{O((K D)^3)} with \eqn{K \approx M(1 + 3/s_{long})} rows — fine to
#' \eqn{M} of a few hundred; the banded structure is not yet exploited.
#' Requires a uniform grid.
#'
#' @param U Numeric matrix (M x D) of noisy observations.
#' @param f_,J_u,dF_dt_,d2F_dt2_,d3F_dt3_ Callable RHS, state Jacobian and
#'   total time-derivative evaluators built from the symbolic engine.
#' @param tt Numeric vector (length M) of UNIFORM time points.
#' @param p Numeric parameter vector \eqn{\hat\theta} (held fixed).
#' @param sigma Data-noise SD, scalar or length-D; required.
#' @param param_cov Optional J x J covariance of \eqn{\hat p}; adds the
#'   parameter channel \eqn{S_p \hat C S_p^T} with
#'   \eqn{S_p = \partial \hat z / \partial p} via the IFT projection, exactly
#'   as in \code{\link{estimate_IC}}.
#' @param L_long,s_long,degs Long-row family: window length (samples), start
#'   stride, and Legendre degrees (degree 0 is dropped with a warning when
#'   \code{use_short}; see Details).
#' @param use_short Logical; include the adjacent \eqn{\psi = 1} chain.
#' @param max_iter,tol EM fixed-point controls (best-iterate guard as in
#'   \code{\link{estimate_IC}}).
#' @param debias Logical (default \code{TRUE}). Subtract the analytic
#'   \eqn{O(\sigma^2)} bias (both channels) from \eqn{\hat z}; gated at twice
#'   the noise-channel SE per entry.
#' @return Named list with \code{U_star} (M x D), \code{P_smooth}
#'   (M x D x D per-time marginal covariance blocks, noise channel plus
#'   parameter channel when supplied), \code{se} (M x D), \code{cov_z}
#'   (MD x MD joint covariance, state-fast layout \eqn{(k-1)D + d}),
#'   \code{cov_method}, \code{bias_o2} (M x D analytic bias estimate or
#'   \code{NULL}), \code{debiased} (logical; whether it was subtracted),
#'   \code{iters}, \code{converged}, and \code{families} (row-family counts).
#' @export
wendy_coupled_state <- function(U, f_, J_u, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p,
                                sigma,
                                param_cov = NULL,
                                L_long    = 21L,
                                s_long    = 2L,
                                degs      = 1:3,
                                use_short = TRUE,
                                max_iter  = 10L,
                                tol       = 1e-8,
                                debias    = TRUE) {
  M <- nrow(U); D <- ncol(U); J <- length(p)
  tt_vec <- as.vector(tt)
  dt <- mean(diff(tt_vec))
  c2 <- dt^2 / 12; c4 <- dt^4 / 720

  if (!(length(sigma) %in% c(1L, D)) || !all(is.finite(sigma)))
    stop("wendy_coupled_state requires a finite sigma (scalar or per-state).",
         call. = FALSE)
  sig_vec <- if (length(sigma) == 1L) rep(sigma, D) else as.numeric(sigma)
  L_long <- as.integer(L_long)
  if (M < L_long + 1L)
    stop(sprintf("wendy_coupled_state needs M >= L_long + 1 (= %d).",
                 L_long + 1L), call. = FALSE)
  degs <- sort(unique(as.integer(degs)))
  if (use_short && any(degs == 0L)) {
    warning("dropping Legendre degree 0: a P0 long row telescopes exactly ",
            "to the sum of the short rows it spans (singular Omega).",
            call. = FALSE)
    degs <- setdiff(degs, 0L)
  }
  if (any(degs < 0L | degs > 3L))
    stop("degs must be within 0:3.", call. = FALSE)

  # anchor-independent per-sample quantities
  input <- rbind(matrix(rep(p, M), nrow = J), t(U), matrix(tt_vec, nrow = 1L))
  F_all <- f_(input)
  Ju_all <- array(0, c(M, D, D))
  for (m in seq_len(M))
    Ju_all[m, , ] <- matrix(as.vector(J_u(c(p, U[m, ], tt_vec[m]))), D, D)

  # ---- row groups: each couples unknowns z_a, z_b (D scalar rows) ----------
  groups <- list()
  add_group <- function(a, b, psi_cf) {
    idx <- a:b; w <- length(idx)
    Wt  <- tt_vec[b] - tt_vec[a]
    xi  <- if (w > 1L) 2 * (tt_vec[idx] - tt_vec[a]) / Wt - 1 else 0
    cfs <- vector("list", 5L); cfs[[1L]] <- psi_cf
    for (j in 2:5) cfs[[j]] <- .poly_deriv(cfs[[j - 1L]])
    sc  <- (2 / Wt)^(0:4)
    phi_a <- vapply(1:5, function(j) sc[j] * .poly_eval(cfs[[j]], -1), numeric(1))
    phi_b <- vapply(1:5, function(j) sc[j] * .poly_eval(cfs[[j]],  1), numeric(1))
    psi_m  <- .poly_eval(cfs[[1L]], xi)
    dpsi_m <- sc[2] * .poly_eval(cfs[[2L]], xi)
    wts <- rep(1, w); wts[1L] <- wts[w] <- 0.5
    groups[[length(groups) + 1L]] <<- list(a = a, b = b, idx = idx,
                                           psi = psi_m * wts,
                                           dpsi = dpsi_m * wts,
                                           phi_a = phi_a, phi_b = phi_b)
  }
  if (use_short) for (k in seq_len(M - 1L)) add_group(k, k + 1L, c(1))
  n_short <- length(groups)
  starts <- unique(pmin(seq(1L, M, by = max(1L, as.integer(s_long))),
                        M - L_long + 1L))
  for (a in starts) for (g in degs + 1L) add_group(a, a + L_long - 1L,
                                                   .coupled_leg[[g]])

  nG <- length(groups); KD <- nG * D; MD <- M * D
  zi <- function(k) ((k - 1L) * D + 1L):(k * D)        # unknown/datum block

  # ---- assemble A (design), Tv (data integrals), X (noise sensitivity) -----
  A  <- matrix(0, KD, MD)
  Tv <- numeric(KD)
  X  <- matrix(0, KD, MD)
  for (g in seq_len(nG)) {
    gr <- groups[[g]]; rr <- zi(g)
    A[cbind(rr, zi(gr$b))] <- A[cbind(rr, zi(gr$b))] + gr$phi_b[1L]
    A[cbind(rr, zi(gr$a))] <- A[cbind(rr, zi(gr$a))] - gr$phi_a[1L]
    Uw <- U[gr$idx, , drop = FALSE]
    Tv[rr] <- dt * as.numeric(crossprod(F_all[gr$idx, , drop = FALSE], gr$psi) +
                              crossprod(Uw, gr$dpsi))
    for (i in seq_along(gr$idx)) {
      m <- gr$idx[i]
      blk <- dt * (gr$psi[i] * Ju_all[m, , ] + gr$dpsi[i] * diag(D))
      X[rr, zi(m)] <- X[rr, zi(m)] + blk
    }
  }

  s2 <- rep(sig_vec^2, times = M)                      # datum layout (m-1)D + c
  Omega <- X %*% (s2 * t(X))

  # SECOND-ORDER variance floor: rows whose first-order sensitivity vanishes
  # (psi-channel with J_u ~ 0) are NOT noiseless — their error is the
  # second-order (1/2) eps' f'' eps term the delta method drops. Without the
  # floor GLS enforces them as exact constraints and the chain accumulates
  # unmodeled O(sigma^2) error. Huu is reused below by the mean debias.
  Huu <- array(0, c(M, D, D, D))           # [m, d, c, c'] = d2 f_d / du_c du_c'
  for (m in seq_len(M)) {
    h0 <- 1e-5 * max(1, sqrt(sum(U[m, ]^2)))
    for (cp in seq_len(D)) {
      up <- U[m, ]; up[cp] <- up[cp] + h0
      dn <- U[m, ]; dn[cp] <- dn[cp] - h0
      Huu[m, , , cp] <- (matrix(as.vector(J_u(c(p, up, tt_vec[m]))), D, D) -
                         matrix(as.vector(J_u(c(p, dn, tt_vec[m]))), D, D)) /
                        (2 * h0)
    }
  }
  v2 <- numeric(KD)
  for (g in seq_len(nG)) {
    gr <- groups[[g]]
    for (d in seq_len(D)) {
      acc <- 0
      for (i in seq_along(gr$idx)) {
        m <- gr$idx[i]
        Hd <- matrix(Huu[m, d, , ], D, D)
        acc <- acc + (dt * gr$psi[i])^2 * 0.5 *
               sum((Hd^2) * outer(sig_vec^2, sig_vec^2))
      }
      v2[(g - 1L) * D + d] <- acc
    }
  }
  Omega <- Omega + diag(v2, nrow = KD)

  # eigen-floored inverse (robust to residual near-collinearity across rows)
  eg  <- eigen(Omega, symmetric = TRUE)
  lam <- pmax(eg$values, 1e-10 * max(eg$values))
  Wm  <- eg$vectors %*% (t(eg$vectors) / lam)
  AtW  <- crossprod(A, Wm)
  AtWA <- AtW %*% A
  proj <- function(rhs) as.numeric(solve(AtWA, AtW %*% rhs))

  # ---- EM via per-endpoint coefficient matrices C_k (5 x D), cached --------
  Cof <- function(zk, tk, p_use = p) {
    inp <- matrix(c(p_use, zk, tk), ncol = 1L)
    fd  <- list(as.vector(f_(inp)), as.vector(dF_dt_(inp)),
                as.vector(d2F_dt2_(inp)), as.vector(d3F_dt3_(inp)))
    Cm <- matrix(0, 5L, D)
    for (j in 1:5) {
      ph <- as.list(numeric(5)); ph[[j]] <- 1
      Cm[j, ] <- c2 * g_deriv_at_endpoint(ph, fd, zk, 1L) -
                 c4 * g_deriv_at_endpoint(ph, fd, zk, 3L)
    }
    Cm
  }
  ep_idx <- sort(unique(unlist(lapply(groups, function(g) c(g$a, g$b)))))
  em_fun <- function(z, p_use = p) {
    Zm <- matrix(z, M, D, byrow = TRUE)
    Ck <- vector("list", M)
    for (k in ep_idx) Ck[[k]] <- Cof(Zm[k, ], tt_vec[k], p_use)
    em <- numeric(KD)
    for (g in seq_len(nG)) {
      gr <- groups[[g]]
      em[zi(g)] <- as.numeric(gr$phi_b %*% Ck[[gr$b]] - gr$phi_a %*% Ck[[gr$a]])
    }
    em
  }

  # ---- fixed point with best-iterate guard ----------------------------------
  z <- as.vector(t(U))
  best <- z; best_res <- Inf
  iters <- 0L; converged <- FALSE
  for (it in seq_len(max_iter)) {
    iters <- it
    EMc <- em_fun(z)
    e   <- as.numeric(A %*% z) - (Tv - EMc)
    rn  <- sqrt(sum(e * (Wm %*% e)))
    if (is.finite(rn) && rn < best_res) { best_res <- rn; best <- z }
    zn <- proj(Tv - EMc)
    if (!all(is.finite(zn))) break
    delta <- sqrt(sum((zn - z)^2) / M)
    z <- zn
    if (is.finite(delta) && delta < tol) { converged <- TRUE; break }
  }
  EMc <- em_fun(z); e <- as.numeric(A %*% z) - (Tv - EMc)
  rn <- sqrt(sum(e * (Wm %*% e)))
  if (!(is.finite(rn) && rn <= best_res)) z <- best

  # ---- IFT covariance with EM feedback --------------------------------------
  Zm <- matrix(z, M, D, byrow = TRUE)
  dC <- vector("list", M)
  for (k in ep_idx) {
    h <- 1e-6 * max(1, sqrt(sum(Zm[k, ]^2)))
    dC[[k]] <- lapply(seq_len(D), function(cc) {
      up <- Zm[k, ]; up[cc] <- up[cc] + h
      dn <- Zm[k, ]; dn[cc] <- dn[cc] - h
      (Cof(up, tt_vec[k]) - Cof(dn, tt_vec[k])) / (2 * h)
    })
  }
  EMp <- matrix(0, KD, MD)
  for (g in seq_len(nG)) {
    gr <- groups[[g]]; rr <- zi(g)
    for (cc in seq_len(D)) {
      EMp[rr, (gr$b - 1L) * D + cc] <- EMp[rr, (gr$b - 1L) * D + cc] +
        as.numeric(gr$phi_b %*% dC[[gr$b]][[cc]])
      EMp[rr, (gr$a - 1L) * D + cc] <- EMp[rr, (gr$a - 1L) * D + cc] -
        as.numeric(gr$phi_a %*% dC[[gr$a]][[cc]])
    }
  }
  P <- solve(AtWA + AtW %*% EMp, AtW)
  G <- P %*% X
  covz <- G %*% (s2 * t(G))
  se_noise <- sqrt(pmax(diag(covz), 0))

  # Parameter channel (law of total variance, as in estimate_IC): S_p =
  # P d(Tv - EM)/dp by central differences; cross term with the noise channel
  # ignored as elsewhere.
  cov_param <- if (!is.null(param_cov)) tryCatch({
    rhs_of_p <- function(p_use) {
      inputp <- rbind(matrix(rep(p_use, M), nrow = J), t(U),
                      matrix(tt_vec, nrow = 1L))
      Fp <- f_(inputp)
      Tvp <- numeric(KD)
      for (g in seq_len(nG)) {
        gr <- groups[[g]]
        Tvp[zi(g)] <- dt * as.numeric(
          crossprod(Fp[gr$idx, , drop = FALSE], gr$psi) +
          crossprod(U[gr$idx, , drop = FALSE], gr$dpsi))
      }
      Tvp - em_fun(z, p_use)
    }
    S_p <- matrix(0, MD, J)
    for (j in seq_len(J)) {
      hj <- 1e-6 * max(1, abs(p[j]))
      pj_up <- p; pj_up[j] <- pj_up[j] + hj
      pj_dn <- p; pj_dn[j] <- pj_dn[j] - hj
      S_p[, j] <- P %*% ((rhs_of_p(pj_up) - rhs_of_p(pj_dn)) / (2 * hj))
    }
    S_p %*% param_cov %*% t(S_p)
  }, error = function(err) NULL) else NULL

  cov_total <- if (!is.null(cov_param)) covz + cov_param else covz

  # ---- O(sigma^2) debias: BOTH channels (see build_ic_bias_o2) --------------
  #   b1 = P E[dTv]  (f'' mean; E[dTv]_(g,d) = dt sum_i psi_i 0.5 sum_c
  #                   Huu[m_i,d,c,c] sig_c^2)
  #   b2 = -P c      (weight feedback; c = sum_n sig_n^2 (T_n S X' + X S T_n')
  #                   (W R X)[, n], T_n = dX/dU_n = dt psi_i Huu)
  # The v2 floor's own U-dependence is O(sigma^4) and ignored. Gated per
  # entry at 2x the noise-channel SE, as in the other estimators.
  bias_z   <- NULL
  debiased <- FALSE
  if (isTRUE(debias)) {
    bias_try <- tryCatch({
      touch_g <- vector("list", M)
      for (g in seq_len(nG)) {
        gr <- groups[[g]]
        for (i in seq_along(gr$idx)) {
          if (gr$psi[i] == 0) next
          m <- gr$idx[i]
          touch_g[[m]] <- c(touch_g[[m]], list(c(g, gr$psi[i])))
        }
      }
      Aw <- Wm %*% (X - (A + EMp) %*% G)               # W (I - (A+EMp)P) X
      Eq <- numeric(KD); cvec <- numeric(KD)
      for (m in seq_len(M)) {
        tg <- touch_g[[m]]
        if (!length(tg)) next
        for (e in tg) {                                # channel 1
          g <- e[1]; psi_i <- e[2]
          for (d in seq_len(D)) {
            tr_d <- 0
            for (cc in seq_len(D))
              tr_d <- tr_d + Huu[m, d, cc, cc] * sig_vec[cc]^2
            Eq[(g - 1L) * D + d] <- Eq[(g - 1L) * D + d] +
              dt * psi_i * 0.5 * tr_d
          }
        }
        for (cp in seq_len(D)) {                       # channel 2
          a_n <- Aw[, (m - 1L) * D + cp]
          w_c <- numeric(D)
          for (cc in seq_len(D)) {
            col <- (m - 1L) * D + cc
            w_c[cc] <- s2[col] * sum(X[, col] * a_n)
          }
          Tw <- numeric(KD); ycc <- numeric(D)
          for (e in tg) {
            g <- e[1]; psi_i <- e[2]
            rr0 <- (g - 1L) * D
            for (d in seq_len(D)) {
              coef <- 0
              for (cc in seq_len(D)) coef <- coef + w_c[cc] * Huu[m, d, cc, cp]
              Tw[rr0 + d] <- Tw[rr0 + d] + dt * psi_i * coef
              for (cc in seq_len(D))
                ycc[cc] <- ycc[cc] + dt * psi_i * Huu[m, d, cc, cp] * a_n[rr0 + d]
            }
          }
          XSy <- numeric(KD)
          for (cc in seq_len(D)) {
            col <- (m - 1L) * D + cc
            XSy <- XSy + X[, col] * (s2[col] * ycc[cc])
          }
          cvec <- cvec + sig_vec[cp]^2 * (Tw + XSy)
        }
      }
      as.numeric(P %*% (Eq - cvec))
    }, error = function(err) NULL)
    if (!is.null(bias_try) && all(is.finite(bias_try)) &&
        all(abs(bias_try) <= 2 * se_noise)) {
      z        <- z - bias_try
      bias_z   <- bias_try
      debiased <- TRUE
    } else if (!is.null(bias_try)) {
      bias_z <- bias_try                               # report, don't apply
    }
  }

  P_arr <- array(0, c(M, D, D))
  for (k in seq_len(M)) P_arr[k, , ] <- cov_total[zi(k), zi(k)]
  list(
    U_star     = matrix(z, M, D, byrow = TRUE),
    P_smooth   = P_arr,
    se         = matrix(sqrt(pmax(diag(cov_total), 0)), M, D, byrow = TRUE),
    cov_z      = cov_total,
    cov_method = if (!is.null(cov_param)) "noise_propagation+param"
                 else "noise_propagation",
    bias_o2    = if (!is.null(bias_z)) matrix(bias_z, M, D, byrow = TRUE)
                 else NULL,
    debiased   = debiased,
    iters      = iters,
    converged  = converged,
    families   = c(short = n_short, long = nG - n_short)
  )
}

#' Estimate the state using the RTS smoother
#'
#' Extended Kalman filter forward pass followed by a Rauch-Tung-Striebel
#' backward smoother on the parameter-conditioned dynamics.
#'
#' The reported posterior covariance follows the law of total variance,
#' \deqn{\mathrm{Cov}(u_k^\star) = \mathrm{Cov}(u_k^\star \mid \hat p)
#'        + S_k \, \hat C \, S_k^T,}
#' splitting trajectory uncertainty into the conditional smoother posterior
#' (data noise + model/discretization error) and the contribution of parameter
#' uncertainty. The conditional pass uses a small model-error process noise
#' \eqn{Q = (0.1\,\sigma)^2 I_D}; parameter uncertainty is folded in separately
#' and \emph{coherently} via the sensitivity
#' \eqn{S_k = \partial u_k^\star / \partial \hat p}, obtained by central
#' differences of the smoothed mean over re-runs at \eqn{\hat p \pm h e_j}, with
#' \eqn{\hat C} the WENDy parameter covariance. This replaces the earlier
#' heuristic of injecting \eqn{(\Delta t_k)^2 \nabla_p f \,\hat C\, \nabla_p f^T}
#' as predict-step process noise, which mis-modelled \eqn{\delta p} as an
#' independent per-step draw (a random walk) rather than a single
#' fixed-but-uncertain value affecting the whole trajectory coherently. Folding
#' requires \code{param_cov} and \code{J_p}; when either is absent (or
#' \code{fold_param_uncertainty = FALSE}) the conditional posterior is returned
#' alone.
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
#' @param fold_param_uncertainty Logical (default \code{TRUE}). When \code{TRUE}
#'   and \code{param_cov}/\code{J_p} are supplied, the parameter-sensitivity term
#'   \eqn{S_k \hat C S_k^T} is added to \code{P_smooth}; otherwise only the
#'   conditional posterior is returned.
#' @return Named list with \code{U_star} (smoothed state), \code{P_smooth}
#'   (total posterior covariance), \code{P_smooth_cond} (the conditional-on-p̂
#'   posterior), \code{P_smooth_param} (the parameter-uncertainty contribution,
#'   or \code{NULL} when not folded), and the intermediate filter/predictor
#'   states.
#' @export
wendy_erts <- function(U, f_, J_u, tt, p, test_function_params,
                       sigma = NULL,
                       u0_init = NULL,
                       P0_init = NULL,
                       param_cov = NULL,
                       J_p = NULL,
                       fold_param_uncertainty = TRUE) {
  tt   <- as.vector(tt)
  mp1  <- nrow(U)
  D    <- ncol(U)
  J    <- length(p)
  I_D  <- diag(D)

  noise_sd <- as.numeric(if (!is.null(sigma)) sigma else estimate_std(U, k = 6))
  noise_sd <- mean(noise_sd)
  R_obs    <- noise_sd^2 * I_D

  # Model/discretization process noise for the conditional pass, modelled as a
  # continuous-time diffusion: the per-step covariance is q_c * dt_k (a rate,
  # state^2 per unit time), so the total injected over [0,T] is q_c*T and is
  # INVARIANT to the sampling density. A constant per-step Q instead grows the
  # total budget with the number of steps, letting a finer grid chase the noise
  # (high-frequency wiggle that does not shrink as data is added). Parameter
  # uncertainty is folded in coherently afterwards (law of total variance), not
  # injected here as per-step process noise.
  q_c <- (noise_sd * 0.1)^2

  u0 <- if (!is.null(u0_init)) as.vector(u0_init) else as.vector(U[1, ])
  P0 <- if (!is.null(P0_init)) P0_init
        else if (!is.null(u0_init)) 0.1 * noise_sd^2 * I_D
        else noise_sd^2 * I_D

  # One EKF forward pass + RTS backward pass at a fixed parameter vector. The
  # smoothed mean (U_star) is always returned; the posterior covariance arrays
  # are built only when want_cov = TRUE (skipped for the sensitivity re-runs).
  erts_pass <- function(p_use, want_cov) {
    S0          <- P0 + R_obs
    K0          <- P0 %*% solve(S0)
    u_filt      <- matrix(0, mp1, D)
    u_filt[1, ] <- u0 + K0 %*% (U[1, ] - u0)
    P_filt      <- array(0, c(mp1, D, D))
    P_filt[1,,] <- (I_D - K0) %*% P0 %*% t(I_D - K0) + K0 %*% R_obs %*% t(K0)

    u_pred  <- matrix(0, mp1, D)
    P_pred  <- array(0, c(mp1, D, D))
    F_store <- array(0, c(mp1 - 1L, D, D))

    for (k in seq_len(mp1 - 1L)) {
      dt_k <- tt[k + 1L] - tt[k]
      uk   <- u_filt[k, ]

      k1 <- as.vector(f_(matrix(c(p_use, uk,                tt[k]          ), ncol = 1)))
      k2 <- as.vector(f_(matrix(c(p_use, uk + 0.5*dt_k*k1, tt[k]+0.5*dt_k ), ncol = 1)))
      k3 <- as.vector(f_(matrix(c(p_use, uk + 0.5*dt_k*k2, tt[k]+0.5*dt_k ), ncol = 1)))
      k4 <- as.vector(f_(matrix(c(p_use, uk +     dt_k*k3, tt[k]+    dt_k  ), ncol = 1)))
      u_pred[k + 1L, ] <- uk + (dt_k / 6) * (k1 + 2*k2 + 2*k3 + k4)

      Ju_k <- matrix(as.vector(J_u(c(p_use, uk, tt[k]))), D, D)  # J[a, b] = df_a/du_b
      Fk   <- I_D + dt_k * Ju_k
      F_store[k,,] <- Fk

      Pk_pred <- Fk %*% P_filt[k,,] %*% t(Fk) + (q_c * dt_k) * I_D
      P_pred[k + 1L,,] <- Pk_pred

      Sk               <- Pk_pred + R_obs
      Kk               <- Pk_pred %*% solve(Sk)
      innov            <- U[k + 1L, ] - u_pred[k + 1L, ]
      u_filt[k + 1L, ] <- u_pred[k + 1L, ] + Kk %*% innov
      P_filt[k + 1L,,] <- (I_D - Kk) %*% Pk_pred %*% t(I_D - Kk) + Kk %*% R_obs %*% t(Kk)
    }

    u_smooth        <- matrix(0, mp1, D)
    u_smooth[mp1, ] <- u_filt[mp1, ]
    P_smooth <- if (want_cov) array(0, c(mp1, D, D)) else NULL
    if (want_cov) P_smooth[mp1,,] <- P_filt[mp1,,]

    for (k in seq(mp1 - 1L, 1L)) {
      Pk <- P_filt[k,,]
      Pp <- P_pred[k + 1L,,]
      Fk <- F_store[k,,]

      Gk <- Pk %*% t(Fk) %*% solve(Pp + 1e-10 * I_D)
      u_smooth[k, ] <- u_filt[k, ] + Gk %*% (u_smooth[k + 1L, ] - u_pred[k + 1L, ])
      if (want_cov)
        P_smooth[k,,] <- Pk + Gk %*% (P_smooth[k + 1L,,] - Pp) %*% t(Gk)
    }

    list(U_star = u_smooth, P_smooth = P_smooth,
         u_filt = u_filt, P_filt = P_filt, u_pred = u_pred, P_pred = P_pred)
  }

  # Conditional-on-p̂ pass: smoothed mean + conditional posterior covariance.
  base     <- erts_pass(p, want_cov = TRUE)
  u_smooth <- base$U_star
  P_cond   <- base$P_smooth

  # Fold parameter uncertainty in coherently via the trajectory sensitivity
  # S_k = ∂u*_k/∂p̂ (central differences over re-runs of the full smoother),
  # adding S_k Ĉ S_k^T to the conditional posterior at every time step.
  fold <- isTRUE(fold_param_uncertainty) && !is.null(param_cov) && !is.null(J_p)
  P_param  <- NULL
  P_smooth <- P_cond
  if (fold) {
    sens <- array(0, c(mp1, D, J))  # ∂u*_k/∂p̂_j
    for (j in seq_len(J)) {
      hj <- 1e-5 * max(1, abs(p[j]))
      pp <- p; pp[j] <- pp[j] + hj
      pm <- p; pm[j] <- pm[j] - hj
      sens[, , j] <- (erts_pass(pp, FALSE)$U_star -
                      erts_pass(pm, FALSE)$U_star) / (2 * hj)
    }
    P_param <- array(0, c(mp1, D, D))
    for (k in seq_len(mp1)) {
      Sk            <- matrix(sens[k, , ], D, J)
      P_param[k,,]  <- Sk %*% param_cov %*% t(Sk)
      P_smooth[k,,] <- P_cond[k,,] + P_param[k,,]
    }
  }

  list(
    U_star         = u_smooth,
    P_smooth       = P_smooth,
    P_smooth_cond  = P_cond,
    P_smooth_param = P_param,
    u_filt         = base$u_filt,
    P_filt         = base$P_filt,
    u_pred         = base$u_pred,
    P_pred         = base$P_pred
  )
}