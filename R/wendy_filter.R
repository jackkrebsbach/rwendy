# Leibniz expansion of g^(n)(t*) where g(t) = phi(t) F(p,u(t),t) + phi'(t) u(t)
# and trajectory derivatives F^(m) are passed in (precomputed via dF_dt_ etc.).
#
# phi_scalars: list of length (order+2) with phi^(0..order+1) at the endpoint.
# f_derivs:    list of length (order+1) with F^(0..order) at the endpoint.
# u_vec:       state at the endpoint (numeric vector of length D).
# order:       derivative order n.
# Binomial coefficients choose(n, 0:n) are constant for the (small, repeated)
# orders this is called with (1 and 3); cache them so the hot loop avoids both
# the per-iteration choose() calls and the seq() dispatch (seq.default was ~14%
# of the IC design-sweep self-time). Identical arithmetic to choose(n, k).
.gderiv_choose <- new.env(parent = emptyenv())
g_choose <- function(n) {
  key <- as.character(n)
  v   <- .gderiv_choose[[key]]
  if (is.null(v)) { v <- choose(n, 0:n); .gderiv_choose[[key]] <- v }
  v
}

g_deriv_at_endpoint <- function(phi_scalars, f_derivs, u_vec, order) {
  n      <- order
  cf     <- g_choose(n)                     # choose(n, 0:n)
  result <- phi_scalars[[n + 2L]] * u_vec   # φ^(n+1)·u
  for (k in 0:n) {
    result <- result + cf[k + 1L] * phi_scalars[[k + 1L]] * f_derivs[[n - k + 1L]]
  }
  if (n >= 1L) {
    for (k in 0:(n - 1L)) {
      result <- result + cf[k + 1L] * phi_scalars[[k + 2L]] * f_derivs[[n - k]]
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

# Boundary-layer IC system helpers

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

# A-priori (r_c, n_bl) selection by the calibrated u0-MSE proxy. For each
# candidate design one EM(4) solve (which also returns the EM(2) fixed point on
# the same built system, via return_em2_u0) gives
#   crit = sum_d Var_d / sigma_d^2  [variance]  + sum_d ( (u0_EM2 - u0_EM4)_d + bias_o2_d )^2 / sigma_d^2  [bias^2]
# (sum over states relative to the data noise, scale-free). Var = diag(cov_u0)
# folds the EM Jacobian AND the parameter channel S_p Chat S_p^T (law of total
# variance); (u0_EM2 - u0_EM4) is the oracle-free Euler-Maclaurin order-
# difference estimate of the truncation defect; bias_o2 is the analytic
# O(sigma^2) statistical bias. Validated against the true Monte-Carlo u0-MSE
# (examples/validation/ic_mse_proxy.R, ic_calibrated_criterion.R,
# ic_realC_overshoot.R): the proxy tracks the oracle argmin across systems and
# noise 0.05-0.20, and -- fed the real Fisher Chat (~sigma^2) via param_cov --
# does NOT over-shoot to the corner on smooth systems the way a variance-only
# criterion does. It self-protects at tiny r_c (the EM Jacobian inflates Var and
# the EM-difference term penalizes truncation), so the grid may roam below the
# pipe radius. Each candidate is a single EM solve (the EM(2) fixed point used
# for the order-difference is recovered from the same build, not a second
# estimate_IC call; for em_order == 2 there is no lower order to difference); the
# selected design is solved once more by the caller. Replaces the earlier
# noise-only no-EM variance objective, whose monotonicity in n_bl/window pinned
# the corner.
select_ic_design <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u,
                             sigma, sig_vec, param_cov, em_order,
                             r_c_grid, n_bl_grid, rc_cap) {
  r_cs <- sort(unique(pmin(pmax(as.integer(r_c_grid), 2L), rc_cap)))
  nbls <- sort(unique(pmax(as.integer(n_bl_grid), 1L)))
  D    <- ncol(U)
  s2   <- sig_vec^2
  rows <- vector("list", length(r_cs) * length(nbls))
  i <- 0L
  for (r_c in r_cs) for (n_bl in nbls) {
    res <- tryCatch({
      # EM(em_order) solve: deployed covariance (noise + param), statbias. With
      # return_em2_u0 the same call also returns the EM(2) fixed point (u0hat_em2)
      # from the same built system, so the EM-order-difference truncation estimate
      # needs no second estimate_IC solve (bit-identical to the former e_lo).
      e_hi <- estimate_IC(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, r_c,
                          J_u = J_u, sigma = sigma, param_cov = param_cov,
                          n_bl = n_bl, combine = "gls", debias = TRUE,
                          em_order = em_order, return_em2_u0 = em_order > 2L,
                          lean = TRUE)
      if (is.null(e_hi$cov_u0) || isTRUE(e_hi$diverged))
        stop("no covariance", call. = FALSE)
      statb <- if (!is.null(e_hi$bias_o2)) e_hi$bias_o2 else rep(0, D)
      u0_hi <- e_hi$u0hat + (if (isTRUE(e_hi$debias_applied)) statb else 0)
      # Truncation defect via the EM-order difference (oracle-free). Needs a
      # lower EM order, so it vanishes for em_order == 2 (criterion reduces to
      # variance + statbias^2 there).
      emb <- if (em_order > 2L) {
        if (isTRUE(e_hi$em2_diverged)) stop("lower-order diverged", call. = FALSE)
        e_hi$u0hat_em2 - u0_hi
      } else rep(0, D)
      vobj <- sum(diag(e_hi$cov_u0) / s2)
      c(vobj + sum((emb + statb)^2 / s2), vobj)
    }, error = function(err) c(NA_real_, NA_real_))
    i <- i + 1L
    rows[[i]] <- data.frame(r_c = r_c, n_bl = n_bl,
                            obj = res[1], var_obj = res[2])
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
#' sweeping a small grid and minimizing a calibrated \eqn{u_0}-MSE proxy
#' \deqn{\textstyle\sum_d \mathrm{Var}_d/\sigma_d^2
#'       + \sum_d ((u_0^{EM2} - u_0^{EM4})_d + b_d)^2/\sigma_d^2,}
#' where \eqn{\mathrm{Var} = \mathrm{diag}(\mathrm{cov\_u0})} folds the EM
#' Jacobian and the parameter channel \eqn{S_p \hat C S_p^T},
#' \eqn{u_0^{EM2} - u_0^{EM4}} is the oracle-free Euler-Maclaurin order-
#' difference estimate of the truncation defect, and \eqn{b} is the analytic
#' \eqn{O(\sigma^2)} statistical bias. Each candidate costs a single EM solve
#' (the lower-order \eqn{u_0^{EM2}} fixed point is recovered from the same build
#' via \code{return_em2_u0}, not a second solve). This replaces the earlier
#' noise-only variance objective, which was
#' monotone in window/count and so pinned the corner; the MSE proxy turns up
#' where the true MSE turns up (validated examples/validation/ic_mse_proxy.R,
#' ic_calibrated_criterion.R, ic_realC_overshoot.R). \code{combine = "ols"}
#' restores the legacy unweighted combine and its \code{max(3, ceiling(r_c/8))}
#' count heuristic.
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
#'   \code{n_bl} are selected by the a-priori MSE-proxy design sweep (see
#'   Details); when \code{NULL} under \code{combine = "ols"} the legacy heuristic
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
#' @param return_em2_u0 Logical (default \code{FALSE}, internal). When
#'   \code{TRUE} and \code{em_order == 4}, additionally solve the EM(2) fixed
#'   point on the same built system and return it as \code{u0hat_em2} (with
#'   \code{em2_diverged}). Used by the design sweep to form the EM-order-
#'   difference truncation estimate without a second \code{estimate_IC} call.
#' @param lean Logical (default \code{FALSE}, internal). When \code{TRUE}, skip
#'   the LS-residual fallback covariance (\code{cov_u0_resid}) whenever the
#'   noise-propagation channel is available -- it is only ever used as a
#'   fallback and the design sweep never reads it. The public default computes
#'   it for reference as before.
#' @param r_c_grid,n_bl_grid Optional integer vectors of design candidates for
#'   the a-priori sweep (defaults \code{c(1, 2, 3, 4) * r_c} capped at
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
#'   columns \code{r_c}, \code{n_bl}, \code{obj} (the MSE proxy), and
#'   \code{var_obj} (its variance part), or \code{NULL} when no sweep ran),
#'   \code{bias_o2} (the analytic
#'   \eqn{O(\sigma^2)} bias estimate of the UNcorrected solve, length-D, or
#'   \code{NULL} when not computed), \code{debias_applied} (whether
#'   \code{bias_o2} was subtracted from \code{u0hat}), \code{u0hat_em2} (the
#'   EM(2) fixed point, or \code{NULL} unless \code{return_em2_u0} and
#'   \code{em_order == 4}), \code{em2_diverged}, \code{iters},
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
                        debias         = TRUE,
                        return_em2_u0  = FALSE,
                        lean           = FALSE) {
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
  # sweep (r_c, n_bl) candidates and keep the design minimizing the calibrated
  # u0-MSE proxy (select_ic_design). The passed r_c (pipeline radius) seeds the
  # grid. Under the OLS combine the legacy heuristic below applies.
  design_table <- NULL
  if (combine == "gls" && is.null(n_bl)) {
    if (is.null(r_c_grid))  r_c_grid  <- c(2L, 5L, 6L, 8L, 10L, 12L, 16L, r_c)
    if (is.null(n_bl_grid)) n_bl_grid <- c(8L)
    sel <- tryCatch(
      select_ic_design(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u,
                       sigma, sig_vec, param_cov, em_order,
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

  # c4_use lets a single built system evaluate the EM defect at either order:
  # c4_use = c4 (default) is the deployed em_order correction; c4_use = 0 drops
  # the h^4 term (and skips the order-3 work), reproducing the EM(2) defect. The
  # design sweep uses this to get the EM(2) fixed point from the EM(4) build
  # without a second estimate_IC call (see return_em2_u0 below).
  em_correction <- function(u0_curr, p_use = p, c4_use = c4) {
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
      g3 <- if (c4_use != 0)
              g_deriv_at_endpoint(phi_t1_k, fd_t1, u_t1, order = 3L)
            else rep(0, D)
      EM[k, ] <- c2 * g1 - c4_use * g3
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

  # Residual over the K_bl x D system, in the same metric the projection
  # minimizes: W-weighted for GLS, Frobenius for OLS.
  residual_norm <- function(u0_curr, c4_use = c4) {
    EM <- em_correction(u0_curr, c4_use = c4_use)
    e  <- outer(B, as.numeric(u0_curr)) - (r_trap - EM)   # K_bl x D
    if (!is.null(gls)) {
      ev <- as.vector(e)
      sqrt(sum(ev * (gls$Wm %*% ev)))
    } else {
      sqrt(sum(e * e))
    }
  }

  # Fixed-point iteration u0 <- proj(r_trap - EM(u0)), tracking the best iterate
  # by residual norm. The iteration contracts only when the EM-correction map
  # has spectral radius < 1 (kappa = O(h^2)); for under-sampled stiff systems
  # the contraction can fail and iterates blow up, so returning the best-seen
  # iterate keeps us no worse than the uncorrected LS.
  run_fixed_point <- function(c4_use = c4) {
    u0       <- proj(r_trap)
    u0_hist  <- list(u0)
    best_u0  <- u0
    best_res <- if (all(is.finite(u0))) residual_norm(u0, c4_use = c4_use) else Inf
    iters    <- 0L
    converged <- FALSE
    diverged  <- FALSE

    for (it in seq_len(max_iter)) {
      iters  <- it
      EM     <- em_correction(u0, c4_use = c4_use)
      rhs    <- r_trap - EM
      u0_new <- proj(rhs)
      u0_hist[[it + 1L]] <- u0_new

      if (!all(is.finite(u0_new))) { diverged <- TRUE; break }

      res_new <- residual_norm(u0_new, c4_use = c4_use)
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

    list(u0 = best_u0, iters = iters, converged = converged,
         diverged = diverged, u0_hist = u0_hist)
  }

  fit       <- run_fixed_point()
  u0        <- fit$u0
  iters     <- fit$iters
  converged <- fit$converged
  diverged  <- fit$diverged
  u0_hist   <- fit$u0_hist

  # Optional EM(2) fixed point on the SAME built system (design sweep only). The
  # EM-order-difference truncation estimate u0_EM2 - u0_EM4 needs the lower-order
  # solution; computing it here (c4_use = 0) reuses this build, GLS weights, and
  # r_trap instead of a second estimate_IC call that would rebuild all of them
  # and redundantly compute a full covariance the sweep never reads. Bit-identical
  # to a standalone em_order = 2, debias = FALSE solve.
  u0hat_em2    <- NULL
  em2_diverged <- FALSE
  if (isTRUE(return_em2_u0) && c4 != 0) {
    fit2         <- run_fixed_point(c4_use = 0)
    u0hat_em2    <- fit2$u0
    em2_diverged <- fit2$diverged
  }

  # Covariance of u0hat. Combines up to three pieces and also returns the
  # implicit-function-theorem sensitivity (P, EMp) that the O(sigma^2) debias
  # reuses:
  #   cov_u0_resid - residual-based closed-form LS variance (fallback)
  #   cov_u0_noise - delta-method propagation of data noise (preferred)
  #   cov_u0_param - law-of-total-variance contribution from Cov(phat)
  compute_covariance <- function(u0) {
    # Residual-based closed-form LS variance (fallback): for the K_bl x 1
    # per-state regression B u0_d = rhs_d, rhs_d = (r_trap - EM(u0))[, d], the
    # per-state residual variance is s_d^2 = ||B u0_d - rhs_d||^2 / (K_bl - 1)
    # and Var(u0_d) = s_d^2 / B^T B. This OVER-states Var(u0hat): the residual
    # norm absorbs the DETERMINISTIC Euler-Maclaurin / trapezoidal truncation
    # mismatch on top of the propagated noise (empirically ~1.8x too wide in
    # SE, ~100% coverage of a nominal-95% interval on the logistic problem).
    # It is only ever used as the fallback when the noise channel is
    # unavailable, so its (em_correction-costing) computation is deferred until
    # after the noise channel below and skipped entirely in `lean` mode (the
    # design sweep) when the noise channel succeeded.
    compute_resid <- function() tryCatch({
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
    # EMp is taken by central differences. For OLS this reduces to the
    # (B^T B I + A)^{-1} form exactly: C Bbold = BtB I_D and C EMp = A.
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
    # win_cols). Calibrated ~95% coverage; for GLS the (near-)minimum-variance
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

    # Compute the residual fallback unless we're in lean mode AND the noise
    # channel already succeeded (lean callers never read cov_u0_resid).
    cov_u0_resid <- if (!lean || is.null(cov_u0_noise)) compute_resid() else NULL

    cov_u0     <- if (!is.null(cov_u0_noise)) cov_u0_noise else cov_u0_resid
    cov_method <- if (!is.null(cov_u0_noise)) "noise_propagation" else "ls_residual"
    if (!is.null(cov_u0) && !is.null(cov_u0_param)) {
      cov_u0     <- cov_u0 + cov_u0_param
      cov_method <- paste0(cov_method, "+param")
    }

    list(cov_u0 = cov_u0, cov_method = cov_method, cov_u0_resid = cov_u0_resid,
         cov_u0_param = cov_u0_param, cov_u0_noise = cov_u0_noise, P = P, EMp = EMp)
  }

  cov_res      <- compute_covariance(u0)
  cov_u0       <- cov_res$cov_u0
  cov_method   <- cov_res$cov_method
  cov_u0_resid <- cov_res$cov_u0_resid
  cov_u0_param <- cov_res$cov_u0_param
  cov_u0_noise <- cov_res$cov_u0_noise
  P            <- cov_res$P
  EMp          <- cov_res$EMp

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

  # O(sigma^2) debias (GLS path)
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
    u0hat_em2      = u0hat_em2,
    em2_diverged   = em2_diverged,
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
#' uncertainty. The conditional pass injects the RK4 \emph{discretization} error
#' as predict-step process noise, estimated per step by step-doubling
#' (Richardson): with \eqn{y_1} the single full RK4 step and \eqn{y_2} two half
#' steps, the local truncation error of \eqn{y_1} is
#' \eqn{e_k = (16/15)(y_1 - y_2)} and \eqn{Q_k = \mathrm{diag}(e_k^2)}. This is
#' \eqn{\sigma}-independent (discretization error depends on \eqn{\Delta t} and
#' the dynamics, not the measurement noise), vanishes as \eqn{O(T\,\Delta t^4)}
#' under grid refinement, and is active only when the grid under-resolves
#' \eqn{f}; it replaces the earlier ad-hoc \eqn{Q = (0.1\,\sigma)^2 I_D}.
#' Parameter uncertainty is folded in separately
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

  # Conditional-pass process noise = the RK4 DISCRETIZATION error of the predict
  # step, estimated per step by step-doubling (Richardson): with y1 the single
  # full RK4 step and y2 two half steps, the local truncation error of y1 is
  # e_k = (16/15)(y1 - y2) (order p = 4, factor 2^p/(2^p - 1)), so Q_k =
  # diag(e_k^2). Unlike the earlier (0.1 sigma)^2 heuristic this is
  # sigma-INDEPENDENT -- discretization error depends on dt and the dynamics,
  # not the measurement noise -- vanishes as O(T dt^4) under grid refinement (so
  # a finer grid cannot chase the noise), and is active only when the grid
  # under-resolves f. A tiny floor keeps Q_k strictly PD. Parameter uncertainty
  # is folded in coherently afterwards (law of total variance), not here.
  q_floor <- 1e-12

  # One RK4 step of f(p_use, ., .) over h; shared by the predict mean and the
  # step-doubling local-error estimate (called three times per step: one full,
  # two half).
  rk4_step <- function(p_use, uk, tk, h) {
    k1 <- as.vector(f_(matrix(c(p_use, uk,            tk        ), ncol = 1)))
    k2 <- as.vector(f_(matrix(c(p_use, uk + .5*h*k1, tk + .5*h ), ncol = 1)))
    k3 <- as.vector(f_(matrix(c(p_use, uk + .5*h*k2, tk + .5*h ), ncol = 1)))
    k4 <- as.vector(f_(matrix(c(p_use, uk +    h*k3, tk +    h ), ncol = 1)))
    uk + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)
  }

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

      y1 <- rk4_step(p_use, uk, tt[k], dt_k)                  # full step = predict mean
      yh <- rk4_step(p_use, uk, tt[k], dt_k / 2)              # step-doubling: two half
      y2 <- rk4_step(p_use, yh, tt[k] + dt_k / 2, dt_k / 2)   #   steps for the error est.
      u_pred[k + 1L, ] <- y1
      e_k <- (16 / 15) * (y1 - y2)                            # RK4 local truncation error

      Ju_k <- matrix(as.vector(J_u(c(p_use, uk, tt[k]))), D, D)  # J[a, b] = df_a/du_b
      Fk   <- I_D + dt_k * Ju_k
      F_store[k,,] <- Fk

      Pk_pred <- Fk %*% P_filt[k,,] %*% t(Fk) + diag(e_k^2 + q_floor, D)
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