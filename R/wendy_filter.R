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

# Coefficient matrix A ((order+2) x D) with g^(order) = sum_m phi^(m) A[m+1, ].
g_coeffs <- function(f_derivs, u_vec, order) {
  n  <- order
  cf <- g_choose(n)
  A  <- matrix(0, n + 2L, length(u_vec))
  A[n + 2L, ] <- A[n + 2L, ] + u_vec
  for (k in 0:n)
    A[k + 1L, ] <- A[k + 1L, ] + cf[k + 1L] * f_derivs[[n - k + 1L]]
  if (n >= 1L)
    for (k in 0:(n - 1L))
      A[k + 2L, ] <- A[k + 2L, ] + cf[k + 1L] * f_derivs[[n - k]]
  A
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

    em_coeffs <- function(f_derivs, u_pt, s) {
      A <- matrix(0, 5L, D)
      A[1:3, ] <- s * -c2 * g_coeffs(f_derivs, u_pt, 1L)
      A + s * c4 * g_coeffs(f_derivs, u_pt, 3L)
    }

    bl_phi_tM %*% em_coeffs(f_derivs_tM, u_tM, 1) +
      bl_phi_t1 %*% em_coeffs(f_derivs_t1, u_t1, -1)
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

    # 5 x (D*J), column-major with d fastest to match jac[k, , j].
    jac_coeffs <- function(Ax, s) {
      out <- matrix(0, 5L, D * J)
      for (j in seq_len(J)) {
        fd <- list(Ax[[1]][, j], Ax[[2]][, j], Ax[[3]][, j], Ax[[4]][, j])
        A  <- matrix(0, 5L, D)
        A[1:3, ] <- s * -c2 * g_coeffs(fd, zeroD, 1L)
        out[, ((j - 1L) * D + 1L):(j * D)] <- A + s * c4 * g_coeffs(fd, zeroD, 3L)
      }
      out
    }

    array(bl_phi_tM %*% jac_coeffs(AM, 1) + bl_phi_t1 %*% jac_coeffs(A1, -1),
          c(K_bl, D, J))
  }
}

# Boundary-layer IC system helpers

# Build the left boundary-layer system for one (r_bl, n_bl) design: trap-weighted
# order-0/1 rows, the boundary vector B = psi_k(t_1), the raw endpoint
# derivatives phi^(orders) at t_1, and the data window the rows touch.
# orders = 0:1 suffices for the design-stage (no-EM) covariance; the estimator
# itself needs 0:4 for the Euler-Maclaurin correction.
#
# include_interior = TRUE additionally stacks the DENSE interior test-function
# block (every admissible center at the same radius, build_test_function_matrix)
# under the BL rows. Interior supports never touch column 1, so phi^(0..4)(t_1)
# is EXACTLY zero for those rows: they load neither the boundary vector B nor
# the Euler-Maclaurin defect. They act purely as GLS control variates -- their
# weak residuals are mean-zero but built from the SAME noisy samples as the BL
# rows, so the BL<->interior cross-covariance block of Omega lets the GLS
# combine subtract the estimable part of the boundary noise (Schur complement
# S_bb - S_bi S_ii^-1 S_ib; zeroing S_bi reverts the variance to BL-only
# exactly). Under OLS the B-only projection ignores them, hence GLS-only.
# Validated in examples/validation/ic_interior_{sanity,mechanism,why}.R and
# ic_bl_split_knob*.R: u0 MSE 7-22x at known p, 1.5-2.2x under phat, coverage
# held with the parameter channel folded.
#
# interior_stride keeps only every s-th admissible center; the default 1 keeps
# them all. The block must keep spanning the WHOLE trajectory -- interior rows
# act through the interior<->interior chain, not as a local control variate, and
# 94-95% of the GLS row weight sits on rows whose support does not touch the BL
# window at all, so truncating the block to a span near t_1 costs 7.5-9.6x in
# a-priori SE (examples/validation/ic_audit_structure.R).
#
# At one center per sample the block can be over-complete enough to let the GLS
# claim cancellation that is not there, and thinning is the lever for that. The
# effect is confined to SMALL radii: over the a-priori grid at p_hat (60 reps
# per cell, examples/validation/ic_phat_design.R) stride 1 vs 2 differ only at
# r_bl = 8 on logistic at 20% noise (u0 NRMSE 0.0263 vs 0.0118), and are
# indistinguishable at every radius the design sweep actually selects (16-20
# there), across logistic / Lotka-Volterra / Lorenz at 5% and 20% noise: picked
# NRMSE ratio 0.973-1.10, coverage 0.95-1.00 either way. So stride stays an
# option rather than a default, and build_ic_noise_quad -- which addresses the
# same over-completeness from the covariance side -- carries it.
#
# The returned K_bl counts ALL stacked rows (downstream helpers use it as the
# equation count); the last K_int of them are interior. em_rows indexes the
# rows with nonzero phi at t_1 (the true BL rows) so the EM loops can skip the
# exactly-zero interior rows.
build_ic_bl_system <- function(tt_vec, r_bl, n_bl, orders = 0:4,
                               include_interior = FALSE,
                               interior_stride = 1L) {
  M <- length(tt_vec)
  blocks <- lapply(orders, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_bl, order = ord,
                               side = "left", n_bl = n_bl))
  K_int <- 0L

  # Interior block: only orders 0/1 carry data (V / Vp rows); orders >= 2 are
  # consumed solely through their t_1 column, which is exactly zero on interior
  # support, so zero rows are exact (and cheap).
  if (include_interior && (2L * r_bl + 1L) <= (M - 2L)) {
    V_int <- lapply(0:1, function(ord)
      build_test_function_matrix(psi, tt_vec, r_bl, order = ord))
    stride <- max(1L, as.integer(interior_stride))
    if (stride > 1L) {
      keep  <- seq(1L, nrow(V_int[[1]]), by = stride)
      V_int <- lapply(V_int, function(V) V[keep, , drop = FALSE])
    }
    K_int <- nrow(V_int[[1]])
    blocks <- lapply(seq_along(orders), function(i) {
      ord <- orders[i]
      rbind(blocks[[i]],
            if (ord <= 1L) V_int[[ord + 1L]] else matrix(0, K_int, M))
    })
  }
  K <- nrow(blocks[[1]])

  apply_trap <- function(V) {
    V[, 1] <- V[, 1] * 0.5
    V[, M] <- V[, M] * 0.5
    V
  }

  bl_phi_t1 <- matrix(0, nrow = K, ncol = length(orders))
  for (i in seq_along(orders)) bl_phi_t1[, i] <- blocks[[i]][, 1]

  win_cols <- which(colSums(abs(blocks[[1]]) + abs(blocks[[2]])) > 0)

  B <- bl_phi_t1[, 1]
  list(V_BL = apply_trap(blocks[[1]]), Vp_BL = apply_trap(blocks[[2]]),
       B = B, BtB = sum(B * B), K_bl = K, K_int = K_int,
       em_rows = which(rowSums(abs(bl_phi_t1)) > 0),
       bl_phi_t1 = bl_phi_t1, win_cols = win_cols)
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

  # J_u is vectorised: ONE (J + D + 1) x nw call returns nw x D^2, flattened
  # column-major, so Jflat[i, (cc - 1) * D + d] is Ju(win[i])[d, cc] -- bitwise
  # identical to the per-sample call. Evaluating one window sample at a time
  # cost 145-167x more here and was 30-38% of the whole design sweep
  # (examples/validation/, audit of 2026-08-14).
  input <- rbind(matrix(rep(p, nw), nrow = length(p)),
                 t(U[win, , drop = FALSE]),
                 matrix(tt_vec[win], nrow = 1L))
  Jflat <- matrix(as.vector(J_u(input)), nrow = nw, ncol = D * D)

  # Block (d, cc) of X is the windowed test-function matrix scaled column-wise
  # by Ju(.)[d, cc], plus Vp on the diagonal block: D^2 matrix operations
  # instead of nw * D^2 single-column writes.
  Vw  <- bl$V_BL[,  win, drop = FALSE]
  Vpw <- bl$Vp_BL[, win, drop = FALSE]

  X  <- matrix(0, KD, nw * D)
  s2 <- numeric(nw * D)
  for (cc in seq_len(D)) {
    cols     <- ((cc - 1L) * nw + 1L):(cc * nw)
    s2[cols] <- sig_vec[cc]^2
    for (d in seq_len(D)) {
      blk <- Vw * rep(Jflat[, (cc - 1L) * D + d], each = K_bl)
      if (d == cc) blk <- blk + Vpw
      X[rowblk(d), cols] <- -dt * blk
    }
  }
  list(X = X, s2 = s2, Bbold = Bbold)
}

# Covariance of the QUADRATIC-in-noise part of the weak residual: the O(sigma^4)
# block that the delta-method Omega = X diag(s2) X^T omits.
#
# With U = u + eta and eta_{m,c} ~ N(0, sigma_c^2) independent,
#   r_{k,d} = -dt sum_m [ V_BL[k,m] f_d(U_m) + Vp_BL[k,m] U_{m,d} ]
#           = r_{k,d}(u) + (X eta)_{k,d} + q_{k,d} + O(eta^3),
#   q_{k,d} = -dt sum_m V_BL[k,m] * 0.5 eta_m' J_uu(m)[d,,] eta_m
# (the Vp term is linear in U, so it contributes nothing here). Isserlis gives
# Cov(0.5 eta'A eta, 0.5 eta'B eta) = 0.5 tr(A S B S) with S = diag(sigma^2), so
#   Omega2[(k,d),(k',d')] = dt^2 sum_m V_BL[k,m] V_BL[k',m] G_m[d,d'],
#   G_m[d,d'] = 0.5 sum_{c,c'} J_uu(m)[d,c,c'] J_uu(m)[d',c,c'] sigma_c^2 sigma_c'^2,
# and there is NO cross-covariance with the linear part because E[eta eta eta]=0.
#
# Why it matters even though it is tiny: tr(Omega2)/tr(Omega1) is only 0.03-3.2%
# on the validated systems, but Omega1 is near-singular (smallest eigenvalues
# ~1e-19) and Omega2 is O(sigma^4) in exactly those directions -- which are the
# directions the GLS combine loads. Omitting it lets the combine claim unbounded
# cancellation where the first-order variance is spuriously zero: an SE that is
# too small AND real excess MSE. It is also what the ad-hoc 1e-10 ridge in
# build_ic_gls_weights was standing in for. MC-validated against the empirical
# covariance of vec(r) (examples/validation/ic_audit_omega2_check.R): relative
# trace error -0.0370 -> -0.0062 on Lorenz at 20% noise, -0.0044 -> -0.0020 on
# LV, shrinking with sigma as O(sigma^4)/O(sigma^2). Deployed effect at 20%
# noise (examples/validation/ic_audit_combined.R): u0 MSE x0.38 (Lorenz M=512,
# coverage 0.69 -> 0.99, z 2.45 -> 0.97), x0.64 (LV M=512), neutral on logistic.
#
# EXACTNESS CAVEAT: the only other O(sigma^4) contribution is Cov(linear, cubic),
# which needs the third derivative of f. That vanishes identically for the
# validated systems (logistic, Lotka-Volterra, Lorenz all have f''' == 0), so
# Omega1 + Omega2 is exact to O(sigma^6) there; for an f with nonzero f''' this
# is an improvement on Omega1 alone but not the complete O(sigma^4) covariance.
#
# Cost is a factor D BELOW forming Omega1: D(D+1)/2 blocks of K x nw times
# nw x K, against Omega1's KD x nwD times nwD x KD. J_uu comes from the shared
# build_ic_hessian_cache. Returns NULL when f is linear in u (J_uu == 0).
# The per-sample curvature array of the O(sigma^4) block,
#   G_m[d, d'] = 0.5 sum_{c,c'} J_uu(m)[d,c,c'] J_uu(m)[d',c,c'] sigma_c^2 sigma_c'^2,
# on an ARBITRARY set of grid columns. It depends on (U[m, ], p, sigma) alone, not
# on the boundary-layer design, so one evaluation on the union window serves every
# radius of a pool -- including the CROSS-radius blocks, which need the same array
# indexed against two different windows.
build_ic_quad_G <- function(U, tt_vec, p, J_u, sig_vec, cols, hess_cache = NULL) {
  D   <- ncol(U)
  nc  <- length(cols)
  juu <- if (!is.null(hess_cache)) hess_cache
         else build_ic_hessian_cache(U, tt_vec, p, J_u, D)
  ss  <- outer(sig_vec^2, sig_vec^2)                     # sigma_c^2 sigma_c'^2
  G   <- array(0, c(nc, D, D))
  for (i in seq_len(nc)) {
    H <- juu(cols[i])
    for (d in seq_len(D)) for (dp in d:D) {
      v <- 0.5 * sum(ss * (matrix(H[d, , ], D, D) * matrix(H[dp, , ], D, D)))
      G[i, d, dp] <- v
      G[i, dp, d] <- v
    }
  }
  G
}

build_ic_noise_quad <- function(bl, U, tt_vec, p, J_u, sig_vec, dt, hess_cache = NULL) {
  D   <- ncol(U)
  K   <- bl$K_bl
  win <- bl$win_cols
  nw  <- length(win)

  Vw  <- bl$V_BL[, win, drop = FALSE]                    # K x nw
  # Columns with V_BL == 0 enter only through Vp, i.e. linearly: no curvature.
  act <- which(colSums(abs(Vw)) > 0)

  G <- array(0, c(nw, D, D))
  if (length(act)) {
    Ga <- build_ic_quad_G(U, tt_vec, p, J_u, sig_vec, win[act], hess_cache)
    for (d in seq_len(D)) for (dp in seq_len(D)) G[act, d, dp] <- Ga[, d, dp]
  }
  if (max(abs(G)) <= 0) return(NULL)                     # linear f: Omega2 == 0

  O2 <- matrix(0, K * D, K * D)
  for (d in seq_len(D)) for (dp in d:D) {
    blk <- dt^2 * (Vw %*% (G[, d, dp] * t(Vw)))
    O2[((d - 1L) * K + 1L):(d * K), ((dp - 1L) * K + 1L):(dp * K)] <- blk
    if (dp != d)
      O2[((dp - 1L) * K + 1L):(dp * K), ((d - 1L) * K + 1L):(d * K)] <- t(blk)
  }
  O2
}

# GLS (BLUE) weights for the stacked BL system. Omega is the covariance of
# vec(r_trap): the errors of the K_bl equations are strongly correlated because
# they integrate the SAME noisy samples, which the unweighted (OLS) combine
# ignores. The linear (delta-method) part is X diag(s2) X^T; sens$Omega2, when
# present, adds the O(sigma^4) quadratic-noise block (see build_ic_noise_quad),
# without which the combine over-trusts the near-null directions of the linear
# part. W = Omega^{-1} (tiny ridge for near-duplicate rows).
# cov_design = (Bbold^T W Bbold)^{-1} is the no-EM a-priori GLS covariance of
# u0hat: it depends only on the design (r_bl, n_bl), sigma, and J_u along the
# observed trajectory -- not on the solved u0 -- so it doubles as the
# design-selection criterion.
build_ic_gls_weights <- function(sens) {
  KD    <- nrow(sens$X)
  Omega <- sens$X %*% (sens$s2 * t(sens$X))
  if (!is.null(sens$Omega2)) Omega <- Omega + sens$Omega2
  # Wm    <- solve(Omega + 1e-9 * mean(diag(Omega)) * diag(KD))
  # Ridge required: rank(Omega) <= nw*D < KD once r_bl exceeds ~12.
  Wm    <- solve(Omega + 1e-11 * mean(diag(Omega)) * diag(KD))
  BtW   <- crossprod(sens$Bbold, Wm)        # D x KD
  BtWB  <- BtW %*% sens$Bbold               # D x D
  list(Wm = Wm, BtW = BtW, BtWB = BtWB, cov_design = solve(BtWB))
}

# Memoizing state-Hessian cache, Juu[d, c, c'] = d2 f_d / du_c du_c' at an
# observed sample, by central differences of the state Jacobian J_u.
#
# The tensors depend on (U[m, ], p, t_m) alone -- not on the boundary-layer
# design -- so one cache serves every candidate of the a-priori sweep as well as
# the final solve.
#
# Filled in CHUNKS rather than one sample at a time. The symbolic callables are
# vectorised: J_u accepts a (J + D + 1) x n input matrix and returns n x D^2
# (one flattened Jacobian per row, column-major, so matrix(out[i, ], D, D)
# reproduces the single-point call bitwise). Evaluating one sample at a time
# therefore costs 2 * D * nw scalar symbolic calls per solve -- 3066 of them on
# Lorenz at M = 512 -- which profiling put at ~30% of a solve once
# build_ic_noise_quad started sharing this cache. Chunking collapses that to
# 2 * D calls per CHUNK_SIZE samples. Values are unchanged: same h0, same
# central differences, same evaluator.
#
# Chunking keeps the lazy contract (consumers sweep window columns in increasing
# order, so a narrow window still touches only the chunks it needs) without
# needing to know the window in advance.
.IC_HESS_CHUNK <- 64L

build_ic_hessian_cache <- function(U, tt_vec, p, J_u, D) {
  store <- new.env(parent = emptyenv())
  M     <- nrow(U)
  J     <- length(p)

  fill <- function(ms) {
    n  <- length(ms)
    Um <- U[ms, , drop = FALSE]
    h0 <- 1e-5 * pmax(1, sqrt(rowSums(Um^2)))       # per-sample FD step
    pm <- matrix(rep(p, n), nrow = J)
    tm <- matrix(tt_vec[ms], nrow = 1L)
    out <- vector("list", D)
    for (cp in seq_len(D)) {
      Uu <- Um; Uu[, cp] <- Uu[, cp] + h0
      Ud <- Um; Ud[, cp] <- Ud[, cp] - h0
      Ju <- J_u(rbind(pm, t(Uu), tm))               # n x D^2
      Jd <- J_u(rbind(pm, t(Ud), tm))
      out[[cp]] <- (matrix(as.vector(Ju), n, D * D) -
                    matrix(as.vector(Jd), n, D * D)) / (2 * h0)
    }
    for (i in seq_len(n)) {
      v <- array(0, c(D, D, D))
      for (cp in seq_len(D)) v[, , cp] <- matrix(out[[cp]][i, ], D, D)
      store[[as.character(ms[i])]] <- v
    }
  }

  function(m) {
    key <- as.character(m)
    v   <- store[[key]]
    if (is.null(v)) {
      ms <- seq.int(m, min(M, m + .IC_HESS_CHUNK - 1L))
      ms <- ms[vapply(ms, function(k) is.null(store[[as.character(k)]]), TRUE)]
      fill(ms)
      v <- store[[key]]
    }
    v
  }
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
# Both channels contract against a single window column at a time: T_n = dX/dU_n
# is supported on the one column that carries n, so of S X' a_n only the D
# entries sharing that column can survive, and likewise only the V_BL-weighted
# row-block sums of a_n. Those are the block diagonals
#   Z[i, cc, cp]  = X[, (cc-1)nw+i]' A[, (cp-1)nw+i]
#   Va[i, d,  cp] = sum_k V_BL[k, win[i]] A[(d-1)K+k, (cp-1)nw+i]
# of X'A and V_BL'A, so the whole accumulation reduces to D^2 column-sum passes
# plus one matvec per channel -- O(KD nw D^2), not O(KD nw^2 D^2).
#
# P is the IFT projection (Mbar^{-1} Bb' W with Mbar = Bb' W (Bb + EMp)),
# shared with the covariance channels. hess_cache is a build_ic_hessian_cache
# closure; one is built here when the caller has none to share.
# Returns list(b1, b2, b = b1 + b2).
build_ic_bias_o2 <- function(bl, sens, gls, P, EMp, U, tt_vec, p, J_u, sig_vec, dt,
                             hess_cache = NULL) {
  D    <- length(sig_vec)
  K    <- bl$K_bl
  win  <- bl$win_cols
  nw   <- length(win)
  X    <- sens$X
  s2   <- sens$s2
  s2v  <- sig_vec^2
  Yvec <- sens$Bbold + EMp
  A    <- gls$Wm %*% (X - Yvec %*% (P %*% X))   # W (I - (Bb+EMp) P) X
  juu  <- if (!is.null(hess_cache)) hess_cache
          else build_ic_hessian_cache(U, tt_vec, p, J_u, D)

  cblk <- function(c_) ((c_ - 1L) * nw + 1L):(c_ * nw)   # X columns of state c_
  rblk <- function(d_) ((d_ - 1L) * K  + 1L):(d_ * K)    # vec rows of state d_

  # Columns with V_BL == 0 carry u only through Vp, i.e. linearly: no f''
  # curvature and no weight feedback.
  Vw     <- bl$V_BL[, win, drop = FALSE]                 # K x nw
  active <- which(colSums(abs(Vw)) > 0)

  Z  <- array(0, c(nw, D, D))
  Va <- array(0, c(nw, D, D))
  for (cp in seq_len(D)) {
    Acp <- A[, cblk(cp), drop = FALSE]                   # KD x nw
    for (cc in seq_len(D))
      Z[, cc, cp] <- colSums(X[, cblk(cc), drop = FALSE] * Acp)
    for (d in seq_len(D))
      Va[, d, cp] <- colSums(Vw * Acp[rblk(d), , drop = FALSE])
  }

  # Per-window-column coefficients: q the f'' mean channel (b1), g the
  # T_n S X' half and cf the X S T_n' half of dOmega/dU_n (both b2).
  q  <- matrix(0, nw, D)
  g  <- matrix(0, D, nw)
  cf <- numeric(nw * D)
  for (i in active) {
    Juu <- juu(win[i])
    for (d in seq_len(D)) {
      tr_d <- 0
      for (cc in seq_len(D)) tr_d <- tr_d + Juu[d, cc, cc] * s2v[cc]
      q[i, d] <- -dt * 0.5 * tr_d
    }
    for (cp in seq_len(D)) {
      w <- s2[(seq_len(D) - 1L) * nw + i] * Z[i, , cp]   # (S X' a_n) at column i
      for (d in seq_len(D))
        g[d, i] <- g[d, i] - dt * s2v[cp] * sum(w * Juu[d, , cp])
      for (cc in seq_len(D)) {
        col     <- (cc - 1L) * nw + i
        cf[col] <- cf[col] -
                   dt * s2v[cp] * s2[col] * sum(Juu[, cc, cp] * Va[i, , cp])
      }
    }
  }

  cvec <- as.vector(X %*% cf)
  for (d in seq_len(D)) cvec[rblk(d)] <- cvec[rblk(d)] + as.vector(Vw %*% g[d, ])

  b1 <- as.numeric(P %*% as.vector(Vw %*% q))
  b2 <- -as.numeric(P %*% cvec)
  list(b1 = b1, b2 = b2, b = b1 + b2)
}

# A-PRIORI design selection: scores every candidate WITHOUT solving.
#
#   obj = sum_d Var_d/sig_d^2 + sum_d (emb_d + statb_d)^2/sig_d^2
#
# is the same MSE proxy as before, but the fixed point enters it only through
# the point at which EMp, P and the h^4 Euler-Maclaurin term are evaluated --
# the build (V_BL, X, Omega, W, Bbold) does not involve u0 at all. So the raw
# first observation U[1, ] is plugged in there instead of the converged u0hat,
# and no candidate runs a Picard iteration.
#
# The EM order-difference term was the one piece defined by two solves. It is
# linearised instead: going EM(2) -> EM(4) adds Delta = -c4 * phi_t1 %*% g3 to
# the EM block, and P maps a residual perturbation to a u0 shift, so
#   u0_EM2 - u0_EM4 = P vec(Delta) + O(||Delta||^2),
# one evaluation rather than a second fixed point.
select_ic_design <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u,
                             sig_vec, param_cov, em_order,
                             r_bl_grid, rc_cap, include_interior = TRUE,
                             hess_cache = NULL, quad_cov = TRUE,
                             interior_stride = 1L, include_bias_o2 = TRUE) {
  r_bls  <- sort(unique(pmin(pmax(as.integer(r_bl_grid), 2L), rc_cap)))
  D      <- ncol(U)
  M      <- nrow(U)
  J      <- length(p)
  s2     <- sig_vec^2
  tt_vec <- as.vector(tt)
  dt     <- mean(diff(tt_vec))
  c2     <- dt^2 / 12
  c4     <- if (em_order >= 4L) dt^4 / 720 else 0
  u0p    <- as.numeric(U[1, ])            # the a-priori plug-in

  # EM correction at an arbitrary u0, plus its h^4 block alone (Delta4).
  em_parts <- function(bl, u0, p_use = p) {
    u   <- as.vector(u0)
    inp <- matrix(c(p_use, u, tt_vec[1]), ncol = 1L)
    fd  <- list(as.vector(f_(inp)),       as.vector(dF_dt_(inp)),
                as.vector(d2F_dt2_(inp)), as.vector(d3F_dt3_(inp)))
    A <- matrix(0, 5L, D)
    A[1:3, ] <- c2 * g_coeffs(fd, u, 1L)
    A4 <- if (c4 != 0) -c4 * g_coeffs(fd, u, 3L) else matrix(0, 5L, D)
    er  <- bl$em_rows
    phi <- bl$bl_phi_t1[er, , drop = FALSE]
    EM <- matrix(0, bl$K_bl, D); EM[er, ] <- phi %*% (A + A4)
    D4 <- matrix(0, bl$K_bl, D); D4[er, ] <- phi %*% A4
    list(EM = EM, Delta4 = D4)
  }

  r_trap_of <- function(bl, p_use) {
    input <- rbind(matrix(rep(p_use, M), nrow = J), t(U), matrix(tt_vec, nrow = 1L))
    -dt * (bl$V_BL %*% f_(input)) - dt * (bl$Vp_BL %*% U)
  }

  rows <- vector("list", length(r_bls))
  best_obj <- Inf
  best_sys <- NULL
  i <- 0L
  for (r_bl in r_bls) {
    n_bl <- as.integer(min(r_bl, rc_cap))
    res <- tryCatch({
      bl   <- build_ic_bl_system(tt_vec, r_bl, n_bl, orders = 0:4,
                                 include_interior = include_interior,
                                 interior_stride = interior_stride)
      sens <- build_ic_noise_sensitivity(bl, U, tt_vec, p, J_u, sig_vec, dt)
      if (isTRUE(quad_cov))
        sens$Omega2 <- tryCatch(
          build_ic_noise_quad(bl, U, tt_vec, p, J_u, sig_vec, dt,
                              hess_cache = hess_cache),
          error = function(err) NULL)
      gls <- build_ic_gls_weights(sens)

      # IFT sensitivities at the plug-in
      h   <- 1e-6 * max(1, sqrt(sum(u0p^2)))
      EMp <- matrix(0, bl$K_bl * D, D)
      for (d_ in seq_len(D)) {
        up <- u0p; up[d_] <- up[d_] + h
        dn <- u0p; dn[d_] <- dn[d_] - h
        EMp[, d_] <- as.vector((em_parts(bl, up)$EM -
                                em_parts(bl, dn)$EM) / (2 * h))
      }
      P <- solve(gls$BtWB + gls$BtW %*% EMp, gls$BtW)

      G   <- P %*% sens$X
      cov <- G %*% (sens$s2 * t(G))
      if (!is.null(sens$Omega2)) cov <- cov + (P %*% sens$Omega2) %*% t(P)

      # Parameter channel S_p C S_p' (law of total variance), as in cov_u0
      if (!is.null(param_cov)) {
        S_p <- matrix(0, D, J)
        for (j in seq_len(J)) {
          hj <- 1e-6 * max(1, abs(p[j]))
          pu <- p; pu[j] <- pu[j] + hj
          pd <- p; pd[j] <- pd[j] - hj
          du <- as.vector(r_trap_of(bl, pu) - em_parts(bl, u0p, pu)$EM)
          dd <- as.vector(r_trap_of(bl, pd) - em_parts(bl, u0p, pd)$EM)
          S_p[, j] <- P %*% ((du - dd) / (2 * hj))
        }
        cov <- cov + S_p %*% param_cov %*% t(S_p)
      }

      statb <- if (isTRUE(include_bias_o2)) {
        b <- tryCatch(
          build_ic_bias_o2(bl, sens, gls, P, EMp, U, tt_vec, p, J_u,
                           sig_vec, dt, hess_cache = hess_cache)$b,
          error = function(err) NULL)
        if (is.null(b) || !all(is.finite(b))) rep(0, D) else b
      } else rep(0, D)

      emb <- if (c4 != 0)
        as.numeric(P %*% as.vector(em_parts(bl, u0p)$Delta4)) else rep(0, D)

      vobj <- sum(diag(cov) / s2)
      if (!is.finite(vobj)) stop("non-finite objective", call. = FALSE)
      list(obj = vobj + sum((emb + statb)^2 / s2), var_obj = vobj,
           sys = list(bl = bl, sens = sens, gls = gls))
    }, error = function(err) list(obj = NA_real_, var_obj = NA_real_, sys = NULL))
    i <- i + 1L
    rows[[i]] <- data.frame(r_bl = r_bl, n_bl = n_bl,
                            obj = res$obj, var_obj = res$var_obj)

    if (is.finite(res$obj) && res$obj < best_obj) {
      best_obj <- res$obj
      best_sys <- res$sys
    }
  }
  tab <- do.call(rbind, rows)
  if (!any(is.finite(tab$obj))) return(NULL)
  best <- which.min(tab$obj)
  list(table = tab, r_bl = tab$r_bl[best], n_bl = tab$n_bl[best], sys = best_sys)
}

# MSG-convention radius grid for the pool, base * 2^(0:4) clamped to rc_cap.
#
# The base is `min_radius`, i.e. find_min_radius_int_error's pick, which lands in
# [4, 8] in 98.8% of validated fits and is M-independent -- the right semantics for
# an ANCHOR, since a dyadic grid fixes both ends and only spans four octaves.
# GUARD: `min_radius_int_error` is filled by compute_r_c_hat's CHANGEPOINT radius
# (~20) whenever control$test_fun != "phi" (R/test_functions.R:336-347), and an
# anchor that high loses rungs to pmin(., rc_cap) -- {19,38,63} at M=128 -- which
# measured WORST of the three candidate bases. Halving until at least `min_rungs`
# distinct radii survive repairs that path and leaves every validated (Phi) cell
# untouched.
.ic_msg_grid <- function(base, rc_cap, min_rungs = 4L) {
  b <- max(2, as.numeric(base))
  repeat {
    g <- sort(unique(pmin(pmax(as.integer(round(b * 2^(0:4))), 2L), rc_cap)))
    if (length(g) >= min_rungs || b <= 2) return(g)
    b <- b / 2
  }
}

# Nearest PSD matrix: symmetrise, then clip negative eigenvalues. The pooled
# covariance is INVERTED (not just read off the diagonal as cov_u0 is), so a
# rounding-level negative eigenvalue is fatal rather than cosmetic.
.ic_psd <- function(A) {
  A <- (A + t(A)) / 2
  ee <- eigen(A, symmetric = TRUE)
  if (all(ee$values >= 0)) return(A)
  ee$vectors %*% (pmax(ee$values, 0) * t(ee$vectors))
}

# POOL the per-radius weak-form estimates instead of picking one.
#
# Each radius yields its own estimate u0hat_r; their EXACT joint error covariance is
#   Sig[a,b] = P_a ( X_a S X_b' + Omega2_ab ) P_b'  +  S_p,a Chat S_p,b'
# every piece of which estimate_IC already builds. Pooling by GLS over the L
# estimates keeps the information the argmin used to discard: cross-radius error
# correlations are only 0.37-0.92 at nr <= 0.4.
#
# Three details are load-bearing (each 1.1-1.6x when wrong; examples/validation/):
#   * pool the ESTIMATES, not the rows. Row-stacking makes Omega rank-deficient
#     (rank <= M*D << K*D) so the ridge, not the data, picks the answer.
#   * `bias_floor`: diag(Sig) += (trunc + bias_o2)^2. MANDATORY -- GLS weights by
#     variance alone, and a small radius has SMALL variance with a HUGE bias, so it
#     hijacks the pool. The floor's value is its RELATIVE structure across members,
#     not its scale (a fitted scalar heterogeneity is much worse).
#   * weights from the correlation SHRUNK toward I, but the covariance REPORTED
#     under the unshrunk Sig, or coverage drops to 0.79-0.87.
# A member is admitted only if its fixed point converged.
pool_ic_radii <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, sig_vec,
                          param_cov = NULL, em_order = 4L, r_bl_grid, rc_cap,
                          include_interior = TRUE, interior_stride = 1L,
                          quad_cov = TRUE, hess_cache = NULL, debias = TRUE,
                          shrink = 0.5, max_iter = 100L, tol = 1e-9) {
  M <- nrow(U); D <- ncol(U); J <- length(p)
  tt_vec <- as.vector(tt); dt <- mean(diff(tt_vec))
  c2 <- dt^2 / 12
  c4 <- if (em_order >= 4L) dt^4 / 720 else 0
  r_bls <- sort(unique(pmin(pmax(as.integer(r_bl_grid), 2L), rc_cap)))
  if (!length(r_bls)) return(NULL)
  if (is.null(hess_cache)) hess_cache <- build_ic_hessian_cache(U, tt_vec, p, J_u, D)
  if (!is.null(param_cov)) param_cov <- .ic_psd(as.matrix(param_cov))
  u0p <- as.numeric(U[1, ])

  # EM correction and its h^4 block alone, at an arbitrary u0
  em_parts <- function(bl, u0, p_use = p) {
    u   <- as.vector(u0)
    inp <- matrix(c(p_use, u, tt_vec[1]), ncol = 1L)
    fd  <- list(as.vector(f_(inp)),       as.vector(dF_dt_(inp)),
                as.vector(d2F_dt2_(inp)), as.vector(d3F_dt3_(inp)))
    A  <- matrix(0, 5L, D)
    A[1:3, ] <- c2 * g_coeffs(fd, u, 1L)
    A4 <- if (c4 != 0) -c4 * g_coeffs(fd, u, 3L) else matrix(0, 5L, D)
    er  <- bl$em_rows
    phi <- bl$bl_phi_t1[er, , drop = FALSE]
    EM <- matrix(0, bl$K_bl, D); EM[er, ] <- phi %*% (A + A4)
    D4 <- matrix(0, bl$K_bl, D); D4[er, ] <- phi %*% A4
    list(EM = EM, Delta4 = D4)
  }
  r_trap_of <- function(bl, p_use) {
    input <- rbind(matrix(rep(p_use, M), nrow = J), t(U), matrix(tt_vec, nrow = 1L))
    -dt * (bl$V_BL %*% f_(input)) - dt * (bl$Vp_BL %*% U)
  }

  fit_one <- function(r_bl) {
    n_bl <- as.integer(min(r_bl, rc_cap))
    bl   <- build_ic_bl_system(tt_vec, r_bl, n_bl, orders = 0:4,
                               include_interior = include_interior,
                               interior_stride = interior_stride)
    if (!is.finite(bl$BtB) || bl$BtB < .Machine$double.eps) return(NULL)
    sens <- build_ic_noise_sensitivity(bl, U, tt_vec, p, J_u, sig_vec, dt)
    if (isTRUE(quad_cov))
      sens$Omega2 <- tryCatch(
        build_ic_noise_quad(bl, U, tt_vec, p, J_u, sig_vec, dt, hess_cache = hess_cache),
        error = function(err) NULL)
    gls <- build_ic_gls_weights(sens)

    r_trap <- r_trap_of(bl, p)
    proj   <- function(rhs) as.numeric(solve(gls$BtWB, gls$BtW %*% as.vector(rhs)))
    resn   <- function(u0) {
      e <- outer(bl$B, as.numeric(u0)) - (r_trap - em_parts(bl, u0)$EM)
      m <- as.numeric(gls$BtW %*% as.vector(e)); sqrt(sum(m * m))
    }
    u0 <- proj(r_trap)
    if (!all(is.finite(u0))) return(NULL)
    best <- u0; bres <- resn(u0); conv <- FALSE
    for (it in seq_len(max_iter)) {
      un <- proj(r_trap - em_parts(bl, u0)$EM)
      if (!all(is.finite(un))) return(NULL)
      rn <- resn(un)
      if (is.finite(rn) && rn < bres) { bres <- rn; best <- un }
      dl <- sqrt(sum((un - u0)^2)); u0 <- un
      if (is.finite(dl) && dl < tol * max(1, sqrt(sum(un^2)))) { conv <- TRUE; break }
    }
    if (!conv) return(NULL)
    u0 <- best

    KD  <- bl$K_bl * D
    EMp <- matrix(0, KD, D); h <- 1e-6 * max(1, sqrt(sum(u0^2)))
    for (e_i in seq_len(D)) {
      up <- u0; up[e_i] <- up[e_i] + h
      dn <- u0; dn[e_i] <- dn[e_i] - h
      EMp[, e_i] <- as.vector((em_parts(bl, up)$EM - em_parts(bl, dn)$EM) / (2 * h))
    }
    P <- tryCatch(solve(gls$BtWB + gls$BtW %*% EMp, gls$BtW), error = function(err) NULL)
    if (is.null(P)) return(NULL)

    G  <- P %*% sens$X
    cn <- G %*% (sens$s2 * t(G))
    if (!is.null(sens$Omega2)) cn <- cn + (P %*% sens$Omega2) %*% t(P)

    S_p <- matrix(0, D, J)
    if (!is.null(param_cov)) for (j in seq_len(J)) {
      hj <- 1e-6 * max(1, abs(p[j])); pu <- p; pd <- p
      pu[j] <- pu[j] + hj; pd[j] <- pd[j] - hj
      du <- as.vector(r_trap_of(bl, pu) - em_parts(bl, u0, pu)$EM)
      dd <- as.vector(r_trap_of(bl, pd) - em_parts(bl, u0, pd)$EM)
      S_p[, j] <- P %*% ((du - dd) / (2 * hj))
    }

    # O(sigma^2) debias, gated exactly as the single-design path gates it
    b <- if (isTRUE(debias)) tryCatch(
      build_ic_bias_o2(bl, sens, gls, P, EMp, U, tt_vec, p, J_u, sig_vec, dt,
                       hess_cache = hess_cache)$b, error = function(err) NULL) else NULL
    if (is.null(b) || !all(is.finite(b))) b <- rep(0, D)
    applied <- FALSE
    if (isTRUE(debias) && any(b != 0) && all(abs(b) <= 2 * sqrt(pmax(diag(cn), 0)))) {
      u0 <- u0 - b; applied <- TRUE
    }
    # linearised EM order-difference (one evaluation, no second fixed point)
    trunc <- if (c4 != 0) as.numeric(P %*% as.vector(em_parts(bl, u0p)$Delta4)) else rep(0, D)

    list(r_bl = r_bl, n_bl = n_bl, u0 = u0, P = P, X = sens$X, S_p = S_p,
         V = bl$V_BL, win = bl$win_cols, K = bl$K_bl, bias = b, trunc = trunc,
         debias_applied = applied, cov_noise = cn)
  }

  ms <- lapply(r_bls, function(r) tryCatch(fit_one(r), error = function(err) NULL))
  ok <- !vapply(ms, is.null, logical(1))
  if (!any(ok)) return(NULL)
  ms <- ms[ok]; L <- length(ms); n <- L * D

  # one curvature array on the union window; cross blocks index it per design
  wu <- sort(unique(unlist(lapply(ms, `[[`, "win"))))
  Gq <- if (isTRUE(quad_cov)) tryCatch(
    build_ic_quad_G(U, tt_vec, p, J_u, sig_vec, wu, hess_cache = hess_cache),
    error = function(err) NULL) else NULL
  if (!is.null(Gq) && max(abs(Gq)) <= 0) Gq <- NULL
  s2c <- sig_vec^2

  Sig <- matrix(0, n, n)
  for (a in seq_len(L)) for (b in a:L) {
    ma <- ms[[a]]; mb <- ms[[b]]
    cols <- intersect(ma$win, mb$win)
    ia <- match(cols, ma$win); ib <- match(cols, mb$win)
    nwa <- length(ma$win);     nwb <- length(mb$win)
    Oab <- matrix(0, ma$K * D, mb$K * D)
    for (cc in seq_len(D))
      Oab <- Oab + s2c[cc] * (ma$X[, (cc - 1L) * nwa + ia, drop = FALSE] %*%
                              t(mb$X[, (cc - 1L) * nwb + ib, drop = FALSE]))
    if (!is.null(Gq)) {
      iu <- match(cols, wu)
      Va <- ma$V[, cols, drop = FALSE]; Vb <- mb$V[, cols, drop = FALSE]
      for (d in seq_len(D)) for (dp in seq_len(D)) {
        ra <- ((d  - 1L) * ma$K + 1L):(d  * ma$K)
        rb <- ((dp - 1L) * mb$K + 1L):(dp * mb$K)
        Oab[ra, rb] <- Oab[ra, rb] + dt^2 * (Va %*% (Gq[iu, d, dp] * t(Vb)))
      }
    }
    blk <- ma$P %*% Oab %*% t(mb$P)
    if (!is.null(param_cov)) blk <- blk + ma$S_p %*% param_cov %*% t(mb$S_p)
    ja <- ((a - 1L) * D + 1L):(a * D); jb <- ((b - 1L) * D + 1L):(b * D)
    Sig[ja, jb] <- blk
    if (b != a) Sig[jb, ja] <- t(blk)
  }
  Sig <- .ic_psd(Sig)

  uvec <- as.vector(vapply(ms, function(m) as.numeric(m$u0), numeric(D)))
  flr  <- as.vector(vapply(ms, function(m) m$trunc + m$bias, numeric(D)))
  Sig0 <- Sig; diag(Sig0) <- diag(Sig0) + flr^2      # model error IS part of cov_u0
  Sw   <- Sig0
  if (shrink > 0) {
    sd_ <- sqrt(pmax(diag(Sw), 0))
    Cr  <- Sw / outer(sd_, sd_)
    Sw  <- ((1 - shrink) * Cr + shrink * diag(n)) * outer(sd_, sd_)
  }
  Ones <- do.call(rbind, replicate(L, diag(D), simplify = FALSE))
  Si <- tryCatch(solve(Sw + 1e-10 * mean(diag(Sw)) * diag(n)), error = function(err) NULL)
  if (is.null(Si)) return(NULL)
  A  <- crossprod(Ones, Si)
  Cw <- tryCatch(solve(A %*% Ones), error = function(err) NULL)
  if (is.null(Cw)) return(NULL)
  W  <- Cw %*% A                                     # D x LD pooling weights

  wn <- vapply(seq_len(L), function(a)
    sqrt(sum(W[, ((a - 1L) * D + 1L):(a * D), drop = FALSE]^2)), numeric(1))
  list(u0hat  = as.numeric(W %*% uvec),
       cov_u0 = W %*% Sig0 %*% t(W),
       w = W, w_norm = wn, Sigma = Sig0,
       r_bl = vapply(ms, `[[`, numeric(1), "r_bl"),
       n_bl = vapply(ms, `[[`, numeric(1), "n_bl"),
       u0_per_radius = uvec,
       bias_o2 = as.vector(vapply(ms, function(m) m$bias, numeric(D))),
       debias_applied = vapply(ms, `[[`, logical(1), "debias_applied"),
       L = L)
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
#' is \code{NULL}, the design \code{(r_bl, n_bl)} is selected a priori by
#' sweeping a small grid and minimizing a calibrated \eqn{u_0}-MSE proxy
#' \deqn{\textstyle\sum_d \mathrm{Var}_d/\sigma_d^2
#'       + \sum_d ((u_0^{EM2} - u_0^{EM4})_d + b_d)^2/\sigma_d^2,}
#' where \eqn{\mathrm{Var} = \mathrm{diag}(\mathrm{cov\_u0})} folds the EM
#' Jacobian and the parameter channel \eqn{S_p \hat C S_p^T},
#' \eqn{u_0^{EM2} - u_0^{EM4}} is the oracle-free Euler-Maclaurin order-
#' difference estimate of the truncation defect, and \eqn{b} is the analytic
#' \eqn{O(\sigma^2)} statistical bias. The selection is genuinely A PRIORI: no
#' candidate is solved. Every term reaches the fixed point only through the
#' point at which \eqn{\partial \mathrm{EM}/\partial u_0}, the IFT projection
#' \eqn{P} and the \eqn{h^4} endpoint term are evaluated, so the raw first
#' observation \code{U[1, ]} is plugged in there, and the order-difference term
#' is linearised as \eqn{u_0^{EM2} - u_0^{EM4} = P\,\mathrm{vec}(\Delta)} with
#' \eqn{\Delta = -c_4 \phi(t_1) g^{(3)}} (one evaluation, not two solves).
#' Validated in \code{examples/validation/ic_apriori_proxy.R}: within-rep rank
#' correlation 0.96-1.00 against the former solve-per-candidate objective, same
#' median radius selected in every cell. This replaces the earlier
#' noise-only variance objective, which was
#' monotone in window/count and so pinned the corner; the MSE proxy turns up
#' where the true MSE turns up (validated examples/validation/ic_mse_proxy.R,
#' ic_calibrated_criterion.R, ic_realC_overshoot.R). \code{combine = "ols"}
#' restores the legacy unweighted combine and its \code{max(3, ceiling(r_bl/8))}
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
#' \eqn{\Omega} carries two pieces. The linear (delta-method) part
#' \eqn{X \mathrm{diag}(\sigma^2) X^T} is the covariance of \eqn{X\eta}; the
#' residual is also quadratic in the noise through \eqn{f''}, contributing an
#' \eqn{O(\sigma^4)} block \eqn{\Omega_2} (see \code{build_ic_noise_quad}) that
#' \code{quad_cov = TRUE} (default) adds to both the weights and
#' \code{cov_u0}. \eqn{\Omega_2} is only 0.03-3\% of \eqn{\mathrm{tr}(\Omega)},
#' but \eqn{\Omega_1} is near-singular and \eqn{\Omega_2} dominates in exactly
#' the near-null directions the GLS combine loads, so omitting it makes the
#' combine over-trust cancellation that is not there. Deployed effect at 20\%
#' noise: \eqn{u_0} MSE \eqn{\times}0.38 on Lorenz (coverage 0.69 -> 0.99,
#' RMS z 2.45 -> 0.97) and \eqn{\times}0.64 on Lotka-Volterra, neutral on
#' logistic (examples/validation/ic_audit_omega2_check.R,
#' ic_audit_combined.R).
#'
#' @param U Numeric matrix (M x D) of observed states.
#' @param f_,dF_dt_,d2F_dt2_,d3F_dt3_ Callable RHS and total time-derivative
#'   evaluators built from the symbolic engine.
#' @param tt Numeric vector (length M) of time points.
#' @param p Numeric parameter vector \eqn{\hat\theta} (held fixed).
#' @param n_bl Optional integer; number of left (peak-inside) BL test functions.
#'   When \code{NULL} (default) and \code{combine = "gls"}, the BL radius
#'   \code{r_bl} is selected from \code{r_bl_grid} by the a-priori MSE proxy and
#'   the count is maxed, \code{n_bl = r_bl} (one BL function per sample out to the
#'   radius; see Details). When \code{NULL} under \code{combine = "ols"} the
#'   legacy heuristic \code{max(3, ceiling(r_bl/8))} applies (small count, wide
#'   placement -- under the unweighted combine, clustered test functions inflate
#'   Var(u0hat) by ~7-19\%; see examples/validation/). An explicit value skips
#'   the selection and uses \code{r_bl}.
#' @param r_bl Optional integer; the boundary-layer radius for the paths that do
#'   NOT sweep -- an explicit \code{n_bl}, \code{combine = "ols"}, or a sweep
#'   that failed. Ignored when the a-priori selection runs, since that picks
#'   \code{r_bl} from \code{r_bl_grid}. Defaults to \code{min(16, (M-1)/2)}.
#' @param max_iter,tol Fixed-point iteration controls. \code{tol} is compared
#'   against the step norm scaled by \code{max(1, ||u_0||)}, i.e. it is a
#'   relative tolerance for states of size \eqn{\ge 1}. The default \code{1e-10}
#'   is deliberately above the achievable roundoff floor so that
#'   \code{converged} is a usable health flag rather than always \code{FALSE}
#'   on large-amplitude states -- callers gate the divergence fallback on it
#'   (see the divergence guard in \code{solveWendy}).
#' @param em_order Either 2 or 4.
#' @param combine \code{"gls"} (default) for the minimum-variance GLS combine
#'   of the BL equations, or \code{"ols"} for the legacy unweighted combine.
#'   GLS requires a valid \code{sigma} (and \code{J_u}) to build
#'   \eqn{\Omega}; it degrades to \code{"ols"} otherwise.
#' @param interior_stride Integer (default \code{1}, i.e. keep every centre).
#'   Retains only every s-th interior control-variate centre. At one centre per
#'   sample the interior block can be over-complete enough for the GLS combine
#'   to claim cancellation that is not there, and thinning is the lever for
#'   that; the effect is confined to small radii and does not reach the radii
#'   the a-priori sweep selects, so it is offered as an option rather than
#'   imposed (measured in \code{examples/validation/ic_phat_design.R}; strides
#'   above 2 break down). See \code{build_ic_bl_system}.
#' @param quad_cov Logical (default \code{TRUE}). Include the
#'   \eqn{O(\sigma^4)} quadratic-noise block \eqn{\Omega_2} in the GLS weights
#'   and in \code{cov_u0} (see \code{build_ic_noise_quad} and Details). GLS path
#'   only.
#' @param debias Logical (default \code{TRUE}). On the GLS path, subtract the
#'   analytic \eqn{O(\sigma^2)} bias (both channels; see Details). Ignored on
#'   the OLS path and on diverged solves.
#' @param return_em2_u0 Logical (default \code{FALSE}, internal). When
#'   \code{TRUE} and \code{em_order == 4}, additionally solve the EM(2) fixed
#'   point on the same built system and return it as \code{u0hat_em2} (with
#'   \code{em2_diverged}). Used by the design sweep to form the EM-order-
#'   difference truncation estimate without a second \code{estimate_IC} call.
#' @param hess_cache Optional memoizing state-Hessian closure (internal, from
#'   \code{build_ic_hessian_cache}). The \eqn{f''} tensors the
#'   \eqn{O(\sigma^2)} debias needs depend on \code{U}, \code{p} and \code{tt}
#'   alone, not on the boundary-layer design, so one cache is built here and
#'   shared across every candidate of the a-priori sweep and the final solve.
#' @param r_bl_grid Optional integer vector of candidate BL radii for the
#'   a-priori selection (default the absolute grid
#'   \code{c(4,8,12,16,20,24,32,40,48)} capped at \code{floor((M-1)/2)}); at each
#'   candidate the count is maxed (\code{n_bl = r_bl}). Only used when
#'   \code{combine = "gls"} and \code{n_bl} is \code{NULL}.
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
#'   \code{"diverged"} / \code{"not_converged"}), plus \code{combine} (the combine actually used,
#'   after any degradation), \code{design} (the design-selection table with
#'   columns \code{r_bl}, \code{n_bl}, \code{obj} (the MSE proxy), and
#'   \code{var_obj} (its variance part), or \code{NULL} when no selection ran),
#'   \code{bias_o2} (the analytic
#'   \eqn{O(\sigma^2)} bias estimate of the UNcorrected solve, length-D, or
#'   \code{NULL} when not computed), \code{debias_applied} (whether
#'   \code{bias_o2} was subtracted from \code{u0hat}), \code{fallback}
#'   (\code{TRUE} only when the design was degenerate and \code{u0hat} is the
#'   raw observation; the divergence guard is applied by the caller, see
#'   \code{solveWendy}), \code{u0hat_em2} (the
#'   EM(2) fixed point, or \code{NULL} unless \code{return_em2_u0} and
#'   \code{em_order == 4}), \code{em2_diverged}, \code{iters},
#'   \code{converged}, \code{diverged}, \code{u0_history},
#'   \code{r_bl} (the BL radius actually used),
#'   \code{n_bl}, \code{K_bl}, \code{em_order}.
#' @export
estimate_IC <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, sigma,
                        param_cov      = NULL,
                        n_bl           = NULL,
                        r_bl           = NULL,
                        max_iter       = 100L,
                        tol            = 1e-9,
                        em_order       = c(4L, 2L),
                        combine        = c("gls", "ols"),
                        include_interior = TRUE,
                        interior_stride  = 1L,
                        quad_cov       = TRUE,
                        r_bl_grid      = NULL,
                        pool_radii     = TRUE,
                        debias         = TRUE,
                        return_em2_u0  = FALSE,
                        hess_cache     = NULL) {
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
  r_bl_in <- r_bl                 # NULL unless the caller pinned a radius
  # BL radius for the paths that do NOT sweep (explicit n_bl, OLS, or a failed
  # sweep). The a-priori sweep selects r_bl from an absolute grid, so no
  # integration-error radius is needed here any more.
  r_bl_fixed <- if (!is.null(r_bl)) min(as.integer(r_bl), rc_cap)
                else                min(16L, rc_cap)
  r_bl <- r_bl_fixed           # BL radius actually used (reset below if swept)

  use_noiseprop <- length(sigma) %in% c(1L, D) && all(is.finite(sigma))
  sig_vec <- if (use_noiseprop) {
    if (length(sigma) == 1L) rep(sigma, D) else as.numeric(sigma)
  } else NULL

  if (combine == "gls" && !use_noiseprop) combine <- "ols"

  # Shared by the sweep candidates and the final debias (design-independent)
  if (is.null(hess_cache) && use_noiseprop){
    hess_cache <- build_ic_hessian_cache(U, tt_vec, p, J_u, D)
  }

  # POOL over radii -- the default. Combines the per-radius estimates by GLS using
  # their exact cross-radius covariance instead of keeping the argmin of an MSE
  # proxy, and costs 0.47-0.69x the sweep it replaces (5 solves vs 13 candidates
  # plus a final solve). Falls through to the sweep if it cannot be formed.
  if (isTRUE(pool_radii) && combine == "gls" && is.null(n_bl) && is.null(r_bl_in)) {
    pgrid <- if (!is.null(r_bl_grid)) r_bl_grid
             else c(4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 80, 100)
    pl <- tryCatch(
      pool_ic_radii(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, sig_vec,
                    param_cov = param_cov, em_order = em_order,
                    r_bl_grid = pgrid, rc_cap = rc_cap,
                    include_interior = include_interior,
                    interior_stride = interior_stride, quad_cov = quad_cov,
                    hess_cache = hess_cache, debias = debias,
                    max_iter = max_iter, tol = tol),
      error = function(err) NULL)
    if (!is.null(pl) && all(is.finite(pl$u0hat))) {
      best  <- which.max(pl$w_norm)
      U_hat <- U; U_hat[1, ] <- pl$u0hat
      return(list(
        U_hat = U_hat, u0hat = pl$u0hat, cov_u0 = pl$cov_u0,
        cov_u0_resid = NULL, cov_u0_param = NULL, cov_method = "pool",
        combine = combine,
        design = data.frame(r_bl = pl$r_bl, n_bl = pl$n_bl, w_norm = pl$w_norm,
                            debias_applied = pl$debias_applied),
        bias_o2 = pl$bias_o2, debias_applied = any(pl$debias_applied),
        fallback = FALSE, u0hat_em2 = NULL, em2_diverged = FALSE,
        iters = NA_integer_, converged = TRUE, diverged = FALSE,
        u0_history = matrix(pl$u0hat, nrow = 1L),
        r_bl = pl$r_bl[best], n_bl = pl$n_bl[best],
        K_bl = NA_integer_, K_int = NA_integer_, em_order = em_order,
        pool_r_bl = pl$r_bl, pool_w = pl$w, pool_L = pl$L,
        u0_per_radius = pl$u0_per_radius))
    }
  }

  design_table <- NULL
  sel_sys      <- NULL          # the winning candidate's already-built system
  if (combine == "gls" && is.null(n_bl)) {
    if (is.null(r_bl_grid)){
      r_bl_grid <- c(4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 80, 100)
    }
    sel <- tryCatch(
      select_ic_design(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u,
                       sig_vec, param_cov, em_order, r_bl_grid, rc_cap,
                       include_interior = include_interior, include_bias_o2 = TRUE,
                       hess_cache = hess_cache, quad_cov = quad_cov,
                       interior_stride = interior_stride),
      error = function(err) NULL)
    if (!is.null(sel)) {
      design_table <- sel$table
      r_bl         <- sel$r_bl
      n_bl         <- sel$n_bl
      sel_sys      <- sel$sys
    }
  }

  n_bl <- if (!is.null(n_bl)) max(1L, as.integer(n_bl))
          else  as.integer(min(r_bl, rc_cap))

  build_system <- function(r_bl_use, n_bl_use) {
    bl <- build_ic_bl_system(tt_vec, r_bl_use, n_bl_use, orders = 0:4,
                             include_interior = isTRUE(include_interior) &&
                                                combine == "gls",
                             interior_stride = interior_stride)
    sens <- if (use_noiseprop) tryCatch(
      build_ic_noise_sensitivity(bl, U, tt_vec, p, J_u, sig_vec, dt),
      error = function(err) NULL) else NULL

    if (isTRUE(quad_cov) && combine == "gls" && !is.null(sens))
      sens$Omega2 <- tryCatch(
        build_ic_noise_quad(bl, U, tt_vec, p, J_u, sig_vec, dt, hess_cache = hess_cache),
        error = function(err) NULL)
    gls <- if (combine == "gls" && !is.null(sens)) tryCatch(
      build_ic_gls_weights(sens),
      error = function(err) NULL) else NULL
    list(bl = bl, sens = sens, gls = gls)
  }

  # The sweep already built this exact design (same bl / sens / Omega2 / gls
  # arguments), so reuse it rather than paying for one more candidate.
  sys <- if (!is.null(sel_sys)) sel_sys else build_system(r_bl, n_bl)
  if (combine == "gls" && is.null(sys$gls)) {
    combine <- "ols"
    if (!is.null(design_table)) {
      design_table <- NULL
      r_bl <- r_bl_fixed
      n_bl <- max(3L, as.integer(ceiling(r_bl / 8)))
    }
    sys <- build_system(r_bl, n_bl)   # rebuilt without the interior block
  }

  bl        <- sys$bl
  sens      <- sys$sens
  gls       <- sys$gls
  V_BL      <- bl$V_BL
  Vp_BL     <- bl$Vp_BL
  bl_phi_t1 <- bl$bl_phi_t1
  B         <- bl$B
  BtB       <- bl$BtB
  K_bl      <- bl$K_bl     # total stacked equations (BL + interior)
  K_int     <- bl$K_int
  em_rows   <- bl$em_rows  # rows with nonzero phi(t_1); interior EM is exactly 0

  if (!is.finite(BtB) || BtB < .Machine$double.eps) {
    u0_obs <- as.numeric(U[1, ])
    U_hat <- U; U_hat[1, ] <- u0_obs
    return(list(U_hat = U_hat, u0hat = u0_obs, cov_u0 = NULL,
                iters = 0L, converged = FALSE, diverged = FALSE,
                u0_history = matrix(u0_obs, nrow = 1),
                r_bl = r_bl, n_bl = n_bl, K_bl = K_bl, K_int = K_int,
                em_order = em_order,
                combine = combine, design = design_table,
                bias_o2 = NULL, debias_applied = FALSE, fallback = TRUE
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

  em_correction <- function(u0_curr, p_use = p, c4_use = c4) {
    u_t1   <- as.vector(u0_curr)
    inp_t1 <- matrix(c(p_use, u_t1, tt_vec[1]), ncol = 1L)
    fd_t1  <- list(
      as.vector(f_(inp_t1)),
      as.vector(dF_dt_(inp_t1)),
      as.vector(d2F_dt2_(inp_t1)),
      as.vector(d3F_dt3_(inp_t1))
    )
    A <- matrix(0, 5L, D)
    A[1:3, ] <- c2 * g_coeffs(fd_t1, u_t1, 1L)
    if (c4_use != 0) A <- A - c4_use * g_coeffs(fd_t1, u_t1, 3L)
    EM <- matrix(0, nrow = K_bl, ncol = D)   # interior rows stay exactly 0
    EM[em_rows, ] <- bl_phi_t1[em_rows, , drop = FALSE] %*% A
    EM
  }

  proj <- if (!is.null(gls)) {
    function(rhs) as.numeric(solve(gls$BtWB, gls$BtW %*% as.vector(rhs)))
  } else {
    function(rhs) as.numeric(crossprod(B, rhs) / BtB)
  }

  residual_norm <- function(u0_curr, c4_use = c4) {
    EM <- em_correction(u0_curr, c4_use = c4_use)
    e  <- outer(B, as.numeric(u0_curr)) - (r_trap - EM)   # K_bl x D
    m  <- if (!is.null(gls)) as.numeric(gls$BtW %*% as.vector(e))   # D
          else               as.numeric(crossprod(B, e))            # D
    sqrt(sum(m * m))
  }

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
      if (is.finite(delta) && delta < tol * max(1, sqrt(sum(u0_new^2)))) {
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
    compute_resid <- function() tryCatch({
      EM_final <- em_correction(u0)
      rhs      <- r_trap - EM_final            # K_bl x D
      e        <- outer(B, u0) - rhs            # K_bl x D residuals
      df       <- max(K_bl - 1L, 1L)
      s2       <- colSums(e * e) / df           # length-D
      diag(s2 / BtB, nrow = D, ncol = D)
    }, error = function(err) NULL)

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

    # Noise channel: P Omega P^T. The linear part is (P X) diag(s2) (P X)^T; the
    # quadratic block (build_ic_noise_quad) has to be carried here as well as in
    # the weights, otherwise the reported SE keeps understating the variance in
    # exactly the near-null directions the combine loads (LV at 20% noise:
    # empirical SD / reported SE was 2.9-4.0 before, ~1.0 after).
    cov_u0_noise <- if (!is.null(sens) && !is.null(P)) tryCatch({
      G <- P %*% sens$X
      C <- G %*% (sens$s2 * t(G))
      if (!is.null(sens$Omega2)) C <- C + (P %*% sens$Omega2) %*% t(P)
      C
    }, error = function(err) NULL) else NULL

    #
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

    cov_u0_resid <- compute_resid()

    cov_u0     <- if (!is.null(cov_u0_noise)) cov_u0_noise else cov_u0_resid
    cov_method <- if (!is.null(cov_u0_noise)) "noise_propagation" else "ls_residual"
    if (!is.null(cov_u0) && !is.null(cov_u0_param)) {
      cov_u0     <- cov_u0 + cov_u0_param
      cov_method <- paste0(cov_method, "+param")
    }

    list(cov_u0 = cov_u0, cov_method = cov_method, cov_u0_resid = cov_u0_resid,
         cov_u0_param = cov_u0_param, cov_u0_noise = cov_u0_noise, P = P, EMp = EMp)
  }

  cov_res <- compute_covariance(u0)
  cov_u0 <- cov_res$cov_u0
  cov_method   <- cov_res$cov_method
  cov_u0_resid <- cov_res$cov_u0_resid
  cov_u0_param <- cov_res$cov_u0_param
  cov_u0_noise <- cov_res$cov_u0_noise
  P <- cov_res$P
  EMp <- cov_res$EMp

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
                       sig_vec, dt, hess_cache = hess_cache)$b,
      error = function(err) NULL)
    if (!is.null(bias_o2) && all(is.finite(bias_o2))) {
      se_gate <- sqrt(pmax(diag(cov_u0_noise), 0))
      if (all(abs(bias_o2) <= 2 * se_gate)) {
        u0             <- u0 - bias_o2
        debias_applied <- TRUE
      }
    }
  }

  # The divergence / non-convergence guard is NOT applied here: the raw iterate
  # and its health flags (converged, diverged) are returned as-is, and the
  # caller decides whether to fall back to U[1, ] (see solveWendy in wendy.R).
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
    fallback       = FALSE,
    u0hat_em2      = u0hat_em2,
    em2_diverged   = em2_diverged,
    iters          = iters,
    converged      = converged,
    diverged       = diverged,
    u0_history     = do.call(rbind, u0_hist),
    r_bl           = r_bl,
    n_bl           = n_bl,
    K_bl           = K_bl,
    K_int          = K_int,
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
#'   \code{u0_init} is supplied, defaults to \code{sigma^2 I_D}. When neither is
#'   supplied the prior is \emph{diffuse} (\code{1e6 * sigma^2 I_D}), so the
#'   \eqn{k = 1} update returns \code{U[1, ]} with covariance \code{sigma^2 I_D}
#'   rather than double-counting \code{U[1, ]} as both prior mean and
#'   measurement (which would report \code{sigma^2 / 2}).
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
wendy_erts <- function(U, f_, J_u, tt, p, test_function_params, sigma = NULL, u0_init = NULL, P0_init = NULL, param_cov = NULL, J_p = NULL,
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
  # e_k = (16/15)(y1 - y2) (order p = 4, factor 2^p/(2^p - 1)), so Q_k = diag(e_k^2)
  q_floor <- 1e-12

  # One RK4 step of the process f(p, ., .)
  rk4_step <- function(p_use, uk, tk, h) {
    k1 <- as.vector(f_(matrix(c(p_use, uk,           tk        ), ncol = 1)))
    k2 <- as.vector(f_(matrix(c(p_use, uk + .5*h*k1, tk + .5*h ), ncol = 1)))
    k3 <- as.vector(f_(matrix(c(p_use, uk + .5*h*k2, tk + .5*h ), ncol = 1)))
    k4 <- as.vector(f_(matrix(c(p_use, uk +    h*k3, tk +    h ), ncol = 1)))
    uk + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)
  }

  # With no prior supplied the filter has no information about u(0) independent
  # of the data, so the honest encoding is a DIFFUSE prior. Taking
  # P0 = noise_sd^2 I instead would double-count U[1, ] -- it would be both the
  # prior mean and the k = 1 measurement, two sources the Kalman update treats
  # as independent when they are the identical number. That gives K0 = I/2 and
  # P_filt[1] = noise_sd^2 / 2, half the variance one noisy sample supports.
  # With P0 = diffuse_scale * R_obs the k = 1 update returns exactly
  # (U[1, ], R_obs) up to O(1/diffuse_scale). Validated in
  # examples/validation/ic_kf_diffuse_prior.R: P_filt[1] doubles as predicted;
  # the RTS pass dilutes it to nothing where later data pins u(0) (logistic/LV
  # SE ratio ~1.00) but on chaos it lifts u(0) coverage .61-.69 -> .75-.78.
  diffuse_scale <- 1e6
  u0 <- if (!is.null(u0_init)) as.vector(u0_init) else as.vector(U[1, ])
  P0 <- if (!is.null(P0_init)) P0_init
        else if (!is.null(u0_init)) noise_sd^2 * I_D
        else diffuse_scale * noise_sd^2 * I_D

  # One EKF forward pass + RTS backward pass at a fixed p̂
  erts_pass <- function(p_use, want_cov) {
    S0 <- P0 + R_obs
    K0 <- P0 %*% solve(S0)

    u_filt <- matrix(0, mp1, D)
    u_filt[1, ] <- u0 + K0 %*% (U[1, ] - u0)
    P_filt <- array(0, c(mp1, D, D))
    P_filt[1,,] <- (I_D - K0) %*% P0 %*% t(I_D - K0) + K0 %*% R_obs %*% t(K0)

    u_pred  <- matrix(0, mp1, D)
    P_pred  <- array(0, c(mp1, D, D))
    F_store <- array(0, c(mp1 - 1L, D, D))

    for (k in seq_len(mp1 - 1L)) {

      dt_k <- tt[k + 1L] - tt[k]  # dt
      uk <- u_filt[k, ] # uk current time step

      # Predict
      y1 <- rk4_step(p_use, uk, tt[k], dt_k)                  # full step = predict mean
      yh <- rk4_step(p_use, uk, tt[k], dt_k / 2)              # step-doubling: two half
      y2 <- rk4_step(p_use, yh, tt[k] + dt_k / 2, dt_k / 2)   #   steps for the error est.

      u_pred[k + 1L, ] <- y1 # Predicted step from model + control

      e_k <- (16 / 15) * (y1 - y2)                            # RK4 local truncation error
      Ju_k <- matrix(as.vector(J_u(c(p_use, uk, tt[k]))), D, D)  # J[a, b] = df_a/du_b
      Fk <- I_D + dt_k * Ju_k
      F_store[k,,] <- Fk

      Pk_pred <- Fk %*% P_filt[k,,] %*% t(Fk) + diag(e_k^2 + q_floor, D) # Predicted covariance
      P_pred[k + 1L,,] <- Pk_pred

      # Update step
      innov <- U[k + 1L, ] - u_pred[k + 1L, ]
      Sk <- Pk_pred + R_obs
      Kk <- Pk_pred %*% solve(Sk)

      u_filt[k + 1L, ] <- u_pred[k + 1L, ] + Kk %*% innov
      P_filt[k + 1L,,] <- (I_D - Kk) %*% Pk_pred %*% t(I_D - Kk) + Kk %*% R_obs %*% t(Kk)
    }

    u_smooth <- matrix(0, mp1, D)
    u_smooth[mp1, ] <- u_filt[mp1, ]
    P_smooth <- if (want_cov) array(0, c(mp1, D, D)) else NULL
    if (want_cov) P_smooth[mp1,,] <- P_filt[mp1,,]

    # RTS Smoother
    for (k in seq(mp1 - 1L, 1L)) {
      Pk <- P_filt[k,,]
      Pp <- P_pred[k + 1L,,]
      Fk <- F_store[k,,]

      Ck <- Pk %*% t(Fk) %*% solve(Pp + 1e-10 * I_D)
      u_smooth[k, ] <- u_filt[k, ] + Ck %*% (u_smooth[k + 1L, ] - u_pred[k + 1L, ])
      if (want_cov){
        P_smooth[k,,] <- Pk + Ck %*% (P_smooth[k + 1L,,] - Pp) %*% t(Ck)
      }
    }

    list(U_star = u_smooth, P_smooth = P_smooth, u_filt = u_filt, P_filt = P_filt, u_pred = u_pred, P_pred = P_pred)
  }

  # Conditional p̂ pass: smoothed mean + conditional posterior covariance
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
    sens <- array(0, c(mp1, D, J))  # ∂u_k/∂p̂_j
    for (j in seq_len(J)) {
      hj <- 1e-5 * max(1, abs(p[j]))
      pp <- p
      pp[j] <- pp[j] + hj
      pm <- p
      pm[j] <- pm[j] - hj
      sens[, , j] <- (erts_pass(pp, FALSE)$U_star - erts_pass(pm, FALSE)$U_star) / (2 * hj)
    }
    P_param <- array(0, c(mp1, D, D))
    for (k in seq_len(mp1)) {
      Sk <- matrix(sens[k, , ], D, J)
      P_param[k,,]  <- Sk %*% param_cov %*% t(Sk)
      # Law of total variance to explain the total uncertainty in the estimate
      # total var = unexplained + explained
      P_smooth[k,,] <- P_cond[k,,] + P_param[k,,]
    }
  }

  list(
    U_star = u_smooth,
    P_smooth = P_smooth,
    P_smooth_cond  = P_cond,
    P_smooth_param = P_param,
    u_filt = base$u_filt,
    P_filt = base$P_filt,
    u_pred = base$u_pred,
    P_pred = base$P_pred
  )
}