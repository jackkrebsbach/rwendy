# Joint state + parameter estimation in a test-function basis.
#
# Loss: L(c, p) = lambda * || Phi.dot * uhat + Phi * F(uhat, p) ||^2     (weak residual)
#               +         || uhat - U_obs ||^2_{Sigma_eps^{-1}}          (data fit)
#               + rho * || uhat^(q) [- ODE^(q)(uhat,p)] ||^2             (optional deriv penalty)
# with the state reconstructed from a basis expansion uhat = B^T c. The basis
# rank is selected adaptively (data modes above the noise floor) by default.
#
# Two test-function systems are used (both on the same grid):
#   * the weak residual (Phi = ctx$V, Phi.dot = ctx$V_prime) and the IRLS warm
#     start use MSG test functions with no boundary layer (forced in solveWendy);
#   * the state basis B = t(svd(V_ssl)$v) are the right singular vectors of a
#     separately built SSL + boundary-layer matrix V_ssl. B has orthonormal rows
#     (B B^T = I_Kb), so B^T B projects onto V_ssl's row space.
#
# Layout (same column-major convention as weak_residual.R). Let K = nrow(V) be
# the number of test functions (weak-residual rows) and Kb = nrow(B) =
# min(K, mp1) the basis dimension (B drops null singular directions, so Kb can
# be smaller than K, e.g. for SSL with boundary layers where K > mp1):
#   vec(c) : d-blocks of length Kb, so c is (Kb x D), col c[ , d].
#   uhat   : (mp1 x D), uhat = t(B) %*% c.
#   weak residual r_weak = vec_{(K,D)}( V %*% F(uhat,p) + Vp %*% uhat ) in R^{K*D}.
#   data residual r_data = vec_{(mp1,D)}( (uhat - U_obs) / col sigma ) in R^{mp1*D}.
# The optimisation variable is theta = c(vec(c), p) in R^{Kb*D + J}.

# State basis B and matching derivative bases, from one SVD of the SSL + BL
# test-function matrix V_ssl = U_s diag(d) W^T (built on ctx$U_sys / ctx$tt_sys;
# test_fun_type and boundary layer forced on regardless of the MSG residual).
# The right singular vectors v_k = (1/d_k) V_ssl^T u_k are smooth grid functions;
# their analytic q-th time-derivatives are the SAME recombination of the q-th
# derivative test functions, v_k^(q) = (1/d_k) V^(q)_ssl^T u_k (mirrors the MSG
# V_prime build). So uhat = B^T c, uhat' = Bp^T c, uhat'' = Bpp^T c, with
#   B   = W^T                  = t(svd(V_ssl)$v)       (Kb x mp1, orthonormal rows)
#   Bp  = diag(1/d) U_s^T V'_ssl                       (Kb x mp1)
#   Bpp = diag(1/d) U_s^T V''_ssl                      (Kb x mp1).
.joint_bases <- function(ctx) {
  ssl_ctl <- modifyList(ctx$control,
                        list(test_fun_type = "SSL", include_boundary_layer = TRUE,
                             p = ctx$control$joint_basis_p %||% 16))
  tf  <- build_full_test_function_matrices_ssl(ctx$U_sys, as.vector(ctx$tt_sys), ssl_ctl)
  sv  <- svd(tf$V)
  inv_d_Ut <- (1 / sv$d) * t(sv$u)                # row k scaled by 1/d_k
  B   <- t(sv$v)                                  # Kb x mp1
  Bp  <- inv_d_Ut %*% tf$V_prime                  # Kb x mp1 (first derivative)
  Bpp <- inv_d_Ut %*% tf$V_pp                     # Kb x mp1 (second derivative)
  list(B = B, Bp = Bp, Bpp = Bpp)
}

# Data-driven rank for the state basis. Project the data onto the full SSL+BL
# basis (orthonormal rows, so each coefficient carries noise variance sigma^2);
# count the modes whose energy sum_d c_kd^2 exceeds tau times the noise floor
# sum_d sigma_d^2, and keep that many of the LEADING (smoothest) modes. The
# energy threshold sets the low-pass cutoff adaptively: at large n the trailing
# noise-only modes are excluded (killing the null-space overfitting), while
# under-sampled / complex signals retain more modes.
.joint_select_rank <- function(Bfull, U, sig, tau = 1) {
  D <- ncol(U)
  if (length(sig) == 1L) sig <- rep(sig, D)
  e <- rowSums((Bfull %*% U)^2)
  max(1L, sum(e > tau * sum(sig^2)))
}

# d r_weak / d vec(uhat)  in R^{K*D x mp1*D}, evaluated at (uhat, p).
# This is build_L() with sigma set to 1: the diagonal (a==a) blocks carry the
# Phi.dot term d(Vp uhat)/d uhat = Vp, and every (a,b) block carries the
# reaction term V %*% diag(d f_a / d u_b).
.joint_state_jac <- function(uhat, p, tt, V, Vp, J_u, K, D, mp1, J) {
  tU    <- t(uhat)
  ttt   <- matrix(tt, nrow = 1L)
  p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
  input <- rbind(p_mat, tU, ttt)
  J_F   <- array(J_u(input), c(mp1, D, D))   # J_F[m, a, b] = d f_a / d u_b @ m
  out   <- matrix(0, K * D, mp1 * D)
  for (a in seq_len(D)) for (b in seq_len(D)) {
    rs  <- (a - 1L) * K   + seq_len(K)
    cs  <- (b - 1L) * mp1 + seq_len(mp1)
    blk <- sweep(V, 2, J_F[, a, b], "*")
    if (a == b) blk <- blk + Vp
    out[rs, cs] <- blk
  }
  out
}

# Build the stacked residual + analytic Jacobian closures for nls.lm.
# Bbig  = kron(I_D, t(B))            (mp1*D x Kb*D), maps vec(c) -> vec(uhat).
# B2c   = kron(diag(1/sigma), t(B))  (mp1*D x Kb*D), d r_data  / d vec(c) (const).
# Bdbig = kron(I_D, t(Bd))           (mp1*D x Kb*D), d vec(uhat^(q))/d vec(c) (const),
#   where Bd is the order-q (1 or 2) derivative basis for the smoothness penalty.
.joint_make_fns <- function(ctx, B, Bd, ode_target = NULL) {
  V    <- ctx$V
  Vp   <- ctx$V_prime
  f_   <- ctx$f_
  J_u  <- ctx$J_u
  J_p  <- ctx$J_p
  J    <- ctx$J
  tt   <- as.vector(ctx$tt_sys)
  Uobs <- ctx$U_sys
  D    <- ncol(Uobs)
  mp1  <- nrow(Uobs)
  K    <- nrow(V)   # weak-residual rows
  Kb   <- nrow(B)   # basis dimension

  sig <- ctx$sig
  if (length(sig) == 1L) sig <- rep(sig, D)

  lambda <- ctx$control$joint_lambda %||% 1
  sl     <- sqrt(lambda)

  rho <- ctx$control$joint_deriv_pen %||% 0   # weight on ||uhat^(q)||^2
  sr  <- sqrt(rho)

  Bt   <- t(B)                                    # mp1 x Kb
  Bbig <- kronecker(diag(D), Bt)                  # mp1*D x Kb*D
  B2c  <- kronecker(diag(1 / sig, nrow = D), Bt)  # mp1*D x Kb*D  (Sigma^{-1/2} Bbig)
  Bdt  <- t(Bd)                                   # mp1 x Kb
  Bdbig <- if (sr > 0 && is.null(ode_target)) kronecker(diag(D), Bdt) else NULL
  KbD  <- Kb * D    # coefficient block length
  KD   <- K * D     # weak-residual block length

  unpack <- function(theta) {
    cmat <- matrix(theta[seq_len(KbD)], Kb, D)
    p    <- theta[KbD + seq_len(J)]
    uhat <- crossprod(B, cmat)                  # t(B) %*% cmat = mp1 x D
    list(cmat = cmat, p = p, uhat = uhat)
  }

  make_input <- function(uhat, p) {
    rbind(matrix(rep(p, mp1), ncol = mp1, nrow = J), t(uhat), matrix(tt, nrow = 1L))
  }
  weak_F <- function(uhat, p) f_(make_input(uhat, p))     # mp1 x D

  # Derivative-penalty residual (length mp1*D, unscaled). Raw smoothness uses
  # target 0 -> ||uhat^(q)||^2; ODE-consistency uses target = ode_target(uhat,p)
  # (= f for q=1, dF/dt for q=2) -> ||uhat^(q) - ODE^(q)(uhat,p)||^2, which
  # penalises only the part of the derivative the ODE does not explain.
  pen_resid <- function(st) {
    d_q <- crossprod(Bd, st$cmat)                         # uhat^(q), mp1 x D
    if (is.null(ode_target)) as.vector(d_q)
    else as.vector(d_q - ode_target(make_input(st$uhat, st$p)))
  }

  residual <- function(theta) {
    st     <- unpack(theta)
    r_weak <- as.vector(V %*% weak_F(st$uhat, st$p) + Vp %*% st$uhat)  # K*D
    r_data <- as.vector(sweep(st$uhat - Uobs, 2, sig, "/"))           # mp1*D
    out    <- c(sl * r_weak, r_data)
    if (sr > 0) out <- c(out, sr * pen_resid(st))                     # mp1*D deriv penalty
    out
  }

  jacobian <- function(theta) {
    st  <- unpack(theta)
    p   <- st$p
    uhat <- st$uhat

    # d r_weak / d vec(c) = L1(uhat,p) %*% Bbig  ;  d r_weak / d p = V %*% J_p
    L1     <- .joint_state_jac(uhat, p, tt, V, Vp, J_u, K, D, mp1, J)   # K*D x mp1*D
    dweak_dc <- L1 %*% Bbig                                             # K*D x Kb*D

    p_mat   <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input   <- rbind(p_mat, t(uhat), matrix(tt, nrow = 1L))
    Jp_eval <- J_p(input)                                              # mp1 x (D*J)
    dweak_dp <- matrix(V %*% Jp_eval, KD, J)                           # K*D x J

    top    <- cbind(sl * dweak_dc, sl * dweak_dp)                      # K*D   x (Kb*D+J)
    bottom <- cbind(B2c, matrix(0, mp1 * D, J))                        # mp1*D x (Kb*D+J)
    Jmat   <- rbind(top, bottom)
    if (sr > 0) {
      if (is.null(ode_target)) {
        deriv <- cbind(sr * Bdbig, matrix(0, mp1 * D, J))             # analytic raw block
      } else {
        # ODE-consistent target depends nonlinearly on (c, p); prototype with a
        # numeric Jacobian on the (cheap) penalty block only.
        deriv <- sr * numDeriv::jacobian(function(th) pen_resid(unpack(th)),
                                         theta, method = "simple")
      }
      Jmat <- rbind(Jmat, deriv)
    }
    Jmat
  }

  list(residual = residual, jacobian = jacobian, unpack = unpack,
       K = K, D = D, J = J)
}

# JOINT method runner. Warm-starts p0 with WENDy-IRLS when not supplied and
# c0 with the projection of the noisy data into the basis, B %*% U_obs.
.run_joint <- function(ctx, p0, G_cont, b_cont) {
  if (is.null(p0)) {
    p0_for_irls <- if (ctx$lip) NULL else .ee_init_p0(ctx)
    p0 <- .run_irls(ctx, p0_for_irls, G_cont, b_cont)$p
  }

  # Resolve an "auto" weak-term weight from the weak-residual variance scale:
  # lambda = B / mean(diag(S(p0))). mean(diag(S)) ~ sigma^2, so lambda ~ 1/sigma^2
  # -- more regularisation at low noise -- and the dimensionless B ~ 300 lands in
  # the state-reconstruction plateau across systems and noise levels.
  lambda <- ctx$control$joint_lambda
  if (is.character(lambda) && identical(lambda, "auto")) {
    Bsc <- ctx$control$joint_lambda_scale %||% 300
    mdS <- tryCatch(mean(diag(ctx$system$S(p0))), error = function(e) NA_real_)
    lambda <- if (is.finite(mdS) && mdS > 0) Bsc / mdS else 1e3
  }
  ctx$control$joint_lambda <- lambda   # .joint_make_fns reads the numeric value

  # Build the SSL + BL state basis B and its derivative basis Bp, then select the
  # rank. A full-rank basis (Kb = mp1) leaves a large null space of the weak
  # operator that the data fills with noise -- the overfitting that worsens as
  # mp1 grows. "auto" keeps only the data modes above the noise floor (adaptive);
  # "full" keeps all; an integer keeps the top-r smooth modes.
  bb     <- .joint_bases(ctx)
  Bfull  <- bb$B
  rk     <- ctx$control$joint_basis_rank
  modes  <- if (is.null(rk) || identical(rk, "full")) {
    seq_len(nrow(Bfull))
  } else if (identical(rk, "auto")) {
    r <- .joint_select_rank(Bfull, ctx$U_sys, ctx$sig,
                            tau = ctx$control$joint_rank_tau %||% 1)
    seq_len(min(r, nrow(Bfull)))
  } else {
    seq_len(min(as.integer(rk), nrow(Bfull)))
  }
  B     <- bb$B[modes, , drop = FALSE]
  # Derivative-penalty basis: order 2 (curvature, default) or order 1.
  order <- ctx$control$joint_deriv_order %||% 2
  Bd    <- if (order == 1) bb$Bp[modes, , drop = FALSE] else bb$Bpp[modes, , drop = FALSE]
  # ODE-consistency target: NULL = raw ||uhat^(q)||^2; else penalise the
  # deviation from the ODE's q-th derivative (f for q=1, dF/dt for q=2), which
  # does not fight the trajectory's genuine structure.
  ode_target <- if (!isTRUE(ctx$control$joint_deriv_ode %||% FALSE)) NULL
                else if (order == 1) ctx$f_ else ctx$dF_dt_
  fns <- .joint_make_fns(ctx, B, Bd, ode_target)

  Kb  <- nrow(B)                                  # basis dimension
  D   <- fns$D
  J   <- fns$J
  KbD <- Kb * D

  # c0 = projection of the data into the basis. Unobserved columns have no data,
  # so initialise them at zero; the state warm start / weak residual fills them.
  U_init <- ctx$U_sys
  obs    <- ctx$control$joint_observed %||% seq_len(ncol(U_init))
  if (length(obs) < ncol(U_init)) U_init[, setdiff(seq_len(ncol(U_init)), obs)] <- 0
  c0     <- B %*% U_init                          # Kb x D
  theta0 <- c(as.vector(c0), p0)

  lm_ctl <- minpack.lm::nls.lm.control(maxiter = ctx$control$max_iterates,
                                       maxfev  = 100 * (length(theta0) + 1))

  # Optional state-only warm start: hold p = p0 and optimise c alone, so the
  # coefficients settle onto the ODE manifold (at p0) before the joint solve.
  # This re-initialises the state from the noisy-data projection to an
  # ODE-consistent one, which can land the joint optimiser in a better basin.
  ws_its <- 0L
  if (isTRUE(ctx$control$joint_state_warmstart %||% FALSE)) {
    res_c <- function(cv) fns$residual(c(cv, p0))
    jac_c <- function(cv) fns$jacobian(c(cv, p0))[, seq_len(KbD), drop = FALSE]
    fit_c <- nls.lm(theta0[seq_len(KbD)], fn = res_c, jac = jac_c, control = lm_ctl)
    theta0 <- c(fit_c$par, p0)
    ws_its <- fit_c$niter
  }

  fit       <- nls.lm(theta0, fn = fns$residual, jac = fns$jacobian, control = lm_ctl)
  theta_hat <- fit$par
  info      <- fit$info
  dev       <- fit$deviance
  its       <- fit$niter

  st <- fns$unpack(theta_hat)

  list(
    p           = st$p,
    c           = st$cmat,
    uhat        = st$uhat,
    converged   = info %in% c(1L, 2L, 3L),
    iterations  = its,
    warmstart_iters = ws_its,
    lambda      = lambda,
    basis_rank  = Kb,
    n_coef      = KbD,
    info        = info,
    deviance    = dev
  )
}
