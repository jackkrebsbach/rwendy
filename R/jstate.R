# JSTATE: joint state + parameter estimation, optimising uhat directly (no
# basis expansion) with a soft-rank-cap projection penalty for smoothing.
#
## Loss: L(uhat, p) = lambda * || V_msg %*% F(uhat, p) + Vp_msg %*% uhat ||^2  (weak residual)
##                  + mu    * || (uhat - U_obs) / sigma ||^2                    (data fit)
##                  + rho   * || (I - B^T B) %*% uhat ||^2                      (out-of-basis penalty)
##
# The MSG matrices come from the ctx-level build, same as IRLS. The smoothness
# basis B is built from a *separate* SSL + boundary-layer test-function matrix
# V_ssl (same as JOINT's .joint_bases): take the SVD, keep the top-K right
# singular vectors (the K smoothest modes well-represented by the test
# functions), so B = t(svd(V_ssl)$v)[1:K, ] has orthonormal rows (B B^T = I_K).
# Boundary-layer augmentation is REQUIRED for non-zero-boundary trajectories
# (e.g. logistic at t=0 with u != 0): without BL, the compactly supported psi
# vanish at the endpoints and span(B^T) cannot represent the true trajectory,
# making the weak residual at proj_B U_true large enough that the trivial
# weak-residual solution (uhat = 0, p = 0, valid when F(0, p) = 0) becomes the
# global optimum.
#
# (I - B^T B) projects onto the orthogonal complement of span(B^T) -- the
# "high-frequency" directions the smooth test functions cannot represent. The
# penalty drives uhat's energy in those directions to zero. As rho -> infty,
# uhat is forced into span(B^T) (the JOINT rank-cap subspace, but as a soft
# constraint). As rho -> 0, the penalty disappears.
#
# Implementation: instead of the mp1-row identity form r_pen = (I - B^T B) uhat
# (rank only mp1-K), we use the orthonormal complement basis B_perp from the
# FULL SVD of V_ssl: B_perp has (mp1 - K) rows, B_perp B_perp^T = I, and
# ||(I - B^T B) uhat||^2 = ||B_perp uhat||^2. (mp1-K)*D residual rows, smaller
# Jacobian, same objective.
# Optimization variable: theta = c(vec(uhat), p) in R^{mp1*D + J}, column-major.

# Build the orthonormal complement basis B_perp (rows span the orthogonal
# complement of the top-K smooth singular vectors of an SSL + BL test-function
# matrix V_ssl). The BL augmentation is required so span(B^T) can represent
# non-zero-boundary trajectories -- without it the basis lops off endpoint
# values and the trivial weak-residual solution (uhat=0, p=0) wins globally.
# Returns list(B_perp, K_keep, K_total).
.jstate_basis_perp <- function(ctx) {
  ssl_ctl <- modifyList(ctx$control,
                        list(test_fun_type = "SSL",
                             include_boundary_layer = TRUE,
                             p = ctx$control$jstate_basis_p %||% 16))
  radius <- ctx$control$jstate_basis_radius
  if (!is.null(radius) && !identical(radius, "auto")) {
    ssl_ctl$fixed_radius <- as.integer(radius)
  } else {
    ssl_ctl$fixed_radius <- NULL                # SSL picks rc from data
  }
  V_ssl <- build_full_test_function_matrices_ssl(ctx$U_sys,
                                                 as.vector(ctx$tt_sys),
                                                 ssl_ctl)$V                # K_tot x mp1
  mp1 <- ncol(V_ssl)
  # FULL right-singular-vector matrix (mp1 x mp1): the first K_eff columns
  # (sorted by descending singular value) span the *row space* of V_ssl --
  # these are the smooth modes the test functions can actually represent. The
  # remaining mp1 - K_eff columns span V_ssl's NULL space and are arbitrary
  # orthonormal vectors with no smoothness guarantee (often wiggly). It is
  # essential to cap K at K_eff: putting null-space directions in B_keep lets
  # the optimizer fit noise spikes through those wiggly modes (visible as
  # local jumps in uhat), defeating the smoothing penalty.
  sv <- svd(V_ssl, nu = 0L, nv = mp1)
  W  <- sv$v                                                                # mp1 x mp1
  # Well-conditioned rank: drop modes whose singular value is below
  # sv_max / jstate_basis_cond_max. Modes below the cut are near-null-space
  # directions of V_ssl (very small sv) -- they're orthonormal only by accident
  # of SVD ordering and tend to be wiggly. Mirrors the MSG k_max /
  # max_test_fun_condition_number machinery (R/test_functions.R:380). The
  # default cond_max = 1e6 falls in the gap below the smooth modes (sv$d[1] ~
  # O(1), the smooth tail decays to ~1e-6) for logistic/LV/Lorenz at mp1
  # 100-1000. True log(sv) elbow / cumulative-energy cutoffs (K ~ 16-35 here)
  # are too aggressive -- they cut into the modes JSTATE actually needs.
  cond_max <- ctx$control$jstate_basis_cond_max %||% 1e6
  k_max    <- ctx$control$jstate_basis_rank_max                   # NULL = no cap
  K_eff    <- sum(sv$d > sv$d[1] / cond_max)
  K_auto   <- min(K_eff, k_max %||% .Machine$integer.max, mp1 - 1L)

  K_user <- ctx$control$jstate_basis_rank %||% "auto"
  K_keep <- if (identical(K_user, "auto")) K_auto
            else max(1L, min(as.integer(K_user), K_eff, mp1 - 1L))
  B_perp <- t(W[, (K_keep + 1L):mp1, drop = FALSE])                         # (mp1 - K) x mp1
  list(B_perp = B_perp, K_keep = K_keep, K_total = mp1, K_eff = K_eff)
}

# Residual + analytic Jacobian closures for nls.lm. theta = c(vec(uhat), p).
# Reuses .joint_state_jac for d r_weak / d vec(uhat).
.jstate_make_fns <- function(ctx, B_perp) {
  V    <- ctx$V                          # MSG weak residual (K x mp1)
  Vp   <- ctx$V_prime
  f_   <- ctx$f_
  J_u  <- ctx$J_u
  J_p  <- ctx$J_p
  J    <- ctx$J
  tt   <- as.vector(ctx$tt_sys)
  Uobs <- ctx$U_sys
  D    <- ncol(Uobs)
  mp1  <- nrow(Uobs)
  K    <- nrow(V)
  Kperp <- nrow(B_perp)

  sig <- ctx$sig
  if (length(sig) == 1L) sig <- rep(sig, D)

  lambda <- ctx$control$jstate_lambda %||% 1
  sl     <- sqrt(lambda)
  rho    <- ctx$control$jstate_pen %||% 0
  sr     <- sqrt(rho)
  mu     <- ctx$control$jstate_data_pen %||% 1     # weight on data-fit term
  sm     <- sqrt(mu)

  obs <- ctx$control$jstate_observed %||% seq_len(D)

  KD     <- K     * D
  mD     <- mp1   * D
  KperpD <- Kperp * D

  # Per-component data-fit weight: sqrt(mu)/sig_d for observed columns, 0 for
  # unobserved (those rows then contribute nothing -- the column is shaped
  # only by the weak residual + penalty). mu is the user-tunable weight on the
  # data fit; sig_d is the noise-variance whitening already baked in.
  data_w_col <- ifelse(seq_len(D) %in% obs, sm / sig, 0)
  data_scale <- rep(data_w_col, each = mp1)            # length mD

  # Penalty Jacobian d r_pen / d vec(uhat) = block-diag(B_perp), built once.
  pen_jac_u <- if (sr > 0) kronecker(diag(D), B_perp) else NULL

  unpack <- function(theta) {
    uhat <- matrix(theta[seq_len(mD)], mp1, D)
    p    <- theta[mD + seq_len(J)]
    list(uhat = uhat, p = p)
  }

  make_input <- function(uhat, p) {
    rbind(matrix(rep(p, mp1), ncol = mp1, nrow = J), t(uhat), matrix(tt, nrow = 1L))
  }
  weak_F <- function(uhat, p) f_(make_input(uhat, p))

  residual <- function(theta) {
    st     <- unpack(theta)
    r_weak <- as.vector(V %*% weak_F(st$uhat, st$p) + Vp %*% st$uhat)   # K*D
    r_data <- data_scale * as.vector(st$uhat - Uobs)                    # mp1*D
    out    <- c(sl * r_weak, r_data)
    if (sr > 0) out <- c(out, sr * as.vector(B_perp %*% st$uhat))       # Kperp*D
    out
  }

  jacobian <- function(theta) {
    st   <- unpack(theta)
    uhat <- st$uhat
    p    <- st$p

    L1     <- .joint_state_jac(uhat, p, tt, V, Vp, J_u, K, D, mp1, J)   # K*D x mp1*D
    p_mat   <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input   <- rbind(p_mat, t(uhat), matrix(tt, nrow = 1L))
    Jp_eval <- J_p(input)                                               # mp1 x (D*J)
    dweak_dp <- matrix(V %*% Jp_eval, KD, J)                            # K*D x J

    top  <- cbind(sl * L1, sl * dweak_dp)                               # K*D    x (mD + J)
    mid  <- cbind(diag(data_scale, nrow = mD), matrix(0, mD, J))        # mp1*D  x (mD + J)
    Jmat <- rbind(top, mid)
    if (sr > 0) {
      bot  <- cbind(sr * pen_jac_u, matrix(0, KperpD, J))               # Kperp*D x (mD + J)
      Jmat <- rbind(Jmat, bot)
    }
    Jmat
  }

  list(residual = residual, jacobian = jacobian, unpack = unpack,
       K = K, Kperp = Kperp, D = D, J = J, mD = mD)
}

# JSTATE runner. Warm-starts p0 with WENDy-IRLS when not supplied and
# uhat0 = U_obs. Optional state-only warm start (uhat alone with p held at p0).
.run_jstate <- function(ctx, p0, G_cont, b_cont) {
  if (is.null(p0)) {
    p0_for_irls <- if (ctx$lip) NULL else .ee_init_p0(ctx)
    p0 <- .run_irls(ctx, p0_for_irls, G_cont, b_cont)$p
  }

  # Auto-tune the weak-term weight: lambda = scale / mean(diag(S(p0))). Mirrors
  # the JOINT heuristic -- mean(diag(S)) ~ sigma^2, so lambda ~ 1/sigma^2.
  lambda <- ctx$control$jstate_lambda
  if (is.character(lambda) && identical(lambda, "auto")) {
    Bsc <- ctx$control$jstate_lambda_scale %||% 300
    mdS <- tryCatch(mean(diag(ctx$system$S(p0))), error = function(e) NA_real_)
    lambda <- if (is.finite(mdS) && mdS > 0) Bsc / mdS else 1e3
  }
  ctx$control$jstate_lambda <- lambda

  bb       <- .jstate_basis_perp(ctx)
  fns      <- .jstate_make_fns(ctx, bb$B_perp)

  D  <- fns$D
  J  <- fns$J
  mD <- fns$mD

  # uhat0 = smooth projection of U_obs onto span(B^T) when the penalty is on,
  # else raw U_obs. Init matters: for ODEs with F(0,p) = 0 (logistic, LV, ...)
  # the weak residual has a TRIVIAL solution at uhat = 0, and at large rho the
  # optimizer can fall into that basin from a U_obs init (the penalty's pull
  # toward the smooth subspace combines with lambda*weak's pull toward zero).
  # Starting in span(B^T) keeps the penalty term ~0 at init and the data term's
  # gradient (toward proj_B U_obs) keeps uhat away from the trivial basin.
  # Unobserved columns are zeroed (weak residual + penalty fill them).
  U_init <- ctx$U_sys
  obs    <- ctx$control$jstate_observed %||% seq_len(ncol(U_init))
  if (length(obs) < ncol(U_init)) U_init[, setdiff(seq_len(ncol(U_init)), obs)] <- 0
  if ((ctx$control$jstate_pen %||% 0) > 0) {
    U_init <- U_init - crossprod(bb$B_perp, bb$B_perp %*% U_init)            # proj onto span(B^T)
  }
  theta0 <- c(as.vector(U_init), p0)

  lm_ctl <- minpack.lm::nls.lm.control(maxiter = ctx$control$max_iterates,
                                       maxfev  = 100 * (length(theta0) + 1))

  # Optional state-only warm start: hold p = p0 and optimise uhat alone. Lands
  # the joint solve in an ODE-consistent basin at p0; off by default since
  # uhat0 = U_obs is usually fine.
  ws_its <- 0L
  if (isTRUE(ctx$control$jstate_state_warmstart %||% FALSE)) {
    res_u <- function(uv) fns$residual(c(uv, p0))
    jac_u <- function(uv) fns$jacobian(c(uv, p0))[, seq_len(mD), drop = FALSE]
    fit_u <- nls.lm(theta0[seq_len(mD)], fn = res_u, jac = jac_u, control = lm_ctl)
    theta0 <- c(fit_u$par, p0)
    ws_its <- fit_u$niter
  }

  fit       <- nls.lm(theta0, fn = fns$residual, jac = fns$jacobian, control = lm_ctl)
  theta_hat <- fit$par
  st        <- fns$unpack(theta_hat)

  list(
    p           = st$p,
    uhat        = st$uhat,
    converged   = fit$info %in% c(1L, 2L, 3L),
    iterations  = fit$niter,
    warmstart_iters = ws_its,
    lambda      = lambda,
    rho         = ctx$control$jstate_pen %||% 0,
    basis_rank  = bb$K_keep,
    basis_rank_eff = bb$K_eff,                    # well-conditioned rank of V_ssl
    n_perp_rows = fns$Kperp,
    info        = fit$info,
    deviance    = fit$deviance
  )
}