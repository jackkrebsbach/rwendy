# Solve (A + lam*I) x = b with ridge bumping; fall back to pseudoinverse if all
# bumps fail. Shared by the Gauss-Newton variants below.
solve_pd <- function(A, b, lam) {
  n <- nrow(A)
  for (i in seq_len(10L)) {
    result <- tryCatch(solve(A + lam * diag(n), b), error = function(e) NULL)
    if (!is.null(result)) return(result)
    lam <- lam * 10
  }
  MASS::ginv(A) %*% b
}

# Broyden-update Gauss-Newton loop with optional adaptive Tikhonov penalty on a
# subset of theta. target_fn(theta) returns a residual vector to be driven to 0.
#
# penalty_mask: length-length(theta) 0/1 vector marking entries that receive the
#   (theta_i - theta_ref_i)^2 penalty. NULL means all ones.
# theta_ref: penalty pulls theta toward this; defaults to theta0.
# alpha = 0 with sigma given means "adaptive": alpha_eff is recomputed each
#   iteration from the current residual norm so the penalty decays naturally.
# Returns the final theta vector.
gauss_newton_iterate <- function(target_fn, theta0,
                                  penalty_mask = NULL, theta_ref = NULL,
                                  max_iter = 100, tol = 1e-8, lambda = 1e-8,
                                  alpha = 0, sigma = NULL) {
  n        <- length(theta0)
  if (is.null(penalty_mask)) penalty_mask <- rep(1, n)
  if (is.null(theta_ref))    theta_ref    <- theta0
  n_mask   <- sum(penalty_mask)
  sig_mean <- if (!is.null(sigma)) mean(as.numeric(sigma)) else NULL

  theta_w <- theta0
  G_w     <- target_fn(theta_w)
  J_w     <- jacobian(target_fn, theta_w)

  for (iter in seq_len(max_iter)) {
    # Recompute alpha each iteration from the current residual so the penalty
    # decays naturally as we converge (large early, small at solution).
    alpha_eff <- if (alpha == 0 && !is.null(sig_mean)) {
      sum(G_w^2) / (n_mask * sig_mean^2)
    } else alpha

    grad_pen <- alpha_eff * penalty_mask * (theta_w - theta_ref)
    R_reg    <- diag(alpha_eff * penalty_mask, nrow = n)

    JtJ   <- t(J_w) %*% J_w + R_reg
    rhs   <- t(J_w) %*% G_w + grad_pen
    delta <- solve_pd(JtJ, rhs, lambda)
    theta_new <- theta_w - delta

    if (anyNA(delta) || !all(is.finite(delta))) break

    G_new <- target_fn(theta_new)
    s <- -delta; y <- G_new - G_w
    ss <- as.numeric(t(s) %*% s)
    if (ss > 1e-14)
      J_w <- J_w + ((y - J_w %*% s) %*% t(s)) / ss

    theta_w <- theta_new
    G_w     <- G_new

    if (max(abs(delta)) < tol) break
  }
  theta_w
}

# Gauss-Newton joint state and parameter estimate
gn <- function(p, U, tt, f_, V, V_prime, max_iter = 100, tol = 1e-10,
               lambda = 1e-8, alpha = 0, sigma = NULL){
  D  <- ncol(U)
  Vp <- as.array(V_prime$contiguous())
  V  <- as.array(V$contiguous())

  u0_vec <- as.vector(U)
  nu <- length(u0_vec)
  np <- length(p)

  F_ <- function(U, p){
    mp1   <- nrow(U)
    J     <- length(p)
    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    f_(rbind(p_mat, tU, ttt))
  }

  joint_target <- function(theta) {
    U_cur <- matrix(theta[seq_len(nu)], ncol = D)
    p_cur <- theta[nu + seq_len(np)]
    as.vector(Vp %*% U_cur + V %*% F_(U_cur, p_cur))
  }

  theta0  <- c(u0_vec, p)
  mask    <- c(rep(1, nu), rep(0, np))
  theta_w <- gauss_newton_iterate(joint_target, theta0,
                                   penalty_mask = mask, theta_ref = theta0,
                                   max_iter = max_iter, tol = tol,
                                   lambda = lambda, alpha = alpha, sigma = sigma)

  list(Uhat = matrix(theta_w[seq_len(nu)], ncol = D),
       p    = theta_w[nu + seq_len(np)])
}

# GLS-Gauss-Newton joint state and parameter estimate.
# Variant of gn(): weights the normal equations by S(p)^{-1}, evaluated once
# at the initial p (held fixed across iterations as a preconditioner). No
# Tikhonov penalty on the state, no Broyden safety floor — intentionally
# minimal so it's easy to compare against gn().
gls_gn <- function(p, U, tt, f_, V, V_prime, S,
                   max_iter = 100, tol = 1e-8, lambda = 1e-8){
  D  <- ncol(U)
  Vp <- as.array(V_prime$contiguous())
  V  <- as.array(V$contiguous())

  u0_vec <- as.vector(U)
  nu <- length(u0_vec)
  np <- length(p)

  F_ <- function(U, p){
    mp1   <- nrow(U)
    J     <- length(p)
    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    f_(rbind(p_mat, tU, ttt))
  }

  joint_target <- function(theta) {
    U_cur <- matrix(theta[seq_len(nu)], ncol = D)
    p_cur <- theta[nu + seq_len(np)]
    as.vector(-Vp %*% U_cur - V %*% F_(U_cur, p_cur))
  }

  theta_w <- c(u0_vec, p)
  G_w     <- joint_target(theta_w)
  J_w     <- jacobian(joint_target, theta_w)

  S_inv <- solve(as.array(S(p)$contiguous()))

  for (iter in seq_len(max_iter)) {
    JtSiJ <- t(J_w) %*% S_inv %*% J_w
    delta <- solve_pd(JtSiJ, t(J_w) %*% S_inv %*% G_w, lambda)
    theta_new <- theta_w - delta

    G_new <- joint_target(theta_new)
    s <- -delta; y <- G_new - G_w
    J_w <- J_w + ((y - J_w %*% s) %*% t(s)) / as.numeric(t(s) %*% s)

    theta_w <- theta_new
    G_w     <- G_new

    if (max(abs(delta)) < tol) break
  }

  list(Uhat = matrix(theta_w[seq_len(nu)], ncol = D),
       p    = theta_w[nu + seq_len(np)])
}

# Gauss-Newton boundary refinement: hold interior U and p fixed, optimize only
# the first and last r_c rows of U against an SSL+BL system built lazily.
# Returns list(Uhat, p, g, build_compute_b) — the latter two reusable by gn_bl.
gn_boundary <- function(p, U, tt,
                        f_, J_u, J_uu, J_up, J_p, J_pp, J_upp,
                        dF_dt_ = NULL, d2F_dt2_ = NULL, d3F_dt3_ = NULL,
                        J_num, lip, sig, control,
                        U_ref = NULL, fixed_radius = NULL,
                        max_iter = 100, tol = 1e-10, lambda = 1e-8,
                        alpha = 0, sigma = NULL) {
  mp1 <- nrow(U)
  D   <- ncol(U)
  tt  <- as.vector(tt)
  device <- if (is.character(control$device)) torch::torch_device(control$device) else control$device

  # Build SSL+BL system lazily.
  # Use U_ref (original noisy data) for r_c selection — noisy data yields a
  # larger change-point radius than the smooth stage-1 estimate, giving the
  # boundary window enough coverage.
  ctrl_ssl_bl <- modifyList(control, list(
    test_fun_type          = "SSL",
    include_boundary_layer = TRUE,
    fixed_radius           = fixed_radius %||% control$fixed_radius
  ))
  wendy_data_bl <- list(U = U_ref %||% U, tt = tt, var = NULL)
  prob <- build_wendy_problem_joint(
    wendy_data_bl,
    f_, J_u, J_uu, J_up, J_p, J_pp, J_upp,
    dF_dt_ = dF_dt_, d2F_dt2_ = d2F_dt2_, d3F_dt3_ = d3F_dt3_,
    J = J_num, lip = lip, sig = sig, device = device, control = ctrl_ssl_bl
  )
  sys <- build_wendy_system_joint(
    list(prob), lip, control$diag_reg, control$use_interp_uncertainty, device
  )
  g               <- sys$g
  build_compute_b <- sys$build_compute_b

  # r_c is an integer grid-step count; clamp so left/right ranges never overlap
  r_c  <-   min(prob$rc, floor((mp1 - 1L) / 2L))
  l_idx   <- seq_len(r_c)
  r_idx   <- seq(mp1 - r_c + 1L, mp1)
  int_idx <- seq(r_c + 1L, mp1 - r_c)

  U_interior <- U[int_idx, , drop = FALSE]
  U_left0    <- U[l_idx,   , drop = FALSE]
  U_right0   <- U[r_idx,   , drop = FALSE]

  n_l     <- length(l_idx) * D
  n_r     <- length(r_idx) * D
  nu_bdry <- n_l + n_r
  u0_vec  <- c(as.vector(U_left0), as.vector(U_right0))

  reconstruct_U <- function(theta) {
    U_l <- matrix(theta[seq_len(n_l)],          nrow = length(l_idx), ncol = D)
    U_r <- matrix(theta[n_l + seq_len(n_r)],    nrow = length(r_idx), ncol = D)
    rbind(U_l, U_interior, U_r)
  }

  joint_target <- function(theta) {
    U_full <- reconstruct_U(theta)
    g_val  <- as.numeric(g(U_full, p, tt)$contiguous())
    b_val  <- as.numeric(build_compute_b(U_full, tt)$contiguous())
    g_val - b_val
  }

  theta_w <- gauss_newton_iterate(joint_target, u0_vec,
                                   max_iter = max_iter, tol = tol,
                                   lambda = lambda, alpha = alpha, sigma = sigma)

  list(
    Uhat            = reconstruct_U(theta_w),
    p               = p,
    g               = g,
    build_compute_b = build_compute_b
  )
}

# Gauss-Newton joint state and parameter estimate with boundary-layer test functions
gn_bl <- function(p, U, tt, build_compute_b, g, max_iter = 100, tol = 1e-10,
                  lambda = 1e-8, alpha = 0, sigma = NULL){
  D <- ncol(U)

  u0_vec <- as.vector(U)
  nu <- length(u0_vec)
  np <- length(p)

  joint_target <- function(theta) {
    U_cur <- matrix(theta[seq_len(nu)], ncol = D)
    p_cur <- theta[nu + seq_len(np)]
    g_val <- as.numeric(g(U_cur, p_cur, tt)$contiguous())
    b_val <- as.numeric(build_compute_b(U_cur, tt)$contiguous())
    g_val - b_val
  }

  theta0  <- c(u0_vec, p)
  mask    <- c(rep(1, nu), rep(0, np))
  theta_w <- gauss_newton_iterate(joint_target, theta0,
                                   penalty_mask = mask, theta_ref = theta0,
                                   max_iter = max_iter, tol = tol,
                                   lambda = lambda, alpha = alpha, sigma = sigma)

  list(Uhat = matrix(theta_w[seq_len(nu)], ncol = D),
       p    = theta_w[nu + seq_len(np)])
}

# Iterative (weak) re-weighted least squares
irls <- function(G, b, L, p0 = NULL, W = NULL, reg = 10e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
  dm <- nrow(G)
  alphaIdm <- reg * diag(rep(1, dm))
  W_mat <- if (!is.null(W)) as.array(W) else NULL
  p <- if(is.null(p0)) lm.fit(G, b)$coefficients else p0
  n <- 0
  SW <- Inf

  sw_pvalues <- numeric(max_its)

  while(n < max_its){
    pn1 <- p
    n <- n + 1

    # We solve the weighted least squares problem
    # https://en.wikipedia.org/wiki/Weighted_least_squares
    # GᵀS⁻¹G = GᵀS⁻¹b, S can be factored using the Cholesky decomposition of S = RᵀR
    # Gᵀ(RᵀR)⁻¹G = Gᵀ(RᵀR)⁻¹b which reduces to (GᵀR⁻¹)R⁻ᵀG =(GᵀR⁻¹)R⁻ᵀb
    # Thus we solve the least squares problem R⁻ᵀG = R⁻ᵀb where R is upper triangular matrix
    Ln <- as.array(L(p)$contiguous())
    Sn <- if (!is.null(W_mat)) (1 - reg) * Ln %*% W_mat %*% t(Ln) + alphaIdm
          else                 (1 - reg) * Ln %*% t(Ln) + alphaIdm
    RT <- t(chol(Sn)) # S = RᵀR R is upper triangular
    G_ <- forwardsolve(RT, G)
    b_ <- forwardsolve(RT, b)
    p <- lm.fit(G_, b_)$coefficients

    denom <- sqrt(sum(pn1^2))
    relative_change <- if (denom == 0 || !is.finite(denom)) Inf else sqrt(sum((p - pn1)^2)) / denom
    if (!is.finite(relative_change)) relative_change <- Inf

    residuals <- b_ - G_ %*% p
    residuals_clean <- residuals[is.finite(residuals)]

    if (length(residuals_clean) < 3) {
      p_val <- 0
    } else {
      sw_test <- shapiro.test(residuals_clean)
      p_val <- sw_test$p.value
      if (is.na(p_val)) p_val <- 0
    }
    sw_pvalues[n] <- p_val

    if(n >= n0){
      SW <- p_val
    } else {
      SW <- 1
    }

    if(relative_change > tau_FP && n < max_its && SW > tau_SW){
      next
    } else {
      break
    }
  }

  sw_pvalues <- sw_pvalues[1:n]

  return(list(
    p = p,
    iterations = n,
    converged = (relative_change <= tau_FP || SW <= tau_SW) && n < max_its,
    relative_change_n = relative_change,
    sw_pvalues = sw_pvalues,
    final_sw_pvalue = tail(sw_pvalues, n=1)
  ))
}

# Nonlinear iterative (weak) re-weighted least squares
nirls <- function(g, b, L, Jp_r, p0, W = NULL, reg = 10e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
  dm <- length(b)
  alphaIdm <- reg * diag(rep(1, dm))
  W_mat <- if (!is.null(W)) as.array(W) else NULL
  p <- p0
  n <- 0
  SW <- Inf

  sw_pvalues <- numeric(max_its)

  while(n < max_its){
    pn1 <- p
    n <- n + 1

    # We solve the nonlinear weighted least squares problem
    # we have the R⁻ᵀb = R⁻ᵀg(p)
    # r = R⁻ᵀ(g(p) - b) and want to minimize  1⁄2 ||r||²

    weighted_residual <- function(p, RT){
      forwardsolve(RT, as.array(g(p)$contiguous()) - b)
    }

    weighted_residual_jacobian <- function(p, RT){
      forwardsolve(RT, as.array(Jp_r(p)$contiguous()))
    }

    Ln <- as.array(L(p)$contiguous())
    Sn <- if (!is.null(W_mat)) (1 - reg) * Ln %*% W_mat %*% t(Ln) + alphaIdm
          else                 (1 - reg) * Ln %*% t(Ln) + alphaIdm
    RT <- t(chol(Sn)) # S = RᵀR R is upper triangular

    p <- nls.lm(p, lower = NULL, upper = NULL, function(p){weighted_residual(p, RT)}, function(p){weighted_residual_jacobian(p, RT)})$par

    denom <- sqrt(sum(pn1^2))
    relative_change <- if (denom == 0 || !is.finite(denom)) Inf else sqrt(sum((p - pn1)^2)) / denom
    if (!is.finite(relative_change)) relative_change <- Inf

    residuals <- weighted_residual(p, RT)

    sw_test <- shapiro.test(residuals)
    p_val <- sw_test$p.value
    sw_pvalues[n] <- p_val

    if(n >= n0){
      SW <- p_val
    } else {
      SW <- 1
    }

    if(relative_change > tau_FP && n < max_its && SW > tau_SW){
      next
    } else {
      break
    }
  }

  sw_pvalues <- sw_pvalues[1:n]

  return(list(
    p = p,
    iterations = n,
    converged = (relative_change <= tau_FP || SW <= tau_SW) && n < max_its,
    relative_change_n = relative_change,
    sw_pvalues = sw_pvalues,
    final_sw_pvalue = tail(sw_pvalues, n=1)
  ))
}

# Weak ordinary least squares
ols <- function(G, b, L, reg = 10e-10){
  J <- ncol(G)
  p <- solve(crossprod(G) + reg * diag(J), crossprod(G, b))
  return(list(p = as.vector(p)))
}

# Weak ordinary nonlinear least squares
nols <- function(g, b, L, Jp_r, p0, reg = 10e-10){
  residual <- function(p){
    as.array(g(p)$contiguous() - b)
  }
  residual_jacobian <- function(p){
    as.array(Jp_r(p)$contiguous())
  }
  p <- nls.lm(p0, lower = NULL, upper = NULL, residual, residual_jacobian)$par
  return(list(p = p))
}

# Maximum likelihood estimation for r(p) ~ N(0, S(p))
mle <- function(p0, wnll, J_wnll, H_wnll, S, Jp_r, control){

  objfun <- function(p) {
      f <- wnll(p)
      g <- J_wnll(p)
      h <- H_wnll(p)
    list(value = f, gradient = g, hessian = h)
  }

  data <-  trust::trust(objfun, p0, rinit = 25, rmax = 200, blather = FALSE) 
  phat <- as.vector(data$argument)

  data$p <- phat
  
  return(data)
}

# Output Error forward based nonlinear least squares (NLS)
output_error <- function(f, U, tt, p0, lower = NULL, upper = NULL) {
  D  <- ncol(U) 
  J  <- length(p0)
  u0_init <- as.vector(U[1, ])

  obs <- data.frame(time = tt)
  for (d in seq_len(D)) obs[[paste0("x", d)]] <- U[, d]

  p_lower <- if (!is.null(lower)) lower else rep(-Inf, J)
  p_upper <- if (!is.null(upper)) upper else rep(Inf, J)

  u0_lower <- rep(-Inf, D)
  u0_upper <- rep(Inf, D)

  theta0 <- c(p0, u0_init)
  lower  <- c(p_lower, u0_lower)
  upper  <- c(p_upper, u0_upper)

  modelODE <- function(t, state, parms) list(as.vector(f(state, parms, t)))

  modelRun <- function(theta) {
    parms <- theta[1:J]
    u0    <-  theta[(J + 1):(J + D)]
    
    sol <- tryCatch({
      capture.output(
        result <- suppressWarnings(
          deSolve::ode(y = u0, times = tt, func = modelODE, parms = parms)
        ),
        type = "output"
      )
      result
    }, error = function(e) NULL)

    if (is.null(sol) || nrow(sol) != length(tt)) return(NULL)
      out <- as.data.frame(sol)
      colnames(out) <- c("time", paste0("x", seq_len(D)))
      out
  }

  costFn <- function(theta) {
    out <- modelRun(theta)
    if (is.null(out)) {
      bad_model <- obs
      bad_model[, -1] <- 1e6
      return(FME::modCost(model = bad_model, obs = obs))  
    }
    FME::modCost(model = out, obs = obs)         
  }

  fit <- tryCatch({
      FME::modFit(f = costFn, p = theta0, method = "Marq", lower = lower, upper = upper)
  }, error = function(e) e)

  n_theta <- J + D

  if (inherits(fit, "error")) {
    return(list(
      p          = rep(NA, J),
      u0         = rep(NA, D),
      cov        = matrix(NA, nrow = n_theta, ncol = n_theta),
      converged  = NA,
      iterations = NA
    ))
  }

  fit_summary <- suppressWarnings(summary(fit))
  cov_scaled  <- fit_summary$cov.scaled
  cov <- if (!is.null(cov_scaled) && !anyNA(cov_scaled)) cov_scaled else NULL
  
  return(list(
    p          = fit$par[1:J],
    u0         = fit$par[(J + 1):(J + D)],
    cov        = cov,
    converged  = fit$info %in% c(1L, 2L, 3L),
    iterations = fit_summary$niter,
    ssr        = fit$ssr
  ))
}

# Shared Equation-Error grid setup: optionally interpolate to max_points,
# smooth columns with splines, return midpoint state, time, and du/dt.
prepare_ee_grid <- function(U, tt, sigma = NULL, max_points = 256, poly_degree = 3) {
  tt <- as.vector(tt)
  if (!is.null(sigma) && nrow(U) < max_points) {
    nsr     <- min(sigma / apply(U, 2, sd))
    tt_fine <- seq(min(tt), max(tt), length.out = max_points)
    method  <- if (nsr <= 0.1) "linear" else paste0("poly_ls_", poly_degree)
    U  <- do.call(cbind, lapply(seq_len(ncol(U)), function(d) {
      fit_col(U[, d], tt, tt_fine, method, sigma = sigma)$fit
    }))
    tt <- tt_fine
  }
  U_smooth <- apply(U, 2, function(col) smooth.spline(tt, col)$y)
  u_mid    <- (U_smooth[-nrow(U_smooth), , drop = FALSE] + U_smooth[-1, , drop = FALSE]) / 2
  dudt     <- (U_smooth[-1, , drop = FALSE] - U_smooth[-nrow(U_smooth), , drop = FALSE]) / diff(tt)
  t_mid    <- (tt[-length(tt)] + tt[-1]) / 2
  list(u_mid = u_mid, t_mid = t_mid, dudt = dudt, n_mid = nrow(u_mid))
}

# Equation Error initial guess for linear-in-parameters ODEs.
# Smooths U with splines, approximates du/dt at midpoints, solves J_p(u,t) p = du/dt.
ee_linear <- function(U, tt, f_, J_p, J, D, sigma = NULL, max_points = 256, poly_degree = 3) {
  grid <- prepare_ee_grid(U, tt, sigma, max_points, poly_degree)
  n_mid <- grid$n_mid
  input <- rbind(matrix(0, nrow = J, ncol = n_mid), t(grid$u_mid), matrix(grid$t_mid, nrow = 1))
  # Affine offset: b(u,t) = f(u, 0, t). For purely linear problems this is zero.
  f0     <- f_(input)
  jp_raw <- J_p(input)  # n_mid x (D*J), D-outer J-inner per row
  A_mat  <- do.call(rbind, lapply(seq_len(n_mid), function(i) {
    matrix(jp_raw[i, ], nrow = D, ncol = J, byrow = TRUE)
  }))
  lm.fit(A_mat, as.vector(t(grid$dudt - f0)))$coefficients
}

# Equation Error initial guess for nonlinear-in-parameters ODEs.
# Minimizes ||f(u_mid, p, t_mid) - du/dt||^2 via Levenberg-Marquardt.
ee_nonlinear <- function(U, tt, f_, J_p, J, D, sigma = NULL, max_points = 256, poly_degree = 3) {
  grid  <- prepare_ee_grid(U, tt, sigma, max_points, poly_degree)
  u_mid <- grid$u_mid; t_mid <- grid$t_mid; dudt <- grid$dudt; n_mid <- grid$n_mid

  build_input <- function(p) {
    p_mat <- matrix(rep(p, n_mid), nrow = J)
    rbind(p_mat, t(u_mid), matrix(t_mid, nrow = 1))
  }
  dm_residual <- function(p) as.vector(t(f_(build_input(p)) - dudt))
  dm_jacobian <- function(p) {
    jp_raw <- J_p(build_input(p))  # n_mid x (D*J), D-outer J-inner per row
    do.call(rbind, lapply(seq_len(n_mid), function(i) {
      matrix(jp_raw[i, ], nrow = D, ncol = J, byrow = TRUE)
    }))
  }
  nls.lm(rep(1, J), fn = dm_residual, jac = dm_jacobian)$par
}

# Equation-error initial p0 used by every method that needs one. The
# linear/nonlinear branch matches f's structure in p; lognormal uses the
# log-transformed data so EE coordinates match f_.
.ee_init_p0 <- function(ctx) {
  ee_degree <- detect_max_state_order(ctx$f_orig_expr, ctx$u_expr) +
               if (ctx$noise_dist == "lognormal") 0L else 1L
  ee_fn <- if (ctx$lip) ee_linear else ee_nonlinear
  ee_fn(ctx$U_processed, ctx$tt_processed, ctx$f_, ctx$J_p, ctx$J, ctx$D,
        sigma = ctx$estimated_sd, poly_degree = ee_degree)
}

# OLS fallback for linear-in-p methods that need a warm-start p0 (MLE, IRLS,
# HYBRID, JOINT stage 1). NIRLS/NOLS use user-supplied or EE p0 directly.
.ols_warmstart <- function(ctx, G_cont, b_cont) {
  ols(G_cont, b_cont, ctx$system$L)$p
}

.run_ols <- function(ctx, p0, G_cont, b_cont) {
  sys <- ctx$system
  if (ctx$lip) ols(G_cont, b_cont, sys$L)
  else         nols(sys$g, b_cont, sys$L, sys$Jp_r, p0 = p0, reg = 10e-10)
}

.run_irls <- function(ctx, p0, G_cont, b_cont) {
  sys <- ctx$system
  if (ctx$lip) {
    if (is.null(p0)) p0 <- .ols_warmstart(ctx, G_cont, b_cont)
    irls(G_cont, b_cont, sys$L, p0 = p0, W = sys$W, max_its = ctx$control$max_iterates)
  } else {
    nirls(sys$g, b_cont, sys$L, sys$Jp_r, p0 = p0, W = sys$W, max_its = ctx$control$max_iterates)
  }
}

.run_mle <- function(ctx, p0, G_cont, b_cont) {
  sys <- ctx$system
  if (ctx$lip && is.null(p0)) p0 <- .ols_warmstart(ctx, G_cont, b_cont)
  mle(p0, sys$wnll, sys$J_wnll, sys$H_wnll, sys$S, sys$Jp_r, ctx$control)
}

# For lognormal noise the user-supplied f sees observed (positive) data, but
# the WENDy system fits in log space; OE drives dz/dt for z = log(u) and
# compares to log-observations.
.f_for_oe <- function(ctx) {
  if (ctx$noise_dist == "lognormal") {
    f_user <- ctx$f
    function(z, p, t) f_user(exp(z), p, t) / exp(z)
  } else {
    ctx$f
  }
}

.run_oe <- function(ctx, p0, G_cont, b_cont) {
  if (ctx$noise_dist == "lognormal") {
    output_error(.f_for_oe(ctx), ctx$U_processed, ctx$tt_processed, p0)
  } else {
    output_error(ctx$f, ctx$U_orig, ctx$tt_orig, p0)
  }
}

.run_hybrid <- function(ctx, p0, G_cont, b_cont) {
  mle_result <- .run_mle(ctx, p0, G_cont, b_cont)
  if (ctx$noise_dist == "lognormal") {
    output_error(.f_for_oe(ctx), ctx$U_processed, ctx$tt_processed, mle_result$p)
  } else {
    output_error(ctx$f, ctx$U_orig, ctx$tt_orig, mle_result$p)
  }
}

# JOINT: MLE-warm-started Gauss-Newton on the joint (U, p) residual, with an
# optional boundary-refinement second pass.
#   stage 1: MLE on p (with OLS warm-start for linear-in-p) + state smoother
#   stage 2: gn() / gls_gn() for MSG or gn_bl() for SSL+BL — joint U and p
#   stage 3 (optional, two_stage = TRUE): gn_boundary refines corner U rows
# control$gn_method ("std" default, or "gls") selects the MSG stage-2 solver.
.run_joint <- function(ctx, p0, G_cont, b_cont) {
  control <- ctx$control
  sys     <- ctx$system
  is_ssl  <- identical(control$test_fun_type, "SSL")
  use_gls <- identical(control$gn_method, "gls")

  # Use the grid that V/V_prime were built on. For lognormal this is already
  # log-space (wendy_data is built post preprocess_data); for sparse data it is
  # the densified interpolated grid.
  U_gn  <- ctx$U_sys
  tt_gn <- ctx$tt_sys

  if (ctx$lip && is.null(p0)) p0 <- .ols_warmstart(ctx, G_cont, b_cont)
  p <- mle(p0, sys$wnll, sys$J_wnll, sys$H_wnll, sys$S, sys$Jp_r, control)$p

  smooth <- if (identical(control$smoother, "erts")) {
    wendy_erts(U_gn, ctx$f_, ctx$J_u, tt_gn, p, control, sigma = ctx$sig)
  } else {
    gp_smooth(U_gn, tt_gn, sigma2_n = as.numeric(ctx$sig)^2)
  }

  res1 <- if (is_ssl) {
    if (use_gls) warning("gn_method = 'gls' has no SSL variant; falling back to gn_bl.", call. = FALSE)
    sj <- ctx$system_joint
    gn_bl(p, smooth$U_star, tt_gn, sj$build_compute_b, sj$g,
          alpha = control$gn_alpha %||% 0, sigma = ctx$sig)
  } else if (use_gls) {
    gls_gn(p, smooth$U_star, tt_gn, ctx$f_, ctx$V, ctx$V_prime, sys$S)
  } else {
    gn(p, smooth$U_star, tt_gn, ctx$f_, ctx$V, ctx$V_prime,
       alpha = control$gn_alpha %||% 0, sigma = ctx$sig)
  }

  if (!isTRUE(control$two_stage)) return(res1)

  # MSG pins the boundary refinement to the stage-1 radius; SSL re-detects
  # its own change-point radius from U_ref (noisy data).
  fixed_radius_arg <- if (is_ssl) NULL else ctx$min_radius
  gn_boundary(
    res1$p, res1$Uhat, tt_gn,
    ctx$f_, ctx$J_u, ctx$J_uu, ctx$J_up, ctx$J_p, ctx$J_pp, ctx$J_upp,
    dF_dt_ = ctx$dF_dt_, d2F_dt2_ = ctx$d2F_dt2_, d3F_dt3_ = ctx$d3F_dt3_,
    J_num = ctx$J, lip = ctx$lip, sig = ctx$sig, control = control,
    U_ref = U_gn, fixed_radius = fixed_radius_arg,
    alpha = control$gn_alpha %||% 0, sigma = ctx$sig
  )
}

.method_runners <- list(
  JOINT  = .run_joint,
  OLS    = .run_ols,
  IRLS   = .run_irls,
  MLE    = .run_mle,
  OE     = .run_oe,
  HYBRID = .run_hybrid
)

# Run a single optimization pass given a pre-built wendy system.
# Returns list(p, data, objective) where objective = wnll(phat) (or SSR for
# OE/HYBRID), used to rank multistart results.
.run_wendy_optimization <- function(ctx, p0) {
  method <- ctx$method
  runner <- .method_runners[[method]]
  if (is.null(runner)) stop("Unknown method: ", method, call. = FALSE)

  # Auto-init p0 via equation error when none is supplied. Linear OLS/IRLS
  # can solve in closed form; everyone else needs a starting p.
  if (is.null(p0) && (!ctx$lip || method == "OE")) {
    p0 <- .ee_init_p0(ctx)
  }

  # OE drives an ODE; the other methods consume base-R arrays from the
  # weak-residual system.
  G_cont <- if (method != "OE") as.array(ctx$system$G$contiguous()) else NULL
  b_cont <- if (method != "OE") as.array(ctx$system$b$contiguous()) else NULL

  data <- runner(ctx, p0, G_cont, b_cont)

  phat <- data$p
  objective <- if (method %in% c("OE", "HYBRID")) {
    data$ssr %||% Inf
  } else {
    tryCatch(ctx$system$wnll(phat), error = function(e) Inf)
  }
  list(p = phat, data = data, objective = objective)
}

#' Optimize a Pre-Built WENDy System
#'
#' Run parameter optimization on a \code{wendy} object that was built with
#' \code{control = list(optimize = FALSE)}.  The expensive setup (symbolics,
#' test functions, covariance assembly) is reused; only the optimizer runs.
#'
#' @param wendy_obj A \code{wendy} object returned by \code{solveWendy} with
#'   \code{control$optimize = FALSE}.
#' @param p0 Numeric vector (single start) or matrix (multistart, one row per
#'   starting point).  \code{NULL} triggers automatic equation-error
#'   initialisation (same as \code{solveWendy}).
#' @param method Optimization method; defaults to the method stored in
#'   \code{wendy_obj}.
#' @param control Optional list of control parameters to override those stored
#'   in \code{wendy_obj}.
#'
#' @return A \code{wendy} object with \code{$phat} and \code{$data} populated.
#'   Multistart runs additionally attach \code{$multistart_results} and
#'   \code{$multistart_objectives}.
#' @export
optimizeWendy <- function(wendy_obj, p0 = NULL,
                          method  = attr(wendy_obj, "method"),
                          control = list()) {
  stored_control <- wendy_obj$opt_ctx$control %||% list()
  merged_control <- modifyList(stored_control, control)

  ctx <- wendy_obj$opt_ctx
  if (is.null(ctx)) {
    stop("wendy_obj does not contain an optimization context. ",
         "Re-run solveWendy() with control = list(optimize = FALSE).",
         call. = FALSE)
  }
  ctx$method  <- method
  ctx$control <- merged_control

  res <- wendy_obj

  if (is.matrix(p0) && nrow(p0) > 1) {
    apply_fn <- merged_control$apply_fn %||% lapply

    all_results <- apply_fn(seq_len(nrow(p0)), function(i) {
      .run_wendy_optimization(ctx, p0[i, ])
    })

    objectives <- sapply(all_results, `[[`, "objective")
    best <- all_results[[which.min(objectives)]]

    res$phat <- best$p
    res$data <- best$data
    res$multistart_results <- all_results
    res$multistart_objectives <- objectives
  } else {
    p0_vec <- if (is.matrix(p0)) p0[1, ] else p0
    result  <- .run_wendy_optimization(ctx, p0_vec)
    res$phat <- result$p
    res$data <- result$data
  }

  attr(res, "method") <- method
  return(res)
}