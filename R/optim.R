# Robust Shapiro-Wilk p-value function
.sw_pvalue <- function(residuals) {
  rc <- residuals[is.finite(residuals)]
  if (length(rc) < 3L) return(0)
  if (length(rc) > 5000L) rc <- rc[round(seq.int(1, length(rc), length.out = 5000L))]
  pv <- shapiro.test(rc)$p.value
  if (is.na(pv)) 0 else pv
}

# Iterative (weak) re-weighted least squares
irls <- function(G, b, L, p0 = NULL, W = NULL, reg = 10e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
  p <- if(is.null(p0)) lm.fit(G, b)$coefficients else p0
  n <- 0
  SW <- Inf

  sw_pvalues <- numeric(max_its)

  while(n < max_its){
    pn1 <- p
    n <- n + 1

    # Weighted least squares via Cholesky: with S = R^T R upper-triangular R,
    # solve R^{-T} G p = R^{-T} b in unweighted form.
    # https://en.wikipedia.org/wiki/Weighted_least_squares
    Ln <- L(p)
    Sn <- if (!is.null(W)) (1 - reg) * Ln %*% (W * t(Ln))
          else             (1 - reg) * tcrossprod(Ln)
    diag(Sn) <- diag(Sn) + reg            # + reg * I, without a dense dm×dm allocation
    RT <- t(chol(Sn)) # S = R^T R, R upper triangular
    Gb_ <- forwardsolve(RT, cbind(G, b))  # single triangular solve for [G | b]
    G_  <- Gb_[, seq_len(ncol(G)), drop = FALSE]
    b_  <- Gb_[, ncol(G) + 1L]
    p <- lm.fit(G_, b_)$coefficients

    denom <- sqrt(sum(pn1^2))
    relative_change <- if (denom == 0 || !is.finite(denom)) Inf else sqrt(sum((p - pn1)^2)) / denom
    if (!is.finite(relative_change)) relative_change <- Inf

    residuals <- b_ - G_ %*% p
    p_val <- .sw_pvalue(residuals)
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
  p <- p0
  n <- 0
  SW <- Inf

  sw_pvalues <- numeric(max_its)

  while(n < max_its){
    pn1 <- p
    n <- n + 1

    # Nonlinear weighted least squares: minimize 0.5 ||R^{-T} (g(p) - b)||^2
    weighted_residual <- function(p, RT){
      forwardsolve(RT, g(p) - b)
    }
    weighted_residual_jacobian <- function(p, RT){
      forwardsolve(RT, Jp_r(p))
    }

    Ln <- L(p)
    Sn <- if (!is.null(W)) (1 - reg) * Ln %*% (W * t(Ln))
          else             (1 - reg) * tcrossprod(Ln)
    diag(Sn) <- diag(Sn) + reg            # + reg * I, without a dense dm×dm allocation
    RT <- t(chol(Sn)) # S = R^T R, R upper triangular

    p <- nls.lm(p, lower = NULL, upper = NULL, function(p){weighted_residual(p, RT)}, function(p){weighted_residual_jacobian(p, RT)})$par

    denom <- sqrt(sum(pn1^2))
    relative_change <- if (denom == 0 || !is.finite(denom)) Inf else sqrt(sum((p - pn1)^2)) / denom
    if (!is.finite(relative_change)) relative_change <- Inf

    residuals <- weighted_residual(p, RT)
    p_val <- .sw_pvalue(residuals)
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
  residual          <- function(p) g(p) - b
  residual_jacobian <- function(p) Jp_r(p)
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

  # Residual vector (model - observed), stacked over all state columns. This is
  # exactly what FME::modCost produced with its default unweighted/unscaled
  # settings, but handed straight to nls.lm without the FME wrapper.
  #
  # On ODE failure, fall back to an Euler-style first-order expansion of the
  # trajectory at t[1]: u(t) ≈ u0 + (t - t[1]) * f(u0, p, t[1]). This makes
  # the residuals a smooth function of every theta component, and (unlike a flat
  # 1e6 or any g(theta)*h(t) penalty) the resulting Jacobian has row x column
  # entanglement so it is full rank and not parallel to fvec. The 1e6 offset
  # keeps failed evaluations dominant over any feasible solution.
  resFn <- function(theta) {
    out <- modelRun(theta)
    if (is.null(out)) {
      parms <- theta[1:J]
      u0    <- theta[(J + 1):(J + D)]
      f0    <- tryCatch(as.vector(f(u0, parms, tt[1])),
                        error = function(e) rep(0, D))
      bad_model <- matrix(0, nrow = length(tt), ncol = D)
      for (d in seq_len(D)) {
        bad_model[, d] <- 1e6 + u0[d] + (tt - tt[1]) * f0[d]
      }
      return(as.vector(bad_model - U))
    }
    model <- as.matrix(out[, paste0("x", seq_len(D)), drop = FALSE])
    as.vector(model - U)
  }

  n_theta <- J + D

  # Only pass box constraints when the user actually supplied finite bounds;
  # otherwise leave them NULL so nls.lm runs unconstrained.
  use_bounds <- any(is.finite(lower)) || any(is.finite(upper))
  fit <- tryCatch({
    if (use_bounds) {
      nls.lm(par = theta0, lower = lower, upper = upper, fn = resFn)
    } else {
      nls.lm(par = theta0, fn = resFn)
    }
  }, error = function(e) e)

  if (inherits(fit, "error")) {
    return(list(
      p          = rep(NA, J),
      u0         = rep(NA, D),
      cov        = matrix(NA, nrow = n_theta, ncol = n_theta),
      converged  = NA,
      iterations = NA
    ))
  }

  # Scaled nonlinear-LS covariance (SSR / (N - npar)) * (JᵀJ)⁻¹, reproducing
  # FME::modFit's summary$cov.scaled. fit$hessian is the Gauss-Newton JᵀJ.
  n_res <- length(fit$fvec)
  cov <- tryCatch({
    resvar <- fit$deviance / (n_res - n_theta)
    cov_scaled <- resvar * solve(fit$hessian)
    if (!anyNA(cov_scaled)) cov_scaled else NULL
  }, error = function(e) NULL)

  return(list(
    p          = fit$par[1:J],
    u0         = fit$par[(J + 1):(J + D)],
    cov        = cov,
    converged  = fit$info %in% c(1L, 2L, 3L),
    iterations = fit$niter,
    ssr        = fit$deviance
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
  jp_raw <- J_p(input)  # n_mid x (D*J), column-major (d, j) d-fast per row
  A_mat  <- do.call(rbind, lapply(seq_len(n_mid), function(i) {
    matrix(jp_raw[i, ], nrow = D, ncol = J)  # d-fast: A[d, j] = df_d/dp_j
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
    jp_raw <- J_p(build_input(p))  # n_mid x (D*J), column-major (d, j) d-fast per row
    do.call(rbind, lapply(seq_len(n_mid), function(i) {
      matrix(jp_raw[i, ], nrow = D, ncol = J)
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

.method_runners <- list(
  OLS    = .run_ols,
  IRLS   = .run_irls,
  MLE    = .run_mle,
  OE     = .run_oe,
  HYBRID = .run_hybrid
)

.run_wendy_optimization <- function(ctx, p0) {
  method <- ctx$method
  runner <- .method_runners[[method]]
  if (is.null(runner)) stop("Unknown method: ", method, call. = FALSE)

  if (is.null(p0) && (!ctx$lip || method == "OE")) {
    warning("p0 not supplied with a nonlinear-in-parameters system; ",
            "initializing p0 with strong-form Equation Error.",
            call. = FALSE)
    p0 <- .ee_init_p0(ctx)
  }

  G_cont <- if (method != "OE") ctx$system$G else NULL
  b_cont <- if (method != "OE") ctx$system$b else NULL

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