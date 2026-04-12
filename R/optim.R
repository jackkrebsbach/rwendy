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

    relative_change <- sqrt(sum((p - pn1)^2)) / sqrt(sum(pn1^2))

    residuals <- b_ - G_ %*% p

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

    relative_change <- sqrt(sum((p - pn1)^2)) / sqrt(sum(pn1^2))

    residuals <- weighted_residual(p, RT)

    print("here")

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
  p <- lm.fit(G, b)$coefficients
  residuals <- b - G %*% p
  sw_test <- shapiro.test(residuals)
  sw_p_value <- sw_test$p.value
  return(list(p = p, sw_p_value = sw_p_value))
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
  residuals <- b - as.array(g(p)$contiguous())
  sw_test <- shapiro.test(residuals)
  sw_p_value <- sw_test$p.value
  return(list(p = p, sw_p_value = sw_p_value))
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

# Output Error forward based nonlinear least squares
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
      bad_model      <- obs
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

# Equation Error initial guess for linear-in-parameters ODEs.
# Smooths U with splines, approximates du/dt at midpoints, solves J_p(u,t) p = du/dt.
ee_linear <- function(U, tt, f_, J_p, J, D, sigma = NULL, max_points = 256, poly_degree = 3) {
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
  t_mid    <- (tt[-length(tt)] + tt[-1]) / 2
  U_smooth <- apply(U, 2, function(col) smooth.spline(tt, col)$y)
  u_mid    <- (U_smooth[-nrow(U_smooth), , drop = FALSE] + U_smooth[-1, , drop = FALSE]) / 2
  dudt     <- (U_smooth[-1, , drop = FALSE] - U_smooth[-nrow(U_smooth), , drop = FALSE]) / diff(tt)
  n_mid    <- nrow(u_mid)
  input    <- rbind(matrix(0, nrow = J, ncol = n_mid), t(u_mid), matrix(t_mid, nrow = 1))
  # Affine offset: b(u,t) = f(u, 0, t). For purely linear problems this is zero.
  f0       <- f_(input)  # n_mid x D
  jp_raw   <- J_p(input)  # n_mid x (D*J), D-outer J-inner per row
  A_mat    <- do.call(rbind, lapply(seq_len(n_mid), function(i) {
    matrix(jp_raw[i, ], nrow = D, ncol = J, byrow = TRUE)
  }))
  lm.fit(A_mat, as.vector(t(dudt - f0)))$coefficients
}

# Equation Error initial guess for nonlinear-in-parameters ODEs.
# Minimizes ||f(u_mid, p, t_mid) - du/dt||^2 via Levenberg-Marquardt.
ee_nonlinear <- function(U, tt, f_, J_p, J, D, sigma = NULL, max_points = 256, poly_degree = 3) {
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
  t_mid    <- (tt[-length(tt)] + tt[-1]) / 2
  U_smooth <- apply(U, 2, function(col) smooth.spline(tt, col)$y)
  u_mid    <- (U_smooth[-nrow(U_smooth), , drop = FALSE] + U_smooth[-1, , drop = FALSE]) / 2
  dudt     <- (U_smooth[-1, , drop = FALSE] - U_smooth[-nrow(U_smooth), , drop = FALSE]) / diff(tt)
  n_mid    <- nrow(u_mid)
  dm_residual <- function(p) {
    p_mat <- matrix(rep(p, n_mid), nrow = J)
    input <- rbind(p_mat, t(u_mid), matrix(t_mid, nrow = 1))
    as.vector(t(f_(input) - dudt))
  }
  dm_jacobian <- function(p) {
    p_mat  <- matrix(rep(p, n_mid), nrow = J)
    input  <- rbind(p_mat, t(u_mid), matrix(t_mid, nrow = 1))
    jp_raw <- J_p(input)  # n_mid x (D*J), D-outer J-inner per row
    do.call(rbind, lapply(seq_len(n_mid), function(i) {
      matrix(jp_raw[i, ], nrow = D, ncol = J, byrow = TRUE)
    }))
  }
  nls.lm(rep(1, J), fn = dm_residual, jac = dm_jacobian)$par
}

# Run a single optimization pass given a pre-built wendy system.
# Returns list(p, data, objective) where objective = wnll(phat) for comparing starts.
.run_wendy_optimization <- function(ctx, p0) {
  system       <- ctx$system
  method       <- ctx$method
  f            <- ctx$f
  U            <- ctx$U_orig
  tt           <- ctx$tt_orig
  lip          <- ctx$lip
  f_orig_expr  <- ctx$f_orig_expr
  u_expr       <- ctx$u_expr
  estimated_sd <- ctx$estimated_sd
  f_           <- ctx$f_
  J_p          <- ctx$J_p
  J            <- ctx$J
  D            <- ctx$D
  noise_dist   <- ctx$noise_dist
  control      <- ctx$control

  g      <- system$g
  b      <- system$b
  G      <- system$G
  L      <- system$L
  W      <- system$W
  Jp_r   <- system$Jp_r
  S      <- system$S
  wnll   <- system$wnll
  J_wnll <- system$J_wnll
  H_wnll <- system$H_wnll

  # Compute p0 via equation error when none is supplied.
  # Needed for: all nonlinear methods, and OE (no p0 fallback).
  # Linear OLS/IRLS don't need p0; MLE+lip and HYBRID+lip have a fallbacks.
  if (is.null(p0) && (!lip || method == "OE")) {
    ee_degree <- detect_max_state_order(f_orig_expr, u_expr) + if (noise_dist == "lognormal") 0L else 1L
    p0 <- if (lip) {
      ee_linear(U, tt, f_, J_p, J, D, sigma = estimated_sd, poly_degree = ee_degree)
    } else {
      ee_nonlinear(U, tt, f_, J_p, J, D, sigma = estimated_sd, poly_degree = ee_degree)
    }
  }
  
  # Convert from torch tensors to R arrays
  b_cont <- as.array(b$contiguous())
  G_cont <-  as.array(G$contiguous())

  data <- switch(method,
    # Ordinary Least Squares on the weak residual
    OLS = if (lip) {
      ols(G_cont, b_cont, L)
    } else {
      nols(g, b_cont, L, Jp_r, p0 = p0, reg = 10e-10)
    },
    # Iterative Reweighted Least Squares on the weak residual
    IRLS = if (lip) {
      if (is.null(p0)) {
        p0 <- ols(G_cont, b_cont, L)$p
      }
      irls(G_cont, b_cont, L, p0 = p0, W = W, max_its = control$max_iterates)
    } else {
      nirls(g, b_cont, L, Jp_r, p0 = p0, W = W, max_its = control$max_iterates)
    },
    # Maximum Likelihood Estimation
    MLE = {
      if (lip && is.null(p0)) {
        p0 <- ols(G_cont, b_cont, L)$p
      }
      mle(p0, wnll, J_wnll, H_wnll, S, Jp_r, control)
    },
    # Output Error
    OE = output_error(f, U, tt, p0),
    # WENDy-MLE + OE
    HYBRID = {
      if (lip && is.null(p0)) {
        p0 <- ols(G_cont, b_cont, L)$p
      }
      mle_result <- mle(p0, wnll, J_wnll, H_wnll, S, Jp_r, control)
      output_error(f, U, tt, mle_result$p)
    }
  )

  phat <- data$p
  objective <- if (method %in% c("OE", "HYBRID")) {
    data$ssr %||% Inf
  } else {
    tryCatch(wnll(phat), error = function(e) Inf)
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