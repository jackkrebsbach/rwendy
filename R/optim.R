# Iterative (weak) re-weighted least squares
irls <- function(G, b, L, W = NULL, reg = 10e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
  dm <- nrow(G)
  alphaIdm <- reg * diag(rep(1, dm))
  W_mat <- if (!is.null(W)) as.array(W) else NULL
  p <- lm.fit(G, b)$coefficients
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
    converged = (relative_change <= tau_FP || SW <= tau_SW),
    relative_change_n = relative_change,
    sw_pvalues = sw_pvalues,
    final_sw_pvalue = tail(sw_pvalues, n=1)
  ))
}

# Nonlinear iterative (weak) re-weighted least squares
nirls <- function(g, b, L, Jp_r, p0, W = NULL, reg = 1e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
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
    converged = (relative_change <= tau_FP || SW <= tau_SW),
    relative_change_n = relative_change,
    sw_pvalues = sw_pvalues,
    final_sw_pvalue = tail(sw_pvalues, n=1)
  ))
}

# Weak ordinary least squares
ols <- function(G, b, L, reg = 1e-10){
  p <- lm.fit(G, b)$coefficients
  residuals <- b - G %*% p
  sw_test <- shapiro.test(residuals)
  sw_p_value <- sw_test$p.value
  return(list(p = p, sw_p_value = sw_p_value))
}

# Weak ordinary nonlinear least squares
nols <- function(g, b, L, Jp_r, p0, reg = 1e-10){
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

  # Extracting estimated parameters from various packages
  # trust.optim -> data$solution = phat
  # trust::trust -> data$argument = phat
  # trust.optim(p0, wnll, J_wnll, method = "BFGS") 
}


# Output Error
output_error <- function(f, U, tt, p0, use_bounds = FALSE, u0_range_factor = 10.0, lower = NULL, upper = NULL) {
  D  <- ncol(U) 
  J  <- length(p0)
  u0_init <- as.vector(U[1, ])

  obs <- data.frame(time = tt)
  for (d in seq_len(D)) obs[[paste0("x", d)]] <- U[, d]

  p_lower <- if (!is.null(lower)) lower else rep(-Inf, J)
  p_upper <- if (!is.null(upper)) upper else rep(Inf, J)

  if (use_bounds) {
      U_mat     <- as.matrix(U)
      u0_scale  <- apply(U_mat, 2, function(col) diff(range(col)) + 1e-6)
      u0_lower  <- u0_init - u0_range_factor * u0_scale
      u0_upper  <- u0_init + u0_range_factor * u0_scale
    } else {
      u0_lower <- rep(-Inf, D)
      u0_upper <- rep(Inf, D)
  }
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

  data <- fit

  return(list(p = data$par[1:J], u0 = data$par[(J+1):(J+D)], data = data))
}


# Equation Error initial guess for linear-in-parameters ODEs.
# Smooths U with splines, approximates du/dt at midpoints, solves J_p(u,t) p = du/dt.
ee_linear <- function(U, tt, f_, J_p, J, D) {
  tt       <- as.vector(tt)
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
ee_nonlinear <- function(U, tt, f_, J_p, J, D) {
  tt       <- as.vector(tt)
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