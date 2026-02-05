# Iterative (weak) re-weighted least squares
irls <- function(G, b, L, reg = 10e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
  dm <- nrow(G)
  alphaIdm <- reg * diag(rep(1, dm))
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
    Sn <- (1 - reg) * Ln %*% t(Ln) + alphaIdm
    RT <- t(chol(Sn)) # S = RᵀR R is upper triangular
    G_ <- forwardsolve(RT, G)
    b_ <- forwardsolve(RT, b)
    p <- lm.fit(G_, b_)$coefficients

    relative_change <- sqrt(sum((p - pn1)^2)) / sqrt(sum(pn1^2))

    residuals <- b - G %*% p

    if(n >= n0){
      sw_test <- shapiro.test(residuals)
      sw_pvalues[n] <- sw_test$p.value
      SW <- sw_test$p.value
    } else {
      sw_pvalues[n] <- 1
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
    final_sw_pvalue = SW
  ))
}

# Nonlinear iterative (weak) re-weighted least squares
nirls <- function(g, b, L, Jp_r, p0, reg = 10e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
  dm <- length(b)
  alphaIdm <- reg * diag(rep(1, dm))
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
    Sn <- (1 - reg) * Ln %*% t(Ln) + alphaIdm
    RT <- t(chol(Sn)) # S = RᵀR R is upper triangular

    p <- nls.lm(p, lower = NULL, upper = NULL, function(p){weighted_residual(p, RT)}, function(p){weighted_residual_jacobian(p, RT)})$par

    relative_change <- sqrt(sum((p - pn1)^2)) / sqrt(sum(pn1^2))

    residuals <- b - g(p) 

    if(n >= n0){
      sw_test <- shapiro.test(residuals)
      sw_pvalues[n] <- sw_test$p.value
      SW <- sw_test$p.value
    } else {
      sw_pvalues[n] <- 1
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
    final_sw_pvalue = SW
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
# Extracting estimated parameters from various packages
# trust.optim -> data$solution = phat
# trust::trust -> data$argument = phat
# trust.optim(p0, wnll, J_wnll, method = "BFGS") 