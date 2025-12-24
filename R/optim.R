
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
    # GᵀS⁻¹G = GᵀS⁻¹b, S can be factored using the Cholesky decomposition of S = RRᵀ
    # Gᵀ(RRᵀ)⁻¹G = Gᵀ(RRᵀ)⁻¹b which reduces to (GᵀR⁻ᵀ)R⁻¹G =(GᵀR⁻ᵀ)R⁻¹b
    # Thus we solve the least squares problem R⁻¹G = R⁻¹b where R is lower triangular matrix
    Ln <- as.array(L(p)$contiguous())
    Sn <- (1 - reg) * Ln %*% t(Ln) + alphaIdm
    RT <- chol(Sn) # Upper triangular matrix S = RRᵀ
    G_ <- forwardsolve(t(RT), G)
    b_ <- forwardsolve(t(RT), b)
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
  
  # The parameter covariance is estimated from the approxiamte distribution
  # b ~ N(Gŵ,σ̂^2 Ŝ) and ŵ = (GᵀG)^-1 Gᵀb 
  # so cov(ŵ) = σ̂^2 (GᵀG)^-1 GᵀŜG(GᵀG)^-1 
  GTG <- t(G) %*% G
  quad <- t(G) %*% Sn %*% G 
  param_covariance <- solve(GTG, t(solve(GTG, t(quad))))

  return(list(
    p = p,
    covp = param_covariance,
    iterations = n,
    converged = (relative_change <= tau_FP || SW <= tau_SW),
    relative_change_n = relative_change,
    sw_pvalues = sw_pvalues,
    final_sw_pvalue = SW
  ))
}

ols <- function(G, b, L, reg = 10e-10){
  dm <- nrow(G)
  alphaIdm <- reg * diag(rep(1, dm))
  p <- lm.fit(G, b)$coefficients
  residuals <- b - G %*% p
  sw_test <- shapiro.test(residuals)
  sw_p_value <- sw_test$p.value
  return(list(p = p, sw_p_value = sw_p_value))
}
