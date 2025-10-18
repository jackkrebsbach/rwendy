
irls <- function(G, b1, L, reg = 10e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
  dm <- nrow(G)
  alphaIdm <- reg * diag(rep(1, dm))
  p <- lm(b1 ~ G + 0)$coefficients
  n <- 0
  SW <- Inf

  # Track Shapiro-Wilk p-values
  sw_pvalues <- numeric(max_its)

  while(n < max_its){
    pn1 <- p
    n <- n + 1

    Ln <- L(p)
    Sn <- (1 - reg) * Ln %*% t(Ln) + alphaIdm
    cholSn <- chol(Sn)
    S_invb <- solve(cholSn, solve(t(cholSn), b1))
    S_invG <- solve(cholSn, solve(t(cholSn), G))
    p <- solve(t(G) %*% S_invG, t(G) %*% S_invb)

    relative_change <- sqrt(sum((p - pn1)^2)) / sqrt(sum(pn1^2))

    residuals <- b1 - G %*% p

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


