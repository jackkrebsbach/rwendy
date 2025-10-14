# Iterative re-weighted least squares (IRLS)
irls <- function(G, b1, L, reg = 10e-10, tau_FP = 1e-6, tau_SW = 1e-4, n0 = 10, max_its = 100){
  dm <- nrow(G)
  alphaIdm <- reg * diag(rep(1, dm))
  p <- lm(b1 ~ G + 0)$coefficients

  n <- 0
  SW <- Inf

  while(n < max_its){
    pn1 <- p
    n <- n + 1

    # IRLS update
    Ln <- L(p)
    Sn <- (1 - reg) * Ln %*% t(Ln) + alphaIdm
    cholSn <- chol(Sn)
    S_invb <- solve(cholSn, solve(t(cholSn), b1))
    S_invG <- solve(cholSn, solve(t(cholSn), G))
    p <- solve(t(G) %*% S_invG, t(G) %*% S_invb)

    relative_change <- sqrt(sum((p - pn1)^2)) / sqrt(sum(pn1^2))

    # Update SW statistic after n0 iterations
    if(n >= n0){
      SW <- max(SW, relative_change)
    }

    if(relative_change > tau_FP && n < max_its && SW > tau_SW){
      next
    } else {
      break
    }
  }

  return(list(p = p, iterations = n, converged = (relative_change <= tau_FP || SW <= tau_SW)))
}
