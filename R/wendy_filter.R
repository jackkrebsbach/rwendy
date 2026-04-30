
# Estimate the initial condition using Euler-MacLaurin Trapezoid correction on the function g(t) := φ(t)f(p,u,t) + φ′(t)u(t)
estimate_u0 <- function(U, f_, J_u, J_t, tt, p, test_function_params){
  # Dimension of system
  mp1 <- nrow(U)
  D <- ncol(U)
  J <- length(p)

  # Build boundary layer test functions 
  data <- compute_r_c_hat(U, tt, test_function_params$S, test_function_params$p)
  radius_c <- data$rc
  diameter <- 2 * radius_c + 1
  dt <- mean(diff(tt))
  r <- dt * radius_c
  lin <- seq(-r, r, length.out = diameter)
  xx <- lin[2:(diameter-1)] # Don't include the endpoints (outside of support for 𝚽 or zero for 𝚿) 

  ph <- test_function_derivative(psi, radius_c, dt, order = 0)
  php <- test_function_derivative(psi, radius_c, dt, order = 1)
  phpp <- test_function_derivative(psi, radius_c, dt, order = 2)

  phxx <- ph(xx)
  phpxx <- php(xx)
  phppxx <- phpp(xx)
  scale <- sqrt(sum(phxx^2))

  p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
  input <- rbind(p_mat, t(U), matrix(tt, nrow = 1L))
  Fp <- f_(input)
  
  n_bl <- max(1L, floor(radius_c / 5))

  v_row <- c(0, phxx / scale, 0)
  vp_row <- c(0, phpxx / scale, 0)
  vpp_row <- c(0, phppxx / scale, 0)

  step <- max(1L, floor(radius_c / n_bl) - 2)

  V <- matrix(0, nrow = n_bl, ncol = mp1)
  Vp <- matrix(0, nrow = n_bl, ncol = mp1)
  Vpp <- matrix(0, nrow = n_bl, ncol = mp1)

  for (i in seq(nrow(V))) {
    shift <- (i - 1) * step + radius_c - 2 * step + 1
    
    l <- length(v_row)

    v_trimmed <- v_row[shift:l]
    vp_trimmed <- vp_row[shift:l]
    vpp_trimmed <- vpp_row[shift:l]

    s <- length(v_trimmed)

    if (s > 0) {
      V[i, 1:s] <- v_trimmed
      Vp[i, 1:s] <- vp_trimmed
      Vpp[i, 1:s] <- vpp_trimmed
    }
  }

  # Solve uncorrected system -ϕ(0)u(0) = ∫g(t)dt for u0 and k boundary layer functions
  B <- V[, 1] # Boundary layer test functions
  Q <- diag(c(0.5, rep(1, mp1 - 2), 0.5) * dt) # Quadrature matrix for Trapezoidal Rule
  residual <- -(V %*% Q %*% Fp  + Vp %*% Q %*% U) # -(Weak residual)
  u00 <- solve(crossprod(B), t(B) %*% residual) # Least squares solution to over determined system

  # Iterate with the Euler-Maclaurin formula
  un0 <- u00
  t0 <- tt[1]

  for (iter in seq_len(5)) {
    input  <- c(p, un0, t0)
    f_val  <- matrix(f_(input), ncol = 1)
    Ju_val <- matrix(as.vector(J_u(input)), nrow = D, ncol = D)
    Jt_val <- as.vector(J_t(input))
    
    # g(t) := φ(t)f(p,u,t) + φ′(t)u(t)
    # g'(t) = 2φ′(t)f(p,u,t) + φ(t)( ∇ᵤf(p,u,t) ⋅ f(p,u,t) + ∇ₜf(p,u,t)) + φ′′(t)u(t) 
    gprime0 <- matrix(t(vapply(seq(length(B)), function(row){
          return(2 * Vp[row, 1] * f_val + V[row, 1] * (as.vector(Ju_val %*% f_val) + Jt_val) + Vpp[row,1] * as.vector(un0))
    }, numeric(D))), nrow = length(B), ncol = D)

    # EM Correction is h^2 /12 (g'(t0) - g'(tf)) but g'(tf) = 0 to estimate initial condition
    un0 <- solve(crossprod(B), t(B) %*% (residual - dt^2 / 12  * gprime0)) # Least squares solution to over determined system
  }

  return(as.numeric(un0))
}
