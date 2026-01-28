torch::torch_set_default_dtype(torch::torch_float64())

# F(p) in -𝚽'U = 𝚽F(p,U,t)
build_F <- function(U, tt, f_, J, device = torch::torch_device("cpu")) {
  mp1 <- nrow(U)
  D   <- ncol(U)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    Fp <- f_(input)
    return(torch::torch_tensor(Fp, dtype = torch::torch_float64(), device = device))
  }
}

# g(p) = 𝚽F(p,u,t)
build_g <- function(V, F_) {
  function(p) {
    torch::torch_matmul(V, F_(p))$reshape(c(-1))
  }
}

# G matrix in the linear system Gp = b - g0 : G ∈ ℝ^{K * D x J}
build_G_matrix <- function(V, U, tt, F_, J, device = torch::torch_device("cpu")){
  K <- nrow(V)
  mp1 <- nrow(U)
  D <- ncol(U)
  G <- matrix(0, nrow = K*D, ncol = J)
  g0 <- F_(rep(0,J))
  for(j in seq(1,J)){
    e_j <- rep(0, J)
    e_j[j] <- 1
    F_e <- F_(e_j) - g0
    result <- torch::torch_matmul(V, F_e)
    g_j <- torch::torch_reshape(result, c(-1))
    G[,j] <- as.numeric(g_j)
  }
  return(torch::torch_tensor(G, dtype = torch::torch_float64(), device = device))
}

# g(p) = Gp in system Gp = b - g0
build_g_linear <- function(G, device = torch::torch_device("cpu")) {
  function(p) {
    p <- torch::torch_tensor(p, dtype = torch::torch_float64(), device = device)
    torch::torch_matmul(G, p)$reshape(c(-1))
  }
}

# ∇ₚr(p) ∈ ℝ^(K*D × J) Jacobian of the weak residual
build_Jp_r <- function(J_p, K, D, J, mp1, V, U, tt, device = torch::torch_device("cpu")){
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_p_eval <- torch::torch_tensor(J_p(input), dtype = torch::torch_float64(), device = device)
    Jp_F <- torch::torch_reshape(J_p_eval, c(mp1, D, J))
    Jp_r <- torch::torch_einsum("km,mdj->kdj", list(V, Jp_F))
    Jp_r <- torch::torch_reshape(Jp_r, c(K * D, J))
  }
}

# ∇ₚr(p) ∈ ℝ^(K*D × J) Jacobian of the weak residual: r(p) = Gp + go - b → ∇ₚr(p) = G 
build_Jp_r_linear <- function(G){
  function(p){G}
}

# ∇ₚ∇ₚr(p) ∈ ℝ^(K*D × J × J) Hessian of the weak residual
build_Hp_r <- function(H_p, K, D, J, mp1, V, U, tt, device = torch::torch_device("cpu")){
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    H_p_eval <- torch::torch_tensor(H_p(input), dtype = torch::torch_float64(), device = device)
    Hp_F <- torch::torch_reshape(H_p_eval, c(mp1, D, J, J))
    Hp_r <- torch::torch_einsum("km,mdab->kdab", list(V, Hp_F))
    Hp_r <- torch::torch_reshape(Hp_r, c(K * D, J, J))
    return(Hp_r)
  }
}

# L₀ where L(p) = L₁(p) + L₀
build_L0 <- function(K, D, mp1, Vp, sig, device = torch::torch_device("cpu")) {
  sig_diag <- torch::torch_diag(sig)
  L0_ <- torch::torch_einsum('km,ab->kamb', list(Vp, sig_diag))
  L0 <- torch::torch_reshape(L0_, c(K * D, mp1 * D))
  return(L0)
}

# L(p) where L(p)Lᵀ(p) = S(p), the covariance matrix of the weak residual
build_L <-function(U, tt, J_u, K, V, L0, sig, J, device = torch::torch_device("cpu")){
  D <- ncol(U)
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_F <- torch::torch_tensor(J_u(input), dtype = torch::torch_float64(), device = device)
    J_F <- torch::torch_reshape(J_F, c(mp1, D, D))
    L1 <- torch::torch_einsum('mab,km,a->kamb', list(J_F, V, sig))
    L1 <- torch::torch_reshape(L1, c(K * D, mp1 * D))
    L <- L1 + L0
    return(L)
  }
}

# L for linear in parameters where LpLᵀp = S(p), the covariance matrix of the weak residual
build_L_linear <- function(U, tt, J_u, K, V, L0, sig, J, device = torch::torch_device("cpu")){
  D <- ncol(U)
  mp1 <- length(tt)
  L1_ <- torch::torch_zeros(c(K * D, mp1 * D, J), dtype = torch::torch_float64(), device = device)

  p_mat_affine <- matrix(rep(rep(0,J), mp1), ncol = mp1, nrow = J)
  input_affine <- rbind(p_mat_affine, t(U), t(tt))
  J_F_affine <- torch::torch_tensor(J_u(input_affine), dtype = torch::torch_float64(), device = device)
  J_F_affine <- torch::torch_reshape(J_F_affine, c(mp1, D, D))
  J_F_affine <- torch::torch_einsum('mab,km,a->kamb', list(J_F_affine, V, sig))
  L1_affine <- torch::torch_reshape(J_F_affine, c(K * D, mp1 * D))

  for(j in seq_len(J)){
    e_j <- rep(0, J)
    e_j[j] <- 1
    p_mat <- matrix(rep(e_j, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_F <- torch::torch_tensor(J_u(input), dtype = torch::torch_float64(), device = device)
    J_F <- torch::torch_reshape(J_F, c(mp1, D, D))
    Jj_F <- torch::torch_einsum('mab,km,a->kamb', list(J_F, V, sig))
    Jj_F <- torch::torch_reshape(Jj_F, c(K * D, mp1 * D))
    L1_[,,j] <- Jj_F - L1_affine
  }

  function(p){
    p <- torch::torch_tensor(p, dtype = torch::torch_float64(), device = device)
    L1 <- torch::torch_einsum('abj,j->ab', list(L1_, p))
    L <- L1 + L1_affine + L0
    return(L)
  }
}

# ∇ₚL(p) Jacobian of the covariance factor
build_Jp_L <-function(U, tt, J_up, K, J, D, V, sig, device = torch::torch_device("cpu")){
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    H_F <- torch::torch_tensor(J_up(input), dtype = torch::torch_float64(), device = device)
    H_F <- torch::torch_reshape(H_F, c(mp1, D, D, J))
    J_ <- torch::torch_einsum("mabj,km,a->kambj", list(H_F, V, sig))
    J_ <- torch::torch_reshape(J_, c(K * D, mp1 * D, J))
    return(J_)
  }
}

# ∇ₚL(p) Linear Jacobian of the Covariance factor L(p) = L₁p + L_affine + L0 → ∇ₚL(p) = L₁
build_Jp_L_linear <- function(U, tt, J_u, K, V, L0, sig, J, device = torch::torch_device("cpu")){
  D <- ncol(U)
  mp1 <- length(tt)
  L1_ <- torch::torch_zeros(c(K * D, mp1 * D, J), dtype = torch::torch_float64(), device = device)

  p_mat_affine <- matrix(rep(rep(0,J), mp1), ncol = mp1, nrow = J)
  input_affine <- rbind(p_mat_affine, t(U), t(tt))
  J_F_affine <- torch::torch_tensor(J_u(input_affine), dtype = torch::torch_float64(), device = device)
  J_F_affine <- torch::torch_reshape(J_F_affine, c(mp1, D, D))
  J_F_affine <- torch::torch_einsum('mab,km,a->kamb', list(J_F_affine, V, sig))
  L1_affine <- torch::torch_reshape(J_F_affine, c(K * D, mp1 * D))

  for(j in seq_len(J)){
    e_j <- rep(0, J)
    e_j[j] <- 1
    p_mat <- matrix(rep(e_j, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_F <- torch::torch_tensor(J_u(input), dtype = torch::torch_float64(), device = device)
    J_F <- torch::torch_reshape(J_F, c(mp1, D, D))
    Jj_F <- torch::torch_einsum('mab,km,a->kamb', list(J_F, V, sig))
    Jj_F <- torch::torch_reshape(Jj_F, c(K * D, mp1 * D))
    L1_[,,j] <- Jj_F  - L1_affine
  }
  return(function(p){L1_})
}

# ∇ₚ∇ₚL(p) Hessian of the covariance factor where ∇ₚ∇ₚS(p) = ∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ + (∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ)ᵀ
build_Hp_L <-function(U, tt, J_upp, K, J, D, V, sig, device = torch::torch_device("cpu")){
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    T_F <- torch::torch_tensor(J_upp(input), dtype = torch::torch_float64(), device = device)
    T_F <- torch::torch_reshape(T_F, c(mp1, D, D, J, J))
    H_ <- torch::torch_einsum("mabcd,km,a->kambcd", list(T_F, V, sig))
    H_ <- torch::torch_reshape(H_, c(K * D, mp1 * D, J, J))
    return(H_)
  }
}

# S(p) is the covariance matrix of the weak residual
build_S <- function(L, diag_reg = 1e-9) {
  function(p) {
    Lp <- L(p)
    Lpt <- torch::torch_transpose(Lp, 1, 2)
    S_ <- torch::torch_matmul(Lp, Lpt)
    WEIGHT <- 1.0 - diag_reg 
    eye <- diag_reg * diag(nrow(S_))
    S <- WEIGHT * S_ + eye
    return(S)
  }
}

# ∇ₚS(p) Jacobian of the covariance matrix
build_J_S <- function(L, Jp_L, J, K, D){
  function(p){
    Lp <- L(p)
    Jp_Lp <- Jp_L(p)
    Lp_t <- Lp$t()
    prt <- torch::torch_einsum("ijk,jl->ilk", list(Jp_Lp, Lp_t))
    Jp_S <- prt + prt$transpose(1, 2)
    return(Jp_S)
  }
}

# Weak negative log likelihood  (Multivariate Gaussian distribution)
build_wnll <- function(S, g, b, K, D){
  constant_term <- 0.5 * K * D * log(2 * pi)
  function(p){
    Sp <- as.array(S(p))
    r <- as.array(g(p) - b)
    cholF <- chol(Sp)
    log_det <- 2 * sum(log(diag(cholF)))
    S_invr <- backsolve(cholF, forwardsolve(t(cholF), r))
    mdist <- (r %*% S_invr)[1,1]
    return(0.5 * (mdist + log_det) + constant_term)
  }
}

# Jacobian of the weak form negative log likelihood 
build_J_wnll <- function(S, Jp_S, Jp_r, g, b, J){
  function(p){
    Sp <- as.array(S(p))

    J_Sp <- as.array(Jp_S(p)$contiguous())
    J_rp <- as.array(Jp_r(p)$contiguous())

    r <- as.array(g(p) - b)

    cholF <- chol(Sp)
    S_inv_rp <- backsolve(cholF, forwardsolve(t(cholF), r))

    gradient <- numeric(J)

    for(j in seq(J)){
      J_S_j <- J_Sp[, , j]
      J_r_j <- J_rp[, j]

      tmp <- J_S_j %*% S_inv_rp

      prt0 <- 2.0 * (J_r_j %*% S_inv_rp)[1,1]
      prt1 <- -1.0 * (S_inv_rp %*% tmp)[1,1]

      fact <- backsolve(cholF, forwardsolve(t(cholF), J_S_j))
      logDetPart <- sum(diag(fact))

      gradient[j] <- 0.5 * (prt0 + prt1 + logDetPart)
    }

    return(gradient)
  }
}

# Hessian of the weak form negative log likelihood 
build_H_wnll <- function(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J) {
  function(p) {
    r <- as.array(g(p) - b)

    Jp_rp <- as.array(Jp_r(p)$contiguous())
    Hp_rp <- as.array(Hp_r(p)$contiguous())

    Lp <- as.array(L(p))                     # L(p)
    Jp_Lp <- as.array(Jp_L(p)$contiguous())  # ∇ₚL(p)
    Hp_Lp <- as.array(Hp_L(p)$contiguous())  # ∇ₚ∇ₚL(p)

    Sp <- as.array(S(p)) # S(p)
    Jp_Sp <- as.array(Jp_S(p)$contiguous()) # ∇ₚS(p)

    S_inv_solve <- make_S_inv_solver(Sp)

    S_inv_rp <- S_inv_solve(r)
    H_wnn <- matrix(0, nrow = J, ncol = J)

    for (j in seq_len(J)) {
      Jp_rp_j <- Jp_rp[, j]
      Jp_Sp_j <- Jp_Sp[, , j]
      Jp_Lp_j <- Jp_Lp[, , j]

      shar_ <- S_inv_solve(Jp_Sp_j) # S⁻¹∂ⱼS

      for (i in j:J) {
        Jp_Sp_i <- Jp_Sp[, , i]  # ∂ᵢS(p)
        Jp_Lp_i <- Jp_Lp[, , i]  # ∂ᵢL(p)
        Jp_rp_i <- Jp_rp[, i]    # ∂ᵢr(p)

        term <- Jp_Sp_i %*% shar_ # ∂ᵢSS⁻¹∂ⱼS

        Hp_Lp_ji <- Hp_Lp[, , j, i] # ∂ᵢ∂ⱼS(p)
        p1 <- Jp_Lp_j %*% t(Jp_Lp_i)
        p2 <- Hp_Lp_ji %*% t(Lp)
        Hp_Sp_ji <- p1 + t(p1) + p2 + t(p2)  # ∂ᵢ∂ⱼS(p)

        Hp_rp_ji <- Hp_rp[, j, i]  # ∂ᵢ∂ⱼ r(p)

        prt0 <- as.numeric(t(Hp_rp_ji) %*% S_inv_rp)
        prt1 <- -1.0 * as.numeric(t(S_inv_solve(Jp_rp_j)) %*% (Jp_Sp_i %*% S_inv_rp))

        inv_factor <- S_inv_solve(Jp_rp_i)
        prt2 <- as.numeric(t(Jp_rp_j) %*% inv_factor)

        prt3 <- -2.0 * as.numeric(t(inv_factor) %*% (Jp_Sp_j %*% S_inv_rp))
        prt4 <- -1.0 * as.numeric(t(S_inv_rp) %*% (Hp_Sp_ji %*% S_inv_rp))
        prt5 <- 2 * as.numeric(t(S_inv_rp) %*% (term %*% S_inv_rp))

        logDetTerm <- -1.0 * sum(diag(S_inv_solve(term))) + sum(diag(S_inv_solve(Hp_Sp_ji)))

        Hij <- 0.5 * (2 * (prt0 + prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)

        H_wnn[j, i] <- Hij

        if (i != j) {
          H_wnn[i, j] <- Hij
        }
      }
    }
    return(H_wnn)
  }
}

# Hessian of the weak form negative log likelihood when linear in parameters
build_H_wnll_linear <- function(S, Jp_S, L, Jp_L, Jp_r, g, b, J) {
  function(p) {

    r <- as.array(g(p) - b)
    Jp_rp <- as.array(Jp_r(p)$contiguous())

    Lp <- as.array(L(p))                     # L(p)
    Jp_Lp <- as.array(Jp_L(p)$contiguous())  # ∇ₚL(p)

    Sp <- as.array(S(p)) # S(p)
    Jp_Sp <- as.array(Jp_S(p)$contiguous()) # ∇ₚS(p)

    S_inv_solve <- make_S_inv_solver(Sp)

    S_inv_rp <- S_inv_solve(r)
    H_wnn <- matrix(0, nrow = J, ncol = J)

    for (j in seq_len(J)) {
      Jp_rp_j <- Jp_rp[, j]
      Jp_Sp_j <- Jp_Sp[, , j]
      Jp_Lp_j <- Jp_Lp[, , j]

      shar_ <- S_inv_solve(Jp_Sp_j) # S⁻¹∂ⱼS

      for (i in j:J) {

        Jp_Sp_i <- Jp_Sp[, , i]  # ∂ᵢS(p)
        Jp_Lp_i <- Jp_Lp[, , i]  # ∂ᵢL(p)
        Jp_rp_i <- Jp_rp[, i]    # ∂ᵢr(p)

        term <- Jp_Sp_i %*% shar_ # ∂ᵢSS⁻¹∂ⱼS

        p1 <- Jp_Lp_j %*% t(Jp_Lp_i) # ∂ᵢ∂ⱼS(p)
        Hp_Sp_ji <- p1 + t(p1)  # ∂ᵢ∂ⱼS(p)

        inv_factor <- S_inv_solve(Jp_rp_i)

        prt1 <- -1.0 * as.numeric(t(S_inv_solve(Jp_rp_j)) %*% (Jp_Sp_i %*% S_inv_rp))
        prt2 <- as.numeric(t(Jp_rp_j) %*% inv_factor)
        prt3 <- -2.0 * as.numeric(t(inv_factor) %*% (Jp_Sp_j %*% S_inv_rp))
        prt4 <- -1.0 * as.numeric(t(S_inv_rp) %*% (Hp_Sp_ji %*% S_inv_rp))
        prt5 <- 2 * as.numeric(t(S_inv_rp) %*% (term %*% S_inv_rp))

        logDetTerm <- -1.0 * sum(diag(S_inv_solve(term))) + sum(diag(S_inv_solve(Hp_Sp_ji)))

        Hij <- 0.5 * (2 * (prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)

        H_wnn[j, i] <- Hij

        if (i != j) {
          H_wnn[i, j] <- Hij
        }
      }
    }
    return(H_wnn)
  }
}

# Robust solver for applying the inverse of S(p) to a vector or matrix
make_S_inv_solver <- function(Sp) {
  cholF <- NULL
  qrF <- NULL
  use_regularization <- FALSE
  diag_reg <- NULL

  cholF <- tryCatch(chol(Sp), error = function(e) NULL)

  if (is.null(cholF)) {
    qrF <- tryCatch(qr(Sp), error = function(e) NULL)

    if (is.null(qrF)) {
      use_regularization <- TRUE
      diag_reg <- 1e-12 * diag(nrow(Sp))
    }
  }

  function(x) {
    if (!is.null(cholF)) {
      return(backsolve(cholF, forwardsolve(t(cholF), x)))
    } else if (!is.null(qrF)) {
      return(solve(qrF, x))
    } else {
      return(solve(Sp + diag_reg, x))
    }
  }
}