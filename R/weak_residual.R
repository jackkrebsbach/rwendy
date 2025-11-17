
build_F <- function(U, tt, f_, J) {
  mp1 <- nrow(U)
  D   <- ncol(U)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    Fp <- f_(input)
    return(torch::torch_tensor(Fp, dtype = torch::torch_float64()))
  }
}


build_g <- function(V, F_) {
  function(p) {
    torch::torch_matmul(V, F_(p))$reshape(-1)
  }
}

build_G_matrix <- function(V, U, tt, F_, J){
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
  return(torch::torch_tensor(G))
}

# ∇ₚ∇ₚL(p) Hessian of the Covariance factor where ∇ₚ∇ₚS(p) = ∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ + (∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ)ᵀ
build_Hp_L <-function(U, tt, J_upp, K, J, D, V, sig){
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    T_F <- torch::torch_tensor(J_upp(input), dtype = torch::torch_float64())
    T_F <- torch::torch_reshape(T_F, c(mp1, D, D, J, J))
    H_ <- torch::torch_einsum("mabcd,km,a->kambcd", list(T_F, V, sig))
    H_ <- torch::torch_reshape(H_, c(K * D, mp1 * D, J, J))
    return(H_)
  }
}

# ∇ₚL(p) Jacobian of the Covariance factor
build_Jp_L <-function(U, tt, J_up, K, J, D, V, sig){
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    H_F <- torch::torch_tensor(J_up(input), dtype = torch::torch_float64())
    H_F <- torch::torch_reshape(H_F, c(mp1, D, D, J))
    J_ <- torch::torch_einsum("mabj,km,a->kambj", list(H_F, V, sig))
    J_ <- torch::torch_reshape(J_, c(K * D, mp1 * D, J))
    return(J_)
  }
}

build_L0 <- function(K, D, mp1, Vp, sig) {
  sig_diag <- torch::torch_diag(sig)
  L0_ <- torch::torch_einsum('km,ab->kabm', list(Vp, sig_diag))$permute(c(1,2,4,3))
  L0 <- torch::torch_reshape(L0_, c(K * D, mp1 * D))
  return(L0)
}

build_L <-function(U, tt, J_u, K, V, L0, sig, J){
  D <- ncol(U)
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_F <- torch::torch_tensor(J_u(input), dtype = torch::torch_float64())
    J_F <- torch::torch_reshape(J_F, c(mp1, D, D))
    L1 <- torch::torch_einsum('mab,km,a->kamb', list(J_F, V, sig))
    L1 <- torch::torch_reshape(L1, c(K * D, mp1 * D))
    L <- L1 + L0
    return(L)
  }
}

build_S <- function(L, REG = 1e-9) {
  function(p) {
    Lp <- L(p)
    Lpt <- torch::torch_transpose(Lp, 1, 2)
    S_ <- torch::torch_matmul(Lp, Lpt)
    WEIGHT <- 1.0 - REG
    eye <- REG * diag(nrow(S_))
    S <- WEIGHT * S_ + eye
    return(S)
  }
}

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

# Jₚr(p) ∈ ℝ^(K*D x J)
build_Jp_r <- function(J_p, K, D, J, mp1, V, U, tt){
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_p_eval <- torch::torch_tensor(J_p(input),dtype = torch::torch_float64())
    Jp_F <- torch::torch_reshape(J_p_eval, c(mp1, D, J))
    Jp_r <- torch::torch_einsum("km,mdj->kdj", list(V, Jp_F))
    Jp_r <- torch::torch_reshape(Jp_r, c(K * D, J))
  }
}

# Hₚr(p) ∈ ℝ^(K*D x J x J)
build_Hp_r <- function(H_p, K, D, J, mp1, V, U, tt){
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    H_p_eval <- torch::torch_tensor(H_p(input), dtype = torch::torch_float64())
    Hp_F <- torch::torch_reshape(H_p_eval, c(mp1, D, J, J))
    Hp_r <- torch::torch_einsum("km,mdab->kdab", list(V, Hp_F))
    Hp_r <- torch::torch_reshape(Hp_r, c(K * D, J, J))
    return(Hp_r)
  }
}

build_wnll <- function(S, g, b, K, D){
  constant_term <- 0.5 * K * D * log(2 * pi)
  function(p){
    Sp <- as.array(S(p))
    r <- as.array(g(p) - b)
    cholF <- chol(Sp)
    #svdF <- svd(Sp)
    #log_det <- sum(log(svdF$d))
    log_det <- 2 * sum(log(diag(cholF)))

    S_invr <- solve(cholF, solve(t(cholF), r))

    mdist <- (r %*% S_invr)[1,1]
    return(0.5 * (mdist + log_det) + constant_term)
  }
}

build_J_wnll <- function(S, Jp_S, Jp_r, g, b, J){
  function(p){
    Sp <- as.array(S(p))

    J_Sp <- as.array(Jp_S(p)$contiguous())
    J_rp <- as.array(Jp_r(p)$contiguous())

    r <- as.array(g(p) - b)

    cholF <- chol(Sp)
    S_inv_rp <- solve(cholF, solve(t(cholF), r))

    gradient <- numeric(J)

    for(j in seq(J)){
      J_S_j <- J_Sp[, , j]
      J_r_j <- J_rp[, j]

      tmp <- J_S_j %*% S_inv_rp

      prt0 <- 2.0 * (J_r_j %*% S_inv_rp)[1,1]
      prt1 <- -1.0 * (S_inv_rp %*% tmp)[1,1]

      fact <- solve(cholF, solve(t(cholF), J_S_j))
      logDetPart <- sum(diag(fact))

      gradient[j] <- 0.5 * (prt0 + prt1 + logDetPart)
    }

    return(gradient)
  }
}

build_H_wnll <- function(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J) {
  function(p) {

    r <- as.array(g(p) - b)

    Jp_rp <- as.array(Jp_r(p)$contiguous())
    Hp_rp <- as.array(Hp_r(p)$contiguous())

    Lp <- as.array(L(p))                     # L(p)
    Jp_Lp <- as.array(Jp_L(p)$contiguous())  # ∇ₚL(p)
    Hp_Lp <- as.array(Hp_L(p)$contiguous())  # ∇ₚ∇ₚL(p)

    Sp <- as.array(S(p))
    Jp_Sp <- as.array(Jp_S(p)$contiguous())

    S_inv_solve <- make_S_inv_solver(Sp)

    S_inv_rp <- S_inv_solve(r)
    H_wnn <- matrix(0, nrow = J, ncol = J)

    for (j in seq_len(J)) {
      Jp_rp_j <- Jp_rp[, j]
      Jp_Sp_j <- Jp_Sp[, , j]
      Jp_Lp_j <- Jp_Lp[, , j]

      shar_ <- S_inv_solve(Jp_Sp_j) # S⁻¹∂ⱼS

      for (i in j:J) {
        # ∂ᵢS(p) (Jacobian information)
        Jp_Sp_i <- Jp_Sp[, , i]
        Jp_Lp_i <- Jp_Lp[, , i]  # ∂ᵢL(p)
        Jp_rp_i <- Jp_rp[, i]    # ∂ᵢr(p)

        term <- Jp_Sp_i %*% shar_ # ∂ᵢSS⁻¹∂ⱼS

        # ∂ᵢ∂ⱼS(p)
        Hp_Lp_ji <- Hp_Lp[, , j, i]
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
