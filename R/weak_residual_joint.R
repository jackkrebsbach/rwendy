# F(p) in -𝚽'U = 𝚽F(p,U,t)
build_F_joint <- function(f_, J, device = torch::torch_device("cpu")) {
  dtype <- torch::torch_float64()
  function(U, p, tt) {
    mp1   <- nrow(U)
    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    Fp <- f_(input)
    return(torch::torch_tensor(Fp, dtype = dtype, device = device))
  }
}

# g(p) = 𝚽F(p,u,t)
#
# When bl is supplied, the augmented residual contribution
#   g_aug(U, p, tt) = ΦF + Φ[,1]⊗U[1,] - Φ[,M]⊗U[M,] + EM(U, p, tt)
# is returned instead. Interior rows of bl$phi_t1 / bl$phi_tM are zero so the
# boundary-term outer products only fire on BL rows; EM is zero-padded above.
build_g_joint <- function(V, F_, bl = NULL) {
  if (is.null(bl)) {
    return(function(U, p, tt) {
      torch::torch_matmul(V, F_(U, p, tt))$reshape(c(-1))
    })
  }

  dtype  <- V$dtype
  device <- V$device
  K_int  <- bl$K_interior
  K_bl   <- bl$K_bl

  function(U, p, tt) {
    D     <- ncol(U)
    M     <- nrow(U)
    VF    <- torch::torch_matmul(V, F_(U, p, tt))                  # K × D
    U_t1  <- torch::torch_tensor(as.numeric(U[1L, ]), dtype = dtype, device = device)
    U_tM  <- torch::torch_tensor(as.numeric(U[M,  ]), dtype = dtype, device = device)
    bdry  <- torch::torch_outer(bl$phi_t1, U_t1) -
             torch::torch_outer(bl$phi_tM, U_tM)                   # K × D
    result <- VF + bdry

    if (!is.null(bl$em_correction)) {
      em      <- bl$em_correction(U, p, tt)                        # K_bl × D
      em_full <- torch::torch_cat(
        list(torch::torch_zeros(c(K_int, D), dtype = dtype, device = device), em),
        dim = 1L
      )
      result <- result + em_full
    }

    result$reshape(c(-1))
  }
}

# G matrix in the linear system Gp = b - g0 : G ∈ ℝ^{K * D x J}
build_G_matrix_joint <- function(V, U, tt, F_, J, device = torch::torch_device("cpu")){
  K   <- nrow(V)
  mp1 <- nrow(U)
  D   <- ncol(U)
  G   <- matrix(0, nrow = K*D, ncol = J)
  g0  <- F_(U, rep(0, J), tt)
  for(j in seq(1, J)){
    e_j <- rep(0, J)
    e_j[j] <- 1
    F_e    <- F_(U, e_j, tt) - g0
    result <- torch::torch_matmul(V, F_e)
    g_j    <- torch::torch_reshape(result, c(-1))
    G[,j]  <- as.numeric(g_j)
  }
  return(torch::torch_tensor(G, dtype = torch::torch_float64(), device = device))
}

# g(p) = Gp in system Gp = b - g0
build_g_linear_joint <- function(G, device = torch::torch_device("cpu")) {
  dtype <- torch::torch_float64()
  function(U, p, tt) {
    p <- torch::torch_tensor(p, dtype = dtype, device = device)
    torch::torch_matmul(G, p)$reshape(c(-1))
  }
}

# ∇ₚr(p) ∈ ℝ^(K*D × J) Jacobian of the weak residual
build_Jp_r_joint <- function(J_p, K, J, V, device = torch::torch_device("cpu")){
  dtype <- torch::torch_float64()
  function(U, p, tt){
    mp1   <- nrow(U)
    D     <- ncol(U)
    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    J_p_eval <- torch::torch_tensor(J_p(input), dtype = dtype, device = device)
    Jp_F <- torch::torch_reshape(J_p_eval, c(mp1, D, J))
    Jp_r <- torch::torch_einsum("km,mdj->kdj", list(V, Jp_F))
    Jp_r <- torch::torch_reshape(Jp_r, c(K * D, J))
  }
}

# ∇ₚr(p) ∈ ℝ^(K*D × J) Jacobian of the weak residual: r(p) = Gp + go - b → ∇ₚr(p) = G
build_Jp_r_linear_joint <- function(G){
  function(U, p, tt){G}
}

# ∇ₚ∇ₚr(p) ∈ ℝ^(K*D × J × J) Hessian of the weak residual
build_Hp_r_joint <- function(H_p, K, J, V, device = torch::torch_device("cpu")){
  dtype <- torch::torch_float64()
  function(U, p, tt){
    mp1   <- nrow(U)
    D     <- ncol(U)
    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    H_p_eval <- torch::torch_tensor(H_p(input), dtype = dtype, device = device)
    Hp_F <- torch::torch_reshape(H_p_eval, c(mp1, D, J, J))
    Hp_r <- torch::torch_einsum("km,mdab->kdab", list(V, Hp_F))
    Hp_r <- torch::torch_reshape(Hp_r, c(K * D, J, J))
    return(Hp_r)
  }
}

# L₀(sig) where L(p,sig) = L₁(p,sig) + L₀(sig)
# Returns a function(sig) -> tensor rather than a fixed tensor so sig can be optimised.
build_L0_joint <- function(K, D, mp1, Vp, device = torch::torch_device("cpu")) {
  function(sig) {
    sig_diag <- torch::torch_diag(sig)
    L0_ <- torch::torch_einsum('km,ab->kamb', list(Vp, sig_diag))
    torch::torch_reshape(L0_, c(K * D, mp1 * D))
  }
}

# L(p, sig) where L(p,sig)Lᵀ(p,sig) = S(p,sig), the covariance matrix of the weak residual
# L0_fn: function(sig) -> L0 tensor (from build_L0_joint).
build_L_joint <- function(J_u, K, V, L0_fn, J, device = torch::torch_device("cpu")){
  dtype <- torch::torch_float64()
  function(U, p, tt, sig){
    D     <- ncol(U)
    mp1   <- length(tt)
    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    J_F <- torch::torch_tensor(J_u(input), dtype = dtype, device = device)
    J_F <- torch::torch_reshape(J_F, c(mp1, D, D))
    L1 <- torch::torch_einsum('mab,km,a->kamb', list(J_F, V, sig))
    L1 <- torch::torch_reshape(L1, c(K * D, mp1 * D))
    L0 <- L0_fn(sig)
    L <- L1 + L0
    return(L)
  }
}

# L for linear in parameters where LpLᵀp = S(p,sig), the covariance matrix of the weak residual
# L0_fn: function(sig) -> L0 tensor (from build_L0_joint).
build_L_linear_joint <- function(J_u, K, V, L0_fn, J, device = torch::torch_device("cpu")){
  dtype <- torch::torch_float64()
  function(U, p, tt, sig){
    D   <- ncol(U)
    mp1 <- length(tt)
    L1_ <- torch::torch_zeros(c(K * D, mp1 * D, J), dtype = dtype, device = device)

    p_mat_affine <- matrix(rep(rep(0, J), mp1), ncol = mp1, nrow = J)
    input_affine <- rbind(p_mat_affine, t(U), matrix(tt, nrow = 1L))
    J_F_affine <- torch::torch_tensor(J_u(input_affine), dtype = dtype, device = device)
    J_F_affine <- torch::torch_reshape(J_F_affine, c(mp1, D, D))
    J_F_affine <- torch::torch_einsum('mab,km,a->kamb', list(J_F_affine, V, sig))
    L1_affine  <- torch::torch_reshape(J_F_affine, c(K * D, mp1 * D))

    for(j in seq_len(J)){
      e_j   <- rep(0, J)
      e_j[j] <- 1
      p_mat <- matrix(rep(e_j, mp1), ncol = mp1, nrow = J)
      input <- rbind(p_mat, t(U), matrix(tt, nrow = 1L))
      J_F   <- torch::torch_tensor(J_u(input), dtype = dtype, device = device)
      J_F   <- torch::torch_reshape(J_F, c(mp1, D, D))
      Jj_F  <- torch::torch_einsum('mab,km,a->kamb', list(J_F, V, sig))
      Jj_F  <- torch::torch_reshape(Jj_F, c(K * D, mp1 * D))
      L1_[,,j] <- Jj_F - L1_affine
    }

    p_t <- torch::torch_tensor(p, dtype = dtype, device = device)
    L1  <- torch::torch_einsum('abj,j->ab', list(L1_, p_t))
    L0  <- L0_fn(sig)
    L   <- L1 + L1_affine + L0
    return(L)
  }
}

# ∇ₚL(p, sig) Jacobian of the covariance factor w.r.t. p
build_Jp_L_joint <- function(J_up, K, J, V, device = torch::torch_device("cpu")){
  dtype <- torch::torch_float64()
  function(U, p, tt, sig){
    D     <- ncol(U)
    mp1   <- length(tt)
    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    H_F <- torch::torch_tensor(J_up(input), dtype = dtype, device = device)
    H_F <- torch::torch_reshape(H_F, c(mp1, D, D, J))
    J_  <- torch::torch_einsum("mabj,km,a->kambj", list(H_F, V, sig))
    J_  <- torch::torch_reshape(J_, c(K * D, mp1 * D, J))
    return(J_)
  }
}

# ∇ₚL(p, sig) Linear Jacobian of the Covariance factor L(p,sig) = L₁p + L_affine + L0 → ∇ₚL = L₁
build_Jp_L_linear_joint <- function(J_u, K, V, J, device = torch::torch_device("cpu")){
  dtype <- torch::torch_float64()
  function(U, p, tt, sig){
    D   <- ncol(U)
    mp1 <- length(tt)
    L1_ <- torch::torch_zeros(c(K * D, mp1 * D, J), dtype = dtype, device = device)

    p_mat_affine <- matrix(rep(rep(0, J), mp1), ncol = mp1, nrow = J)
    input_affine <- rbind(p_mat_affine, t(U), matrix(tt, nrow = 1L))
    J_F_affine <- torch::torch_tensor(J_u(input_affine), dtype = dtype, device = device)
    J_F_affine <- torch::torch_reshape(J_F_affine, c(mp1, D, D))
    J_F_affine <- torch::torch_einsum('mab,km,a->kamb', list(J_F_affine, V, sig))
    L1_affine  <- torch::torch_reshape(J_F_affine, c(K * D, mp1 * D))

    for(j in seq_len(J)){
      e_j    <- rep(0, J)
      e_j[j] <- 1
      p_mat  <- matrix(rep(e_j, mp1), ncol = mp1, nrow = J)
      input  <- rbind(p_mat, t(U), matrix(tt, nrow = 1L))
      J_F    <- torch::torch_tensor(J_u(input), dtype = dtype, device = device)
      J_F    <- torch::torch_reshape(J_F, c(mp1, D, D))
      Jj_F   <- torch::torch_einsum('mab,km,a->kamb', list(J_F, V, sig))
      Jj_F   <- torch::torch_reshape(Jj_F, c(K * D, mp1 * D))
      L1_[,,j] <- Jj_F - L1_affine
    }
    return(L1_)
  }
}

# ∇ₚ∇ₚL(p, sig) Hessian of the covariance factor
build_Hp_L_joint <- function(J_upp, K, J, V, device = torch::torch_device("cpu")){
  dtype <- torch::torch_float64()
  function(U, p, tt, sig){
    D     <- ncol(U)
    mp1   <- length(tt)
    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    T_F <- torch::torch_tensor(J_upp(input), dtype = dtype, device = device)
    T_F <- torch::torch_reshape(T_F, c(mp1, D, D, J, J))
    H_  <- torch::torch_einsum("mabcd,km,a->kambcd", list(T_F, V, sig))
    H_  <- torch::torch_reshape(H_, c(K * D, mp1 * D, J, J))
    return(H_)
  }
}

# S(p, sig) is the covariance matrix of the weak residual
build_S_joint <- function(L, W = NULL, diag_reg = 10e-10) {
  WEIGHT    <- 1.0 - diag_reg
  eye_cache <- NULL
  function(U, p, tt, sig) {
    Lp  <- L(U, p, tt, sig)
    Lpt <- torch::torch_transpose(Lp, 1, 2)
    S_  <- if (!is.null(W)) torch::torch_matmul(Lp, torch::torch_matmul(W, Lpt))
           else              torch::torch_matmul(Lp, Lpt)
    if (is.null(eye_cache)) {
      eye_cache <<- diag_reg * torch::torch_eye(S_$size(1L), dtype = S_$dtype, device = S_$device)
    }
    return(WEIGHT * S_ + eye_cache)
  }
}

# ∇ₚS(p, sig) Jacobian of the covariance matrix w.r.t. p
build_J_S_joint <- function(L, Jp_L, J, W = NULL, diag_reg = 1e-10){
  WEIGHT <- 1.0 - diag_reg
  function(U, p, tt, sig){
    Lp    <- L(U, p, tt, sig)
    Jp_Lp <- Jp_L(U, p, tt, sig)
    # With W: ∂S/∂pₖ = (∂L/∂pₖ) W Lᵀ + L W (∂L/∂pₖ)ᵀ  →  use W Lᵀ in einsum
    WLp_t <- if (!is.null(W)) torch::torch_mm(W, Lp$t()) else Lp$t()
    prt   <- torch::torch_einsum("ijk,jl->ilk", list(Jp_Lp, WLp_t))
    Jp_S  <- WEIGHT * (prt + prt$transpose(1, 2))
    return(Jp_S)
  }
}

# Weak negative log likelihood  (Multivariate Gaussian distribution)
build_wnll_joint <- function(S, g, K, D){
  constant_term <- 0.5 * K * D * log(2 * pi)
  function(U, p, tt, b, sig){
    Sp    <- S(U, p, tt, sig)
    r     <- g(U, p, tt) - b
    cholL <- torch::linalg_cholesky(Sp)
    log_det <- 2.0 * as.numeric(torch::torch_diag(cholL)$log()$sum())
    S_invr  <- torch::torch_cholesky_solve(r$unsqueeze(-1L), cholL, upper = FALSE)$squeeze(-1L)
    mdist   <- as.numeric(r$dot(S_invr))
    return(0.5 * (mdist + log_det) + constant_term)
  }
}

# Jacobian of the weak form negative log likelihood w.r.t. p
build_J_wnll_joint <- function(S, Jp_S, Jp_r, g, J){
  function(U, p, tt, b, sig){
    Sp   <- S(U, p, tt, sig)
    J_Sp <- Jp_S(U, p, tt, sig)   # (KD, KD, J)
    J_rp <- Jp_r(U, p, tt)        # (KD, J)
    r    <- g(U, p, tt) - b       # (KD,)
    KD   <- Sp$size(1L)

    cholL   <- torch::linalg_cholesky(Sp)
    S_inv_r <- torch::torch_cholesky_solve(r$unsqueeze(-1L), cholL, upper = FALSE)$squeeze(-1L)

    # prt0: 2 * (∂ⱼr)ᵀ S⁻¹r  for all j → (J,)
    prt0_all <- 2.0 * torch::torch_mv(J_rp$t()$contiguous(), S_inv_r)

    # tmp_all: (∂ⱼS)(S⁻¹r) for all j → (KD, J)
    tmp_all  <- torch::torch_einsum("klj,l->kj", list(J_Sp, S_inv_r))
    # prt1: -(S⁻¹r)ᵀ tmp_j for all j → (J,)
    prt1_all <- -1.0 * torch::torch_mv(tmp_all$t()$contiguous(), S_inv_r)

    # tr(S⁻¹ ∂ⱼS) for all j: solve KD × (KD*J) system, reshape, extract diagonals
    J_Sp_2d    <- J_Sp$reshape(c(KD, KD * J))
    J_S_sol_2d <- torch::torch_cholesky_solve(J_Sp_2d, cholL, upper = FALSE)
    logDet_all <- torch::torch_einsum("iij->j", J_S_sol_2d$reshape(c(KD, KD, J)))

    return(as.numeric(0.5 * (prt0_all + prt1_all + logDet_all)))
  }
}

# Hessian of the weak form negative log likelihood w.r.t. p
build_H_wnll_joint <- function(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, J, W = NULL, diag_reg = 1e-10) {
  WEIGHT <- 1.0 - diag_reg
  function(U, p, tt, b, sig) {
    r      <- g(U, p, tt) - b         # (KD,)
    Jp_rp  <- Jp_r(U, p, tt)          # (KD, J)
    Hp_rp  <- Hp_r(U, p, tt)          # (KD, J, J)
    Lp     <- L(U, p, tt, sig)        # (KD, mp1*D)
    Jp_Lp  <- Jp_L(U, p, tt, sig)    # (KD, mp1*D, J)
    Hp_Lp  <- Hp_L(U, p, tt, sig)    # (KD, mp1*D, J, J)
    Sp     <- S(U, p, tt, sig)        # (KD, KD)
    Jp_Sp  <- Jp_S(U, p, tt, sig)    # (KD, KD, J)

    # Precompute W Lᵀ once (used in p2 = ∂ᵢ∂ⱼL · W · Lᵀ)
    LpW_t <- if (!is.null(W)) torch::torch_mm(W, Lp$t()$contiguous()) else Lp$t()$contiguous()

    S_inv_solve <- make_S_inv_solver(Sp)
    S_inv_r     <- S_inv_solve(r)  # (KD,)

    H_vals <- torch::torch_zeros(c(J, J), dtype = Sp$dtype, device = Sp$device)

    for (j in seq_len(J)) {
      Jp_rp_j   <- Jp_rp[, j]
      Jp_Sp_j   <- Jp_Sp[, , j]$contiguous()   # (KD, KD)
      Jp_Lp_j   <- Jp_Lp[, , j]$contiguous()   # (KD, mp1*D)
      # Precompute ∂ⱼL · W once per outer iteration
      Jp_Lp_j_W <- if (!is.null(W)) torch::torch_mm(Jp_Lp_j, W) else Jp_Lp_j

      shar_          <- S_inv_solve(Jp_Sp_j)                        # S⁻¹∂ⱼS: (KD, KD)
      S_inv_rp_j     <- S_inv_solve(Jp_rp_j)                        # S⁻¹∂ⱼr: (KD,)
      Jp_Sp_j_S_invr <- torch::torch_matmul(Jp_Sp_j, S_inv_r)      # (∂ⱼS)(S⁻¹r): (KD,)

      for (i in j:J) {
        Jp_Sp_i  <- Jp_Sp[, , i]$contiguous()   # (KD, KD)
        Jp_Lp_i  <- Jp_Lp[, , i]$contiguous()   # (KD, mp1*D)
        Jp_rp_i  <- Jp_rp[, i]                    # (KD,)

        term     <- torch::torch_mm(Jp_Sp_i, shar_)  # ∂ᵢS S⁻¹∂ⱼS: (KD, KD)

        Hp_Lp_ji <- Hp_Lp[, , j, i]$contiguous()  # (KD, mp1*D)
        # p1 = ∂ⱼL · W · ∂ᵢLᵀ,  p2 = ∂ᵢ∂ⱼL · W · Lᵀ
        p1       <- torch::torch_mm(Jp_Lp_j_W, Jp_Lp_i$t()$contiguous())  # (KD, KD)
        p2       <- torch::torch_mm(Hp_Lp_ji, LpW_t)                        # (KD, KD)
        Hp_Sp_ji <- WEIGHT * (p1 + p1$t() + p2 + p2$t())

        Hp_rp_ji   <- Hp_rp[, j, i]  # (KD,)
        inv_factor <- S_inv_solve(Jp_rp_i)  # (KD,)

        prt0 <- Hp_rp_ji$dot(S_inv_r)
        prt1 <- -1.0 * S_inv_rp_j$dot(torch::torch_matmul(Jp_Sp_i, S_inv_r))
        prt2 <- Jp_rp_j$dot(inv_factor)
        prt3 <- -2.0 * inv_factor$dot(Jp_Sp_j_S_invr)
        prt4 <- -1.0 * S_inv_r$dot(torch::torch_matmul(Hp_Sp_ji, S_inv_r))
        prt5 <- 2.0 * S_inv_r$dot(torch::torch_matmul(term, S_inv_r))

        logDetTerm <- -1.0 * torch::torch_diag(S_inv_solve(term))$sum() +
                             torch::torch_diag(S_inv_solve(Hp_Sp_ji))$sum()

        Hij <- 0.5 * (2 * (prt0 + prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)

        H_vals[j, i] <- Hij
        if (i != j) H_vals[i, j] <- Hij
      }
    }
    return(as.matrix(H_vals$cpu()))
  }
}

# Hessian of the weak form negative log likelihood when linear in parameters
build_H_wnll_linear_joint <- function(S, Jp_S, L, Jp_L, Jp_r, g, J, W = NULL, diag_reg = 1e-10) {
  WEIGHT <- 1.0 - diag_reg
  function(U, p, tt, b, sig) {
    r      <- g(U, p, tt) - b         # (KD,)
    Jp_rp  <- Jp_r(U, p, tt)          # (KD, J)
    Jp_Lp  <- Jp_L(U, p, tt, sig)    # (KD, mp1*D, J)
    Sp     <- S(U, p, tt, sig)        # (KD, KD)
    Jp_Sp  <- Jp_S(U, p, tt, sig)    # (KD, KD, J)

    S_inv_solve <- make_S_inv_solver(Sp)
    S_inv_r     <- S_inv_solve(r)  # (KD,)

    H_vals <- torch::torch_zeros(c(J, J), dtype = Sp$dtype, device = Sp$device)

    for (j in seq_len(J)) {
      Jp_rp_j   <- Jp_rp[, j]
      Jp_Sp_j   <- Jp_Sp[, , j]$contiguous()   # (KD, KD)
      Jp_Lp_j   <- Jp_Lp[, , j]$contiguous()   # (KD, mp1*D)
      Jp_Lp_j_W <- if (!is.null(W)) torch::torch_mm(Jp_Lp_j, W) else Jp_Lp_j

      shar_          <- S_inv_solve(Jp_Sp_j)                        # S⁻¹∂ⱼS: (KD, KD)
      S_inv_rp_j     <- S_inv_solve(Jp_rp_j)                        # S⁻¹∂ⱼr: (KD,)
      Jp_Sp_j_S_invr <- torch::torch_matmul(Jp_Sp_j, S_inv_r)      # (∂ⱼS)(S⁻¹r): (KD,)

      for (i in j:J) {
        Jp_Sp_i <- Jp_Sp[, , i]$contiguous()   # (KD, KD)
        Jp_Lp_i <- Jp_Lp[, , i]$contiguous()   # (KD, mp1*D)
        Jp_rp_i <- Jp_rp[, i]                    # (KD,)

        term <- torch::torch_mm(Jp_Sp_i, shar_)  # ∂ᵢS S⁻¹∂ⱼS: (KD, KD)

        # p1 = ∂ⱼL · W · ∂ᵢLᵀ (Hp_L = 0 in linear case, so no p2 term)
        p1       <- torch::torch_mm(Jp_Lp_j_W, Jp_Lp_i$t()$contiguous())  # (KD, KD)
        Hp_Sp_ji <- WEIGHT * (p1 + p1$t())

        inv_factor <- S_inv_solve(Jp_rp_i)  # (KD,)

        prt1 <- -1.0 * S_inv_rp_j$dot(torch::torch_matmul(Jp_Sp_i, S_inv_r))
        prt2 <- Jp_rp_j$dot(inv_factor)
        prt3 <- -2.0 * inv_factor$dot(Jp_Sp_j_S_invr)
        prt4 <- -1.0 * S_inv_r$dot(torch::torch_matmul(Hp_Sp_ji, S_inv_r))
        prt5 <- 2.0 * S_inv_r$dot(torch::torch_matmul(term, S_inv_r))

        logDetTerm <- -1.0 * torch::torch_diag(S_inv_solve(term))$sum() +
                             torch::torch_diag(S_inv_solve(Hp_Sp_ji))$sum()

        Hij <- 0.5 * (2 * (prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)

        H_vals[j, i] <- Hij
        if (i != j) H_vals[i, j] <- Hij
      }
    }
    return(as.matrix(H_vals$cpu()))
  }
}

# --- Block-diagonal stacking helpers for multiple interpolants ---
# Each L_m(U,p,tt,sig) has shape K_m*D x mp1*D.  The block-diagonal assembly gives
# shape (K_total*D) x (M * mp1*D), so that L_block @ L_block^T is block-
# diagonal with blocks S_1, ..., S_M — no cross-covariance between interpolants.

# Assemble a constant block-diagonal tensor from a list of L0 matrices.
# @keywords internal
build_L0_block_joint <- function(L0_list, K_list, D, mp1, device) {
  M <- length(L0_list); K_total <- sum(K_list)
  out <- torch::torch_zeros(c(K_total * D, M * mp1 * D),
                            dtype = torch::torch_float64(), device = device)
  k_off <- 0L; m_off <- 0L
  for (i in seq_along(L0_list)) {
    k_i <- K_list[i] * D; m_i <- mp1 * D
    out[(k_off + 1):(k_off + k_i), (m_off + 1):(m_off + m_i)] <- L0_list[[i]]
    k_off <- k_off + k_i; m_off <- m_off + m_i
  }
  out
}

# Returns (U, p, tt, sig) -> block-diagonal (K_total*D) x (M*mp1*D) tensor.
# @keywords internal
build_L_block_joint <- function(L_fns, K_list, D, device) {
  M <- length(L_fns); K_total <- sum(K_list)
  dtype <- torch::torch_float64()
  function(U, p, tt, sig) {
    mp1   <- length(tt)
    mats  <- lapply(L_fns, function(lf) lf(U, p, tt, sig))
    out   <- torch::torch_zeros(c(K_total * D, M * mp1 * D),
                                dtype = dtype, device = device)
    k_off <- 0L; m_off <- 0L
    for (i in seq_along(mats)) {
      k_i <- K_list[i] * D; m_i <- mp1 * D
      out[(k_off + 1):(k_off + k_i), (m_off + 1):(m_off + m_i)] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}

# Returns (U, p, tt, sig) -> block-diagonal (K_total*D) x (M*mp1*D) x J tensor.
# @keywords internal
build_Jp_L_block_joint <- function(Jp_L_fns, K_list, D, J, device) {
  M <- length(Jp_L_fns); K_total <- sum(K_list)
  dtype <- torch::torch_float64()
  function(U, p, tt, sig) {
    mp1   <- length(tt)
    mats  <- lapply(Jp_L_fns, function(jf) jf(U, p, tt, sig))
    out   <- torch::torch_zeros(c(K_total * D, M * mp1 * D, J),
                                dtype = dtype, device = device)
    k_off <- 0L; m_off <- 0L
    for (i in seq_along(mats)) {
      k_i <- K_list[i] * D; m_i <- mp1 * D
      out[(k_off + 1):(k_off + k_i), (m_off + 1):(m_off + m_i), ] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}

# Returns (U, p, tt, sig) -> block-diagonal (K_total*D) x (M*mp1*D) x J x J tensor.
# @keywords internal
build_Hp_L_block_joint <- function(Hp_L_fns, K_list, D, J, device) {
  M <- length(Hp_L_fns); K_total <- sum(K_list)
  dtype <- torch::torch_float64()
  function(U, p, tt, sig) {
    mp1   <- length(tt)
    mats  <- lapply(Hp_L_fns, function(hf) hf(U, p, tt, sig))
    out   <- torch::torch_zeros(c(K_total * D, M * mp1 * D, J, J),
                                dtype = dtype, device = device)
    k_off <- 0L; m_off <- 0L
    for (i in seq_along(mats)) {
      k_i <- K_list[i] * D; m_i <- mp1 * D
      out[(k_off + 1):(k_off + k_i), (m_off + 1):(m_off + m_i), , ] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}

# Analytic gradient of RSS w.r.t. U: ∂RSS/∂U = 2*(∂r/∂U)ᵀ r
# Uses the same (∂r/∂U)ᵀ-times-vector pattern as J_u_wnll term1, with r in place of α.
# Returns function(U, p, tt, b) -> (mp1 x D) matrix.
build_J_u_rss_joint <- function(J_u_fn, g, K, J_num, V, Vp, device) {
  dtype <- torch::torch_float64()
  function(U, p, tt, b) {
    D   <- ncol(U)
    mp1 <- length(tt)

    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J_num)
    input <- rbind(p_mat, tU, ttt)

    J_u_val <- torch::torch_tensor(J_u_fn(input), dtype = dtype, device = device)
    J_u_val <- torch::torch_reshape(J_u_val, c(mp1, D, D))

    r    <- g(U, p, tt) - b
    r_2d <- r$reshape(c(K, D))

    # ∂r[k,d]/∂U[m',d'] = V[k,m']*J_u[m',d,d'] + Vp[k,m']*δ(d,d')
    term_F <- torch::torch_einsum("km,mde,kd->me", list(V, J_u_val, r_2d))
    term_b <- torch::torch_einsum("km,kd->md",     list(Vp, r_2d))

    2 * as.matrix((term_F + term_b)$cpu())
  }
}

# Analytic gradient of WNLL w.r.t. U: ∂WNLL/∂U = (∂r/∂U)ᵀ α - αᵀ(∂L/∂U)β + tr(Γᵀ ∂L/∂U)
# where α = S⁻¹r, β = Lᵀα, Γ = S⁻¹L.
# sig is now passed at call time so the gradient can be computed at any sig.
# Returns function(U, p, tt, b, sig) -> (mp1 x D) matrix.
build_J_u_wnll_joint <- function(J_u_fn, J_uu_fn, S, g, L, K, J_num, V, Vp, device) {
  dtype <- torch::torch_float64()
  function(U, p, tt, b, sig) {
    D   <- ncol(U)
    mp1 <- length(tt)

    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J_num)
    input <- rbind(p_mat, tU, ttt)

    J_u_val  <- torch::torch_tensor(J_u_fn(input),  dtype = dtype, device = device)
    J_u_val  <- torch::torch_reshape(J_u_val,  c(mp1, D, D))

    J_uu_val <- torch::torch_tensor(J_uu_fn(input), dtype = dtype, device = device)
    J_uu_val <- torch::torch_reshape(J_uu_val, c(mp1, D, D, D))

    Sp  <- S(U, p, tt, sig)
    r   <- g(U, p, tt) - b
    Lp  <- L(U, p, tt, sig)

    cholL_S  <- torch::linalg_cholesky(Sp)
    alpha    <- torch::torch_cholesky_solve(r$unsqueeze(-1L), cholL_S, upper = FALSE)$squeeze(-1L)
    beta_vec <- torch::torch_mv(Lp$t()$contiguous(), alpha)
    Gamma    <- torch::torch_cholesky_solve(Lp$contiguous(), cholL_S, upper = FALSE)

    alpha_2d <- alpha$reshape(c(K, D))
    beta_2d  <- beta_vec$reshape(c(mp1, D))
    Gamma_4d <- Gamma$reshape(c(K, D, mp1, D))

    # Term 1: (∂r/∂U)ᵀ α
    # ∂r[k,d]/∂U[m',d'] = V[k,m']*J_u[m',d,d'] + Vp[k,m']*δ(d,d')
    term1a <- torch::torch_einsum("km,mde,kd->me", list(V, J_u_val,  alpha_2d))
    term1b <- torch::torch_einsum("km,kd->md",     list(Vp, alpha_2d))
    term1  <- term1a + term1b

    # Term 2: -αᵀ (∂L/∂U) β
    # (∂L/∂U[m',d'])[k*D+a, m*D+b] = V[k,m']*sig[a]*J_uu[m',a,b,d']*δ(m=m')
    term2 <- -torch::torch_einsum("ka,km,a,mabd,mb->md",
                                  list(alpha_2d, V, sig, J_uu_val, beta_2d))

    # Term 3: tr(Γᵀ ∂L/∂U)
    # = Σ_{k,a,b} Γ[k,a,m',b]*V[k,m']*sig[a]*J_uu[m',a,b,d']
    term3 <- torch::torch_einsum("kamb,km,a,mabd->md",
                                 list(Gamma_4d, V, sig, J_uu_val))

    as.matrix((term1 + term2 + term3)$cpu())
  }
}

# Analytic gradient of WNLL w.r.t. sig: ∂WNLL/∂σ_d = ⟨∂L/∂σ_d, Δ⟩_F
# where Δ = (1-ε)(Γ W - α β_Wᵀ), α = S⁻¹r, β_W = W Lᵀα, Γ = S⁻¹L.
# ∂L/∂σ_d has two contributions (L1 Jacobian term and L0 diagonal term):
#   einsum("mab,km,kamb->a", J_u, V, Δ_4d)  +  einsum("km,kama->a", Vp, Δ_4d)
# Returns function(U, p, tt, b, sig) -> length-D numeric vector.
build_J_sig_wnll_joint <- function(J_u_fn, S, g, L, K, J_num, V, Vp, W = NULL, diag_reg = 1e-10, device) {
  WEIGHT <- 1.0 - diag_reg
  dtype  <- torch::torch_float64()
  function(U, p, tt, b, sig) {
    D   <- ncol(U)
    mp1 <- length(tt)

    tU    <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J_num)
    input <- rbind(p_mat, tU, ttt)

    J_u_val <- torch::torch_tensor(J_u_fn(input), dtype = dtype, device = device)
    J_u_val <- torch::torch_reshape(J_u_val, c(mp1, D, D))

    Sp <- S(U, p, tt, sig)
    r  <- g(U, p, tt) - b
    Lp <- L(U, p, tt, sig)

    cholL_S <- torch::linalg_cholesky(Sp)
    alpha   <- torch::torch_cholesky_solve(r$unsqueeze(-1L), cholL_S, upper = FALSE)$squeeze(-1L)
    Gamma   <- torch::torch_cholesky_solve(Lp$contiguous(), cholL_S, upper = FALSE)

    beta    <- torch::torch_mv(Lp$t()$contiguous(), alpha)
    beta_W  <- if (!is.null(W)) torch::torch_mv(W, beta) else beta
    Gamma_W <- if (!is.null(W)) torch::torch_mm(Gamma, W) else Gamma

    Delta    <- WEIGHT * (Gamma_W - torch::torch_outer(alpha, beta_W))
    Delta_4d <- Delta$reshape(c(K, D, mp1, D))

    # L1 contribution: ∂L1/∂σ_d picks the d-th row of J_u (a=d)
    term1 <- torch::torch_einsum("mab,km,kamb->a", list(J_u_val, V, Delta_4d))
    # L0 contribution: ∂L0/∂σ_d is diagonal in (a,b), picks a=b=d
    term2 <- torch::torch_einsum("km,kama->a",     list(Vp, Delta_4d))

    as.numeric((term1 + term2)$cpu())
  }
}

# Robust solver for applying the inverse of S(p) to a vector or matrix (torch)
# Sp: torch tensor (n x n), positive (semi-)definite
# Returns a closure that applies S⁻¹ to a 1D (n,) or 2D (n, m) torch tensor.
make_S_inv_solver <- function(Sp) {
  cholL <- tryCatch(torch::linalg_cholesky(Sp), error = function(e) NULL)

  if (!is.null(cholL)) {
    function(x) {
      squeezed <- (x$dim() == 1L)
      if (squeezed) x <- x$unsqueeze(-1L)
      z <- torch::torch_cholesky_solve(x$contiguous(), cholL, upper = FALSE)
      if (squeezed) z <- z$squeeze(-1L)
      z
    }
  } else {
    n      <- Sp$size(1)
    reg    <- 1e-12 * torch::torch_eye(n, dtype = Sp$dtype, device = Sp$device)
    Sp_reg <- (Sp + reg)$contiguous()
    function(x) {
      squeezed <- (x$dim() == 1L)
      if (squeezed) x <- x$unsqueeze(-1L)
      z <- torch::linalg_solve(Sp_reg, x$contiguous())
      if (squeezed) z <- z$squeeze(-1L)
      z
    }
  }
}
