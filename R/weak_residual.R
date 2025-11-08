
build_F <- function(U, tt, f, J) {
  mp1 <- nrow(U)
  D   <- ncol(U)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    Fp <- f(input)
    return(torch::torch_tensor(Fp))
  }
}

build_g <- function(V, F_) {
  function(p) {
    result <- torch::torch_matmul(V, F_(p))
    result <- torch::torch_transpose(result, dim0 = 1, dim1 = 2)
    torch::torch_reshape(result, c(-1))
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
    g_j <- as.vector(torch::as_array(torch::torch_matmul(V, F_e)))
    G[,j] <- g_j
  }
  return(torch::torch_tensor(G))
}

# ∇ₚ∇ₚL(p) Hessian of the Covariance factor where ∇ₚ∇ₚS(p) = ∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ + (∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ)ᵀ
build_Hp_L <-function(U, tt, J_upp, K, J, D, V, sig){
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    T_F <- torch::torch_reshape(torch::torch_tensor(J_upp(input)), c(mp1, D, D, J, J))
    H_ <- torch::torch_einsum("mabcd,km,a->kambcd", list(T_F, V, sig))$permute(c(2,1,4,3,5,6))
    H_ <- torch::torch_reshape(H_, c(K * D, mp1 * D, J, J))
    return(H_)
  }
}

build_Jp_L <-function(U, tt, J_up, K, J, D, V, sig){
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    H_F <- torch::torch_reshape(torch::torch_tensor(J_up(input)), c(mp1, D, D, J))
    J_ <- torch::torch_einsum("mabj,km,a->kambj", list(H_F, V, sig))$permute(c(2,1,4,3,5))
    J_ <- torch::torch_reshape(J_, c(K * D, mp1 * D, J))
    return(J_)
  }
}

build_L0 <- function(K, D, mp1, Vp, sig) {
  sig_diag <- torch::torch_diag(sig)
  L0_ <- torch::torch_einsum('km,ab->kabm', list(Vp, sig_diag))$permute(c(2,1,3,4))
  L0 <- torch::torch_reshape(L0_, c(K * D, mp1 * D))
  return(L0)
}

build_L <-function(U, tt, J_u, K, V, L0, sig, J){
  D <- ncol(U)
  mp1 <- length(tt)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_F <- torch::torch_reshape(torch::torch_tensor(J_u(input)), c(mp1, D, D))
    L1 <- torch::torch_einsum('mab,km,a->kamb', list(J_F, V, sig))$permute(c(2,1,4,3))
    L1 <- torch::torch_reshape(L1, c(K * D, mp1 * D))
    L <- L1 + L0
    return(L)
  }
}

build_S <- function(L, REG = 1e-10) {
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

build_Jp_r <- function(J_p, K, D, J, mp1, V, U, tt){
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))

    Jp_F <- torch::torch_reshape(torch::torch_tensor(J_p(input)), c(mp1, D, J))
    Jp_r <- torch::torch_einsum("km,mdj->kdj", list(V, Jp_F))

    Jp_r <- Jp_r$permute(c(2,1,3))
    Jp_r <- torch::torch_reshape(Jp_r, c(K * D, J))
  }
}

build_Hp_r <- function(H_p, K, D, J, mp1, V, U, tt){
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))

    Hp_F <- torch::torch_reshape(torch::torch_tensor(H_p(input)), c(mp1, D, J, J))
    Hp_r <- torch::torch_einsum("km,mdab->kdab", list(V, Hp_F))

    Hp_r <- Hp_r$permute(c(2,1,3,4))
    Hp_r <- torch::torch_reshape(Hp_r, c(K * D, J, J))

    return(Hp_r)
  }
}

build_wnll <- function(S, g, b, K, D){
  constant_term <- 0.5 * K * D * log(2 * pi)
  function(p){
    Sp <- S(p)
    r <- g(p) - b
    cholF <- torch::torch_cholesky(Sp)

    log_det <- 2 * torch::torch_sum(torch::torch_log(torch::torch_diag(cholF)))
    S_invr <- torch::torch_cholesky_solve(r$unsqueeze(2), cholF)$squeeze(2)
    mdist <- torch::torch_matmul(r$unsqueeze(1), S_invr)$squeeze()

    result <- 0.5 * (mdist + log_det) + constant_term
    return(result)
  }
}

build_J_wnll <- function(S, Jp_S, Jp_r, g, b, J){
  function(p){
    Sp <- S(p)
    J_Sp <- Jp_S(p)
    J_rp <- Jp_r(p)
    J_rp_t <- J_rp$t()
    r <- g(p) - b
    cholF <- torch::torch_cholesky(Sp)

    # Batched across all J parameters. prt0, prt1, logDetPart are vectors
    S_inv_rp <- torch::torch_cholesky_solve(r$unsqueeze(2), cholF)
    tmp <- torch::torch_einsum("ijk,j->ik", list(J_Sp, S_inv_rp$squeeze(2)))$unsqueeze(2)
    prt0 <- 2 * torch::torch_matmul(J_rp_t, S_inv_rp)$squeeze(2)
    prt1 <- -1.0 * torch::torch_einsum("ij,ijk -> k", list(S_inv_rp, tmp))
    fact <- torch::torch_cholesky_solve(J_Sp, cholF)
    diags <- torch::torch_diagonal(fact, dim1=1, dim2=2)
    logDetPart <- torch::torch_sum(diags, dim=2)

    gradient <- 0.5 * (prt0 + prt1 + logDetPart)
    return(gradient)
  }
}

build_H_wnll <- function(S, J_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J) {
  function(p) {
    r <- g(p) - b
    Jp_rp <- Jp_r(p)
    Hp_rp <- Hp_r(p)
    Lp <- L(p)
    Jp_Lp <- Jp_L(p)
    Hp_Lp <- Hp_L(p)
    Sp <- S(p)
    Jp_Sp <- J_S(p)

    S_inv_solve <- build_inv_solver(Sp)
    S_inv_rp <- S_inv_solve(r)

    H_wnn <- torch::torch_zeros(c(J, J), dtype = Sp$dtype, device = Sp$device)

    for (j in seq_len(J)) {
      Jp_rp_j <- Jp_rp[, j]
      Jp_Sp_j <- Jp_Sp[, , j]
      Jp_Lp_j <- Jp_Lp[, , j]
      shar_ <- S_inv_solve(Jp_Sp_j)

      for (i in j:J) {
        Jp_Sp_i <- Jp_Sp[, , i]
        Jp_Lp_i <- Jp_Lp[, , i]
        Jp_rp_i <- Jp_rp[, i]

        term <- torch::torch_mm(Jp_Sp_i, shar_)

        Hp_Lp_ji <- Hp_Lp[, , j, i]
        p1 <- torch::torch_mm(Jp_Lp_j, Jp_Lp_i$t())
        p2 <- torch::torch_mm(Hp_Lp_ji, Lp$t())
        Hp_Sp_ji <- p1 + p1$t() + p2 + p2$t()

        Hp_rp_ji <- Hp_rp[, j, i]

        prt0 <- torch::torch_dot(Hp_rp_ji, S_inv_rp)$item()
        prt1 <- -1.0 * torch::torch_dot(S_inv_solve(Jp_rp_j), torch::torch_mv(Jp_Sp_i, S_inv_rp))$item()

        inv_factor <- S_inv_solve(Jp_rp_i)
        prt2 <- torch::torch_dot(Jp_rp_j, inv_factor)$item()
        prt3 <- -2.0 * torch::torch_dot(inv_factor, torch::torch_mv(Jp_Sp_j, S_inv_rp))$item()
        prt4 <- -1.0 * torch::torch_dot(S_inv_rp, torch::torch_mv(Hp_Sp_ji, S_inv_rp))$item()
        prt5 <- 2.0 * torch::torch_dot(S_inv_rp, torch::torch_mv(term, S_inv_rp))$item()

        logDetTerm <- -1.0 * torch::torch_trace(S_inv_solve(term))$item() +
          torch::torch_trace(S_inv_solve(Hp_Sp_ji))$item()

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

build_inv_solver <- function(Sp) {
  cholF <- tryCatch(torch::torch_cholesky(Sp), error = function(e) NULL)
  if (!is.null(cholF)) {
    return(function(x) {
      is_vector <- length(x$shape) == 1
      if (is_vector) {
        result <- torch::torch_cholesky_solve(x$unsqueeze(2), cholF)$squeeze(2)
      } else {
        result <- torch::torch_cholesky_solve(x$unsqueeze(-1), cholF)$squeeze(-1)
      }
      return(result)
    })
  }
  tryCatch({
    test_solve <- torch_solve(torch_eye(nrow(Sp), device = Sp$device), Sp)
    return(function(x) {
      is_vector <- length(x$shape) == 1
      if (is_vector) {
        result <- torch::torch_solve(x$unsqueeze(2), Sp)[[1]]$squeeze(2)
      } else {
        result <- torch::torch_solve(x$unsqueeze(-1), Sp)[[1]]$squeeze(-1)
      }
      return(result)
    })
  }, error = function(e) {
    diag_reg <- 1e-12 * torch::torch_eye(Sp$shape[1], device = Sp$device)
    Sp_reg <- Sp + diag_reg
    return(function(x) {
      is_vector <- length(x$shape) == 1
      if (is_vector) {
        result <- torch::torch_solve(x$unsqueeze(2), Sp_reg)[[1]]$squeeze(2)
      } else {
        result <- torch::torch_solve(x$unsqueeze(-1), Sp_reg)[[1]]$squeeze(-1)
      }
      return(result)
    })
  })
}
