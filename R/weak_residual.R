
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
    Sp <- S(p)
    r <- g(p) - b
    cholF <- torch::torch_cholesky(Sp)

    #log_det <- 2 * torch::torch_sum(torch::torch_log(torch::torch_diag(cholF)))
    svd_result <- torch::torch_svd(Sp, some = TRUE, compute_uv = FALSE)
    singular_values <- svd_result[[2]]
    log_det <- torch::torch_sum(torch::torch_log(singular_values))

    S_invr <- torch::torch_cholesky_solve(r$unsqueeze(2), cholF)$squeeze(2)
    mdist <- torch::torch_matmul(r$unsqueeze(1), S_invr)$squeeze()

    result <- 0.5 * (mdist + log_det) + constant_term
    return(as.numeric(result))
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
    return(as.numeric(gradient))
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

    return(as.matrix(H_wnn))
  }
}

build_inv_solver <- function(Sp, max_jitter = 1e-3, verbose = FALSE) {
  device <- Sp$device
  n <- Sp$shape[1]

  jitter_values <- c(0, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, max_jitter)

  for (jitter in jitter_values) {
    cholF <- tryCatch({
      Sp_work <- if (jitter > 0) {
        Sp + torch::torch_eye(n, device = device) * jitter
      } else {
        Sp
      }
      torch::torch_cholesky(Sp_work)
    }, error = function(e) NULL)

    if (!is.null(cholF)) {
      if (verbose && jitter > 0) {
        message(sprintf("Cholesky succeeded with jitter = %.2e", jitter))
      }

      jitter_used <- jitter

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
  }

  if (verbose) {
    message("Cholesky failed, falling back to LU decomposition")
  }

  for (jitter in c(0, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3)) {
    tryCatch({
      Sp_work <- if (jitter > 0) {
        Sp + torch::torch_eye(n, device = device) * jitter
      } else {
        Sp
      }

      test_vec <- torch::torch_randn(n, device = device)
      test_result <- torch::torch_linalg_solve(Sp_work, test_vec)

      residual <- torch::torch_norm(
        torch::torch_matmul(Sp_work, test_result) - test_vec
      )$item()

      if (residual < 1e-3) {  # Solution is acceptable
        if (verbose && jitter > 0) {
          message(sprintf("LU solve succeeded with jitter = %.2e", jitter))
        }

        Sp_stable <- Sp_work

        return(function(x) {
          torch::torch_linalg_solve(Sp_stable, x)
        })
      }
    }, error = function(e) NULL)
  }

  if (verbose) {
    warning("All standard methods failed, using pseudoinverse (may be inaccurate)")
  }

  lambda <- torch::torch_trace(Sp)$item() / n * 1e-3
  Sp_reg <- Sp + torch::torch_eye(n, device = device) * lambda

  tryCatch({
    Sp_inv <- torch::torch_linalg_pinv(Sp_reg, rtol = 1e-6)

    return(function(x) {
      torch::torch_matmul(Sp_inv, x)
    })
  }, error = function(e) {
    stop("Matrix is too ill-conditioned to invert reliably. Consider reformulating the problem.")
  })
}
