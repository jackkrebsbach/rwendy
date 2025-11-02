{
library(deSolve)
library(plotly)
library(trust)
library(symengine)
source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/noise.R")

test_params <- list(
  radius_params = 2^(0:3),
  radius_min_time = 0.1,
  radius_max_time = 5.0,
  k_max = 200,
  max_test_fun_condition_number = 1e4,
  min_test_fun_info_number = 0.95
)

f <- function(u, p, t) {
  du1 <- p[1] * (u[2] - u[1])
  du2 <- u[1] * (p[2] - u[3]) - u[2]
  du3 <- u[1] * u[2] - p[3] * u[3]
  c(du1, du2, du3)
}

noise_sd <- 0.05
p_star <- c(10.0, 28.0, 8.0 / 3.0)
p0 <- c(12.0, 21, 4.0)
u0 <- c(2, 1, 1)
npoints <- 256
t_span <- c(0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

test_fun_matrices <- build_full_test_function_matrices(U, tt, test_params, TRUE)

V <- test_fun_matrices$V
V_torch <- torch::torch_tensor(V)

Vp <- test_fun_matrices$V_prime
Vp_torch <- torch::torch_tensor(Vp)

u <- do.call(c, lapply(1:ncol(U), function(i) symengine::S(paste0("u", i))))
p <- do.call(c, lapply(1:length(p0), function(i) symengine::S(paste0("p", i))))
t <- symengine::S("t")

f_expr <- f(u, p, t)

J_u_sym <- compute_symbolic_jacobian(f_expr, u)
J_up_sym <- compute_symbolic_jacobian(J_u_sym, p)
J_p_sym <- compute_symbolic_jacobian(f_expr, p)
J_pp_sym <- compute_symbolic_jacobian(J_p_sym, p)
J_upp_sym <- compute_symbolic_jacobian(J_up_sym, p)

vars <- c(p, u ,t)

f_ <- build_fn(f_expr, vars)
J_u <- build_fn(J_u_sym, vars)
J_up <- build_fn(J_up_sym, vars)
J_p <- build_fn(J_p_sym, vars)
J_pp <- build_fn(J_pp_sym, vars)
J_upp <- build_fn(J_upp_sym, vars)

mp1 <- nrow(U)
D <- ncol(U)
J <- length(p_star)
K <- nrow(V)
p <- c(12.0, 21, 4.0)

sig <- estimate_std(U, k = 6)
sig_torch <- torch::torch_tensor(sig)

b <- -1 * as.vector(Vp %*% U)
b_torch <- torch::torch_tensor(b)
}

build_H_wnll <- function(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J) {
  function(p) {

    r <- g(p) - b

    Jp_rp <- Jp_r(p)
    Hp_rp <- Hp_r(p)

    Lp <- L(p)        # L(p)
    Jp_Lp <- Jp_L(p)  # ∇ₚL(p)
    Hp_Lp <- Hp_L(p)  # ∇ₚ∇ₚL(p)

    Sp <- S(p)
    Jp_Sp <- Jp_S(p)

    S_inv_solve <- function(x) {
      tryCatch({
        cholF <- chol(Sp)
        return(solve(cholF, solve(t(cholF), x)))
      }, error = function(e1) {
        tryCatch({
          qrF <- qr(Sp)
          return(solve(qrF, x))
        }, error = function(e2) {
          diag_reg <- 1e-12 * diag(nrow(Sp))
          return(solve(Sp + diag_reg, x))
        })
      })
    }

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

# Not working quite yet
build_H_wnll_torch <- function(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J) {
  function(p) {
    r <- g(p) - b

    Jp_rp <- Jp_r(p)
    Hp_rp <- Hp_r(p)

    Lp <- L(p)
    Jp_Lp <- Jp_L(p)
    Hp_Lp <- Hp_L(p)

    Sp <- S(p)
    Jp_Sp <- Jp_S(p)

    cholF <- torch::torch_cholesky(Sp)

    S_inv_rp <- torch::torch_cholesky_solve(r$unsqueeze(2), cholF)

    S_inv_Jp_Sp <- torch::torch_cholesky_solve(Jp_Sp, cholF)

    H_wnn <- torch::torch_zeros(c(J, J), dtype = Sp$dtype, device = Sp$device)

    for (j in seq_len(J)) {
      Jp_rp_j_vec <- Jp_rp[, j]
      Jp_Sp_j <- Jp_Sp[, , j]
      S_inv_Jp_Sp_j <- S_inv_Jp_Sp[, , j]
      Jp_Lp_j <- Jp_Lp[, , j]

      Jp_rp_j_col <- Jp_rp[, j]$reshape(c(-1, 1))
      S_inv_Jp_rp_j_vec <- torch::torch_cholesky_solve(Jp_rp_j_col, cholF)$squeeze(2)

      for (i in j:J) {
        S_inv_Jp_Sp_i <- S_inv_Jp_Sp[, , i]
        Jp_Sp_i <- Jp_Sp[, , i]
        Jp_Lp_i <- Jp_Lp[, , i]

        Jp_rp_i_vec <- Jp_rp[, i]

        Jp_rp_i_col <- Jp_rp[, i]$reshape(c(-1, 1))
        S_inv_Jp_rp_i_vec <- torch::torch_cholesky_solve(Jp_rp_i_col, cholF)$squeeze(2)

        Hp_rp_ji <- Hp_rp[, j, i]

        term <- torch::torch_matmul(Jp_Sp_i, S_inv_Jp_Sp_j)

        Hp_Lp_ji <- Hp_Lp[, , j, i]
        p1 <- torch::torch_matmul(Jp_Lp_j, Jp_Lp_i$t())
        p2 <- torch::torch_matmul(Hp_Lp_ji, Lp$t())
        Hp_Sp_ji <- p1 + p1$t() + p2 + p2$t()

        S_inv_Hp_Sp_ji <- torch::torch_cholesky_solve(Hp_Sp_ji, cholF)
        S_inv_term <- torch::torch_cholesky_solve(term, cholF)

        prt0 <- torch::torch_dot(Hp_rp_ji, S_inv_rp$squeeze(2))

        prt1 <- -1.0 * torch::torch_dot(
          S_inv_Jp_rp_j_vec,
          torch::torch_matmul(Jp_Sp_i, S_inv_rp)$squeeze(2)
        )

        prt2 <- torch::torch_dot(Jp_rp_j_vec, S_inv_Jp_rp_i_vec)

        prt3 <- -2.0 * torch::torch_dot(
          S_inv_Jp_rp_i_vec,
          torch::torch_matmul(Jp_Sp_j, S_inv_rp)$squeeze(2)
        )

        prt4 <- -1.0 * torch::torch_dot(
          S_inv_rp$squeeze(2),
          torch::torch_matmul(Hp_Sp_ji, S_inv_rp)$squeeze(2)
        )

        prt5 <- 2.0 * torch::torch_dot(
          S_inv_rp$squeeze(2),
          torch::torch_matmul(term, S_inv_rp)$squeeze(2)
        )

        logDetTerm <- -1.0 * torch::torch_trace(S_inv_term) + torch::torch_trace(S_inv_Hp_Sp_ji)

        Hij <- 0.5 * (2.0 * (prt0 + prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)

        H_wnn[j, i] <- Hij

        if (i != j) {
          H_wnn[i, j] <- Hij
        }
      }
    }

    return(H_wnn)
  }
}



H_wnll <- build_H_wnll(S, J_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J)
H_wnll_torch <- build_H_wnll_torch(S_torch, J_S_torch, L_torch, Jp_L_torch, Hp_L_torch, Jp_r_torch, Hp_r_torch, g_torch, b_torch, J)

H_wnllp <- H_wnll(p)
H_wnll_torchp <- H_wnll_torch(p)

all.equal(H_wnllp, as.array(H_wnll_torchp))


build_J_wnll <- function(S, Jp_S, Jp_r, g, b, J){
  function(p){
    Sp <- S(p)
    J_Sp <- Jp_S(p)
    J_rp <- Jp_r(p)
    r <- g(p) - b

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

build_J_wnll_torch <- function(S, Jp_S, Jp_r, g, b, J){
  function(p){
    Sp <- S(p)
    J_Sp <- Jp_S(p)
    J_rp <- Jp_r(p)
    J_rp_t <- J_rp$t()
    r <- g(p) - b
    cholF <- torch::torch_cholesky(Sp)

    # Batched across all parameters J
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

J_wnll <- build_J_wnll(S, J_S, Jp_r, g, b, J)
J_wnll_torch <- build_J_wnll_torch(S_torch, J_S_torch, Jp_r_torch, g_torch, b_torch, J)

p <- p_star
J_wnllp <- J_wnll(p)
J_wnll_torchp <- J_wnll_torch(p)

all.equal(J_wnllp, as.array(J_wnll_torchp))

build_wnll <- function(S, g, b, K, D){
  constant_term <- 0.5 * K * D * log(2 * pi)
  function(p){
    Sp <- S(p)
    r <- g(p) - b

    cholF <- chol(Sp)
    log_det <- 2 * sum(log(diag(cholF)))

    S_invr <- solve(cholF, solve(t(cholF), r))

    mdist <- (r %*% S_invr)[1,1]
    return(0.5 * (mdist + log_det) + constant_term)
  }
}

build_wnll_torch <- function(S, g, b, K, D){
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

wnll <- build_wnll(S, g, b, K, D)
wnll_torch <- build_wnll_torch(S_torch, g_torch, b_torch, K, D)

wnllp <- wnll(p)
wnll_torchp <- wnll_torch(p)

all.equal(wnllp, as.array(wnll_torchp))


build_J_S <- function(L, Jp_L, J, K, D){
  function(p){
    Lp <- L(p)
    Jp_Lp <- Jp_L(p)

    Jp_S <- array(0, dim = c(K*D, K*D, J))
    for(j in 1:J){
      Jp_L_j <-  Jp_Lp[,,j]
      prt <- Jp_L_j %*% t(Lp)
      Jp_S[,,j] <- prt + t(prt)
    }
    return(Jp_S)
  }
}

build_J_S_torch <- function(L, Jp_L, J, K, D){
  function(p){
    Lp <- L(p)
    Jp_Lp <- Jp_L(p)
    Lp_t <- Lp$t()
    prt <- torch::torch_einsum("ijk,jl->ilk", list(Jp_Lp, Lp_t))
    Jp_S <- prt + prt$transpose(1, 2)
    return(Jp_S)
  }
}

J_S <- build_J_S(L, Jp_L, J, K, D)
J_S_torch <- build_J_S_torch(L_torch, Jp_L_torch, J, K, D)

p <- c(12.0, 21, 4.0)
J_Sp <- J_S(p)
J_Sp_torch <- J_S_torch(p)

all.equal(J_Sp, as.array(J_Sp_torch))

build_S <- function(L, REG = 1e-10) {
  function(p) {
    Lp <- L(p)
    S_ <- tcrossprod(Lp)
    WEIGHT <- 1.0 - REG
    eye <- REG * diag(nrow(S_))
    S <- WEIGHT * S_ + eye
    return(S)
  }
}

build_S_torch <- function(L, REG = 1e-10) {
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

S <- build_S(L)
S_torch <- build_S_torch(L_torch)

Sp <- S(p)
Sp_torch <- S_torch(p)

all.equal(Sp, as.array(Sp_torch))

# ∇ₚ∇ₚL(p) Hessian of the Covariance factor where ∇ₚ∇ₚS(p) = ∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ + (∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ)ᵀ
build_Hp_L <-function(U, tt, J_upp, K, J, D, V, sig){
  D <- ncol(U)
  mp1 <- length(tt)

  function(p){
    T_F <- array(0, dim = c(mp1, D, D, J, J))

    for (i in 1:mp1) {
      t <- tt[i]
      u <- U[i, ]
      val <- J_upp(c(p, u, t))
      T_F[i, , , ,] <- val
    }

    H_ <- array(0, dim = c(K, D, mp1, D, J, J))

    for (k in 1:K) {
      for (d1 in 1:D) {
        for (m in 1:mp1) {
          for (d2 in 1:D) {
            for (j1 in 1:J) {
              for (j2 in 1:J) {
                H_[k, d1, m, d2, j1, j2] <- T_F[m, d1, d2, j1, j2] * V[k, m] * sig[d1]

              }
            }
          }
        }
      }
    }

    Hp_L <- array(H_, dim = c(K * D, mp1 * D, J, J))
    return(Hp_L)
  }
}

# ∇ₚ∇ₚL(p) Hessian of the Covariance factor where ∇ₚ∇ₚS(p) = ∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ + (∇ₚ∇ₚLLᵀ + ∇ₚL∇ₚLᵀ)ᵀ
build_Hp_L_torch <-function(U, tt, J_upp, K, J, D, V, sig){
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

p <- c(12.0, 21, 4.0)

Hp_L <- build_Hp_L(U, tt, J_upp, K, J, D, V, sig)
Hp_L_torch <- build_Hp_L_torch(U, tt, J_upp, K, J, D, V_torch, sig_torch)

Hp_Lp <- Hp_L(p)
Hp_L_torchp <- Hp_L_torch(p)

all.equal(Hp_Lp, as.array(Hp_L_torchp))


build_Jp_L <-function(U, tt, J_up, K, J, D, V, sig){
  mp1 <- length(tt)
  function(p){
    H_F <- array(0, dim = c(mp1, D, D, J))

    for (i in 1:mp1) {
      t <- tt[i]
      u <- U[i, ]
      val <- J_up(c(p, u, t))
      H_F[i, , ,] <- val
    }

    J_ <- array(0, dim = c(K, D, mp1, D, J))

    for (k in 1:K) {
      for (d1 in 1:D) {
        for (m in 1:mp1) {
          for (d2 in 1:D) {
            for (j in 1:J) {
              J_[k, d1, m, d2, j] <- H_F[m, d1, d2, j] * V[k, m] * sig[d1]
            }
          }
        }
      }
    }
    Jp_L <- array(J_, dim = c(K * D, mp1 * D, J))
    return(Jp_L)
  }
}

build_Jp_L_torch <-function(U, tt, J_up, K, J, D, V, sig){
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

Jp_L_torch <- build_Jp_L_torch(U, tt, J_up, K, J, D, V_torch, sig_torch)
Jp_L <- build_Jp_L(U, tt, J_up, K, J, D, V, sig)

p <- c(12.0, 21, 4.0)
Jp_L_torchp <- Jp_L_torch(p)
Jp_Lp <- Jp_L(p)

all.equal(Jp_Lp, as.array(Jp_L_torchp))


build_L0 <- function(K, D, mp1, Vp, sig) {
  L0_ <- array(0, dim = c(K, D, D, mp1))

  for (k in 1:K) {
    for (d in 1:D) {
      for (m in 1:mp1) {
        L0_[k, d, d, m] <- Vp[k,m] * sig[d]
      }
    }
  }
  L0_ <- aperm(L0_, c(1,2,4,3))
  L0 <- array(L0_, dim = c(K * D, mp1 * D))
  return(L0)
}

build_L0_torch <- function(K, D, mp1, Vp, sig) {
  sig_diag <- torch::torch_diag(sig)
  L0_ <- torch::torch_einsum('km,ab->kabm', list(Vp, sig_diag))$permute(c(2,1,3,4))
  L0 <- torch::torch_reshape(L0_, c(K * D, mp1 * D))
  return(L0)
}

L0 <- build_L0(K, D, mp1, Vp, sig)
L0_torch <- build_L0_torch(K, D, mp1 ,Vp_torch, sig_torch)

all.equal(L0, as.array(L0_torch))

build_L <-function(U, tt, J_u, K, V, Vp, sig){
  D <- ncol(U)
  mp1 <- length(tt)
  L0 <- build_L0(K, D, mp1, Vp, sig)

  function(p){
    J_F <- array(0, dim = c(mp1, D, D))
    for (i in 1:mp1) {
      u <- U[i, ]
      t <- tt[i]
      val <- J_u(c(p, u, t))
      J_F[i, ,] <- val
    }

    L1 <- array(0, dim = c(K, D, mp1, D))
    for (k in 1:K) {
      for (d1 in 1:D) {
        for (m in 1:mp1) {
          for (d2 in 1:D) {
            L1[k, d1, m, d2] <- V[k, m] * J_F[m, d1, d2] * sig[d1]
          }
        }
      }
    }

    L1 <- array(L1, dim = c(K * D, mp1 * D))
    L <- L1 + L0
    return(L)
  }
}

build_L_torch <-function(U, tt, J_u, K, V, L0, sig){
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

L <- build_L(U, tt, J_u, K, V, Vp, sig)
L_torch <- build_L_torch(U, tt, J_u, K, V_torch, L0_torch, sig_torch)
Lp <- L(p)
L_torchp <- L_torch(p)
all.equal(Lp, as.array(L_torchp))


build_Hp_r <- function(H_p, K, D, J, mp1, V, U, tt){
  function(p){
    Hp_F <- array(0, dim = c(mp1, D, J, J))

    for (i in 1:mp1) {
      t <- tt[i]
      u <- U[i, ]
      input_vec <- c(p, u, t)
      val <- H_p(input_vec)
      Hp_F[i, , ,] <- val
    }

    Hp_r <- array(0, dim = c(K, mp1, D, J, J))
    for (k in 1:K) {
      for (m in 1:mp1) {
        for (d1 in 1:D) {
          for (j1 in 1:J) {
            for (j2 in 1:J) {
              Hp_r[k, m, d1, j1, j2] <- V[k, m] * Hp_F[m, d1, j1, j2]
            }
          }
        }
      }
    }
    Hp_r <- apply(Hp_r, c(1, 3, 4, 5), sum)
    Hp_r <- array(Hp_r, dim = c(K * D, J, J))

    return(Hp_r)
  }
}

build_Hp_r_torch <- function(H_p, K, D, J, mp1, V, U, tt){
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))

    Hp_F <- torch::torch_reshape(torch::torch_tensor(H_p(input)), c(mp1, D, J, J))
    Hp_r <- torch::torch_einsum("km,mdab->kdab", list(V, Hp_F))

    Hp_r <- Hp_r$permute(c(3,2,1,4))
    Hp_r <- torch::torch_reshape(Hp_r, c(K * D, J, J))

    return(Hp_r)
  }
}

Hp_r <- build_Hp_r(J_pp, K, D, J, mp1, V, U, tt)
Hp_r_torch <- build_Hp_r_torch(J_pp, K, D, J, mp1, V_torch, U, tt)

Hp_rp <- Hp_r(p)
Hp_r_torchp <- Hp_r_torch(p)

all.equal(Hp_rp, as.array(Hp_r_torchp))


build_Jp_r <- function(J_p, K, D, J, mp1, V, U, tt){
  function(p){
    Jp_F <- array(0, dim = c(mp1, D, J))

    for (i in 1:mp1) {
      t <- tt[i]
      u <- U[i, ]
      val <- J_p(c(p, u, t))
      Jp_F[i, ,] <- val
    }

    Jp_r <- array(0, dim = c(K, mp1, D, J))
    for (k in 1:K) {
      for (m in 1:mp1) {
        for (d1 in 1:D) {
          for (j1 in 1:J) {
            Jp_r[k, m, d1, j1] <- V[k, m] * Jp_F[m, d1, j1]
          }
        }
      }
    }
    Jp_r <- apply(Jp_r, c(1, 3, 4), sum)
    Jp_r <- array(Jp_r, dim = c(K * D, J))
    return(Jp_r)
  }
}

build_Jp_r_torch <- function(J_p, K, D, J, mp1, V, U, tt){
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))

    Jp_F <- torch::torch_reshape(torch::torch_tensor(J_p(input)), c(mp1, D, D))
    Jp_r <- torch::torch_einsum("km,mdj->kdj", list(V, Jp_F))

    Jp_r <- Jp_r$permute(c(2,1,3))
    Jp_r <- torch::torch_reshape(Jp_r, c(K * D, J))
  }
}

Jp_r <- build_Jp_r(J_p, K, D, J, mp1, V, U, tt)
Jp_r_torch <- build_Jp_r_torch(J_p, K, D, J, mp1, V_torch, U, tt)

Jp_rp <- Jp_r(p)
Jp_r_torchp <- Jp_r_torch(p)

all.equal(Jp_rp, as.array(Jp_r_torchp))

build_F <- function(U, tt, f, vars) {
  mp1 <- nrow(U)
  D   <- ncol(U)
  function(p) {
    out <- matrix(0, nrow = mp1, ncol = D)
    for (i in 1:mp1) {
      u <- U[i, ]
      t <- tt[i]
      input_vec <- c(p, u, t)
      udot <- f(input_vec)
      out[i,] <- udot
    }
    return(out)
  }
}

build_F_torch <- function(U, tt, f, vars, J) {
  mp1 <- nrow(U)
  D   <- ncol(U)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    Fp <- f(input)
    return(torch::torch_tensor(Fp))
  }
}

F_ <- build_F(U, tt, f_, vars)
F_torch <- build_F_torch(U, tt, f_, vars, J)

p <- c(12.0, 21, 4.0)
Fp <- F_(p)
F_torchp <- as.array(F_torch(p))

all.equal(Fp, F_torchp)

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
    g_j <- as.vector(V %*% F_e)
    G[,j] <- g_j
  }
  return(G)
}

build_G_matrix_torch <- function(V, U, tt, F_, J){
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

G_torch <- build_G_matrix_torch(V_torch, U, tt, F_torch, J)
G <- build_G_matrix(V, U, tt, F_, J)

all.equal(Fp, F_torchp)


build_g <- function(V, F_) {
  function(p) {
    as.vector(V %*% F_(p))
  }
}

build_g_torch <- function(V, F_) {
  function(p) {
    result <- torch::torch_matmul(V, F_(p))
    result <- torch::torch_transpose(result, dim0 = 1, dim1 = 2)
    torch::torch_reshape(result, c(-1))
  }
}

g <- build_g(V, F_)
g_torch <- build_g_torch(V_torch, F_torch)

p <- c(12.0, 21, 4.0)
gp <- g(p)
g_torchp <- as.array(g_torch(p))

all.equal(gp, g_torchp)









