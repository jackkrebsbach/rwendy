{
library(deSolve)
library(plotly)
library(trust)
library(symengine)
source("./R/symbolics.R")

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

u <- do.call(c, lapply(1:ncol(U), function(i) S(paste0("u", i))))
p <- do.call(c, lapply(1:length(p0), function(i) S(paste0("p", i))))
t <- S("t")

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

}



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









