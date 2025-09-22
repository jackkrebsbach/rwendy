source("./R/noise.R")
source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/weak_residual.R")
source("./R/wendy.R")

set.seed(123)


build_L0 <- function(K, D, mp1, Vp, sig) {
  L0_ <- array(0, dim = c(K, D, D, mp1))

  for (d in 1:D) {
    L0_[, d, d, ] <- Vp * sig[d]
  }

  L0 <- aperm(L0_, c(1, 2, 4, 3))
  dim(L0) <- c(K * D, mp1 * D)
  return(L0)
}


build_L <-function(U, tt, J_u, K, V, Vp, sig){
  D <- ncol(U)
  mp1 <- length(tt)
  L0 <- build_L0(K, D, mp1, Vp, sig)

  function(p){
    J_F <- array(0, dim = c(mp1, D, D))
    for (i in 1:mp1) {
      t <- tt[i]
      u_row <- U[i, ]
      input_vec <- c(p, u_row, t)
      J_F[i, ,] <- as.array(J_u(input_vec))
    }

    L1_ <- array(0, dim = c(K, D, mp1, D))
    for (k in 1:K) {
      for (d1 in 1:D) {
        for (m in 1:mp1) {
          for (d2 in 1:D) {
            L1_[k, d1, m, d2] <- J_F[m, d1, d2] * V[k, m] * sig[d1]
          }
        }
      }
    }

    L1 <- aperm(L1_, c(1, 2, 3, 4))
    dim(L1) <- c(K * D, mp1 * D)
    L <- L1 + L0
    matrix(L, nrow = K*D, ncol = D*mp1)

  }
}

D   <- 3
mp1 <- 4
K   <- 2

U   <- matrix(rnorm(mp1 * D), nrow = mp1, ncol = D)
tt  <- seq_len(mp1)
V   <- matrix(runif(K * mp1), nrow = K, ncol = mp1)
Vp  <- matrix(runif(K * mp1), nrow = K, ncol = mp1)
sig <- runif(D)

J_u <- function(x) {
  matrix(runif(D * D), nrow = D, ncol = D)
}


f1 <- build_L(U, tt, J_u, K, V, Vp, sig)
f2 <- build_L_2(U, tt, J_u, K, V, Vp, sig)

p <- runif(2)

L1 <- f1(p)
L2 <- f2(p)

all.equal(L1, L2)
max(abs(L1 - L2))
