# RHS of ODE system at different data points
build_F <- function(U, tt, f) {
  n_rows <- nrow(U)
  n_f    <- ncol(U)
  function(p) {
    out <- matrix(0, nrow = n_rows, ncol = n_f)
    for (i in 1:n_rows) {
      u <- U[i, ]
      t <- tt[i]
      input_vec <- c(p, u, t)
      udot <- f(input_vec)
      out[i,] <- udot
    }
    return(out)
  }
}
#  <V,F(p)> = VGp when f is Linear in Parameters
build_G_matrix <- function(V, U, tt, f, F_, J){
  K <- nrow(V)
  mp1 <- nrow(U)
  D <- ncol(U)
  G <- matrix(0, nrow = K*D, ncol = J)

  for(j in seq(1,J)){
    e_j <- rep(0, J)
    e_j[j] <- 1
    F_j <- F_(e_j)
    g_j <- as.vector(V %*% F_j)
    G[,j] <- g_j
  }
  return(G)
}

#  <V,F(p)> =  f(p) is nonlinear in parameters
build_g <- function(V, F_) {
  function(p) {
       as.vector(V %*% F_(p))
  }
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
            L1[k, d1, m, d2] <- J_F[m, d1, d2] * V[k, m] * sig[d1]
          }
        }
      }
    }

    dim(L1) <- c(K * D, mp1 * D)
    L <- L1 + L0
    matrix(L, nrow = K * D, ncol = D * mp1)

  }
}

build_S <- function(L, REG = 1e-9) {
  function(p) {
    Lp <- L(p)
    S_ <- tcrossprod(Lp)
    WEIGHT <- 1.0 - REG
    eye <- REG * diag(nrow(S_))
    S <- WEIGHT * S_ + eye
    return(S)
  }
}


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

    dim(J_) <- c(K * D, mp1 * D, J)
    return(J_)
  }
}

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

    dim(H_) <- c(K * D, mp1 * D, J, J)
    return(H_)
  }
}


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
           Jp_r[k, m, d1, j1] <- Jp_F[m, d1, j1] * V[k, m]
          }
        }
      }
    }
    Jp_r <- apply(Jp_r, c(1, 3, 4), sum)
    dim(Jp_r) <- c(K * D, J)
    return(Jp_r)
  }
}

# ∇ᵤr(p) ∈ ℝ^(K*D x J x J)
build_Ju_r <- function(J_u, K, D, J, mp1, V, U, tt){
  function(p){
    Ju_F <- array(0, dim = c(mp1, D, D))

    for (i in 1:mp1) {
      t <- tt[i]
      u <- U[i, ]
      val <- J_u(c(p, u, t))
      Ju_F[i, ,] <- val
    }

    Ju_r <- array(0, dim = c(K, mp1, D, D))
    for (k in 1:K) {
      for (m in 1:mp1) {
        for (d1 in 1:D) {
          for (d2 in 1:D) {
            Ju_r[k, m, d1, d2] <- Ju_F[m, d1, d2] * V[k,m]
          }
        }
      }
    }

    dim(Ju_r) <- c(K * D, mp1 * D)
    return(Ju_r)
  }
}

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
              Hp_r[k, m, d1, j1, j2] <- Hp_F[m, d1, j1, j2] * V[k, m]
          }
         }
        }
      }
    }
    Hp_r <- apply(Hp_r, c(1, 3, 4, 5), sum)
    dim(Hp_r) <- c(K * D, J, J)
    return(Hp_r)
  }
}

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

build_J_wnll <- function(S, Jp_S, Jp_r, g, b, J){
  function(p){
    Sp <- S(p)
    J_Sp <- Jp_S(p)
    J_rp <- Jp_r(p)
    r <- g(p) - b

    cholF <- chol(Sp)
    S_inv_r <- solve(cholF, solve(t(cholF), r))

    gradient <- numeric(J)

    for(j in seq(J)){
      J_S_j <- J_Sp[, , j]
      J_r_j <- J_rp[, j]

      tmp <- J_S_j %*% S_inv_r

      prt0 <- 2.0 * (J_r_j %*% S_inv_r)[1,1]
      prt1 <- -1.0 * (S_inv_r %*% tmp)[1,1]

      fact <- solve(cholF, solve(t(cholF), J_S_j))
      logDetPart <- sum(diag(fact))

      gradient[j] <- 0.5 * (prt0 + prt1 + logDetPart)
    }

    return(gradient)
  }
}

build_H_wnll <- function(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Ju_r, Hp_r, g, b, J) {
  function(p) {
    Sp <- S(p)
    Jp_Sp <- Jp_S(p)
    Jp_rp <- Jp_r(p)
    Ju_rp <- Ju_r(p)  # ∇ᵤr(p) ∈ ℝ^(K*D x D*mp1)
    Hp_rp <- Hp_r(p)

    Lp <- L(p)        # L(p)
    Jp_Lp <- Jp_L(p)  # ∇ₚL(p)
    Hp_Lp <- Hp_L(p)  # ∇ₚ∇ₚL(p)

    r <- g(p) - b

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
