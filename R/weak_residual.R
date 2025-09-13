
build_F <- function(U, tt, f) {
  n_rows <- nrow(U)
  n_f    <- ncol(U)

  function(p) {
    out <- matrix(0, nrow = n_rows, ncol = n_f)

    for (i in 1:n_rows) {
      u_row <- U[i, ]
      t_val <- tt[i]

      input_vec <- c(u_row, p, t_val)
       udot <- f(input_vec)
       out[i,] <- udot
    }

    return(out)
  }
}

build_g <- function(V, F_) {
  function(p) {
       as.vector(V %*% F_(p))
  }
}


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

build_S <- function(L){
  function(p){
    Lp <- L(p)
    S <- tcrossprod(Lp)
  }
}



build_Jp_L <-function(U, tt, J_up, K, J, D, V, sig){
  D <- ncol(U)
  mp1 <- length(tt)
  function(p){
    H_F <- array(0, dim = c(mp1, D, D, J))

    for (i in 1:mp1) {
      t <- tt[i]
      u_row <- U[i, ]
      input_vec <- c(p, u_row, t)
      val <- as.array(J_up(input_vec))
      dim(val) <- c(D, D, J)
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
      u_row <- U[i, ]
      input_vec <- c(p, u_row, t)
      val <- as.array(J_upp(input_vec))
      dim(val) <- c(D, D, J, J)
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


build_Jp_r <- function(J_p, K, D, J, mp1, V, sig){
  function(p){
    Jp_F <- array(0, dim = c(mp1, D, J))

    for (i in 1:mp1) {
      t <- tt[i]
      u_row <- U[i, ]
      input_vec <- c(p, u_row, t)
      val <- as.array(J_p(input_vec))
      dim(val) <- c(D, J)
      Jp_F[i, ,] <- val
    }

    Jp_r <- array(0, dim = c(K, D, J))
    for (k in 1:K) {
      for (d1 in 1:D) {
        for (j1 in 1:J) {
          Jp_r[k, d1, j1] <- sum(Jp_F[, d1, j1] * V[k, ] * sig[d1])
        }
      }
    }

    dim(Jp_r) <- c(K * D, J)
    return(Jp_r)
  }
}


build_wnll <- function(S, g, b){
  function(p){
    Sp <- S(p)
    r <- g(p) - b
    cholF <- chol(Sp)
    log_det <- 2 * sum(log(diag(cholF)))
    S_invr <- solve(cholF, solve(t(cholF), r))
    mdist <- r %*% S_invr
    return(0.5 * (mdist + log_det))
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

      prt0 <- 2.0 * sum(J_r_j * S_inv_r)
      prt1 <- -1.0 * sum(S_inv_r * tmp)

      fact <- solve(cholF, solve(t(cholF), J_S_j))
      logDetPart <- sum(diag(fact))

      gradient[j] <- prt0 + prt1 + logDetPart
    }

    return(gradient)
  }
}

build_H_wnll <- function(S, Jp_S, Jp_r, g, b, J){
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

      prt0 <- 2.0 * sum(J_r_j * S_inv_r)
      prt1 <- -1.0 * sum(S_inv_r * tmp)

      fact <- solve(cholF, solve(t(cholF), J_S_j))
      logDetPart <- sum(diag(fact))

      gradient[j] <- prt0 + prt1 + logDetPart
    }

    return(gradient)
  }
}
