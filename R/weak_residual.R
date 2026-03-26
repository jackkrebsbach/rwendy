# F(p) in -𝚽'U = 𝚽F(p,U,t)
build_F <- function(U, tt, f_, J) {
  mp1 <- nrow(U)
  D   <- ncol(U)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    f_(input)  # returns mp1 x D matrix
  }
}

# g(p) = 𝚽F(p,u,t)
build_g <- function(V, F_) {
  function(p) {
    as.vector(V %*% F_(p))
  }
}

# G matrix in the linear system Gp = b - g0 : G ∈ ℝ^{K * D x J}
build_G_matrix <- function(V, U, tt, F_, J) {
  K   <- nrow(V)
  mp1 <- nrow(U)
  D   <- ncol(U)
  G   <- matrix(0, nrow = K * D, ncol = J)
  g0  <- F_(rep(0, J))
  for (j in seq_len(J)) {
    e_j    <- rep(0, J); e_j[j] <- 1
    F_e    <- F_(e_j) - g0
    G[, j] <- as.vector(V %*% F_e)
  }
  G
}

# g(p) = Gp in system Gp = b - g0
build_g_linear <- function(G) {
  function(p) as.vector(G %*% p)
}

# ∇ₚr(p) ∈ ℝ^(K*D × J) Jacobian of the weak residual
# J_p returns (mp1 x D*J) with column c = ∂f_{c/J+1}/∂p_{c%J+1}
# Reshape as [m, j, d]: Jp_F[m,j,d] = ∂f_d/∂p_j (p-index first, column-major convention)
build_Jp_r <- function(J_p, K, D, J, mp1, V, U, tt) {
  function(p) {
    p_mat  <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input  <- rbind(p_mat, t(U), t(tt))
    Jp_F   <- array(J_p(input), c(mp1, J, D))
    result <- array(V %*% matrix(Jp_F, mp1, J * D), c(K, J, D))
    matrix(aperm(result, c(1L, 3L, 2L)), K * D, J)
  }
}

# ∇ₚr(p) ∈ ℝ^(K*D × J) Jacobian of the weak residual: r(p) = Gp + go - b → ∇ₚr(p) = G
build_Jp_r_linear <- function(G) {
  function(p) G
}

# ∇ₚ∇ₚr(p) ∈ ℝ^(K*D × J × J) Hessian of the weak residual
# J_pp returns (mp1 x D*J*J) with column c = ∂²f_{c/J/J+1}/(∂p_{(c/J)%J+1} ∂p_{c%J+1})
# Reshape as [m, j2, j1, d]: Hp_F[m,j2,j1,d] = ∂²f_d/(∂p_j1 ∂p_j2)
build_Hp_r <- function(H_p, K, D, J, mp1, V, U, tt) {
  function(p) {
    p_mat  <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input  <- rbind(p_mat, t(U), t(tt))
    Hp_F   <- array(H_p(input), c(mp1, J, J, D))
    result <- array(V %*% matrix(Hp_F, mp1, J * J * D), c(K, J, J, D))
    array(aperm(result, c(1L, 4L, 3L, 2L)), c(K * D, J, J))
  }
}

# L₀ where L(p) = L₁(p) + L₀
# einsum 'km,ab->kamb': result[k,a,m,b] = Vp[k,m] * sig_diag[a,b]
build_L0 <- function(K, D, mp1, Vp, sig) {
  sig_diag <- diag(as.numeric(sig), nrow = D)
  L0_      <- array(0, c(K, D, mp1, D))
  for (a in seq_len(D)) {
    for (b in seq_len(D)) {
      L0_[, a, , b] <- sig_diag[a, b] * Vp
    }
  }
  matrix(L0_, K * D, mp1 * D)
}

# Helper: compute L1 matrix (K*D x mp1*D) from J_F (mp1 x D x D), V, sig
# einsum 'mab,km,a->kamb': result[k,a,m,b] = J_F[m,a,b] * V[k,m] * sig[a]
.compute_L1_mat <- function(J_F, V, sig, K, D, mp1) {
  L1_ <- array(0, c(K, D, mp1, D))
  for (a in seq_len(D)) {
    # Use matrix() to guard against dimension dropping when D==1
    sig_JF_a <- sig[a] * matrix(J_F[, a, ], mp1, D)  # mp1 x D
    for (b in seq_len(D)) {
      # L1_[k,a,m,b] = V[k,m] * sig_JF_a[m,b]
      L1_[, a, , b] <- V * matrix(rep(sig_JF_a[, b], each = K), K, mp1)
    }
  }
  matrix(L1_, K * D, mp1 * D)
}

# L(p) where L(p)Lᵀ(p) = S(p), the covariance matrix of the weak residual
build_L <- function(U, tt, J_u, K, V, L0, sig, J) {
  D   <- ncol(U)
  mp1 <- length(tt)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_F   <- array(J_u(input), c(mp1, D, D))
    L1    <- .compute_L1_mat(J_F, V, sig, K, D, mp1)
    L1 + L0
  }
}

# L for linear in parameters where L(p)Lᵀ(p) = S(p)
build_L_linear <- function(U, tt, J_u, K, V, L0, sig, J) {
  D   <- ncol(U)
  mp1 <- length(tt)

  p_mat_affine <- matrix(0, nrow = J, ncol = mp1)
  input_affine <- rbind(p_mat_affine, t(U), t(tt))
  J_F_affine   <- array(J_u(input_affine), c(mp1, D, D))
  L1_affine    <- .compute_L1_mat(J_F_affine, V, sig, K, D, mp1)

  L1_ <- array(0, c(K * D, mp1 * D, J))
  for (j in seq_len(J)) {
    e_j   <- rep(0, J); e_j[j] <- 1
    p_mat <- matrix(rep(e_j, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_F   <- array(J_u(input), c(mp1, D, D))
    Jj_F  <- .compute_L1_mat(J_F, V, sig, K, D, mp1)
    L1_[, , j] <- Jj_F - L1_affine
  }

  function(p) {
    L1 <- matrix(matrix(L1_, K * D * mp1 * D, J) %*% p, K * D, mp1 * D)
    L1 + L1_affine + L0
  }
}

# ∇ₚL(p) Jacobian of the covariance factor
# J_up_sym has shape (J, D, D) in practice: J_up_sym[j,a,b] = ∂²f_a/(∂u_b ∂p_j)
# After evaluation: H_F[m,j,a,b] = ∂²f_a/(∂u_b ∂p_j)
# result[k,a,m,b,j] = H_F[m,j,a,b] * V[k,m] * sig[a]
build_Jp_L <- function(U, tt, J_up, K, J, D, V, sig) {
  mp1 <- length(tt)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    # Reshape as [m, j, a, b]: p-derivative index is first (column-major convention)
    H_F   <- array(J_up(input), c(mp1, J, D, D))
    Jp_L_ <- array(0, c(K, D, mp1, D, J))
    for (a in seq_len(D)) {
      for (b in seq_len(D)) {
        # H_F[m, j, a, b] for all j: use dims 2=j, 3=a, 4=b
        sig_HF_ab <- sig[a] * matrix(H_F[, , a, b], mp1, J)  # mp1 x J
        for (j in seq_len(J)) {
          Jp_L_[, a, , b, j] <- V * matrix(rep(sig_HF_ab[, j], each = K), K, mp1)
        }
      }
    }
    array(Jp_L_, c(K * D, mp1 * D, J))
  }
}

# ∇ₚL(p) Linear Jacobian of the Covariance factor L(p) = L₁p + L_affine + L0 → ∇ₚL(p) = L₁
build_Jp_L_linear <- function(U, tt, J_u, K, V, L0, sig, J) {
  D   <- ncol(U)
  mp1 <- length(tt)

  p_mat_affine <- matrix(0, nrow = J, ncol = mp1)
  input_affine <- rbind(p_mat_affine, t(U), t(tt))
  J_F_affine   <- array(J_u(input_affine), c(mp1, D, D))
  L1_affine    <- .compute_L1_mat(J_F_affine, V, sig, K, D, mp1)

  L1_ <- array(0, c(K * D, mp1 * D, J))
  for (j in seq_len(J)) {
    e_j   <- rep(0, J); e_j[j] <- 1
    p_mat <- matrix(rep(e_j, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    J_F   <- array(J_u(input), c(mp1, D, D))
    Jj_F  <- .compute_L1_mat(J_F, V, sig, K, D, mp1)
    L1_[, , j] <- Jj_F - L1_affine
  }

  function(p) L1_
}

# ∇ₚ∇ₚL(p) Hessian of the covariance factor
# J_upp_sym has shape (J,J,D,D): J_upp_sym[j1,j2,a,b] = ∂³f_a/(∂u_b ∂p_{j2} ∂p_{j1})
# After evaluation: T_F[m,j1,j2,a,b] = ∂³f_a/(∂u_b ∂p_{j2} ∂p_{j1})
# result[k,a,m,b,j1,j2] = T_F[m,j1,j2,a,b] * V[k,m] * sig[a]
build_Hp_L <- function(U, tt, J_upp, K, J, D, V, sig) {
  mp1 <- length(tt)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, t(U), t(tt))
    # Reshape as [m, j1, j2, a, b]: p-derivative indices are first
    T_F   <- array(J_upp(input), c(mp1, J, J, D, D))
    Hp_L_ <- array(0, c(K, D, mp1, D, J, J))
    for (a in seq_len(D)) {
      for (b in seq_len(D)) {
        for (j1 in seq_len(J)) {
          for (j2 in seq_len(J)) {
            sig_TF <- sig[a] * T_F[, j1, j2, a, b]  # mp1 vector
            Hp_L_[, a, , b, j1, j2] <- V * matrix(rep(sig_TF, each = K), K, mp1)
          }
        }
      }
    }
    array(Hp_L_, c(K * D, mp1 * D, J, J))
  }
}

# S(p) is the covariance matrix of the weak residual
build_S <- function(L, W = NULL, diag_reg = 10e-10) {
  function(p) {
    Lp     <- L(p)
    Lpt    <- t(Lp)
    S_     <- if (!is.null(W)) Lp %*% W %*% Lpt else Lp %*% Lpt
    WEIGHT <- 1.0 - diag_reg
    n      <- nrow(S_)
    WEIGHT * S_ + diag_reg * diag(n)
  }
}

# ∇ₚS(p) Jacobian of the covariance matrix
# einsum "ijk,jl->ilk": result[i,l,k] = sum_j Jp_Lp[i,j,k] * WLp_t[j,l]
build_J_S <- function(L, Jp_L, J, K, D, W = NULL, diag_reg = 1e-10) {
  WEIGHT <- 1.0 - diag_reg
  function(p) {
    Lp     <- L(p)
    Jp_Lp  <- Jp_L(p)
    WLp_t  <- if (!is.null(W)) W %*% t(Lp) else t(Lp)
    # For each j: prt[,,j] = Jp_Lp[,,j] %*% WLp_t
    KD     <- nrow(Lp)
    prt    <- array(0, c(KD, KD, J))
    for (k in seq_len(J)) {
      prt[, , k] <- Jp_Lp[, , k] %*% WLp_t
    }
    WEIGHT * (prt + aperm(prt, c(2L, 1L, 3L)))
  }
}

# Weak negative log likelihood  (Multivariate Gaussian distribution)
build_wnll <- function(S, g, b, K, D) {
  constant_term <- 0.5 * K * D * log(2 * pi)
  function(p) {
    Sp    <- S(p)
    r     <- g(p) - b
    cholF <- chol(Sp)
    log_det <- 2 * sum(log(diag(cholF)))
    S_invr  <- backsolve(cholF, forwardsolve(t(cholF), r))
    mdist   <- (r %*% S_invr)[1, 1]
    0.5 * (mdist + log_det) + constant_term
  }
}

# Jacobian of the weak form negative log likelihood
build_J_wnll <- function(S, Jp_S, Jp_r, g, b, J) {
  function(p) {
    Sp     <- S(p)
    J_Sp   <- Jp_S(p)
    J_rp   <- Jp_r(p)
    r      <- g(p) - b
    cholF  <- chol(Sp)
    S_inv_rp <- backsolve(cholF, forwardsolve(t(cholF), r))
    gradient <- numeric(J)
    for (j in seq_len(J)) {
      J_S_j <- J_Sp[, , j]
      J_r_j <- J_rp[, j]
      tmp       <- J_S_j %*% S_inv_rp
      prt0      <- 2.0 * (J_r_j %*% S_inv_rp)[1, 1]
      prt1      <- -1.0 * (S_inv_rp %*% tmp)[1, 1]
      fact      <- backsolve(cholF, forwardsolve(t(cholF), J_S_j))
      logDetPart <- sum(diag(fact))
      gradient[j] <- 0.5 * (prt0 + prt1 + logDetPart)
    }
    gradient
  }
}

# Hessian of the weak form negative log likelihood
build_H_wnll <- function(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J, W = NULL, diag_reg = 1e-10) {
  WEIGHT <- 1.0 - diag_reg
  function(p) {
    r      <- g(p) - b
    Jp_rp  <- Jp_r(p)
    Hp_rp  <- Hp_r(p)
    Lp     <- L(p)
    Jp_Lp  <- Jp_L(p)
    Hp_Lp  <- Hp_L(p)
    W_mat  <- W
    LpW_t  <- if (!is.null(W_mat)) W_mat %*% t(Lp) else t(Lp)
    Sp     <- S(p)
    Jp_Sp  <- Jp_S(p)

    S_inv_solve <- make_S_inv_solver(Sp)
    S_inv_rp    <- S_inv_solve(r)
    H_wnn       <- matrix(0, nrow = J, ncol = J)

    for (j in seq_len(J)) {
      Jp_rp_j    <- Jp_rp[, j]
      Jp_Sp_j    <- Jp_Sp[, , j]
      Jp_Lp_j    <- Jp_Lp[, , j]
      Jp_Lp_j_W  <- if (!is.null(W_mat)) Jp_Lp_j %*% W_mat else Jp_Lp_j
      shar_       <- S_inv_solve(Jp_Sp_j)

      for (i in j:J) {
        Jp_Sp_i  <- Jp_Sp[, , i]
        Jp_Lp_i  <- Jp_Lp[, , i]
        Jp_rp_i  <- Jp_rp[, i]
        term     <- Jp_Sp_i %*% shar_
        Hp_Lp_ji <- Hp_Lp[, , j, i]
        p1       <- Jp_Lp_j_W %*% t(Jp_Lp_i)
        p2       <- Hp_Lp_ji %*% LpW_t
        Hp_Sp_ji <- WEIGHT * (p1 + t(p1) + p2 + t(p2))
        Hp_rp_ji <- Hp_rp[, j, i]

        prt0 <- as.numeric(t(Hp_rp_ji) %*% S_inv_rp)
        prt1 <- -1.0 * as.numeric(t(S_inv_solve(Jp_rp_j)) %*% (Jp_Sp_i %*% S_inv_rp))
        inv_factor <- S_inv_solve(Jp_rp_i)
        prt2 <- as.numeric(t(Jp_rp_j) %*% inv_factor)
        prt3 <- -2.0 * as.numeric(t(inv_factor) %*% (Jp_Sp_j %*% S_inv_rp))
        prt4 <- -1.0 * as.numeric(t(S_inv_rp) %*% (Hp_Sp_ji %*% S_inv_rp))
        prt5 <-  2.0 * as.numeric(t(S_inv_rp) %*% (term %*% S_inv_rp))
        logDetTerm <- -1.0 * sum(diag(S_inv_solve(term))) + sum(diag(S_inv_solve(Hp_Sp_ji)))

        Hij <- 0.5 * (2 * (prt0 + prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)
        H_wnn[j, i] <- Hij
        if (i != j) H_wnn[i, j] <- Hij
      }
    }
    H_wnn
  }
}

# Hessian of the weak form negative log likelihood when linear in parameters
build_H_wnll_linear <- function(S, Jp_S, L, Jp_L, Jp_r, g, b, J, W = NULL, diag_reg = 1e-10) {
  WEIGHT <- 1.0 - diag_reg
  function(p) {
    r      <- g(p) - b
    Jp_rp  <- Jp_r(p)
    Lp     <- L(p)
    Jp_Lp  <- Jp_L(p)
    Sp     <- S(p)
    Jp_Sp  <- Jp_S(p)
    W_mat  <- W

    S_inv_solve <- make_S_inv_solver(Sp)
    S_inv_rp    <- S_inv_solve(r)
    H_wnn       <- matrix(0, nrow = J, ncol = J)

    for (j in seq_len(J)) {
      Jp_rp_j   <- Jp_rp[, j]
      Jp_Sp_j   <- Jp_Sp[, , j]
      Jp_Lp_j   <- Jp_Lp[, , j]
      Jp_Lp_j_W <- if (!is.null(W_mat)) Jp_Lp_j %*% W_mat else Jp_Lp_j
      shar_      <- S_inv_solve(Jp_Sp_j)

      for (i in j:J) {
        Jp_Sp_i  <- Jp_Sp[, , i]
        Jp_Lp_i  <- Jp_Lp[, , i]
        Jp_rp_i  <- Jp_rp[, i]
        term     <- Jp_Sp_i %*% shar_
        p1       <- Jp_Lp_j_W %*% t(Jp_Lp_i)
        Hp_Sp_ji <- WEIGHT * (p1 + t(p1))
        inv_factor <- S_inv_solve(Jp_rp_i)
        prt1 <- -1.0 * as.numeric(t(S_inv_solve(Jp_rp_j)) %*% (Jp_Sp_i %*% S_inv_rp))
        prt2 <-  as.numeric(t(Jp_rp_j) %*% inv_factor)
        prt3 <- -2.0 * as.numeric(t(inv_factor) %*% (Jp_Sp_j %*% S_inv_rp))
        prt4 <- -1.0 * as.numeric(t(S_inv_rp) %*% (Hp_Sp_ji %*% S_inv_rp))
        prt5 <-  2.0 * as.numeric(t(S_inv_rp) %*% (term %*% S_inv_rp))
        logDetTerm <- -1.0 * sum(diag(S_inv_solve(term))) + sum(diag(S_inv_solve(Hp_Sp_ji)))

        Hij <- 0.5 * (2 * (prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)
        H_wnn[j, i] <- Hij
        if (i != j) H_wnn[i, j] <- Hij
      }
    }
    H_wnn
  }
}

# --- Block-diagonal stacking helpers for multiple interpolants ---

# Assemble a constant block-diagonal matrix from a list of L0 matrices.
# @keywords internal
build_L0_block <- function(L0_list, K_list, D, mp1) {
  M       <- length(L0_list)
  K_total <- sum(K_list)
  out     <- matrix(0, K_total * D, M * mp1 * D)
  k_off   <- 0L; m_off <- 0L
  for (i in seq_along(L0_list)) {
    k_i <- K_list[i] * D; m_i <- mp1 * D
    out[(k_off + 1):(k_off + k_i), (m_off + 1):(m_off + m_i)] <- L0_list[[i]]
    k_off <- k_off + k_i; m_off <- m_off + m_i
  }
  out
}

# Returns p -> block-diagonal (K_total*D) x (M*mp1*D) matrix.
# @keywords internal
build_L_block <- function(L_fns, K_list, D, mp1) {
  M       <- length(L_fns)
  K_total <- sum(K_list)
  function(p) {
    mats  <- lapply(L_fns, function(lf) lf(p))
    out   <- matrix(0, K_total * D, M * mp1 * D)
    k_off <- 0L; m_off <- 0L
    for (i in seq_along(mats)) {
      k_i <- K_list[i] * D; m_i <- mp1 * D
      out[(k_off + 1):(k_off + k_i), (m_off + 1):(m_off + m_i)] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}

# Returns p -> block-diagonal (K_total*D) x (M*mp1*D) x J array.
# @keywords internal
build_Jp_L_block <- function(Jp_L_fns, K_list, D, mp1, J) {
  M       <- length(Jp_L_fns)
  K_total <- sum(K_list)
  function(p) {
    mats  <- lapply(Jp_L_fns, function(jf) jf(p))
    out   <- array(0, c(K_total * D, M * mp1 * D, J))
    k_off <- 0L; m_off <- 0L
    for (i in seq_along(mats)) {
      k_i <- K_list[i] * D; m_i <- mp1 * D
      out[(k_off + 1):(k_off + k_i), (m_off + 1):(m_off + m_i), ] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}

# Returns p -> block-diagonal (K_total*D) x (M*mp1*D) x J x J array.
# @keywords internal
build_Hp_L_block <- function(Hp_L_fns, K_list, D, mp1, J) {
  M       <- length(Hp_L_fns)
  K_total <- sum(K_list)
  function(p) {
    mats  <- lapply(Hp_L_fns, function(hf) hf(p))
    out   <- array(0, c(K_total * D, M * mp1 * D, J, J))
    k_off <- 0L; m_off <- 0L
    for (i in seq_along(mats)) {
      k_i <- K_list[i] * D; m_i <- mp1 * D
      out[(k_off + 1):(k_off + k_i), (m_off + 1):(m_off + m_i), , ] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}

# Robust solver for applying the inverse of S(p) to a vector or matrix
make_S_inv_solver <- function(Sp) {
  cholF            <- NULL
  qrF              <- NULL
  use_regularization <- FALSE
  diag_reg_mat     <- NULL

  cholF <- tryCatch(chol(Sp), error = function(e) NULL)

  if (is.null(cholF)) {
    qrF <- tryCatch(qr(Sp), error = function(e) NULL)
    if (is.null(qrF)) {
      use_regularization <- TRUE
      diag_reg_mat       <- 1e-12 * diag(nrow(Sp))
    }
  }

  function(x) {
    if (!is.null(cholF)) {
      return(backsolve(cholF, forwardsolve(t(cholF), x)))
    } else if (!is.null(qrF)) {
      return(solve(qrF, x))
    } else {
      return(solve(Sp + diag_reg_mat, x))
    }
  }
}
