# Weak residual components (column-major convention).
#
# Row layout for K*D-vectors and (K*D × mp1*D)-matrices:
#   row r = k + (d-1)*K, k fast.   (equivalently: as.vector of a (K, D) matrix)
# Column layout for (× mp1*D)-matrices:
#   col c = m + (b-1)*mp1, m fast. (equivalently: as.vector of an (mp1, D) matrix)
# Symbolic outputs (J_p, J_u, J_up, J_upp, J_pp) follow the same column-major
# convention after the loop-order fix in compute_symbolic_jacobian.

# F(p) ∈ ℝ^{mp1 × D}; columns are f_d(p, U, t) evaluated at the mp1 grid points.
build_F <- function(U, tt, f_, J) {
  mp1 <- nrow(U); D <- ncol(U)
  tU  <- t(U)
  ttt <- matrix(tt, nrow = 1L)
  function(p) {
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    f_(input)   # (mp1 × D)
  }
}

# g(p) = vec(V F(p)) ∈ ℝ^{K*D}
build_g <- function(V, F_) {
  function(p) as.vector(V %*% F_(p))
}

# G ∈ ℝ^{K*D × J} for the linear-in-p system Gp = b - g0.
build_G_matrix <- function(V, U, tt, F_, J){
  K   <- nrow(V); mp1 <- nrow(U); D <- ncol(U)
  G   <- matrix(0, nrow = K * D, ncol = J)
  g0  <- as.vector(V %*% F_(rep(0, J)))
  for (j in seq_len(J)) {
    e_j <- rep(0, J); e_j[j] <- 1
    G[, j] <- as.vector(V %*% F_(e_j)) - g0
  }
  G
}

# g(p) = G p
build_g_linear <- function(G) {
  function(p) as.vector(G %*% p)
}

# ∇ₚ r(p) ∈ ℝ^{K*D × J}, Jacobian of the weak residual.
build_Jp_r <- function(J_p, K, D, J, mp1, V, U, tt){
  tU  <- t(U); ttt <- matrix(tt, nrow = 1L)
  function(p){
    p_mat   <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input   <- rbind(p_mat, tU, ttt)
    Jp_eval <- J_p(input)                 # (mp1 × D*J), column-major (d, j) d fast
    matrix(V %*% Jp_eval, K * D, J)       # (K*D × J), rows (k, d) k fast
  }
}

# ∇ₚr(p) = G when linear in p.
build_Jp_r_linear <- function(G) {
  function(p) G
}

# ∇_p ∇_p r(p) ∈ ℝ^{K*D × J × J}, Hessian of the weak residual.
build_Hp_r <- function(H_p, K, D, J, mp1, V, U, tt){
  tU  <- t(U); ttt <- matrix(tt, nrow = 1L)
  function(p){
    p_mat  <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input  <- rbind(p_mat, tU, ttt)
    H_eval <- H_p(input)                                 # (mp1 × D*J*J)
    array(V %*% H_eval, c(K * D, J, J))                  # (K*D, J, J)
  }
}

# L0 ∈ ℝ^{K*D × mp1*D}. Block diagonal: block (a,a) = sig[a] * Vp.
build_L0 <- function(K, D, mp1, Vp, sig) {
  if (length(sig) == 1L) sig <- rep(sig, D)
  kronecker(diag(sig, nrow = D), Vp)
}

# L(p) where L(p)L(p)^T = S(p), covariance of the weak residual.
build_L <- function(U, tt, J_u, K, V, L0, sig, J){
  D   <- ncol(U); mp1 <- length(tt)
  tU  <- t(U); ttt <- matrix(tt, nrow = 1L)
  if (length(sig) == 1L) sig <- rep(sig, D)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    J_F   <- array(J_u(input), c(mp1, D, D))        # J_F[m, a, b] = ∂f_a/∂u_b @ m
    L1    <- matrix(0, K * D, mp1 * D)
    for (a in seq_len(D)) for (b in seq_len(D)) {
      rs <- (a - 1L) * K + seq_len(K)
      cs <- (b - 1L) * mp1 + seq_len(mp1)
      L1[rs, cs] <- sig[a] * sweep(V, 2, J_F[, a, b], "*")
    }
    L1 + L0
  }
}

# L(p) when f is linear in p: L(p) = sum_j p_j L1_j + L_affine + L0.
build_L_linear <- function(U, tt, J_u, K, V, L0, sig, J){
  D   <- ncol(U); mp1 <- length(tt)
  tU  <- t(U); ttt <- matrix(tt, nrow = 1L)
  if (length(sig) == 1L) sig <- rep(sig, D)

  build_L1_at <- function(p_val){
    p_mat <- matrix(rep(p_val, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    J_F   <- array(J_u(input), c(mp1, D, D))
    out   <- matrix(0, K * D, mp1 * D)
    for (a in seq_len(D)) for (b in seq_len(D)) {
      rs <- (a - 1L) * K + seq_len(K)
      cs <- (b - 1L) * mp1 + seq_len(mp1)
      out[rs, cs] <- sig[a] * sweep(V, 2, J_F[, a, b], "*")
    }
    out
  }

  L1_affine <- build_L1_at(rep(0, J))
  L1_basis  <- array(0, c(K * D, mp1 * D, J))
  for (j in seq_len(J)) {
    e_j <- rep(0, J); e_j[j] <- 1
    L1_basis[, , j] <- build_L1_at(e_j) - L1_affine
  }

  function(p) {
    L1 <- matrix(L1_basis, K * D * mp1 * D, J) %*% p
    matrix(L1, K * D, mp1 * D) + L1_affine + L0
  }
}

# ∇ₚ L(p) ∈ ℝ^{K*D × mp1*D × J}
build_Jp_L <- function(U, tt, J_up, K, J, D, V, sig){
  mp1 <- length(tt)
  tU  <- t(U); ttt <- matrix(tt, nrow = 1L)
  if (length(sig) == 1L) sig <- rep(sig, D)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    H_F   <- array(J_up(input), c(mp1, D, D, J))    # H_F[m, a, b, j] = ∂²f_a/∂u_b∂p_j
    out   <- array(0, c(K * D, mp1 * D, J))
    for (j in seq_len(J)) for (a in seq_len(D)) for (b in seq_len(D)) {
      rs <- (a - 1L) * K + seq_len(K)
      cs <- (b - 1L) * mp1 + seq_len(mp1)
      out[rs, cs, j] <- sig[a] * sweep(V, 2, H_F[, a, b, j], "*")
    }
    out
  }
}

# ∇ₚ L(p) = L1_basis (a constant tensor) when f is linear in p.
build_Jp_L_linear <- function(U, tt, J_u, K, V, L0, sig, J){
  D   <- ncol(U); mp1 <- length(tt)
  tU  <- t(U); ttt <- matrix(tt, nrow = 1L)
  if (length(sig) == 1L) sig <- rep(sig, D)

  build_L1_at <- function(p_val){
    p_mat <- matrix(rep(p_val, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    J_F   <- array(J_u(input), c(mp1, D, D))
    out   <- matrix(0, K * D, mp1 * D)
    for (a in seq_len(D)) for (b in seq_len(D)) {
      rs <- (a - 1L) * K + seq_len(K)
      cs <- (b - 1L) * mp1 + seq_len(mp1)
      out[rs, cs] <- sig[a] * sweep(V, 2, J_F[, a, b], "*")
    }
    out
  }

  L1_affine <- build_L1_at(rep(0, J))
  L1_basis  <- array(0, c(K * D, mp1 * D, J))
  for (j in seq_len(J)) {
    e_j <- rep(0, J); e_j[j] <- 1
    L1_basis[, , j] <- build_L1_at(e_j) - L1_affine
  }
  function(p) L1_basis
}

# ∇ₚ∇ₚ L(p) ∈ ℝ^{K*D × mp1*D × J × J}
build_Hp_L <- function(U, tt, J_upp, K, J, D, V, sig){
  mp1 <- length(tt)
  tU  <- t(U); ttt <- matrix(tt, nrow = 1L)
  if (length(sig) == 1L) sig <- rep(sig, D)
  function(p){
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    T_F   <- array(J_upp(input), c(mp1, D, D, J, J))   # T_F[m,a,b,j1,j2] = ∂³f_a/∂u_b∂p_j1∂p_j2
    out   <- array(0, c(K * D, mp1 * D, J, J))
    for (j1 in seq_len(J)) for (j2 in seq_len(J))
      for (a in seq_len(D)) for (b in seq_len(D)) {
        rs <- (a - 1L) * K + seq_len(K)
        cs <- (b - 1L) * mp1 + seq_len(mp1)
        out[rs, cs, j1, j2] <- sig[a] * sweep(V, 2, T_F[, a, b, j1, j2], "*")
      }
    out
  }
}

# S(p) = L(p) W L(p)^T + diag_reg I (W is a diagonal, stored as a length-mp1*D vector
# or NULL if uniform). 
build_S <- function(L, W = NULL, diag_reg = 1e-10) {
  WEIGHT    <- 1 - diag_reg
  eye_cache <- NULL
  function(p) {
    Lp <- L(p)
    Sp <- if (!is.null(W)) Lp %*% (W * t(Lp)) else tcrossprod(Lp)
    if (is.null(eye_cache)) eye_cache <<- diag_reg * diag(nrow(Sp))
    WEIGHT * Sp + eye_cache
  }
}

# ∇ₚ S(p) ∈ ℝ^{K*D × K*D × J}: ∂_j S = ∂_j L · W · L^T + (·)^T
build_J_S <- function(L, Jp_L, J, K, D, W = NULL, diag_reg = 1e-10){
  WEIGHT <- 1 - diag_reg
  function(p){
    Lp    <- L(p)
    Jp_Lp <- Jp_L(p)                                       # (K*D, mp1*D, J)
    WLp_t <- if (!is.null(W)) W * t(Lp) else t(Lp)         # (mp1*D, K*D)
    KD    <- nrow(Lp)
    prt   <- array(0, c(KD, KD, J))
    for (j in seq_len(J)) prt[, , j] <- Jp_Lp[, , j] %*% WLp_t
    WEIGHT * (prt + aperm(prt, c(2, 1, 3)))
  }
}

# Cached Cholesky solver for repeated S^{-1} applications.
# Returns a closure that applies S^{-1} to a vector or matrix.
make_S_inv_solver <- function(Sp) {
  U_ <- tryCatch(chol(Sp), error = function(e) NULL)
  if (!is.null(U_)) {
    function(x) backsolve(U_, forwardsolve(t(U_), x))
  } else {
    Sp_reg <- Sp + 1e-12 * diag(nrow(Sp))
    function(x) solve(Sp_reg, x)
  }
}

# Weak negative log-likelihood (multivariate Gaussian).
build_wnll <- function(S, g, b, K, D){
  const <- 0.5 * K * D * log(2 * pi)
  function(p){
    Sp <- S(p)
    r <- g(p) - b
    U_ <- chol(Sp)
    log_det <- 2 * sum(log(diag(U_)))
    S_invr <- backsolve(U_, forwardsolve(t(U_), r))
    0.5 * (sum(r * S_invr) + log_det) + const
  }
}

# Jacobian of the weak NLL.
build_J_wnll <- function(S, Jp_S, Jp_r, g, b, J){
  function(p){
    Sp   <- S(p)
    J_Sp <- Jp_S(p)
    J_rp <- Jp_r(p)
    r    <- g(p) - b
    KD   <- nrow(Sp)

    U_       <- chol(Sp)
    S_inv_r  <- backsolve(U_, forwardsolve(t(U_), r))

    # 2 (∂_j r)^T S^{-1} r → (J,)
    prt0 <- 2 * as.vector(crossprod(J_rp, S_inv_r))

    # (∂_j S) S^{-1} r for all j → (KD, J); then -(S^{-1} r)^T · that → (J,)
    tmp <- matrix(0, KD, J)
    for (j in seq_len(J)) tmp[, j] <- J_Sp[, , j] %*% S_inv_r
    prt1 <- -as.vector(crossprod(tmp, S_inv_r))

    # tr(S^{-1} ∂_j S) for all j: solve KD × (KD*J) system, reshape, take diagonals.
    J_Sp_2d    <- matrix(J_Sp, KD, KD * J)
    J_S_sol_2d <- backsolve(U_, forwardsolve(t(U_), J_Sp_2d))
    J_S_sol_3d <- array(J_S_sol_2d, c(KD, KD, J))
    logDet     <- vapply(seq_len(J),
                         function(j) sum(diag(J_S_sol_3d[, , j])),
                         numeric(1))

    as.numeric(0.5 * (prt0 + prt1 + logDet))
  }
}

# Hessian of the weak NLL.
build_H_wnll <- function(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J, W = NULL, diag_reg = 1e-10) {
  WEIGHT <- 1 - diag_reg
  function(p) {
    r <- g(p) - b
    Jp_rp <- Jp_r(p) # (KD, J)
    Hp_rp <- Hp_r(p) # (KD, J, J)
    Lp <- L(p) # (KD, mp1*D)
    Jp_Lp <- Jp_L(p) # (KD, mp1*D, J)
    Hp_Lp <- Hp_L(p) # (KD, mp1*D, J, J)
    Sp <- S(p)
    Jp_Sp <- Jp_S(p) # (KD, KD, J)

    LpW_t       <- if (!is.null(W)) W * t(Lp) else t(Lp)
    S_inv_solve <- make_S_inv_solver(Sp)
    S_inv_r     <- S_inv_solve(r)

    H_vals <- matrix(0, J, J)

    for (j in seq_len(J)) {
      Jp_rp_j <- Jp_rp[, j]
      Jp_Sp_j <- Jp_Sp[, , j]
      Jp_Lp_j <- Jp_Lp[, , j]
      Jp_Lp_j_W <- if (!is.null(W)) sweep(Jp_Lp_j, 2, W, "*") else Jp_Lp_j

      shar_          <- S_inv_solve(Jp_Sp_j)
      S_inv_rp_j     <- S_inv_solve(Jp_rp_j)
      Jp_Sp_j_S_invr <- Jp_Sp_j %*% S_inv_r

      for (i in j:J) {
        Jp_Sp_i <- Jp_Sp[, , i]
        Jp_Lp_i <- Jp_Lp[, , i]
        Jp_rp_i <- Jp_rp[, i]

        term     <- Jp_Sp_i %*% shar_
        Hp_Lp_ji <- Hp_Lp[, , j, i]

        p1       <- Jp_Lp_j_W %*% t(Jp_Lp_i)
        p2       <- Hp_Lp_ji  %*% LpW_t
        Hp_Sp_ji <- WEIGHT * (p1 + t(p1) + p2 + t(p2))

        Hp_rp_ji   <- Hp_rp[, j, i]
        inv_factor <- S_inv_solve(Jp_rp_i)

        prt0 <- sum(Hp_rp_ji * S_inv_r)
        prt1 <- -sum(S_inv_rp_j * (Jp_Sp_i %*% S_inv_r))
        prt2 <- sum(Jp_rp_j * inv_factor)
        prt3 <- -2 * sum(inv_factor * Jp_Sp_j_S_invr)
        prt4 <- -sum(S_inv_r * (Hp_Sp_ji %*% S_inv_r))
        prt5 <- 2 * sum(S_inv_r * (term %*% S_inv_r))

        logDetTerm <- -sum(diag(S_inv_solve(term))) + sum(diag(S_inv_solve(Hp_Sp_ji)))

        Hij <- 0.5 * (2 * (prt0 + prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)
        H_vals[j, i] <- Hij
        if (i != j) H_vals[i, j] <- Hij
      }
    }
    H_vals
  }
}

# Hessian of the weak NLL when linear in parameters (Hp_L = 0, Hp_r = 0).
build_H_wnll_linear <- function(S, Jp_S, L, Jp_L, Jp_r, g, b, J, W = NULL, diag_reg = 1e-10) {
  WEIGHT <- 1 - diag_reg
  function(p) {
    r     <- g(p) - b
    Jp_rp <- Jp_r(p)
    Jp_Lp <- Jp_L(p)
    Sp    <- S(p)
    Jp_Sp <- Jp_S(p)

    S_inv_solve <- make_S_inv_solver(Sp)
    S_inv_r     <- S_inv_solve(r)

    H_vals <- matrix(0, J, J)

    for (j in seq_len(J)) {
      Jp_rp_j <- Jp_rp[, j]
      Jp_Sp_j <- Jp_Sp[, , j]
      Jp_Lp_j <- Jp_Lp[, , j]
      Jp_Lp_j_W <- if (!is.null(W)) sweep(Jp_Lp_j, 2, W, "*") else Jp_Lp_j

      shar_          <- S_inv_solve(Jp_Sp_j)
      S_inv_rp_j     <- S_inv_solve(Jp_rp_j)
      Jp_Sp_j_S_invr <- Jp_Sp_j %*% S_inv_r

      for (i in j:J) {
        Jp_Sp_i <- Jp_Sp[, , i]
        Jp_Lp_i <- Jp_Lp[, , i]
        Jp_rp_i <- Jp_rp[, i]

        term     <- Jp_Sp_i %*% shar_
        p1       <- Jp_Lp_j_W %*% t(Jp_Lp_i)
        Hp_Sp_ji <- WEIGHT * (p1 + t(p1))

        inv_factor <- S_inv_solve(Jp_rp_i)

        prt1 <- -sum(S_inv_rp_j * (Jp_Sp_i %*% S_inv_r))
        prt2 <- sum(Jp_rp_j * inv_factor)
        prt3 <- -2 * sum(inv_factor * Jp_Sp_j_S_invr)
        prt4 <- -sum(S_inv_r * (Hp_Sp_ji %*% S_inv_r))
        prt5 <- 2 * sum(S_inv_r * (term %*% S_inv_r))

        logDetTerm <- -sum(diag(S_inv_solve(term))) + sum(diag(S_inv_solve(Hp_Sp_ji)))

        Hij <- 0.5 * (2 * (prt1 + prt2) + prt3 + prt4 + prt5 + logDetTerm)
        H_vals[j, i] <- Hij
        if (i != j) H_vals[i, j] <- Hij
      }
    }
    H_vals
  }
}

# ---- Block-diagonal stacking helpers for multiple interpolants ----
# Each L_m has shape (K_m*D, mp1*D). The block-diagonal assembly is
# (K_total*D, M*mp1*D), so L_block · L_block^T is block-diagonal in S_1,...,S_M.

# Constant block-diagonal L0 from a list of L0 matrices.
# @keywords internal
build_L0_block <- function(L0_list, K_list, D, mp1) {
  M       <- length(L0_list)
  K_total <- sum(K_list)
  out     <- matrix(0, K_total * D, M * mp1 * D)
  k_off   <- 0L; m_off <- 0L
  for (i in seq_along(L0_list)) {
    k_i <- K_list[i] * D; m_i <- mp1 * D
    out[(k_off + 1L):(k_off + k_i), (m_off + 1L):(m_off + m_i)] <- L0_list[[i]]
    k_off <- k_off + k_i; m_off <- m_off + m_i
  }
  out
}

# Returns p -> block-diagonal (K_total*D, M*mp1*D) matrix.
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
      out[(k_off + 1L):(k_off + k_i), (m_off + 1L):(m_off + m_i)] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}

# Returns p -> block-diagonal (K_total*D, M*mp1*D, J) array.
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
      out[(k_off + 1L):(k_off + k_i), (m_off + 1L):(m_off + m_i), ] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}

# Returns p -> block-diagonal (K_total*D, M*mp1*D, J, J) array.
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
      out[(k_off + 1L):(k_off + k_i), (m_off + 1L):(m_off + m_i), , ] <- mats[[i]]
      k_off <- k_off + k_i; m_off <- m_off + m_i
    }
    out
  }
}
