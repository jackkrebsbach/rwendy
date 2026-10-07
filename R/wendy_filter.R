# Leibniz expansion of g^(n)(t*) where g(t) = phi(t) F(p,u(t),t) + phi'(t) u(t)
# and trajectory derivatives F^(m) are passed as arguments

.binomial_cache <- new.env(parent = emptyenv())

g_choose <- function(n) {
  key <- as.character(n)
  binomial <- .binomial_cache[[key]]
  if (is.null(binomial)) {
    binomial <- choose(n, 0:n)
    .binomial_cache[[key]] <- binomial
  }
  binomial
}

# Coefficient matrix A ((order+2) x D) with g^(order) = sum_m phi^(m) A[m+1, ].
g_coeffs <- function(f_derivs, u, order) {
  binomial <- g_choose(order)
  A <- matrix(0, order + 2L, length(u))
  A[order + 2L, ] <- A[order + 2L, ] + u
  for (k in 0:order) {
    A[k + 1L, ] <- A[k + 1L, ] + binomial[k + 1L] * f_derivs[[order - k + 1L]]
  }
  if (order >= 1L) {
    for (k in 0:(order - 1L)) {
      A[k + 2L, ] <- A[k + 2L, ] + binomial[k + 1L] * f_derivs[[order - k]]
    }
  }
  A
}

g_deriv_at_endpoint <- function(phi_derivs, f_derivs, u, order) {
  binomial <- g_choose(order)
  result <- phi_derivs[[order + 2L]] * u
  for (k in 0:order) {
    result <- result + binomial[k + 1L] * phi_derivs[[k + 1L]] * f_derivs[[order - k + 1L]]
  }
  if (order >= 1L) {
    for (k in 0:(order - 1L)) {
      result <- result + binomial[k + 1L] * phi_derivs[[k + 2L]] * f_derivs[[order - k]]
    }
  }
  result
}

block_indices <- function(i, size) ((i - 1L) * size + 1L):(i * size)

# Phi0 = I_D ⊗ φ(t_0), filled block-wise (cheaper than kronecker).
build_Phi0 <- function(B, D) {
  K <- length(B)
  Phi0 <- matrix(0, K * D, D)
  for (d in seq_len(D)) {
    Phi0[block_indices(d, K), d] <- B
  }
  Phi0
}

# Euler-Maclaurin correction for the boundary-layer weak residual.
build_em_correction <- function(bl_phi_t1, bl_phi_tM, f_, dF_dt_, d2F_dt2_, d3F_dt3_, dt, scale = 1.0) {
  if (is.null(bl_phi_t1) || nrow(bl_phi_t1) == 0L) return(NULL)
  c2 <- scale * dt^2 / 12
  c4 <- scale * dt^4 / 720

  function(U, p, tt) {
    D <- ncol(U)
    M <- nrow(U)
    u_first <- as.numeric(U[1L, ])
    u_last  <- as.numeric(U[M, ])

    time_derivs <- function(u, t) {
      input <- matrix(c(p, u, t), ncol = 1L)
      list(
        as.vector(f_(input)),
        as.vector(dF_dt_(input)),
        as.vector(d2F_dt2_(input)),
        as.vector(d3F_dt3_(input))
      )
    }

    em_coeffs <- function(f_derivs, u, endpoint_sign) {
      A <- matrix(0, 5L, D)
      A[1:3, ] <- endpoint_sign * -c2 * g_coeffs(f_derivs, u, 1L)
      A + endpoint_sign * c4 * g_coeffs(f_derivs, u, 3L)
    }

    bl_phi_tM %*% em_coeffs(time_derivs(u_last, tt[M]), u_last, 1) +
      bl_phi_t1 %*% em_coeffs(time_derivs(u_first, tt[1L]), u_first, -1)
  }
}

# Analytic p-Jacobian of the boundary-layer EM correction
build_em_jacobian <- function(bl_phi_t1, bl_phi_tM, J_p, dF_dt_p_, d2F_dt2_p_, d3F_dt3_p_, dt, D, J, scale = 1.0) {
  if (is.null(bl_phi_t1) || nrow(bl_phi_t1) == 0L) return(NULL)
  K_bl <- nrow(bl_phi_t1)
  c2 <- scale * dt^2 / 12
  c4 <- scale * dt^4 / 720
  zero_state <- numeric(D)

  function(U, p, tt) {
    M <- nrow(U)

    parameter_jacobians <- function(u, t) {
      input <- matrix(c(p, u, t), ncol = 1L)
      list(matrix(as.vector(J_p(input)),        D, J),
           matrix(as.vector(dF_dt_p_(input)),   D, J),
           matrix(as.vector(d2F_dt2_p_(input)), D, J),
           matrix(as.vector(d3F_dt3_p_(input)), D, J))
    }

    jacobian_coeffs <- function(jacobians, endpoint_sign) {
      out <- matrix(0, 5L, D * J)
      for (j in seq_len(J)) {
        f_derivs <- lapply(jacobians, function(jacobian) jacobian[, j])
        A <- matrix(0, 5L, D)
        A[1:3, ] <- endpoint_sign * -c2 * g_coeffs(f_derivs, zero_state, 1L)
        out[, block_indices(j, D)] <- A + endpoint_sign * c4 * g_coeffs(f_derivs, zero_state, 3L)
      }
      out
    }

    jacobians_t1 <- parameter_jacobians(as.numeric(U[1L, ]), tt[1L])
    jacobians_tM <- parameter_jacobians(as.numeric(U[M, ]), tt[M])

    array(bl_phi_tM %*% jacobian_coeffs(jacobians_tM, 1) + bl_phi_t1 %*% jacobian_coeffs(jacobians_t1, -1),
          c(K_bl, D, J))
  }
}

build_ic_bl_system <- function(tt, r_bl, n_bl, orders = 0:4, include_interior = FALSE, interior_stride = 1L) {
  M <- length(tt)
  blocks <- lapply(orders, function(order) {
    build_boundary_layer_block(psi, tt, r_bl, order = order, side = "left", n_bl = n_bl)
  })
  n_interior <- 0L

  if (include_interior && (2L * r_bl + 1L) <= (M - 2L)) {
    interior <- lapply(0:1, function(order) build_test_function_matrix(psi, tt, r_bl, order = order))
    stride <- max(1L, as.integer(interior_stride))
    if (stride > 1L) {
      keep <- seq(1L, nrow(interior[[1]]), by = stride)
      interior <- lapply(interior, function(V) V[keep, , drop = FALSE])
    }
    n_interior <- nrow(interior[[1]])
    blocks <- lapply(seq_along(orders), function(i) {
      order <- orders[i]
      rbind(blocks[[i]],
            if (order <= 1L) interior[[order + 1L]] else matrix(0, n_interior, M))
    })
  }
  n_equations <- nrow(blocks[[1]])

  trapezoid <- function(V) {
    V[, 1] <- V[, 1] * 0.5
    V[, M] <- V[, M] * 0.5
    V
  }

  phi_t1 <- matrix(0, nrow = n_equations, ncol = length(orders))
  for (i in seq_along(orders)) {
    phi_t1[, i] <- blocks[[i]][, 1]
  }

  B <- phi_t1[, 1]
  list(V = trapezoid(blocks[[1]]), Vp = trapezoid(blocks[[2]]),
       B = B, BtB = sum(B * B),
       n_equations = n_equations, n_interior = n_interior,
       boundary_rows = which(rowSums(abs(phi_t1)) > 0),
       phi_t1 = phi_t1,
       support = which(colSums(abs(blocks[[1]]) + abs(blocks[[2]])) > 0))
}

build_ic_noise_sensitivity <- function(bl_system, U, tt, p, J_u, noise_sd, dt) {
  D <- ncol(U)
  K <- bl_system$n_equations
  support <- bl_system$support
  n_support <- length(support)

  input <- rbind(matrix(rep(p, n_support), nrow = length(p)),
                 t(U[support, , drop = FALSE]),
                 matrix(tt[support], nrow = 1L))
  jacobian <- matrix(as.vector(J_u(input)), nrow = n_support, ncol = D * D)

  V_support  <- bl_system$V[, support, drop = FALSE]
  Vp_support <- bl_system$Vp[, support, drop = FALSE]

  L <- matrix(0, K * D, n_support * D)
  column_variance <- numeric(n_support * D)
  for (a in seq_len(D)) {
    columns <- block_indices(a, n_support)
    column_variance[columns] <- noise_sd[a]^2
    for (d in seq_len(D)) {
      block <- V_support * rep(jacobian[, (a - 1L) * D + d], each = K)
      if (d == a) block <- block + Vp_support
      L[block_indices(d, K), columns] <- -dt * block
    }
  }
  list(L = L, column_variance = column_variance, Phi0 = build_Phi0(bl_system$B, D))
}

# Covariance of the QUADRATIC-in-noise part of the weak residual O(sigma^4)
build_ic_noise_quad <- function(bl_system, U, tt, p, J_u, noise_sd, dt, hess_cache = NULL) {
  D <- ncol(U)
  K <- bl_system$n_equations
  support <- bl_system$support
  state_hessian <- if (!is.null(hess_cache)) hess_cache else build_ic_hessian_cache(U, tt, p, J_u, D)

  V_support <- bl_system$V[, support, drop = FALSE]
  # Columns with V == 0 enter only through Vp, i.e. linearly: no curvature.
  curved <- which(colSums(abs(V_support)) > 0)
  variance_products <- outer(noise_sd^2, noise_sd^2)

  G <- array(0, c(length(support), D, D))
  for (i in curved) {
    H <- state_hessian(support[i])
    for (d1 in seq_len(D)) {
      for (d2 in d1:D) {
        value <- 0.5 * sum(variance_products * (matrix(H[d1, , ], D, D) * matrix(H[d2, , ], D, D)))
        G[i, d1, d2] <- value
        G[i, d2, d1] <- value
      }
    }
  }
  if (max(abs(G)) <= 0) return(NULL)

  Omega2 <- matrix(0, K * D, K * D)
  for (d1 in seq_len(D)) {
    for (d2 in d1:D) {
      block <- dt^2 * (V_support %*% (G[, d1, d2] * t(V_support)))
      Omega2[block_indices(d1, K), block_indices(d2, K)] <- block
      if (d2 != d1) Omega2[block_indices(d2, K), block_indices(d1, K)] <- t(block)
    }
  }
  Omega2
}

build_ic_gls_weights <- function(sensitivity, ridge) {
  L <- sensitivity$L
  Omega <- L %*% (sensitivity$column_variance * t(L))
  if (!is.null(sensitivity$Omega2)) Omega <- Omega + sensitivity$Omega2
  Omega_chol  <- chol(Omega + ridge * mean(diag(Omega)) * diag(nrow(L)))
  Omega_solve <- function(rhs) backsolve(Omega_chol, backsolve(Omega_chol, rhs, transpose = TRUE))
  Phi0tW      <- t(Omega_solve(sensitivity$Phi0))
  Phi0tWPhi0  <- Phi0tW %*% sensitivity$Phi0
  list(Omega_solve = Omega_solve, Phi0tW = Phi0tW, Phi0tWPhi0 = Phi0tWPhi0, cov_design = solve(Phi0tWPhi0))
}

# Memoizing state-Hessian cache: H[d, a, b] = d2 f_d / du_a du_b at sample m.
.IC_HESSIAN_CHUNK <- 64L

build_ic_hessian_cache <- function(U, tt, p, J_u, D) {
  store <- new.env(parent = emptyenv())
  M <- nrow(U)
  J <- length(p)

  fill <- function(samples) {
    n <- length(samples)
    U_samples <- U[samples, , drop = FALSE]
    h <- 1e-5 * pmax(1, sqrt(rowSums(U_samples^2)))
    p_rows <- matrix(rep(p, n), nrow = J)
    t_row <- matrix(tt[samples], nrow = 1L)
    jacobian_derivs <- vector("list", D)
    for (b in seq_len(D)) {
      U_up <- U_samples
      U_up[, b] <- U_up[, b] + h
      U_down <- U_samples
      U_down[, b] <- U_down[, b] - h
      jacobian_up   <- J_u(rbind(p_rows, t(U_up), t_row))
      jacobian_down <- J_u(rbind(p_rows, t(U_down), t_row))
      jacobian_derivs[[b]] <- (matrix(as.vector(jacobian_up), n, D * D) -
                               matrix(as.vector(jacobian_down), n, D * D)) / (2 * h)
    }
    for (i in seq_len(n)) {
      hessian <- array(0, c(D, D, D))
      for (b in seq_len(D)) {
        hessian[, , b] <- matrix(jacobian_derivs[[b]][i, ], D, D)
      }
      store[[as.character(samples[i])]] <- hessian
    }
  }

  function(m) {
    key <- as.character(m)
    if (is.null(store[[key]])) {
      samples <- seq.int(m, min(M, m + .IC_HESSIAN_CHUNK - 1L))
      samples <- samples[vapply(samples, function(k) is.null(store[[as.character(k)]]), TRUE)]
      fill(samples)
    }
    store[[key]]
  }
}

# Analytic O(sigma^2) bias of the feasible-GLS fixed point (both channels).
build_ic_bias_o2 <- function(bl_system, sensitivity, gls, P, EM_jacobian, U, tt, p, J_u, noise_sd, dt, hess_cache = NULL) {
  D <- length(noise_sd)
  K <- bl_system$n_equations
  support <- bl_system$support
  n_support <- length(support)
  L <- sensitivity$L
  column_variance <- sensitivity$column_variance
  noise_var <- noise_sd^2
  residual_jacobian <- sensitivity$Phi0 + EM_jacobian
  A <- gls$Omega_solve(L - residual_jacobian %*% (P %*% L))
  state_hessian <- if (!is.null(hess_cache)) hess_cache else build_ic_hessian_cache(U, tt, p, J_u, D)

  V_support <- bl_system$V[, support, drop = FALSE]
  curved <- which(colSums(abs(V_support)) > 0)

  Z  <- array(0, c(n_support, D, D))
  VA <- array(0, c(n_support, D, D))
  for (b in seq_len(D)) {
    A_b <- A[, block_indices(b, n_support), drop = FALSE]
    for (a in seq_len(D)) {
      Z[, a, b] <- colSums(L[, block_indices(a, n_support), drop = FALSE] * A_b)
    }
    for (d in seq_len(D)) {
      VA[, d, b] <- colSums(V_support * A_b[block_indices(d, K), , drop = FALSE])
    }
  }

  curvature      <- matrix(0, n_support, D)
  feedback_via_V <- matrix(0, D, n_support)
  feedback_via_L <- numeric(n_support * D)
  for (i in curved) {
    H <- state_hessian(support[i])
    for (d in seq_len(D)) {
      hessian_trace <- 0
      for (a in seq_len(D)) {
        hessian_trace <- hessian_trace + H[d, a, a] * noise_var[a]
      }
      curvature[i, d] <- -dt * 0.5 * hessian_trace
    }
    for (b in seq_len(D)) {
      weighted_Z <- column_variance[(seq_len(D) - 1L) * n_support + i] * Z[i, , b]
      for (d in seq_len(D)) {
        feedback_via_V[d, i] <- feedback_via_V[d, i] - dt * noise_var[b] * sum(weighted_Z * H[d, , b])
      }
      for (a in seq_len(D)) {
        column <- (a - 1L) * n_support + i
        feedback_via_L[column] <- feedback_via_L[column] -
          dt * noise_var[b] * column_variance[column] * sum(H[, a, b] * VA[i, , b])
      }
    }
  }

  feedback <- as.vector(L %*% feedback_via_L)
  for (d in seq_len(D)) {
    rows <- block_indices(d, K)
    feedback[rows] <- feedback[rows] + as.vector(V_support %*% feedback_via_V[d, ])
  }

  nonlinearity_bias <- as.numeric(P %*% as.vector(V_support %*% curvature))
  feedback_bias <- -as.numeric(P %*% feedback)
  list(nonlinearity = nonlinearity_bias, feedback = feedback_bias, total = nonlinearity_bias + feedback_bias)
}

# A-priori design selection: scores every candidate radius without solving.
#   obj = sum_d Var_d / sigma_d^2 + sum_d (truncation_bias_d + statistical_bias_d)^2 / sigma_d^2
select_ic_design <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u,
                             noise_sd, param_cov, em_order,
                             r_bl_grid, r_bl_max, include_interior = TRUE,
                             hess_cache = NULL, quad_cov = TRUE, ridge,
                             interior_stride = 1L, include_bias_o2 = TRUE) {
  radii <- sort(unique(pmin(pmax(as.integer(r_bl_grid), 2L), r_bl_max)))
  D <- ncol(U)
  M <- nrow(U)
  J <- length(p)
  noise_var <- noise_sd^2
  tt <- as.vector(tt)
  dt <- mean(diff(tt))
  c2 <- dt^2 / 12
  c4 <- if (em_order >= 4L) dt^4 / 720 else 0
  u0_observed <- as.numeric(U[1, ])

  em_parts <- function(bl_system, u0, params = p) {
    u0 <- as.vector(u0)
    input <- matrix(c(params, u0, tt[1]), ncol = 1L)
    f_derivs <- list(as.vector(f_(input)),       as.vector(dF_dt_(input)),
                     as.vector(d2F_dt2_(input)), as.vector(d3F_dt3_(input)))
    A2 <- matrix(0, 5L, D)
    A2[1:3, ] <- c2 * g_coeffs(f_derivs, u0, 1L)
    A4 <- if (c4 != 0) -c4 * g_coeffs(f_derivs, u0, 3L) else matrix(0, 5L, D)
    rows <- bl_system$boundary_rows
    phi <- bl_system$phi_t1[rows, , drop = FALSE]
    EM <- matrix(0, bl_system$n_equations, D)
    EM[rows, ] <- phi %*% (A2 + A4)
    Delta4 <- matrix(0, bl_system$n_equations, D)
    Delta4[rows, ] <- phi %*% A4
    list(EM = EM, Delta4 = Delta4)
  }

  trapezoid_residual <- function(bl_system, params) {
    input <- rbind(matrix(rep(params, M), nrow = J), t(U), matrix(tt, nrow = 1L))
    -dt * (bl_system$V %*% f_(input)) - dt * (bl_system$Vp %*% U)
  }

  rows <- vector("list", length(radii))
  best_objective <- Inf
  best_system <- NULL
  for (i in seq_along(radii)) {
    r_bl <- radii[i]
    n_bl <- as.integer(min(r_bl, r_bl_max))
    candidate <- tryCatch({
      bl_system <- build_ic_bl_system(tt, r_bl, n_bl, orders = 0:4,
                                      include_interior = include_interior,
                                      interior_stride = interior_stride)
      sensitivity <- build_ic_noise_sensitivity(bl_system, U, tt, p, J_u, noise_sd, dt)
      if (isTRUE(quad_cov)) {
        sensitivity$Omega2 <- tryCatch(
          build_ic_noise_quad(bl_system, U, tt, p, J_u, noise_sd, dt, hess_cache = hess_cache),
          error = function(err) NULL)
      }
      gls <- build_ic_gls_weights(sensitivity, ridge)

      h <- 1e-6 * max(1, sqrt(sum(u0_observed^2)))
      EM_jacobian <- matrix(0, bl_system$n_equations * D, D)
      for (d in seq_len(D)) {
        u0_up <- u0_observed
        u0_up[d] <- u0_up[d] + h
        u0_down <- u0_observed
        u0_down[d] <- u0_down[d] - h
        EM_jacobian[, d] <- as.vector((em_parts(bl_system, u0_up)$EM -
                                       em_parts(bl_system, u0_down)$EM) / (2 * h))
      }
      P <- solve(gls$Phi0tWPhi0 + gls$Phi0tW %*% EM_jacobian, gls$Phi0tW)

      G <- P %*% sensitivity$L
      cov_u0 <- G %*% (sensitivity$column_variance * t(G))
      if (!is.null(sensitivity$Omega2)) cov_u0 <- cov_u0 + (P %*% sensitivity$Omega2) %*% t(P)

      if (!is.null(param_cov)) {
        S_p <- matrix(0, D, J)
        for (j in seq_len(J)) {
          h_j <- 1e-6 * max(1, abs(p[j]))
          p_up <- p
          p_up[j] <- p_up[j] + h_j
          p_down <- p
          p_down[j] <- p_down[j] - h_j
          rhs_up   <- as.vector(trapezoid_residual(bl_system, p_up) - em_parts(bl_system, u0_observed, p_up)$EM)
          rhs_down <- as.vector(trapezoid_residual(bl_system, p_down) - em_parts(bl_system, u0_observed, p_down)$EM)
          S_p[, j] <- P %*% ((rhs_up - rhs_down) / (2 * h_j))
        }
        cov_u0 <- cov_u0 + S_p %*% param_cov %*% t(S_p)
      }

      statistical_bias <- if (isTRUE(include_bias_o2)) {
        bias <- tryCatch(
          build_ic_bias_o2(bl_system, sensitivity, gls, P, EM_jacobian, U, tt, p, J_u,
                           noise_sd, dt, hess_cache = hess_cache)$total,
          error = function(err) NULL)
        if (is.null(bias) || !all(is.finite(bias))) rep(0, D) else bias
      } else {
        rep(0, D)
      }

      truncation_bias <- if (c4 != 0) {
        as.numeric(P %*% as.vector(em_parts(bl_system, u0_observed)$Delta4))
      } else {
        rep(0, D)
      }

      variance_objective <- sum(diag(cov_u0) / noise_var)
      if (!is.finite(variance_objective)) stop("non-finite objective", call. = FALSE)
      list(obj = variance_objective + sum((truncation_bias + statistical_bias)^2 / noise_var),
           var_obj = variance_objective,
           ic_system = list(bl_system = bl_system, sensitivity = sensitivity, gls = gls))
    }, error = function(err) list(obj = NA_real_, var_obj = NA_real_, ic_system = NULL))

    rows[[i]] <- data.frame(r_bl = r_bl, n_bl = n_bl, obj = candidate$obj, var_obj = candidate$var_obj)

    if (is.finite(candidate$obj) && candidate$obj < best_objective) {
      best_objective <- candidate$obj
      best_system <- candidate$ic_system
    }
  }
  design_table <- do.call(rbind, rows)
  if (!any(is.finite(design_table$obj))) return(NULL)
  best <- which.min(design_table$obj)
  list(table = design_table, r_bl = design_table$r_bl[best], n_bl = design_table$n_bl[best],
       ic_system = best_system)
}

#' Estimate u(0) via iterative defect-correction on left BL test functions
#'
#' Solves the boundary-layer weak-form equations for the initial condition by
#' fixed-point iteration,
#' \deqn{\Phi_0 u_0^{(n+1)} = r_{trap} - \Delta_{EM}(u_0^{(n)}),}
#' where \eqn{r_{trap}} is the trapezoidal weak residual of the observed data
#' and \eqn{\Delta_{EM}} is the Euler-Maclaurin correction (order 2 or 4) at
#' the first time point.
#'
#' @details
#' The equations share noisy samples, so their errors have covariance
#' \eqn{\Omega = L \mathrm{diag}(\sigma^2) L^T + \Omega_2}, where \eqn{L} is the
#' noise sensitivity and \eqn{\Omega_2} the \eqn{O(\sigma^4)} quadratic-noise
#' term (\code{quad_cov}). With \code{inverse = "gls"} the equations are
#' combined by GLS,
#' \deqn{u_0 = (\Phi_0^T \Omega^{-1} \Phi_0)^{-1} \Phi_0^T \Omega^{-1} \mathrm{vec}(r),}
#' and with \code{inverse = "ols"} by unweighted least squares.
#'
#' When \code{n_bl} is \code{NULL} on the GLS path, the BL radius is chosen
#' from \code{r_bl_grid} without solving, by minimizing the \eqn{u_0} MSE proxy
#' \deqn{\sum_d \left(\mathrm{Var}_d + (\delta_d + b_d)^2\right) / \sigma_d^2,}
#' where \eqn{\delta} is the EM(2) minus EM(4) truncation estimate and \eqn{b}
#' the \eqn{O(\sigma^2)} bias, both evaluated at \code{U[1, ]}.
#'
#' With \code{debias = TRUE}, the analytic \eqn{O(\sigma^2)} bias of the GLS
#' estimate is subtracted unless any component exceeds twice its standard
#' error.
#'
#' @param U Numeric matrix (M x D) of observed states.
#' @param f_,dF_dt_,d2F_dt2_,d3F_dt3_ Callable RHS and total time-derivative
#'   evaluators built from the symbolic engine.
#' @param tt Numeric vector (length M) of time points.
#' @param p Numeric parameter vector \eqn{\hat\theta} (held fixed).
#' @param J_u Callable state Jacobian, with entry \eqn{[a, b] = \partial f_a /
#'   \partial u_b}.
#' @param sigma Noise standard deviation, scalar or length D. If invalid, the
#'   GLS path and noise-propagated covariance are unavailable.
#' @param param_cov Optional J x J covariance of \code{p}. When given,
#'   \code{cov_u0} includes the parameter channel \eqn{S_p \hat C S_p^T}.
#' @param n_bl Number of BL test functions. If \code{NULL}, the GLS path
#'   selects the radius from \code{r_bl_grid} with \code{n_bl = r_bl}, and the
#'   OLS path uses \code{max(3, ceiling(r_bl / 8))}.
#' @param r_bl BL radius used when no selection runs. Defaults to
#'   \code{min(16, floor((M - 1) / 2))}.
#' @param max_iter,tol Fixed-point iteration limit and relative step
#'   tolerance; converged when the step is below
#'   \code{tol * max(1, ||u0||)}.
#' @param em_order Euler-Maclaurin order, 2 or 4.
#' @param inverse \code{"gls"} (default) or \code{"ols"} combination of the BL
#'   equations. GLS falls back to OLS when its weights cannot be built.
#' @param include_interior Add interior test functions as GLS control
#'   variates. GLS only.
#' @param interior_stride Keep every s-th interior test function.
#' @param quad_cov Include the quadratic-noise term \eqn{\Omega_2} in the
#'   weights and \code{cov_u0}. GLS only.
#' @param ridge Relative ridge added to \eqn{\Omega} before factorizing,
#'   \code{ridge * mean(diag(Omega))}. \eqn{\Omega} is near-singular, so this
#'   sets the weight of its near-null directions and can shift \code{u0hat}
#'   by a fraction of its SE. GLS only.
#' @param r_bl_grid Candidate radii for the selection. Defaults to
#'   \code{c(4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 80, 100)}, capped at
#'   \code{floor((M - 1) / 2)}.
#' @param debias Subtract the analytic \eqn{O(\sigma^2)} bias (see Details).
#'   GLS only.
#' @param hess_cache Optional state-Hessian cache from
#'   \code{build_ic_hessian_cache}; built internally if \code{NULL}.
#' @return A list with:
#' \describe{
#'   \item{u0hat, U_hat}{The estimated \eqn{u_0}, and \code{U} with its first
#'     row replaced by it.}
#'   \item{cov_u0}{Covariance of \code{u0hat}: noise propagation (or the LS
#'     residual variance when \code{sigma} is invalid), plus the parameter
#'     channel when \code{param_cov} is given. \code{NULL} if the iteration
#'     diverged.}
#'   \item{cov_u0_resid, cov_u0_param}{The LS residual covariance and the
#'     parameter-channel covariance on their own.}
#'   \item{cov_method}{\code{"noise_propagation"} or \code{"ls_residual"},
#'     with \code{"+param"} when the parameter channel is included, or
#'     \code{"diverged"}.}
#'   \item{inverse}{The combination actually used; \code{"gls"} falls back to
#'     \code{"ols"} when the weights cannot be built.}
#'   \item{design}{The radius-selection table (\code{r_bl}, \code{n_bl},
#'     \code{obj}, \code{var_obj}), or \code{NULL} if no selection ran.}
#'   \item{bias_o2, debias_applied}{The \eqn{O(\sigma^2)} bias estimate and
#'     whether it was subtracted from \code{u0hat}.}
#'   \item{fallback}{\code{TRUE} when the design was degenerate and
#'     \code{u0hat} is \code{U[1, ]}.}
#'   \item{iters, converged, diverged, u0_history}{Fixed-point diagnostics.}
#'   \item{r_bl, n_bl, K_bl, K_int, em_order}{The design used: BL radius, BL
#'     test-function count, total and interior equation counts, and EM order.}
#' }
#' @export
estimate_IC <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, sigma,
                        param_cov        = NULL,
                        n_bl             = NULL,
                        r_bl             = NULL,
                        max_iter         = 100L,
                        tol              = 1e-9,
                        em_order         = c(4L, 2L),
                        inverse          = c("gls", "ols"),
                        include_interior = TRUE,
                        interior_stride  = 1L,
                        quad_cov         = TRUE,
                        ridge            = 1e-11,
                        r_bl_grid        = NULL,
                        debias           = TRUE,
                        hess_cache       = NULL) {

  em_order <- as.integer(em_order[1])
  if (!em_order %in% c(2L, 4L)) {
    stop("em_order must be 2 or 4", call. = FALSE)
  }
  inverse <- match.arg(inverse)

  M  <- nrow(U)
  D  <- ncol(U)
  J  <- length(p)
  tt <- as.vector(tt)
  dt <- mean(diff(tt))

  r_bl_max   <- floor((M - 1L) / 2L)
  r_bl_fixed <- if (!is.null(r_bl)) min(as.integer(r_bl), r_bl_max) else min(16L, r_bl_max)
  r_bl       <- r_bl_fixed

  use_noise_propagation <- length(sigma) %in% c(1L, D) && all(is.finite(sigma))
  noise_sd <- if (use_noise_propagation) {
    if (length(sigma) == 1L) rep(sigma, D) else as.numeric(sigma)
  } else {
    NULL
  }

  if (inverse == "gls" && !use_noise_propagation) inverse <- "ols"

  if (is.null(hess_cache) && use_noise_propagation) {
    hess_cache <- build_ic_hessian_cache(U, tt, p, J_u, D)
  }

  design_table    <- NULL
  selected_system <- NULL
  if (inverse == "gls" && is.null(n_bl)) {
    if (is.null(r_bl_grid)) {
      r_bl_grid <- c(4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 80, 100)
    }
    selection <- tryCatch(
      select_ic_design(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u,
                       noise_sd, param_cov, em_order, r_bl_grid, r_bl_max,
                       include_interior = include_interior,
                       hess_cache = hess_cache, quad_cov = quad_cov, ridge = ridge,
                       interior_stride = interior_stride),
      error = function(err) NULL)
    if (!is.null(selection)) {
      design_table <- selection$table
      r_bl <- selection$r_bl
      n_bl <- selection$n_bl
      selected_system <- selection$ic_system
    }
  }

  n_bl <- if (!is.null(n_bl))        max(1L, as.integer(n_bl))
          else if (inverse == "ols") max(3L, as.integer(ceiling(r_bl / 8)))
          else                       as.integer(min(r_bl, r_bl_max))

  build_system <- function(r_bl, n_bl) {
    bl_system <- build_ic_bl_system(tt, r_bl, n_bl, orders = 0:4,
                                    include_interior = isTRUE(include_interior) && inverse == "gls",
                                    interior_stride = interior_stride)
    sensitivity <- if (use_noise_propagation) {
      tryCatch(build_ic_noise_sensitivity(bl_system, U, tt, p, J_u, noise_sd, dt),
               error = function(err) NULL)
    } else {
      NULL
    }
    if (isTRUE(quad_cov) && inverse == "gls" && !is.null(sensitivity)) {
      sensitivity$Omega2 <- tryCatch(
        build_ic_noise_quad(bl_system, U, tt, p, J_u, noise_sd, dt, hess_cache = hess_cache),
        error = function(err) NULL)
    }
    gls <- if (inverse == "gls" && !is.null(sensitivity)) {
      tryCatch(build_ic_gls_weights(sensitivity, ridge), error = function(err) NULL)
    } else {
      NULL
    }
    list(bl_system = bl_system, sensitivity = sensitivity, gls = gls)
  }

  ic_system <- if (!is.null(selected_system)) selected_system else build_system(r_bl, n_bl)
  if (inverse == "gls" && is.null(ic_system$gls)) {
    inverse <- "ols"
    if (!is.null(design_table)) {
      design_table <- NULL
      r_bl <- r_bl_fixed
      n_bl <- max(3L, as.integer(ceiling(r_bl / 8)))
    }
    ic_system <- build_system(r_bl, n_bl)
  }

  bl_system <- ic_system$bl_system
  sensitivity <- ic_system$sensitivity
  gls <- ic_system$gls
  V <- bl_system$V
  Vp <- bl_system$Vp
  phi_t1 <- bl_system$phi_t1
  B <- bl_system$B
  BtB <- bl_system$BtB
  K_bl <- bl_system$n_equations
  K_int <- bl_system$n_interior
  boundary_rows <- bl_system$boundary_rows

  if (!is.finite(BtB) || BtB < .Machine$double.eps) {
    u0_observed <- as.numeric(U[1, ])
    U_hat <- U
    U_hat[1, ] <- u0_observed
    return(list(U_hat = U_hat, u0hat = u0_observed, cov_u0 = NULL,
                iters = 0L, converged = FALSE, diverged = FALSE,
                u0_history = matrix(u0_observed, nrow = 1),
                r_bl = r_bl, n_bl = n_bl, K_bl = K_bl, K_int = K_int,
                em_order = em_order,
                inverse = inverse, design = design_table,
                bias_o2 = NULL, debias_applied = FALSE, fallback = TRUE))
  }

  trapezoid_residual <- function(params = p) {
    input <- rbind(matrix(rep(params, M), nrow = J), t(U), matrix(tt, nrow = 1L))
    -dt * (V %*% f_(input)) - dt * (Vp %*% U)
  }
  r_trap <- trapezoid_residual()

  c2 <- dt^2 / 12
  c4 <- if (em_order >= 4L) dt^4 / 720 else 0

  em_correction <- function(u0, params = p, include_em4 = TRUE) {
    u0 <- as.vector(u0)
    input <- matrix(c(params, u0, tt[1]), ncol = 1L)
    f_derivs <- list(
      as.vector(f_(input)),
      as.vector(dF_dt_(input)),
      as.vector(d2F_dt2_(input)),
      as.vector(d3F_dt3_(input))
    )
    A <- matrix(0, 5L, D)
    A[1:3, ] <- c2 * g_coeffs(f_derivs, u0, 1L)
    if (include_em4 && c4 != 0) A <- A - c4 * g_coeffs(f_derivs, u0, 3L)
    EM <- matrix(0, nrow = K_bl, ncol = D)
    EM[boundary_rows, ] <- phi_t1[boundary_rows, , drop = FALSE] %*% A
    EM
  }

  project <- if (!is.null(gls)) {
    function(rhs) as.numeric(solve(gls$Phi0tWPhi0, gls$Phi0tW %*% as.vector(rhs)))
  } else {
    function(rhs) as.numeric(crossprod(B, rhs) / BtB)
  }

  residual_norm <- function(u0, include_em4 = TRUE) {
    residual <- outer(B, as.numeric(u0)) - (r_trap - em_correction(u0, include_em4 = include_em4))
    weighted <- if (!is.null(gls)) as.numeric(gls$Phi0tW %*% as.vector(residual))
                else               as.numeric(crossprod(B, residual))
    sqrt(sum(weighted * weighted))
  }

  run_fixed_point <- function(include_em4 = TRUE) {
    u0            <- project(r_trap)
    history       <- list(u0)
    best_u0       <- u0
    best_residual <- if (all(is.finite(u0))) residual_norm(u0, include_em4) else Inf
    iters         <- 0L
    converged     <- FALSE
    diverged      <- FALSE

    for (iter in seq_len(max_iter)) {
      iters   <- iter
      u0_next <- project(r_trap - em_correction(u0, include_em4 = include_em4))
      history[[iter + 1L]] <- u0_next

      if (!all(is.finite(u0_next))) {
        diverged <- TRUE
        break
      }

      residual <- residual_norm(u0_next, include_em4)
      if (is.finite(residual) && residual < best_residual) {
        best_residual <- residual
        best_u0       <- u0_next
      }

      step <- sqrt(sum((u0_next - u0)^2))
      u0   <- u0_next
      if (is.finite(step) && step < tol * max(1, sqrt(sum(u0_next^2)))) {
        converged <- TRUE
        break
      }
    }

    list(u0 = best_u0, iters = iters, converged = converged,
         diverged = diverged, history = history)
  }

  fit <- run_fixed_point()
  u0  <- fit$u0

  EM_jacobian <- tryCatch({
    jacobian <- matrix(0, K_bl * D, D)
    h <- 1e-6 * max(1, sqrt(sum(u0^2)))
    for (d in seq_len(D)) {
      u0_up <- u0
      u0_up[d] <- u0_up[d] + h
      u0_down <- u0
      u0_down[d] <- u0_down[d] - h
      jacobian[, d] <- as.vector((em_correction(u0_up) - em_correction(u0_down)) / (2 * h))
    }
    jacobian
  }, error = function(err) NULL)

  P <- if (!is.null(EM_jacobian)) {
    tryCatch({
      if (!is.null(gls)) {
        solve(gls$Phi0tWPhi0 + gls$Phi0tW %*% EM_jacobian, gls$Phi0tW)
      } else {
        Phi0 <- if (!is.null(sensitivity)) sensitivity$Phi0 else build_Phi0(B, D)
        solve(BtB * diag(D) + crossprod(Phi0, EM_jacobian), t(Phi0))
      }
    }, error = function(err) NULL)
  } else {
    NULL
  }

  cov_u0_noise <- if (!is.null(sensitivity) && !is.null(P)) {
    tryCatch({
      G <- P %*% sensitivity$L
      covariance <- G %*% (sensitivity$column_variance * t(G))
      if (!is.null(sensitivity$Omega2)) covariance <- covariance + (P %*% sensitivity$Omega2) %*% t(P)
      covariance
    }, error = function(err) NULL)
  } else {
    NULL
  }

  cov_u0_param <- if (!is.null(param_cov) && !is.null(P)) {
    tryCatch({
      rhs <- function(params) as.vector(trapezoid_residual(params) - em_correction(u0, params))
      S_p <- matrix(0, D, J)
      for (j in seq_len(J)) {
        h <- 1e-6 * max(1, abs(p[j]))
        p_up <- p
        p_up[j] <- p_up[j] + h
        p_down <- p
        p_down[j] <- p_down[j] - h
        S_p[, j] <- P %*% ((rhs(p_up) - rhs(p_down)) / (2 * h))
      }
      S_p %*% param_cov %*% t(S_p)
    }, error = function(err) NULL)
  } else {
    NULL
  }

  cov_u0_resid <- tryCatch({
    residual <- outer(B, u0) - (r_trap - em_correction(u0))
    # Interior rows have B = 0, so they only inflate the residual variance.
    residual <- residual[boundary_rows, , drop = FALSE]
    dof <- max(length(boundary_rows) - 1L, 1L)
    diag(colSums(residual * residual) / dof / BtB, nrow = D, ncol = D)
  }, error = function(err) NULL)

  cov_u0     <- if (!is.null(cov_u0_noise)) cov_u0_noise else cov_u0_resid
  cov_method <- if (!is.null(cov_u0_noise)) "noise_propagation" else "ls_residual"
  if (!is.null(cov_u0) && !is.null(cov_u0_param)) {
    cov_u0     <- cov_u0 + cov_u0_param
    cov_method <- paste0(cov_method, "+param")
  }

  if (fit$diverged) {
    cov_u0       <- NULL
    cov_u0_param <- NULL
    cov_method   <- "diverged"
  }

  bias_o2 <- NULL
  debias_applied <- FALSE
  if (isTRUE(debias) && !fit$diverged && !is.null(gls) && !is.null(sensitivity) &&
      !is.null(P) && !is.null(EM_jacobian) && !is.null(cov_u0_noise)) {
    bias_o2 <- tryCatch(
      build_ic_bias_o2(bl_system, sensitivity, gls, P, EM_jacobian, U, tt, p, J_u,
                       noise_sd, dt, hess_cache = hess_cache)$total,
      error = function(err) NULL)
    if (!is.null(bias_o2) && all(is.finite(bias_o2))) {
      se <- sqrt(pmax(diag(cov_u0_noise), 0))
      if (all(abs(bias_o2) <= 2 * se)) {
        u0 <- u0 - bias_o2
        debias_applied <- TRUE
      }
    }
  }

  U_hat <- U
  U_hat[1, ] <- u0

  list(
    U_hat          = U_hat,
    u0hat          = u0,
    cov_u0         = cov_u0,
    cov_u0_resid   = cov_u0_resid,
    cov_u0_param   = cov_u0_param,
    cov_method     = cov_method,
    inverse        = inverse,
    design         = design_table,
    bias_o2        = bias_o2,
    debias_applied = debias_applied,
    fallback       = FALSE,
    iters          = fit$iters,
    converged      = fit$converged,
    diverged       = fit$diverged,
    u0_history     = do.call(rbind, fit$history),
    r_bl           = r_bl,
    n_bl           = n_bl,
    K_bl           = K_bl,
    K_int          = K_int,
    em_order       = em_order
  )
}

#' Estimate the state using the RTS smoother
#'
#' Extended Kalman filter forward pass followed by a Rauch-Tung-Striebel
#' backward smoother on the parameter-conditioned dynamics.
#'
#' The reported posterior covariance follows the law of total variance,
#' \deqn{\mathrm{Cov}(u_k^\star) = \mathrm{Cov}(u_k^\star \mid \hat p)
#'        + S_k \, \hat C \, S_k^T,}
#' splitting trajectory uncertainty into the conditional smoother posterior
#' (data noise + model/discretization error) and the contribution of parameter
#' uncertainty. The conditional pass injects the RK4 \emph{discretization} error
#' as predict-step process noise, estimated per step by step-doubling
#' (Richardson): with \eqn{y_1} the single full RK4 step and \eqn{y_2} two half
#' steps, the local truncation error of \eqn{y_1} is
#' \eqn{e_k = (16/15)(y_1 - y_2)} and \eqn{Q_k = \mathrm{diag}(e_k^2)}. This is
#' \eqn{\sigma}-independent (discretization error depends on \eqn{\Delta t} and
#' the dynamics, not the measurement noise), vanishes as \eqn{O(T\,\Delta t^4)}
#' under grid refinement, and is active only when the grid under-resolves
#' \eqn{f}; it replaces the earlier ad-hoc \eqn{Q = (0.1\,\sigma)^2 I_D}.
#' Parameter uncertainty is folded in separately
#' and \emph{coherently} via the sensitivity
#' \eqn{S_k = \partial u_k^\star / \partial \hat p}, obtained by central
#' differences of the smoothed mean over re-runs at \eqn{\hat p \pm h e_j}, with
#' \eqn{\hat C} the WENDy parameter covariance. This replaces the earlier
#' heuristic of injecting \eqn{(\Delta t_k)^2 \nabla_p f \,\hat C\, \nabla_p f^T}
#' as predict-step process noise, which mis-modelled \eqn{\delta p} as an
#' independent per-step draw (a random walk) rather than a single
#' fixed-but-uncertain value affecting the whole trajectory coherently. Folding
#' requires \code{param_cov} and \code{J_p}; when either is absent (or
#' \code{fold_param_uncertainty = FALSE}) the conditional posterior is returned
#' alone.
#'
#' @param U Numeric matrix (mp1 x D) of noisy observations.
#' @param f_ Callable f(p, u, t) RHS evaluator built from the symbolic engine.
#' @param J_u Callable Jacobian d f / d u evaluator from the symbolic engine.
#' @param tt Numeric vector of time points (length mp1).
#' @param p Numeric vector of parameter estimates.
#' @param test_function_params Unused; retained for API compatibility.
#' @param sigma Optional scalar or per-state noise SD; if NULL it is estimated
#'   from U via \code{estimate_std}.
#' @param u0_init Optional length-D vector to seed the filter at \code{tt[1]}
#'   (e.g. \code{estimate_IC()$u0hat}). Defaults to \code{U[1, ]}.
#' @param P0_init Optional D x D prior covariance for \code{u0_init}. When
#'   \code{u0_init} is supplied, defaults to \code{sigma^2 I_D}. When neither is
#'   supplied the prior is \emph{diffuse} (\code{1e6 * sigma^2 I_D}), so the
#'   \eqn{k = 1} update returns \code{U[1, ]} with covariance \code{sigma^2 I_D}
#'   rather than double-counting \code{U[1, ]} as both prior mean and
#'   measurement (which would report \code{sigma^2 / 2}).
#' @param param_cov Optional J x J parameter covariance \eqn{\hat{C}}.
#' @param J_p Optional callable returning the D x J Jacobian
#'   \eqn{\nabla_p f(p, u, t)} (flattened d-fast, like the rest of the package).
#' @param fold_param_uncertainty Logical (default \code{TRUE}). When \code{TRUE}
#'   and \code{param_cov}/\code{J_p} are supplied, the parameter-sensitivity term
#'   \eqn{S_k \hat C S_k^T} is added to \code{P_smooth}; otherwise only the
#'   conditional posterior is returned.
#' @return Named list with \code{U_star} (smoothed state), \code{P_smooth}
#'   (total posterior covariance), \code{P_smooth_cond} (the conditional-on-p̂
#'   posterior), \code{P_smooth_param} (the parameter-uncertainty contribution,
#'   or \code{NULL} when not folded), and the intermediate filter/predictor
#'   states.
#' @export
wendy_erts <- function(U, f_, J_u, tt, p, test_function_params, sigma = NULL, u0_init = NULL, P0_init = NULL, param_cov = NULL, J_p = NULL,
                       fold_param_uncertainty = TRUE) {
  tt   <- as.vector(tt)
  mp1  <- nrow(U)
  D    <- ncol(U)
  J    <- length(p)
  I_D  <- diag(D)

  noise_sd <- as.numeric(if (!is.null(sigma)) sigma else estimate_std(U, k = 6))
  noise_sd <- mean(noise_sd)
  R_obs    <- noise_sd^2 * I_D

  # Conditional-pass process noise = the RK4 DISCRETIZATION error of the predict
  # step, estimated per step by step-doubling (Richardson): with y1 the single
  # full RK4 step and y2 two half steps, the local truncation error of y1 is
  # e_k = (16/15)(y1 - y2) (order p = 4, factor 2^p/(2^p - 1)), so Q_k = diag(e_k^2)
  q_floor <- 1e-12

  # One RK4 step of the process f(p, ., .)
  rk4_step <- function(p_use, uk, tk, h) {
    k1 <- as.vector(f_(matrix(c(p_use, uk,           tk        ), ncol = 1)))
    k2 <- as.vector(f_(matrix(c(p_use, uk + .5*h*k1, tk + .5*h ), ncol = 1)))
    k3 <- as.vector(f_(matrix(c(p_use, uk + .5*h*k2, tk + .5*h ), ncol = 1)))
    k4 <- as.vector(f_(matrix(c(p_use, uk +    h*k3, tk +    h ), ncol = 1)))
    uk + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)
  }

  # With no prior supplied the filter has no information about u(0) independent
  # of the data, so the honest encoding is a DIFFUSE prior. Taking
  # P0 = noise_sd^2 I instead would double-count U[1, ] -- it would be both the
  # prior mean and the k = 1 measurement, two sources the Kalman update treats
  # as independent when they are the identical number. That gives K0 = I/2 and
  # P_filt[1] = noise_sd^2 / 2, half the variance one noisy sample supports.
  # With P0 = diffuse_scale * R_obs the k = 1 update returns exactly
  # (U[1, ], R_obs) up to O(1/diffuse_scale). Validated in
  # examples/validation/ic_kf_diffuse_prior.R: P_filt[1] doubles as predicted;
  # the RTS pass dilutes it to nothing where later data pins u(0) (logistic/LV
  # SE ratio ~1.00) but on chaos it lifts u(0) coverage .61-.69 -> .75-.78.
  diffuse_scale <- 1e6
  u0 <- if (!is.null(u0_init)) as.vector(u0_init) else as.vector(U[1, ])
  P0 <- if (!is.null(P0_init)) P0_init
        else if (!is.null(u0_init)) noise_sd^2 * I_D
        else diffuse_scale * noise_sd^2 * I_D

  # One EKF forward pass + RTS backward pass at a fixed p̂
  erts_pass <- function(p_use, want_cov) {
    S0 <- P0 + R_obs
    K0 <- P0 %*% solve(S0)

    u_filt <- matrix(0, mp1, D)
    u_filt[1, ] <- u0 + K0 %*% (U[1, ] - u0)
    P_filt <- array(0, c(mp1, D, D))
    P_filt[1,,] <- (I_D - K0) %*% P0 %*% t(I_D - K0) + K0 %*% R_obs %*% t(K0)

    u_pred  <- matrix(0, mp1, D)
    P_pred  <- array(0, c(mp1, D, D))
    F_store <- array(0, c(mp1 - 1L, D, D))

    for (k in seq_len(mp1 - 1L)) {

      dt_k <- tt[k + 1L] - tt[k]  # dt
      uk <- u_filt[k, ] # uk current time step

      # Predict
      y1 <- rk4_step(p_use, uk, tt[k], dt_k)                  # full step = predict mean
      yh <- rk4_step(p_use, uk, tt[k], dt_k / 2)              # step-doubling: two half
      y2 <- rk4_step(p_use, yh, tt[k] + dt_k / 2, dt_k / 2)   #   steps for the error est.

      u_pred[k + 1L, ] <- y1 # Predicted step from model + control

      e_k <- (16 / 15) * (y1 - y2)                            # RK4 local truncation error
      Ju_k <- matrix(as.vector(J_u(c(p_use, uk, tt[k]))), D, D)  # J[a, b] = df_a/du_b
      Fk <- I_D + dt_k * Ju_k
      F_store[k,,] <- Fk

      Pk_pred <- Fk %*% P_filt[k,,] %*% t(Fk) + diag(e_k^2 + q_floor, D) # Predicted covariance
      P_pred[k + 1L,,] <- Pk_pred

      # Update step
      innov <- U[k + 1L, ] - u_pred[k + 1L, ]
      Sk <- Pk_pred + R_obs
      Kk <- Pk_pred %*% solve(Sk)

      u_filt[k + 1L, ] <- u_pred[k + 1L, ] + Kk %*% innov
      P_filt[k + 1L,,] <- (I_D - Kk) %*% Pk_pred %*% t(I_D - Kk) + Kk %*% R_obs %*% t(Kk)
    }

    u_smooth <- matrix(0, mp1, D)
    u_smooth[mp1, ] <- u_filt[mp1, ]
    P_smooth <- if (want_cov) array(0, c(mp1, D, D)) else NULL
    if (want_cov) P_smooth[mp1,,] <- P_filt[mp1,,]

    # RTS Smoother
    for (k in seq(mp1 - 1L, 1L)) {
      Pk <- P_filt[k,,]
      Pp <- P_pred[k + 1L,,]
      Fk <- F_store[k,,]

      Ck <- Pk %*% t(Fk) %*% solve(Pp + 1e-10 * I_D)
      u_smooth[k, ] <- u_filt[k, ] + Ck %*% (u_smooth[k + 1L, ] - u_pred[k + 1L, ])
      if (want_cov){
        P_smooth[k,,] <- Pk + Ck %*% (P_smooth[k + 1L,,] - Pp) %*% t(Ck)
      }
    }

    list(U_star = u_smooth, P_smooth = P_smooth, u_filt = u_filt, P_filt = P_filt, u_pred = u_pred, P_pred = P_pred)
  }

  # Conditional p̂ pass: smoothed mean + conditional posterior covariance
  base     <- erts_pass(p, want_cov = TRUE)
  u_smooth <- base$U_star
  P_cond   <- base$P_smooth

  # Fold parameter uncertainty in coherently via the trajectory sensitivity
  # S_k = ∂u*_k/∂p̂ (central differences over re-runs of the full smoother),
  # adding S_k Ĉ S_k^T to the conditional posterior at every time step.
  fold <- isTRUE(fold_param_uncertainty) && !is.null(param_cov) && !is.null(J_p)
  P_param  <- NULL
  P_smooth <- P_cond
  if (fold) {
    sens <- array(0, c(mp1, D, J))  # ∂u_k/∂p̂_j
    for (j in seq_len(J)) {
      hj <- 1e-5 * max(1, abs(p[j]))
      pp <- p
      pp[j] <- pp[j] + hj
      pm <- p
      pm[j] <- pm[j] - hj
      sens[, , j] <- (erts_pass(pp, FALSE)$U_star - erts_pass(pm, FALSE)$U_star) / (2 * hj)
    }
    P_param <- array(0, c(mp1, D, D))
    for (k in seq_len(mp1)) {
      Sk <- matrix(sens[k, , ], D, J)
      P_param[k,,]  <- Sk %*% param_cov %*% t(Sk)
      # Law of total variance to explain the total uncertainty in the estimate
      # total var = unexplained + explained
      P_smooth[k,,] <- P_cond[k,,] + P_param[k,,]
    }
  }

  list(
    U_star = u_smooth,
    P_smooth = P_smooth,
    P_smooth_cond  = P_cond,
    P_smooth_param = P_param,
    u_filt = base$u_filt,
    P_filt = base$P_filt,
    u_pred = base$u_pred,
    P_pred = base$P_pred
  )
}