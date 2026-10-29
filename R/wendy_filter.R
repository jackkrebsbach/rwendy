# Endpoint estimation uses boundary-layer weak equations with Euler-Maclaurin
# quadrature corrections. Shared helpers below also support parameter fitting
# and prediction.

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

# Leibniz coefficients for g = phi * f + phi' * u along the ODE trajectory.
# f_derivs contains f and its total time derivatives. The (order + 2) x D
# result satisfies g^(order) = sum_m phi^(m) * A[m + 1, ].
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

# Evaluate the same expansion from endpoint values of the test derivatives.
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

# Stack equations by state component: Phi = kronecker(I_D, B), where B contains
# phi(t_0) for the initial endpoint or -phi(t_m) for the final endpoint.
build_endpoint_matrix <- function(B, D) {
  K <- length(B)
  Phi <- matrix(0, K * D, D)
  for (d in seq_len(D)) {
    Phi[block_indices(d, K), d] <- B
  }
  Phi
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

    bl_phi_tM %*% em_coeffs(time_derivs(u_last, tt[M]), u_last, 1) + bl_phi_t1 %*% em_coeffs(time_derivs(u_first, tt[1L]), u_first, -1)
  }
}

# Parameter Jacobian of the Euler-Maclaurin correction, holding U fixed.
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

# Assemble endpoint and optional interior equations. V and Vp include the
# trapezoidal endpoint weights; phi_endpoint retains unweighted derivatives.
build_ic_bl_system <- function(tt, r_bl, n_bl, orders = 0:4, include_interior = FALSE,
                               interior_stride = 1L, endpoint = c("initial", "final"),
                               interior_rows = NULL) {
  endpoint <- match.arg(endpoint)
  M <- length(tt)
  endpoint_index <- if (endpoint == "initial") 1L else M
  # A signed B lets both endpoints use r_trap = -T[phi * f + phi' * u].
  boundary_sign <- if (endpoint == "initial") 1 else -1
  side <- if (endpoint == "initial") "left" else "right"
  blocks <- lapply(orders, function(order) {
    build_boundary_layer_block(psi, tt, r_bl, order = order, side = side, n_bl = n_bl)
  })
  n_interior <- 0L

  interior <- NULL
  if (include_interior && !is.null(interior_rows)) {
    interior <- list(interior_rows$V, interior_rows$Vp)
  } else if (include_interior && (2L * r_bl + 1L) <= (M - 2L)) {
    interior <- lapply(0:1, function(order) build_test_function_matrix(psi, tt, r_bl, order = order))
    stride <- max(1L, as.integer(interior_stride))
    if (stride > 1L) {
      keep <- seq(1L, nrow(interior[[1]]), by = stride)
      interior <- lapply(interior, function(V) V[keep, , drop = FALSE])
    }
  }

  if (!is.null(interior)) {
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

  phi_endpoint <- matrix(0, nrow = n_equations, ncol = length(orders))
  for (i in seq_along(orders)) {
    phi_endpoint[, i] <- blocks[[i]][, endpoint_index]
  }

  B <- boundary_sign * phi_endpoint[, 1]
  list(V = trapezoid(blocks[[1]]), Vp = trapezoid(blocks[[2]]),
       B = B, BtB = sum(B * B),
       n_equations = n_equations, n_interior = n_interior,
       boundary_rows = which(rowSums(abs(phi_endpoint)) > 0),
       phi_endpoint = phi_endpoint, endpoint_index = endpoint_index,
       boundary_sign = boundary_sign,
       support = which(colSums(abs(blocks[[1]]) + abs(blocks[[2]])) > 0))
}

# Linearize the trapezoidal weak residual at U. Columns of L are ordered by
# time within each state component, restricted to the test-function support.
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
  list(L = L, column_variance = column_variance, Phi = build_endpoint_matrix(bl_system$B, D))
}

# Covariance of the quadratic measurement-noise term for independent Gaussian
# errors. The Hessians are evaluated at U; this is separate from quadrature error.
build_ic_noise_quad <- function(bl_system, U, tt, p, J_uu, noise_sd, dt, hess_cache = NULL) {
  D <- ncol(U)
  K <- bl_system$n_equations
  support <- bl_system$support
  state_hessian <- if (!is.null(hess_cache)) hess_cache else build_ic_hessian_cache(U, tt, p, J_uu, D)

  V_support <- bl_system$V[, support, drop = FALSE]
  # Samples with V = 0 contribute only through the linear Vp * U term.
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

# Factor the weak-residual covariance with a relative diagonal ridge.
# Reuse this factorization for the GLS projection and uncertainty calculations.
build_ic_gls_weights <- function(sensitivity, ridge) {
  L <- sensitivity$L
  if (is.null(sensitivity$state_cov)) {
    L_cov <- sweep(L, 2L, sensitivity$column_variance, `*`)
    Omega <- L %*% (sensitivity$column_variance * t(L))
  } else {
    L_cov <- L %*% sensitivity$state_cov
    Omega <- L_cov %*% t(L)
    Omega <- (Omega + t(Omega)) / 2
  }
  if (!is.null(sensitivity$Omega2)) Omega <- Omega + sensitivity$Omega2
  Omega_chol  <- chol(Omega + ridge * mean(diag(Omega)) * diag(nrow(L)))
  Omega_solve <- function(rhs) backsolve(Omega_chol, backsolve(Omega_chol, rhs, transpose = TRUE))
  PhitW <- t(Omega_solve(sensitivity$Phi))
  PhitWPhi  <- PhitW %*% sensitivity$Phi
  list(Omega_solve = Omega_solve, PhitW = PhitW, PhitWPhi = PhitWPhi,
       cov_design = solve(PhitWPhi), Omega = Omega, L_cov = L_cov, ridge = ridge * mean(diag(Omega)))
}

# Cache H[d, a, b] = d2 f_d / (du_a du_b) at each observed sample.
build_ic_hessian_cache <- function(U, tt, p, J_uu, D) {
  M <- nrow(U)
  input <- rbind(matrix(rep(p, M), nrow = length(p)), t(U), matrix(tt, nrow = 1L))
  hessians <- array(J_uu(input), c(M, D, D, D))
  function(m) array(hessians[m, , , ], c(D, D, D))
}

# Approximate O(sigma^2) bias from vector-field curvature and data-dependent
# GLS weights. This local calculation holds p and noise_sd fixed and omits
# endpoint-correction curvature and derivatives of the covariance ridge.
build_ic_bias_o2 <- function(bl_system, sensitivity, gls, P, EM_jacobian, U, tt, p, J_uu, noise_sd, dt, hess_cache = NULL) {
  D <- length(noise_sd)
  K <- bl_system$n_equations
  support <- bl_system$support
  n_support <- length(support)
  L <- sensitivity$L
  column_variance <- sensitivity$column_variance
  noise_var <- noise_sd^2
  residual_jacobian <- sensitivity$Phi + EM_jacobian
  A <- gls$Omega_solve(L - residual_jacobian %*% (P %*% L))
  state_hessian <- if (!is.null(hess_cache)) hess_cache else build_ic_hessian_cache(U, tt, p, J_uu, D)

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

# Score candidate radii before any endpoint solve, using U and its observed
# endpoint. The MSE proxy is sum_d [Var_d + (delta_d + b_d)^2] / sigma_d^2,
# where delta is the projected EM(2)-minus-EM(4) indicator and b is the
# estimated statistical bias. Failed candidates receive nonfinite scores.
select_ic_design <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, J_uu,
                             noise_sd, param_cov, em_order,
                             r_bl_grid, r_bl_max, include_interior = TRUE,
                             hess_cache = NULL, quad_cov = TRUE, ridge,
                             interior_stride = 1L, include_bias_o2 = TRUE,
                             endpoint = c("initial", "final"), interior_rows = NULL) {
  endpoint <- match.arg(endpoint)
  radii <- sort(unique(pmin(pmax(as.integer(r_bl_grid), 2L), r_bl_max)))
  D <- ncol(U)
  M <- nrow(U)
  J <- length(p)
  noise_var <- noise_sd^2
  tt <- as.vector(tt)
  dt <- mean(diff(tt))
  c2 <- dt^2 / 12
  c4 <- if (em_order >= 4L) dt^4 / 720 else 0
  endpoint_index <- if (endpoint == "initial") 1L else M
  u_endpoint_observed <- as.numeric(U[endpoint_index, ])

  em_parts <- function(bl_system, u_endpoint, params = p) {
    u_endpoint <- as.vector(u_endpoint)
    input <- matrix(c(params, u_endpoint, tt[endpoint_index]), ncol = 1L)
    f_derivs <- list(as.vector(f_(input)),       as.vector(dF_dt_(input)),
                     as.vector(d2F_dt2_(input)), as.vector(d3F_dt3_(input)))
    A2 <- matrix(0, 5L, D)
    A2[1:3, ] <- c2 * g_coeffs(f_derivs, u_endpoint, 1L)
    A4 <- if (c4 != 0) -c4 * g_coeffs(f_derivs, u_endpoint, 3L) else matrix(0, 5L, D)
    rows <- bl_system$boundary_rows
    phi <- bl_system$boundary_sign * bl_system$phi_endpoint[rows, , drop = FALSE]
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
                                      interior_stride = interior_stride, endpoint = endpoint,
                                      interior_rows = interior_rows)
      sensitivity <- build_ic_noise_sensitivity(bl_system, U, tt, p, J_u, noise_sd, dt)
      if (isTRUE(quad_cov)) {
        sensitivity$Omega2 <- tryCatch(
          build_ic_noise_quad(bl_system, U, tt, p, J_uu, noise_sd, dt, hess_cache = hess_cache),
          error = function(err) NULL)
      }
      gls <- build_ic_gls_weights(sensitivity, ridge)

      # Linearize at the observed endpoint so selection remains a priori.
      h <- 1e-6 * max(1, sqrt(sum(u_endpoint_observed^2)))
      EM_jacobian <- matrix(0, bl_system$n_equations * D, D)
      for (d in seq_len(D)) {
        u_endpoint_up <- u_endpoint_observed
        u_endpoint_up[d] <- u_endpoint_up[d] + h
        u_endpoint_down <- u_endpoint_observed
        u_endpoint_down[d] <- u_endpoint_down[d] - h
        EM_jacobian[, d] <- as.vector((em_parts(bl_system, u_endpoint_up)$EM -
                                       em_parts(bl_system, u_endpoint_down)$EM) / (2 * h))
      }
      P <- solve(gls$PhitWPhi + gls$PhitW %*% EM_jacobian, gls$PhitW)

      G <- P %*% sensitivity$L
      cov_u_endpoint <- G %*% (sensitivity$column_variance * t(G))
      variance_linear <- diag(cov_u_endpoint)
      variance_quadratic <- variance_param <- numeric(D)
      if (!is.null(sensitivity$Omega2)) {
        covariance_quadratic <- (P %*% sensitivity$Omega2) %*% t(P)
        variance_quadratic <- diag(covariance_quadratic)
        cov_u_endpoint <- cov_u_endpoint + covariance_quadratic
      }

      if (!is.null(param_cov)) {
        S_p <- matrix(0, D, J)
        for (j in seq_len(J)) {
          h_j <- 1e-6 * max(1, abs(p[j]))
          p_up <- p
          p_up[j] <- p_up[j] + h_j
          p_down <- p
          p_down[j] <- p_down[j] - h_j
          rhs_up   <- as.vector(trapezoid_residual(bl_system, p_up) - em_parts(bl_system, u_endpoint_observed, p_up)$EM)
          rhs_down <- as.vector(trapezoid_residual(bl_system, p_down) - em_parts(bl_system, u_endpoint_observed, p_down)$EM)
          S_p[, j] <- P %*% ((rhs_up - rhs_down) / (2 * h_j))
        }
        covariance_param <- S_p %*% param_cov %*% t(S_p)
        variance_param <- diag(covariance_param)
        cov_u_endpoint <- cov_u_endpoint + covariance_param
      }

      # Score the raw estimator; the later debiasing safeguard does not alter
      # this bias penalty. include_bias_o2 is an internal score option.
      statistical_bias <- if (isTRUE(include_bias_o2)) {
        bias <- tryCatch(
          build_ic_bias_o2(bl_system, sensitivity, gls, P, EM_jacobian, U, tt, p, J_uu,
                           noise_sd, dt, hess_cache = hess_cache)$total,
          error = function(err) NULL)
        if (is.null(bias) || !all(is.finite(bias))) rep(0, D) else bias
      } else {
        rep(0, D)
      }

      # The projected fourth-order term is an error indicator, not the
      # remaining quadrature error after applying the EM(4) correction.
      truncation_bias <- if (c4 != 0) {
        as.numeric(P %*% as.vector(em_parts(bl_system, u_endpoint_observed)$Delta4))
      } else {
        rep(0, D)
      }

      variance_objective <- sum(diag(cov_u_endpoint) / noise_var)
      if (!is.finite(variance_objective)) stop("non-finite objective", call. = FALSE)
      list(obj = variance_objective + sum((truncation_bias + statistical_bias)^2 / noise_var),
           var_obj = variance_objective,
           var_linear_obj = sum(variance_linear / noise_var),
           var_quadratic_obj = sum(variance_quadratic / noise_var),
           var_param_obj = sum(variance_param / noise_var),
           quad_obj = sum(truncation_bias^2 / noise_var),
           stat_bias_obj = sum(statistical_bias^2 / noise_var),
           bias_cross_obj = 2 * sum(truncation_bias * statistical_bias / noise_var),
           quadrature_bias = truncation_bias, statistical_bias = statistical_bias,
           variance = diag(cov_u_endpoint),
           variance_linear = variance_linear, variance_quadratic = variance_quadratic,
           variance_param = variance_param,
           ic_system = list(bl_system = bl_system, sensitivity = sensitivity, gls = gls))
    }, error = function(err) list(obj = NA_real_, var_obj = NA_real_,
                                 var_linear_obj = NA_real_, var_quadratic_obj = NA_real_,
                                 var_param_obj = NA_real_,
                                 quad_obj = NA_real_, stat_bias_obj = NA_real_,
                                 bias_cross_obj = NA_real_,
                                 quadrature_bias = rep(NA_real_, D),
                                 statistical_bias = rep(NA_real_, D),
                                 variance = rep(NA_real_, D),
                                 variance_linear = rep(NA_real_, D),
                                 variance_quadratic = rep(NA_real_, D),
                                 variance_param = rep(NA_real_, D), ic_system = NULL))

    rows[[i]] <- data.frame(r_bl = r_bl, n_bl = n_bl, obj = candidate$obj,
                           var_obj = candidate$var_obj, quad_obj = candidate$quad_obj,
                           var_linear_obj = candidate$var_linear_obj,
                           var_quadratic_obj = candidate$var_quadratic_obj,
                           var_param_obj = candidate$var_param_obj,
                           stat_bias_obj = candidate$stat_bias_obj,
                           bias_cross_obj = candidate$bias_cross_obj)
    for (d in seq_len(D)) {
      rows[[i]][[paste0("quad_bias_", d)]] <- candidate$quadrature_bias[d]
      rows[[i]][[paste0("stat_bias_", d)]] <- candidate$statistical_bias[d]
      rows[[i]][[paste0("variance_", d)]] <- candidate$variance[d]
      rows[[i]][[paste0("variance_linear_", d)]] <- candidate$variance_linear[d]
      rows[[i]][[paste0("variance_quadratic_", d)]] <- candidate$variance_quadratic[d]
      rows[[i]][[paste0("variance_param_", d)]] <- candidate$variance_param[d]
    }

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

#' Estimate an observed endpoint via boundary-layer fixed-point iteration
#'
#' Estimates the state at the final (default) or initial observed time, with
#' the model parameters held fixed. Uses right or left boundary-layer (BL)
#' test functions, respectively. The observed trajectory in the trapezoidal
#' sum stays fixed throughout the iteration.
#'
#' @details
#' For the final endpoint \eqn{t_m}, choose test functions whose values and
#' derivatives vanish at the initial endpoint. With
#' \eqn{g_k = \phi_k f + \phi_k' u}, integration by parts gives
#' \eqn{\int g_k\,dt = \phi_k(t_m) u_m}. Euler-Maclaurin then gives the
#' fixed-point equations
#' \deqn{\phi_k(t_m) u_m^{(n+1)} = T[g_k(U)] -
#'   \frac{\Delta t^2}{12} g_k'(t_m; u_m^{(n)}) +
#'   \frac{\Delta t^4}{720} g_k'''(t_m; u_m^{(n)}),}
#' where \eqn{T} is the composite trapezoidal sum and the fourth-order term
#' is omitted for \code{em_order = 2}. The initial-endpoint formula uses
#' \eqn{-T[g_k(U)]} with the same correction evaluated at the initial time.
#'
#' Internally both endpoints use \eqn{r_{trap} = -T[g(U)]} and the signed
#' design \eqn{\Phi_e = I_D \otimes (s_e \phi(t_e))}, with \eqn{s_e = 1}
#' initially and \eqn{s_e = -1} finally. The signed correction is
#' \eqn{\Delta_e = s_e (\Delta t^2 g'(t_e)/12 - \Delta t^4 g'''(t_e)/720)}.
#' Thus \eqn{\Phi_e u_e^{(n+1)} = \mathrm{vec}(r_{trap} - \Delta_e(u_e^{(n)}))}.
#'
#' The equations share noisy samples. Their covariance is approximated by
#' \eqn{\Omega = L \Sigma L^T + \Omega_2}, where \eqn{L} is the weak-residual
#' Jacobian evaluated at \code{U}. Measurement errors are assumed independent
#' across times and state components, so \eqn{\Sigma} is diagonal, with the
#' squared component noise scales repeated over time. The optional
#' \eqn{O(\sigma^4)} term \eqn{\Omega_2} accounts for quadratic measurement
#' noise under the Gaussian assumption. With \code{inverse = "gls"}, the
#' equations are combined by generalized least squares (GLS),
#' \deqn{u_e = (\Phi_e^T \Omega^{-1} \Phi_e)^{-1} \Phi_e^T \Omega^{-1} \mathrm{vec}(r),}
#' and with \code{inverse = "ols"} by unweighted least squares.
#'
#' The reported covariance holds the selected design and GLS weights fixed
#' and linearizes the fixed-point equation at the returned iterate, including
#' the endpoint derivative \eqn{J_{EM} = \partial\,\mathrm{vec}(\Delta_e)/\partial u_e}.
#' With \eqn{W = \Omega^{-1}} (or \eqn{W = I} for OLS), the response map is
#' \eqn{P = [\Phi_e^T W (\Phi_e + J_{EM})]^{-1}\Phi_e^T W}, giving the
#' noise contribution \eqn{P\Omega P^T}. The GLS weights use the stated ridge.
#' Parameter uncertainty, when supplied, is added through a separate
#' sensitivity term; its cross-covariance with measurement noise is omitted.
#' The covariance does not propagate uncertainty in the noise estimate,
#' radius selection, or the statistical bias correction.
#'
#' When \code{n_bl} is \code{NULL} on the GLS path, the BL radius is chosen
#' from \code{r_bl_grid} before any endpoint solve, by minimizing the MSE proxy
#' \deqn{\sum_d \left(\mathrm{Var}_d + (\delta_d + b_d)^2\right) / \sigma_d^2,}
#' where \eqn{\delta} is the projected EM(2)-minus-EM(4) quadrature indicator
#' and \eqn{b} estimates the \eqn{O(\sigma^2)} statistical bias. The indicator
#' is zero for \code{em_order = 2}; for order 4, it measures the change in
#' correction, not the remaining quadrature error. All terms use the observed
#' trajectory and endpoint. This score targets the raw estimate and includes
#' the statistical bias penalty regardless of \code{debias}.
#'
#' With \code{debias = TRUE}, the GLS bias estimate includes vector-field
#' curvature and feedback from data-dependent weights. It is subtracted only
#' when the iteration has not diverged, the required covariance calculations
#' are available, and every component has magnitude at most twice its
#' measurement-noise standard error. This safeguard applies to the whole
#' correction vector and uses the noise contribution to covariance, excluding
#' parameter uncertainty. It does not estimate the correction's uncertainty.
#'
#' The returned iterate minimizes the norm of the projected weak residual
#' among the finite iterates visited, before any bias correction. Convergence
#' is determined by the relative step tolerance. Direct calls return these
#' diagnostics even when the iteration limit is reached. For initial-condition
#' estimation through \code{\link{solveWendy}}, a nonconverged or diverged weak
#' estimate is replaced by the first observation before any optional hybrid fit.
#'
#' @param U Numeric matrix (M x D) of observed states.
#' @param f_,dF_dt_,d2F_dt2_,d3F_dt3_ Callable RHS and total time-derivative
#'   evaluators built from the symbolic engine.
#' @param tt Numeric vector (length M) of equally spaced, increasing times.
#' @param p Numeric parameter vector \eqn{\hat\theta} (held fixed).
#' @param J_u Callable state Jacobian, with entry \eqn{[a, b] = \partial f_a /
#'   \partial u_b}.
#' @param J_uu Callable state Hessian, with entry \eqn{[a, b, c] = \partial^2 f_a /
#'   \partial u_b \partial u_c}.
#' @param sigma Noise standard deviation, scalar or length D. If invalid, the
#'   GLS path and noise-propagated covariance are unavailable.
#' @param param_cov Optional J x J covariance of \code{p}, where J is the
#'   number of parameters. Adds \eqn{S_p \hat C S_p^T} to \code{cov_u} when
#'   the parameter sensitivity can be computed.
#' @param n_bl Number of BL test functions. If \code{NULL}, the GLS path
#'   selects the radius from \code{r_bl_grid} with \code{n_bl = r_bl}, and the
#'   OLS path uses \code{max(3, ceiling(r_bl / 8))}.
#' @param r_bl BL radius in grid intervals, used when automatic selection
#'   does not supply a design. Defaults to \code{min(16, floor((M - 1) / 2))}.
#'   Supply \code{n_bl} as well to fix the design on the GLS path.
#' @param max_iter,tol Fixed-point iteration limit and relative step
#'   tolerance; converged when the step is below
#'   \code{tol * max(1, ||u||)}.
#' @param em_order Euler-Maclaurin correction order, 2 or 4 (default).
#' @param inverse \code{"gls"} (default) or \code{"ols"} combination of the BL
#'   equations. GLS falls back to OLS when its weights cannot be built.
#' @param include_interior Add interior weak equations to the GLS system.
#'   Their endpoint values are zero, but their noise correlations with the
#'   boundary equations can affect the estimate. Ignored for OLS.
#' @param interior_stride Keep every s-th automatically generated interior
#'   test function.
#' @param interior_rows Optional \code{list(V, Vp)} of interior test functions
#'   (K x M) and their derivatives, replacing the automatically generated
#'   interior equations when \code{include_interior = TRUE}. Rows must vanish
#'   at both endpoints. Supplied rows bypass \code{interior_stride};
#'   boundary-layer test functions are unchanged.
#' @param quad_cov Include the quadratic-noise term \eqn{\Omega_2} in the
#'   weights and \code{cov_u}. This is measurement-noise variance, not
#'   numerical-quadrature error. GLS only.
#' @param ridge Relative diagonal regularization for the GLS weights. Adds
#'   \code{ridge * mean(diag(Omega))} to each diagonal entry before
#'   factorization. GLS only.
#' @param r_bl_grid Candidate radii for the selection. Defaults to
#'   \code{c(4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 80, 100)}, capped at
#'   \code{floor((M - 1) / 2)}.
#' @param debias Apply the estimated \eqn{O(\sigma^2)} statistical bias
#'   correction when its safeguard passes (see Details). GLS only.
#' @param hess_cache Optional state-Hessian cache from
#'   \code{build_ic_hessian_cache}; built internally if \code{NULL}.
#' @param endpoint \code{"final"} (default) for \code{U[M, ]}, or
#'   \code{"initial"} for \code{U[1, ]}.
#' @return A list with:
#' \describe{
#'   \item{u_hat, U_hat}{The estimated endpoint state, and \code{U} with only
#'     the selected endpoint row replaced by it.}
#'   \item{endpoint, endpoint_index, endpoint_time}{The selected endpoint,
#'     its row index in \code{U}, and its time in \code{tt}.}
#'   \item{cov_u}{Approximate covariance before the statistical bias correction:
#'     noise propagation, with an LS residual estimate as a fallback when
#'     propagation is unavailable, plus the parameter contribution when
#'     available. \code{NULL} if the iteration diverged. See Details for the
#'     uncertainty sources included.}
#'   \item{cov_u_resid, cov_u_param}{The LS residual covariance and the
#'     parameter-channel covariance on their own.}
#'   \item{cov_method}{\code{"noise_propagation"} or \code{"ls_residual"},
#'     with \code{"+param"} when the parameter channel is included, or
#'     \code{"diverged"}.}
#'   \item{inverse}{The combination actually used; \code{"gls"} falls back to
#'     \code{"ols"} when the weights cannot be built.}
#'   \item{design}{The radius-selection table (\code{r_bl}, \code{n_bl},
#'     \code{obj}, \code{var_obj}, \code{quad_obj}, \code{stat_bias_obj},
#'     and \code{bias_cross_obj}), or \code{NULL} if selection was unavailable
#'     or skipped. The four objective terms sum to \code{obj}; component
#'     diagnostics \code{quad_bias_d}, \code{stat_bias_d}, and \code{variance_d}
#'     are in physical units (squared units for variance).
#'     \code{var_linear_obj}, \code{var_quadratic_obj}, and \code{var_param_obj}
#'     sum to \code{var_obj}, with component values in
#'     \code{variance_linear_d}, \code{variance_quadratic_d}, and
#'     \code{variance_param_d}. These are the plug-in covariance contributions;
#'     the linear term uses the Jacobian at noisy observations and need not
#'     equal the first-order covariance at the true trajectory.}
#'   \item{bias_o2, debias_applied}{The estimated \eqn{O(\sigma^2)} bias and
#'     whether it was subtracted from \code{u_hat}. \code{bias_o2} is
#'     \code{NULL} when the correction was disabled or could not be computed.}
#'   \item{fallback}{\code{TRUE} when the design was degenerate and
#'     \code{u_hat} is the observed endpoint.}
#'   \item{iters, converged, diverged, u_history}{Iteration count, relative-step
#'     convergence flag, nonfinite-update flag, and the iterates before any
#'     statistical bias correction.}
#'   \item{r_bl, n_bl, K_bl, K_int, em_order}{The design used: BL radius, BL
#'     test-function count, total and interior equation counts, and EM order.}
#' }
#' @seealso \code{\link{estimate_IC}} for the initial-condition interface;
#'   \code{\link{predict_next}} for an unobserved next point.
#' @export
estimate_endpoint <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, J_uu, sigma,
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
                              hess_cache       = NULL,
                              endpoint         = c("final", "initial"),
                              interior_rows    = NULL) {

  endpoint <- match.arg(endpoint)

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
  endpoint_index <- if (endpoint == "initial") 1L else M

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
    hess_cache <- build_ic_hessian_cache(U, tt, p, J_uu, D)
  }

  design_table    <- NULL
  selected_system <- NULL
  # Selection uses the observed endpoint; the fixed-point solve starts only
  # after a radius and its GLS system have been chosen.
  if (inverse == "gls" && is.null(n_bl)) {
    if (is.null(r_bl_grid)) {
      r_bl_grid <- c(4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 80, 100)
    }
    selection <- tryCatch(
      select_ic_design(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, J_uu,
                       noise_sd, param_cov, em_order, r_bl_grid, r_bl_max,
                       include_interior = include_interior,
                       hess_cache = hess_cache, quad_cov = quad_cov, ridge = ridge,
                       interior_stride = interior_stride, endpoint = endpoint,
                       interior_rows = interior_rows),
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
                                    interior_stride = interior_stride, endpoint = endpoint,
                                    interior_rows = interior_rows)
    sensitivity <- if (use_noise_propagation) {
      tryCatch(build_ic_noise_sensitivity(bl_system, U, tt, p, J_u, noise_sd, dt),
               error = function(err) NULL)
    } else {
      NULL
    }
    if (isTRUE(quad_cov) && inverse == "gls" && !is.null(sensitivity)) {
      sensitivity$Omega2 <- tryCatch(
        build_ic_noise_quad(bl_system, U, tt, p, J_uu, noise_sd, dt, hess_cache = hess_cache),
        error = function(err) NULL)
    }
    gls <- if (inverse == "gls" && !is.null(sensitivity)) {
      tryCatch(build_ic_gls_weights(sensitivity, ridge), error = function(err) NULL)
    } else {
      NULL
    }
    list(bl_system = bl_system, sensitivity = sensitivity, gls = gls)
  }

  # Reuse the selected system so scoring and fitting use the same weights.
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
  phi_endpoint <- bl_system$boundary_sign * bl_system$phi_endpoint
  B <- bl_system$B
  BtB <- bl_system$BtB
  K_bl <- bl_system$n_equations
  K_int <- bl_system$n_interior
  boundary_rows <- bl_system$boundary_rows

  if (!is.finite(BtB) || BtB < .Machine$double.eps) {
    u_observed <- as.numeric(U[endpoint_index, ])
    U_hat <- U
    U_hat[endpoint_index, ] <- u_observed
    return(list(U_hat = U_hat, u_hat = u_observed, cov_u = NULL,
                endpoint = endpoint, endpoint_index = endpoint_index,
                endpoint_time = tt[endpoint_index],
                iters = 0L, converged = FALSE, diverged = FALSE,
                u_history = matrix(u_observed, nrow = 1),
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

  em_correction <- function(u, params = p, include_em4 = TRUE) {
    u <- as.vector(u)
    input <- matrix(c(params, u, tt[endpoint_index]), ncol = 1L)
    f_derivs <- list(
      as.vector(f_(input)),
      as.vector(dF_dt_(input)),
      as.vector(d2F_dt2_(input)),
      as.vector(d3F_dt3_(input))
    )
    A <- matrix(0, 5L, D)
    A[1:3, ] <- c2 * g_coeffs(f_derivs, u, 1L)
    if (include_em4 && c4 != 0) A <- A - c4 * g_coeffs(f_derivs, u, 3L)
    EM <- matrix(0, nrow = K_bl, ncol = D)
    EM[boundary_rows, ] <- phi_endpoint[boundary_rows, , drop = FALSE] %*% A
    EM
  }

  project <- if (!is.null(gls)) {
    function(rhs) as.numeric(solve(gls$PhitWPhi, gls$PhitW %*% as.vector(rhs)))
  } else {
    function(rhs) as.numeric(crossprod(B, rhs) / BtB)
  }

  residual_norm <- function(u, include_em4 = TRUE) {
    residual <- outer(B, as.numeric(u)) - (r_trap - em_correction(u, include_em4 = include_em4))
    weighted <- if (!is.null(gls)) as.numeric(gls$PhitW %*% as.vector(residual))
                else               as.numeric(crossprod(B, residual))
    sqrt(sum(weighted * weighted))
  }

  run_fixed_point <- function(include_em4 = TRUE) {
    u             <- project(r_trap)
    history       <- list(u)
    best_u        <- u
    best_residual <- if (all(is.finite(u))) residual_norm(u, include_em4) else Inf
    iters         <- 0L
    converged     <- FALSE
    diverged      <- FALSE

    for (iter in seq_len(max_iter)) {
      iters  <- iter
      u_next <- project(r_trap - em_correction(u, include_em4 = include_em4))
      history[[iter + 1L]] <- u_next

      if (!all(is.finite(u_next))) {
        diverged <- TRUE
        break
      }

      # Retain the best projected residual, including when iteration stops
      # before convergence. The bias correction is applied after this solve.
      residual <- residual_norm(u_next, include_em4)
      if (is.finite(residual) && residual < best_residual) {
        best_residual <- residual
        best_u       <- u_next
      }

      step <- sqrt(sum((u_next - u)^2))
      u    <- u_next
      if (is.finite(step) && step < tol * max(1, sqrt(sum(u_next^2)))) {
        converged <- TRUE
        break
      }
    }

    list(u = best_u, iters = iters, converged = converged,
         diverged = diverged, history = history)
  }

  fit <- run_fixed_point()
  u   <- fit$u

  EM_jacobian <- tryCatch({
    jacobian <- matrix(0, K_bl * D, D)
    h <- 1e-6 * max(1, sqrt(sum(u^2)))
    for (d in seq_len(D)) {
      u_up <- u
      u_up[d] <- u_up[d] + h
      u_down <- u
      u_down[d] <- u_down[d] - h
      jacobian[, d] <- as.vector((em_correction(u_up) - em_correction(u_down)) / (2 * h))
    }
    jacobian
  }, error = function(err) NULL)

  P <- if (!is.null(EM_jacobian)) {
    tryCatch({
      if (!is.null(gls)) {
        solve(gls$PhitWPhi + gls$PhitW %*% EM_jacobian, gls$PhitW)
      } else {
        Phi <- if (!is.null(sensitivity)) sensitivity$Phi else build_endpoint_matrix(B, D)
        solve(BtB * diag(D) + crossprod(Phi, EM_jacobian), t(Phi))
      }
    }, error = function(err) NULL)
  } else {
    NULL
  }

  # Propagate noise at the chosen design with GLS weights held fixed.
  # This covariance is not adjusted for the later statistical bias correction.
  cov_u_noise <- if (!is.null(sensitivity) && !is.null(P)) {
    tryCatch({
      G <- P %*% sensitivity$L
      covariance <- G %*% (sensitivity$column_variance * t(G))
      if (!is.null(sensitivity$Omega2)) covariance <- covariance + (P %*% sensitivity$Omega2) %*% t(P)
      covariance
    }, error = function(err) NULL)
  } else {
    NULL
  }

  cov_u_param <- if (!is.null(param_cov) && !is.null(P)) {
    tryCatch({
      rhs <- function(params) as.vector(trapezoid_residual(params) - em_correction(u, params))
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

  cov_u_resid <- tryCatch({
    residual <- outer(B, u) - (r_trap - em_correction(u))
    # Use boundary rows for the residual-based endpoint variance: interior
    # rows have B = 0 and do not enter the unweighted endpoint projection.
    residual <- residual[boundary_rows, , drop = FALSE]
    dof <- max(length(boundary_rows) - 1L, 1L)
    diag(colSums(residual * residual) / dof / BtB, nrow = D, ncol = D)
  }, error = function(err) NULL)

  cov_u      <- if (!is.null(cov_u_noise)) cov_u_noise else cov_u_resid
  cov_method <- if (!is.null(cov_u_noise)) "noise_propagation" else "ls_residual"
  if (!is.null(cov_u) && !is.null(cov_u_param)) {
    cov_u      <- cov_u + cov_u_param
    cov_method <- paste0(cov_method, "+param")
  }

  if (fit$diverged) {
    cov_u       <- NULL
    cov_u_param <- NULL
    cov_method  <- "diverged"
  }

  bias_o2 <- NULL
  debias_applied <- FALSE
  if (isTRUE(debias) && !fit$diverged && !is.null(gls) && !is.null(sensitivity) &&
      !is.null(P) && !is.null(EM_jacobian) && !is.null(cov_u_noise)) {
    bias_o2 <- tryCatch(
      build_ic_bias_o2(bl_system, sensitivity, gls, P, EM_jacobian, U, tt, p, J_uu,
                       noise_sd, dt, hess_cache = hess_cache)$total,
      error = function(err) NULL)
    if (!is.null(bias_o2) && all(is.finite(bias_o2))) {
      # Apply the whole correction only if every component is within two
      # noise-only standard errors; parameter covariance is excluded here.
      se <- sqrt(pmax(diag(cov_u_noise), 0))
      if (all(abs(bias_o2) <= 2 * se)) {
        u <- u - bias_o2
        debias_applied <- TRUE
      }
    }
  }

  U_hat <- U
  U_hat[endpoint_index, ] <- u

  list(
    endpoint       = endpoint,
    endpoint_index = endpoint_index,
    endpoint_time  = tt[endpoint_index],
    U_hat          = U_hat,
    u_hat          = u,
    cov_u          = cov_u,
    cov_u_resid    = cov_u_resid,
    cov_u_param    = cov_u_param,
    cov_method     = cov_method,
    inverse        = inverse,
    design         = design_table,
    bias_o2        = bias_o2,
    debias_applied = debias_applied,
    fallback       = FALSE,
    iters          = fit$iters,
    converged      = fit$converged,
    diverged       = fit$diverged,
    u_history      = do.call(rbind, fit$history),
    r_bl           = r_bl,
    n_bl           = n_bl,
    K_bl           = K_bl,
    K_int          = K_int,
    em_order       = em_order
  )
}

#' Estimate the initial condition via boundary-layer fixed-point iteration
#'
#' Calls \code{\link{estimate_endpoint}} with \code{endpoint = "initial"},
#' using left boundary-layer weak equations to estimate the first state.
#' Parameters and the observed trajectory are held fixed during iteration.
#'
#' The default GLS path selects a radius before solving when \code{n_bl}
#' is \code{NULL}. Supply both \code{r_bl} and \code{n_bl} to fix the design.
#' See \code{\link{estimate_endpoint}} for the radius score, convergence
#' diagnostics, bias-correction safeguard, and covariance approximations.
#'
#' @inheritParams estimate_endpoint
#' @param param_cov Optional J x J covariance of \code{p}, where J is the
#'   number of parameters. Adds \eqn{S_p \hat C S_p^T} to \code{cov_u0} when
#'   the parameter sensitivity can be computed.
#' @param quad_cov Include the quadratic-noise term \eqn{\Omega_2} in the
#'   weights and \code{cov_u0}. This is measurement-noise variance, not
#'   numerical-quadrature error. GLS only.
#' @param ridge Relative diagonal regularization for the GLS weights. Adds
#'   \code{ridge * mean(diag(Omega))} to each diagonal entry before
#'   factorization. GLS only.
#' @return The result of \code{estimate_endpoint}, with \code{u_hat},
#'   \code{cov_u}, \code{cov_u_resid}, \code{cov_u_param}, and
#'   \code{u_history} named \code{u0hat}, \code{cov_u0},
#'   \code{cov_u0_resid}, \code{cov_u0_param}, and \code{u0_history},
#'   respectively. The endpoint metadata is omitted. \code{U_hat} replaces
#'   only the first row of \code{U}. The covariance describes the estimate
#'   before the statistical bias correction and does not include uncertainty
#'   in that correction.
#' @seealso \code{\link{estimate_endpoint}} for the formulation and for
#'   estimating the final observed state.
#' @export
estimate_IC <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, J_uu, sigma,
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
                        hess_cache       = NULL,
                        interior_rows    = NULL) {
  fit <- estimate_endpoint(
    U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p, J_u, J_uu, sigma,
    param_cov = param_cov, n_bl = n_bl, r_bl = r_bl,
    max_iter = max_iter, tol = tol, em_order = em_order, inverse = inverse,
    include_interior = include_interior, interior_stride = interior_stride,
    quad_cov = quad_cov, ridge = ridge, r_bl_grid = r_bl_grid,
    debias = debias, hess_cache = hess_cache, endpoint = "initial",
    interior_rows = interior_rows)
  aliases <- c(u_hat = "u0hat", cov_u = "cov_u0", cov_u_resid = "cov_u0_resid",
               cov_u_param = "cov_u0_param", u_history = "u0_history")
  rename <- names(fit) %in% names(aliases)
  names(fit)[rename] <- unname(aliases[names(fit)[rename]])
  fit[c("endpoint", "endpoint_index", "endpoint_time")] <- NULL
  fit
}

prediction_endpoint_jacobian <- function(endpoint_rhs, u) {
  D <- length(u)
  columns <- lapply(seq_len(D), function(a) {
    step <- .Machine$double.eps^(1 / 3) * max(1, abs(u[a]))
    up <- down <- u
    up[a] <- up[a] + step
    down[a] <- down[a] - step
    as.vector(endpoint_rhs(up) - endpoint_rhs(down)) / (2 * step)
  })
  do.call(cbind, columns)
}

certify_prediction_root <- function(u, update, endpoint_rhs, projection, tol) {
  residual <- u - update(u)
  defect <- sqrt(sum(residual^2))
  threshold <- tol * max(1, sqrt(sum(u^2)))
  out <- list(converged = FALSE, defect = defect, correction = Inf, rcond = 0)
  if (!is.finite(defect) || defect > threshold) return(out)
  jacobian <- diag(length(u)) - projection %*% prediction_endpoint_jacobian(endpoint_rhs, u)
  if (any(!is.finite(jacobian))) return(out)
  out$rcond <- rcond(jacobian)
  if (!is.finite(out$rcond) || out$rcond <= sqrt(.Machine$double.eps)) return(out)
  correction <- tryCatch(solve(jacobian, residual), error = function(e) Inf)
  out$correction <- sqrt(sum(correction^2))
  out$converged <- is.finite(out$correction) && out$correction <= threshold
  out
}

build_prediction_covariance <- function(U, tt, p, J_uu, noise_sd, dt,
                                         observed_system, sensitivity, gls,
                                         r_observed, endpoint_rhs, u_hat, ridge,
                                         state_cov = NULL, mse_proxy = FALSE,
                                         hess_cache = NULL) {
  M <- nrow(U)
  D <- ncol(U)
  K <- observed_system$n_equations
  support <- observed_system$support
  n_support <- length(support)
  Phi <- sensitivity$Phi
  L <- sensitivity$L
  PhitW <- if (is.null(gls)) t(Phi) else gls$PhitW
  L_cov <- if (!is.null(gls)) gls$L_cov else if (!is.null(sensitivity$state_cov)) {
    L %*% sensitivity$state_cov
  } else {
    sweep(L, 2L, sensitivity$column_variance, `*`)
  }

  endpoint_jacobian <- prediction_endpoint_jacobian(endpoint_rhs, u_hat)
  H <- PhitW %*% (Phi - endpoint_jacobian)
  rhs <- PhitW %*% L
  P <- solve(H, PhitW)
  G_fixed <- P %*% L
  cov_fixed <- P %*% L_cov %*% t(G_fixed)
  noise_curvature <- feedback <- numeric(K * D)
  state_hessian <- if (is.null(hess_cache)) build_ic_hessian_cache(U, tt, p, J_uu, D) else hess_cache
  if (mse_proxy && !is.null(gls)) {
    # Columns of this matrix are W (L - (Phi-b_z) P L) Sigma.
    residual_cov <- gls$Omega_solve(L_cov -
      (Phi - endpoint_jacobian) %*% (P %*% L_cov))
  }

  if (!is.null(gls)) {
    # F(z, X) = Phi' W(X) e(z, X) = 0 gives
    # H dz/dX_j = Phi' W [L_j + dOmega_tilde/dX_j W e].
    # Although W is fixed during iteration, it is data dependent across
    # noise realizations. Its derivative includes the relative ridge.
    residual <- as.vector(Phi %*% u_hat) - as.vector(r_observed + endpoint_rhs(u_hat))
    weighted_residual <- gls$Omega_solve(residual)
    V <- observed_system$V[, support, drop = FALSE]
    curved <- which(colSums(abs(V)) > 0)
    for (i in curved) {
      hessian <- state_hessian(support[i])
      columns <- i + (seq_len(D) - 1L) * n_support
      weighted_L_i <- L_cov[, columns, drop = FALSE]
      score <- as.vector(crossprod(weighted_L_i, weighted_residual))
      if (mse_proxy) {
        sample_cov <- if (is.null(sensitivity$state_cov)) {
          diag(sensitivity$column_variance[columns], D)
        } else sensitivity$state_cov[columns, columns, drop = FALSE]
        for (d in seq_len(D)) {
          rows <- block_indices(d, K)
          noise_curvature[rows] <- noise_curvature[rows] -
            0.5 * dt * V[, i] * sum(hessian[d, , ] * sample_cov)
        }
      }
      for (b in seq_len(D)) {
        # Only D columns of L depend on the state at this sample.
        dL <- matrix(0, K * D, D)
        for (d in seq_len(D)) {
          dL[block_indices(d, K), ] <- -dt * outer(V[, i], hessian[d, , b])
        }
        omega_score <- dL %*% score +
          weighted_L_i %*% crossprod(dL, weighted_residual)
        ridge_derivative <- 2 * ridge * sum(dL * weighted_L_i) / nrow(L)
        column <- columns[b]
        rhs[, column] <- rhs[, column] +
          as.vector(PhitW %*% (omega_score + ridge_derivative * weighted_residual))
        if (mse_proxy) {
          a <- residual_cov[, column]
          feedback <- feedback + as.vector(
            dL %*% crossprod(weighted_L_i, a) +
              weighted_L_i %*% crossprod(dL, a) + ridge_derivative * a)
        }
      }
    }
  }

  jacobian <- matrix(0, D, M * D)
  columns <- as.vector(outer(support, (seq_len(D) - 1L) * M, `+`))
  jacobian[, columns] <- solve(H, rhs)
  # R column-major layout: time varies fastest, followed by state component.
  cross_covariance <- if (is.null(state_cov)) {
    rep(noise_sd^2, each = M) * t(jacobian)
  } else state_cov %*% t(jacobian)
  covariance <- jacobian %*% cross_covariance
  covariance <- (covariance + t(covariance)) / 2
  # Interleave the appended time within each state block, as in vec(U_hat).
  old <- as.vector(outer(seq_len(M), (seq_len(D) - 1L) * (M + 1L), `+`))
  next_rows <- seq_len(D) * (M + 1L)
  cov_state <- matrix(0, (M + 1L) * D, (M + 1L) * D)
  cov_state[old, old] <- if (is.null(state_cov)) diag(rep(noise_sd^2, each = M)) else state_cov
  cov_state[old, next_rows] <- cross_covariance
  cov_state[next_rows, old] <- t(cross_covariance)
  cov_state[next_rows, next_rows] <- covariance
  out <- list(prediction_jacobian = jacobian, cov_past_next = cross_covariance,
              cov_u = covariance, cov_state = cov_state, mse_proxy = NULL)
  if (mse_proxy) {
    # Curvature of the missing endpoint contribution also feeds back through
    # the fixed point. Contract its Hessian with the frozen-weight variance.
    eig <- eigen((cov_fixed + t(cov_fixed)) / 2, symmetric = TRUE)
    endpoint_curvature <- numeric(K * D)
    for (a in which(eig$values > 0)) {
      direction <- eig$vectors[, a] * sqrt(eig$values[a])
      step <- 1e-3 * max(1, sqrt(sum(u_hat^2))) / sqrt(sum(direction^2))
      endpoint_curvature <- endpoint_curvature + as.vector(
        endpoint_rhs(u_hat + step * direction) - 2 * endpoint_rhs(u_hat) +
          endpoint_rhs(u_hat - step * direction)) / (2 * step^2)
    }
    statistical_bias <- as.vector(P %*% (noise_curvature - feedback + endpoint_curvature))
    truncation_bias <- as.vector(P %*% as.vector(
      endpoint_rhs(u_hat, 2L) - endpoint_rhs(u_hat, 4L)))
    var_obj <- sum(diag(covariance) / noise_sd^2)
    bias_obj <- sum((truncation_bias + statistical_bias)^2 / noise_sd^2)
    out$mse_proxy <- list(obj = var_obj + bias_obj, var_obj = var_obj, bias_obj = bias_obj,
                          truncation_bias = truncation_bias, statistical_bias = statistical_bias)
  }
  out
}

#' Predict the next trajectory point via boundary-layer fixed-point iteration
#'
#' Extends the time grid by one step and solves the right boundary-layer
#' equations for the unobserved next state, with supplied parameters held
#' fixed. The input can contain previous predictions with their full joint
#' covariance. The initial guess is either the last supplied state or its forward
#' Euler step. No parameter estimation is performed.
#'
#' @details
#' Write \eqn{h = \Delta t}, \eqn{t_+ = t_m + h}, and
#' \eqn{g_k = \phi_k f + \phi_k' u}. On the extended grid, separate the
#' trapezoidal sum into its observed contribution \eqn{Q_{k,obs}} and the
#' missing endpoint contribution. The fixed-point equations are
#' \deqn{\phi_k(t_+) z^{(n+1)} = Q_{k,obs}
#'   + \frac{h}{2}\left[\phi_k(t_+) f(p,z^{(n)},t_+)
#'                       + \phi_k'(t_+) z^{(n)}\right]
#'   - \frac{h^2}{12} g_k'(t_+;z^{(n)})
#'   + \frac{h^4}{720} g_k'''(t_+;z^{(n)}).}
#' The last observed sample is now an interior quadrature point with full
#' weight. Both the missing trapezoidal term and the Euler-Maclaurin
#' correction are recomputed at each iterate. The observed samples stay fixed.
#' The fourth-order term is omitted for \code{em_order = 2}.
#'
#' By default, GLS combines the right boundary-layer equations with interior
#' equations spanning the observed history. Interior equations have zero
#' endpoint design values, but contribute through their covariance with the
#' boundary equations. GLS propagates the full covariance of the supplied
#' trajectory through the observed weak residual, \eqn{\Omega=L\Sigma L^T},
#' retaining correlations across times and components. Weights are fixed
#' across iterations. The unobserved state is not treated as
#' an independent noisy observation. Both initializations therefore solve
#' the same equations with the same weights for any fixed candidate design.
#' This prototype does not debias the prediction. Optional OLS uses only
#' the boundary-layer equations, since interior equations do not affect its
#' unweighted endpoint projection.
#'
#' When \code{n_bl = NULL}, GLS selects a radius from \code{r_bl_grid},
#' using \code{n_bl = r_bl} as in \code{estimate_endpoint}. Each candidate
#' is solved because the future endpoint is unobserved. Among converged
#' candidates, selection minimizes the MSE proxy
#' \deqn{\sum_d [\mathrm{Var}_d+(\delta_d+b_d)^2]/s_d^2.}
#' Here \eqn{\mathrm{Var}} comes from the full supplied covariance;
#' \eqn{\delta} is the linearized EM(2)-minus-EM(4) truncation indicator;
#' and \eqn{b} estimates the leading noise bias from weak-residual curvature,
#' GLS weight feedback, and missing-endpoint curvature. This is a proxy,
#' not an exact MSE or a higher-order covariance calculation. No future
#' observations or true future states enter the score. The scales \eqn{s_d}
#' are \code{sigma}, or RMS marginal SDs inferred from \code{state_cov}.
#' Candidate radii are capped to retain interior equations when possible.
#' Supplying \code{n_bl} fixes the design and disables selection.
#'
#' The default \code{radius_selection = "guarded_mse"} restricts this search
#' to radii \eqn{r \ge k+4}, where
#' \eqn{k = M - \mathtt{n\_observed} + 1} is the forecast step beyond the
#' original observations. This keeps at least four original observations
#' inside the endpoint-centered boundary test function while feasible.
#' The score itself is unchanged. If the largest capped candidate is below
#' \eqn{k+4}, only that largest candidate is considered and
#' \code{radius_guard_limited} is \code{TRUE}. The guard is a support
#' heuristic, not a covariance calibration. Use
#' \code{radius_selection = "mse"} for the original unguarded search.
#'
#' For covariance, write \eqn{X = \mathrm{vec}(U)}, with Gaussian observation
#' noise of known covariance \eqn{\Sigma}; the supplied parameters, design,
#' and noise SD are held fixed. Write the converged prediction as
#' \eqn{z = h(X)} and the signed weak residual as
#' \eqn{e = \Phi z - r_{obs}(X) - b(z)}, where \eqn{b} contains the missing
#' quadrature and Euler-Maclaurin terms. The fixed point satisfies
#' \eqn{\Phi^T W(X)e = 0}, with
#' \eqn{W = \widetilde\Omega^{-1}},
#' \eqn{\Omega = L\Sigma L^T}, \eqn{L = D_X r_{obs}}, and
#' \eqn{\widetilde\Omega = \Omega + \lambda\,\mathrm{mean}(\mathrm{diag}(\Omega))I}.
#' Implicit differentiation gives, column by column,
#' \deqn{G_j = \frac{\partial h}{\partial X_j}
#'   = [\Phi^T W(\Phi-D_z b)]^{-1}\Phi^T W
#'     [L_j + (\partial_{X_j}\widetilde\Omega)W e].}
#' Here
#' \eqn{\partial_j\Omega = (\partial_j L)\Sigma L^T + L\Sigma(\partial_j L)^T};
#' the relative ridge is differentiated too. Weights are constant during
#' one solve but vary across hypothetical observation-noise realizations.
#' Their derivative is therefore included. Endpoint derivatives and state
#' Hessians are evaluated by central differences. OLS uses \eqn{W=I} and
#' omits the weight derivative.
#'
#' On a smooth converged branch with the usual Gaussian integrability
#' conditions, Stein's identity gives the exact cross-covariance
#' \deqn{\mathrm{Cov}(X,h(X)) = \Sigma\,\mathrm{E}[G(X)^T].}
#' The returned \code{cov_past_next} is \eqn{\Sigma G(X)^T} evaluated at
#' the observed data, a single-realization Stein estimate (up to numerical
#' differentiation), rather than the expectation itself. The returned
#' \code{cov_u} is the first-order approximation \eqn{G\Sigma G^T}; it is
#' not an exact nonlinear variance or a mean squared prediction error.
#' These form the blocks of the local augmented covariance
#' \deqn{\begin{pmatrix}\Sigma & \Sigma G^T \\
#'                       G\Sigma & G\Sigma G^T\end{pmatrix}.}
#' Without \code{state_cov}, errors are independent over times and components,
#' with component SDs supplied in \code{sigma}. Covariance refers to repeated
#' measurements of the same fixed trajectory; it excludes parameter
#' uncertainty and deterministic discretization bias. No covariance is
#' returned for an unconverged solve. Appending a prediction does not add
#' independent observation noise.
#'
#' For recursive prediction, pass \code{U_hat}, \code{tt_hat}, and
#' \code{cov_state} to the next call as \code{U}, \code{tt}, and
#' \code{state_cov}, and carry forward \code{n_observed} unchanged. Without
#' this count, a new call assumes every supplied row has been observed.
#' The returned covariance is ordered as
#' \code{as.vector(U_hat)}, with time varying fastest within each component.
#' It is generally singular because the appended prediction introduces no
#' new independent randomness; a positive semidefinite input is accepted.
#' The numerical Jacobian holds the selected radius and incoming covariance
#' fixed, while differentiating the weak residual and its GLS weights.
#' With automatic radius selection or previously predicted inputs, the
#' returned blocks are local covariance approximations: they do not include
#' radius-switching effects or derivatives of an upstream covariance estimate.
#' The fixed-design Gaussian Stein identity above is not an exact identity
#' for the resulting adaptive recursive estimator.
#'
#' For subsequent first-order propagation through a weak residual, let
#' \eqn{A} and \eqn{B} be its partial Jacobians with respect to the past
#' observations and the appended prediction. Its updated covariance is
#' \deqn{A\Sigma A^T + A C B^T + B C^T A^T + B P B^T,}
#' where \eqn{C=\mathrm{Cov}(X,z)} and \eqn{P=\mathrm{Cov}(z)}.
#' Substituting the local estimates gives
#' \eqn{(A+BG)\Sigma(A+BG)^T}. The cross terms are required because the
#' prediction reuses the observations. Passing \code{cov_state} into the
#' next call performs this propagation when constructing its GLS weights.
#'
#' @inheritParams estimate_endpoint
#' @param U Numeric matrix (M x D) of supplied states, including any previously
#'   predicted points. Supply their joint covariance via \code{state_cov}.
#' @param p Known parameter vector, held fixed. Use the true parameters for
#'   simulation comparisons.
#' @param J_u Callable state Jacobian; required for GLS or covariance.
#' @param J_uu Callable state Hessian; required for GLS or covariance.
#' @param sigma Observation noise standard deviation, positive scalar or
#'   length D, treated as known and fixed. Required for GLS or covariance
#'   unless \code{state_cov} is supplied. With \code{state_cov}, this only
#'   provides component scales for the MSE proxy; it adds no extra noise.
#' @param initial_guess \code{"last"} (default) starts from \code{U[M, ]};
#'   \code{"euler"} starts from \eqn{U_m + h f(p,U_m,t_m)}.
#' @param n_bl Number of BL test functions. If \code{NULL}, GLS selects the
#'   radius with \code{n_bl = r_bl}. OLS (and \code{max_iter = 0}) uses
#'   \code{max(3, ceiling(r_bl / 8))}.
#' @param r_bl Fixed BL radius when \code{n_bl} is supplied, capped at
#'   \code{floor(M / 2)}. Defaults to \code{min(16, floor(M / 2))}.
#' @param r_bl_grid Candidate radii for automatic GLS selection. Defaults to
#'   \code{c(4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 80, 100)}, capped at
#'   \code{floor((M-2)/2)} when retaining interior equations is possible,
#'   otherwise \code{floor(M/2)}. Ignored when \code{n_bl} is supplied.
#' @param radius_selection \code{"guarded_mse"} (default) applies the
#'   original-data support guard before minimizing the MSE proxy;
#'   \code{"mse"} searches all capped candidates. Only used for automatic
#'   GLS selection; explicit \code{n_bl} keeps the requested fixed design.
#' @param radius_search \code{"full"} (default) evaluates every eligible
#'   radius. \code{"adaptive"} starts at the nearest candidate to
#'   \code{radius_hint} and its two neighbors, with a full-grid fallback on
#'   candidate failure, an outward local minimum, poor root conditioning,
#'   or an MSE proxy increase greater than \code{radius_score_ratio}.
#'   A local search is a heuristic and need not find the global grid minimum.
#' @param radius_hint,radius_score Previous selected radius and MSE proxy,
#'   or \code{NULL}. With no hint the full grid is searched. To reproduce a
#'   continued adaptive forecast, pass the previous result's \code{r_bl}
#'   and \code{mse_proxy$obj}, together with its trajectory and covariance.
#' @param radius_refresh Positive integer full-grid refresh interval for
#'   adaptive search, default 5. Refreshes occur at forecast steps
#'   \code{1, 1 + radius_refresh, ...} relative to \code{n_observed}.
#' @param radius_score_ratio Factor triggering full-grid recovery when the
#'   best local proxy exceeds the previous score; default 4, must exceed 1.
#' @param n_observed Number of rows in the original contiguous observed
#'   prefix of \code{U}, a positive integer no greater than \code{nrow(U)}.
#'   Defaults to \code{nrow(U)} for a new forecast. Later rows are previous
#'   predictions. Pass the returned value unchanged when continuing a
#'   forecast, together with its full \code{state_cov}. Observation updates
#'   can extend this prefix, as described in \code{update_prediction}.
#' @param max_iter,tol Fixed-point iteration limit and relative tolerance.
#'   Convergence requires the freshly evaluated fixed-point defect and
#'   Newton correction norm to be at most \code{tol * max(1, ||u||)}, with
#'   reciprocal root-Jacobian condition number above the square root of
#'   machine precision. The best iterate is checked independently of the
#'   preceding step. Set \code{max_iter = 0} to inspect the initial guess.
#' @param inverse \code{"gls"} (default) or \code{"ols"}. GLS requires
#'   \code{J_u} and \code{sigma}; a failed factorization raises an error.
#' @param include_interior Add interior equations as GLS control variates.
#'   Defaults to \code{TRUE}; ignored for OLS.
#' @param ridge Relative ridge added to the observed weak-residual covariance
#'   before factorizing, \code{ridge * mean(diag(Omega))}. GLS only.
#' @param compute_covariance Compute the prediction Jacobian, past-prediction
#'   cross-covariance and linearized prediction variance after convergence.
#'   Defaults to \code{TRUE} for GLS and \code{FALSE} for OLS. Selection uses
#'   covariance internally even when these outputs are disabled.
#' @param compute_mse_proxy Also compute the MSE proxy at a fixed design.
#'   GLS only. Automatic radius selection always computes the proxy.
#' @param state_cov Optional full (M D) by (M D) covariance in
#'   \code{as.vector(U)} order. Must be finite, symmetric, and positive
#'   semidefinite, with positive total variance. Includes all cross-time and
#'   cross-component correlations. Overrides the independent-noise covariance
#'   implied by \code{sigma}; held fixed within this prediction's linearization.
#' @return A list with:
#' \describe{
#'   \item{u_hat, U_hat, tt_hat}{The predicted state, the observations with
#'     that state appended, and the extended time grid. If unconverged,
#'     \code{u_hat} is the finite iterate with the smallest fixed-point defect.}
#'   \item{u_initial, initial_guess}{The actual initial state and the
#'     initialization method.}
#'   \item{prediction_time}{The next time \eqn{t_m + \Delta t}.}
#'   \item{n_observed, forecast_step}{Original observed-prefix length and
#'     prediction's step beyond that prefix. Carry \code{n_observed} forward
#'     when continuing with the returned trajectory and covariance.}
#'   \item{radius_selection, r_bl_min, radius_guard_limited}{Selection rule
#'     actually used (\code{"fixed"} for an explicit design or OLS), guarded
#'     minimum radius (\code{NA} when unused), and whether the available
#'     radius cap prevented the requested original-data support.}
#'   \item{u_history, defect_history}{All iterates, starting with the initial
#'     guess, and the Euclidean norms of their fixed-point defects.}
#'   \item{iters, converged, diverged, defect}{Iteration count, convergence
#'     flag, nonfinite-update flag, and the defect of the returned prediction.}
#'   \item{root_correction, root_rcond}{Newton correction norm and reciprocal
#'     condition number used to certify the returned root. A correction of
#'     \code{Inf} means the residual or conditioning check did not pass.}
#'   \item{radius_search, radius_search_reason, radius_candidates}{Actual
#'     search (\code{"full"}, \code{"adaptive"}, or \code{"fixed"}), why it
#'     was used, and number of candidates evaluated. Full searches include
#'     refreshes and recovery. No candidate is evaluated twice in one solve.}
#'   \item{radius_score_jump}{Whether the selected proxy still exceeds the
#'     previous score by \code{radius_score_ratio}. This is diagnostic;
#'     it does not impose an absolute statistical accuracy threshold.}
#'   \item{r_bl, n_bl, K_bl, K_int, em_order, inverse}{The fixed design and
#'     equation combination used.}
#'   \item{prediction_jacobian}{D by (M D) Jacobian of the converged prediction
#'     with respect to \code{as.vector(U)}, including GLS weight dependence
#'     at the fixed selected design and fixed incoming covariance.}
#'   \item{cov_past_next}{(M D) by D cross-covariance estimate. Rows
#'     \code{n + (0:(D-1)) * M} give the D by D block for past sample n.}
#'   \item{cov_u}{D by D linearized covariance of \code{u_hat}.}
#'   \item{cov_state}{Full augmented covariance ordered as
#'     \code{as.vector(U_hat)}. Pass to the next call as \code{state_cov}.}
#'   \item{weak_residual_cov, weak_residual_ridge}{The propagated
#'     \eqn{L\Sigma L^T} and scalar ridge added to its diagonal for GLS.
#'     Both are \code{NULL} for OLS.}
#'   \item{design, mse_proxy}{Candidate table with \code{r_bl}, \code{n_bl},
#'     \code{obj}, \code{var_obj}, \code{bias_obj}, and \code{status}, and
#'     selected score with its component bias vectors. The design table is
#'     \code{NULL} at a fixed design; its proxy is computed only when requested.
#'     Failed or unconverged candidates are excluded.}
#'   \item{covariance_status}{\code{"computed"}, \code{"disabled"}, or
#'     \code{"not_converged"}. The four state covariance/Jacobian outputs are
#'     \code{NULL} unless computed.}
#' }
#' @references Gaussian integration by parts (Stein's identity):
#'   \url{https://www.math.ucla.edu/~biskup/PIMS/PDFs/lecture5.pdf}, Lemma 5.2.
#' @seealso \code{\link{estimate_endpoint}} for an observed endpoint;
#'   \code{\link{predict_steps}} for recursive forecasts.
#' @export
predict_next <- function(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, p,
                         J_u = NULL, J_uu = NULL, sigma = NULL,
                         initial_guess = c("last", "euler"),
                         n_bl = NULL, r_bl = NULL,
                         max_iter = 100L, tol = 1e-9,
                         em_order = c(4L, 2L), inverse = c("gls", "ols"),
                         include_interior = TRUE, interior_stride = 1L,
                         ridge = 1e-11, compute_covariance = inverse == "gls",
                         state_cov = NULL, r_bl_grid = NULL,
                         compute_mse_proxy = FALSE,
                         radius_selection = c("guarded_mse", "mse"),
                         n_observed = nrow(U),
                         radius_search = c("full", "adaptive"),
                         radius_hint = NULL, radius_score = NULL,
                         radius_refresh = 5L, radius_score_ratio = 4) {
  initial_guess <- match.arg(initial_guess)
  inverse <- match.arg(inverse)
  radius_selection <- match.arg(radius_selection)
  radius_search <- match.arg(radius_search)
  if (!is.numeric(radius_refresh) || length(radius_refresh) != 1L ||
      !is.finite(radius_refresh) || radius_refresh < 1 || radius_refresh != floor(radius_refresh)) {
    stop("radius_refresh must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(radius_score_ratio) || length(radius_score_ratio) != 1L ||
      !is.finite(radius_score_ratio) || radius_score_ratio <= 1) {
    stop("radius_score_ratio must be finite and greater than 1", call. = FALSE)
  }
  if (!is.null(radius_hint) && (!is.numeric(radius_hint) || length(radius_hint) != 1L ||
      !is.finite(radius_hint) || radius_hint < 2 || radius_hint != floor(radius_hint))) {
    stop("radius_hint must be NULL or an integer of at least 2", call. = FALSE)
  }
  if (!is.null(radius_score) && (!is.numeric(radius_score) || length(radius_score) != 1L ||
      !is.finite(radius_score) || radius_score < 0)) {
    stop("radius_score must be NULL or a finite nonnegative previous MSE proxy", call. = FALSE)
  }
  if (!is.logical(compute_covariance) || length(compute_covariance) != 1L ||
      is.na(compute_covariance)) {
    stop("compute_covariance must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(compute_mse_proxy) || length(compute_mse_proxy) != 1L ||
      is.na(compute_mse_proxy)) {
    stop("compute_mse_proxy must be TRUE or FALSE", call. = FALSE)
  }
  if (compute_mse_proxy && inverse != "gls") {
    stop("the MSE proxy requires inverse = 'gls'", call. = FALSE)
  }
  em_order <- as.integer(em_order[1L])
  if (!em_order %in% c(2L, 4L)) stop("em_order must be 2 or 4", call. = FALSE)
  if (!is.matrix(U) || !is.numeric(U) || nrow(U) < 4L || ncol(U) < 1L ||
      any(!is.finite(U))) {
    stop("U must be a finite numeric matrix with at least four rows", call. = FALSE)
  }
  M <- nrow(U)
  D <- ncol(U)
  if (!is.numeric(n_observed) || length(n_observed) != 1L || !is.finite(n_observed) ||
      n_observed < 1 || n_observed > M || n_observed != floor(n_observed)) {
    stop("n_observed must be a positive integer no greater than nrow(U)", call. = FALSE)
  }
  n_observed <- as.integer(n_observed)
  forecast_step <- M - n_observed + 1L
  tt <- as.vector(tt)
  if (!is.numeric(tt) || length(tt) != M || any(!is.finite(tt))) {
    stop("tt must contain one finite time per row of U", call. = FALSE)
  }
  dt <- mean(diff(tt))
  if (dt <= 0 || any(abs(diff(tt) - dt) > 1e-8 * dt)) {
    stop("tt must be equally spaced and increasing", call. = FALSE)
  }
  if (!is.numeric(p) || any(!is.finite(p))) stop("p must be finite numeric parameters", call. = FALSE)
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 0 || max_iter != as.integer(max_iter)) {
    stop("max_iter must be a nonnegative integer", call. = FALSE)
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be positive", call. = FALSE)
  select_radius <- inverse == "gls" && is.null(n_bl) && max_iter > 0L
  r_bl_min <- NA_integer_
  radius_guard_limited <- FALSE
  r_bl_max <- floor(M / 2L)
  if (is.null(r_bl)) r_bl <- min(16L, r_bl_max)
  if (length(r_bl) != 1L || !is.finite(r_bl) || r_bl < 2L || r_bl != as.integer(r_bl)) {
    stop("r_bl must be an integer of at least 2", call. = FALSE)
  }
  r_bl <- min(as.integer(r_bl), r_bl_max)
  if (is.null(n_bl)) n_bl <- max(3L, as.integer(ceiling(r_bl / 8)))
  if (length(n_bl) != 1L || !is.finite(n_bl) || n_bl < 1L || n_bl != as.integer(n_bl)) {
    stop("n_bl must be a positive integer", call. = FALSE)
  }
  n_bl <- as.integer(n_bl)

  # A recursively appended state is correlated with the existing trajectory;
  # do not replace this matrix by marginal variances or add measurement noise.
  if (!is.null(state_cov)) {
    if (!is.matrix(state_cov) || !is.numeric(state_cov) ||
        !identical(dim(state_cov), c(M * D, M * D)) || any(!is.finite(state_cov))) {
      stop("state_cov must be a finite (M D) by (M D) numeric matrix in as.vector(U) order", call. = FALSE)
    }
    scale <- max(abs(state_cov))
    if (scale == 0 || max(abs(state_cov - t(state_cov))) > 1e-10 * scale) {
      stop("state_cov must be symmetric with positive total variance", call. = FALSE)
    }
    state_cov <- (state_cov + t(state_cov)) / 2
    values <- eigen(state_cov, symmetric = TRUE, only.values = TRUE)$values
    if (min(values) < -1e-10 * max(abs(values))) {
      stop("state_cov must be positive semidefinite (singular is allowed)", call. = FALSE)
    }
  }
  noise_sd <- NULL
  if (inverse == "gls" || compute_covariance) {
    if (!is.function(J_u) || !is.function(J_uu)) stop("GLS or covariance requires J_u and J_uu", call. = FALSE)
    if (is.null(sigma) && !is.null(state_cov)) {
      marginal_variance <- colMeans(matrix(diag(state_cov), M, D))
      noise_sd <- sqrt(pmax(marginal_variance, .Machine$double.eps * max(marginal_variance)))
    } else {
      if (!length(sigma) %in% c(1L, D) || any(!is.finite(sigma)) || any(sigma <= 0)) {
        stop("supply state_cov or positive finite sigma of length 1 or D", call. = FALSE)
      }
      noise_sd <- rep_len(sigma, D)
    }
  }
  if (inverse == "gls" && (length(ridge) != 1L || !is.finite(ridge) || ridge < 0)) {
    stop("ridge must be nonnegative", call. = FALSE)
  }

  prediction_time <- tt[M] + dt
  tt_hat <- c(tt, prediction_time)
  observed <- seq_len(M)
  input <- rbind(matrix(rep(p, M), nrow = length(p)), t(U), matrix(tt, nrow = 1L))
  f_observed <- f_(input)
  hess_cache <- if (inverse == "gls" || compute_covariance) {
    build_ic_hessian_cache(U, tt, p, J_uu, D)
  } else NULL

  fit_radius <- function(r_bl, n_bl) {
    bl_system <- build_ic_bl_system(tt_hat, r_bl, n_bl,
                                    include_interior = isTRUE(include_interior) && inverse == "gls",
                                    interior_stride = interior_stride, endpoint = "final")
    r_observed <- -dt * (bl_system$V[, observed, drop = FALSE] %*% f_observed +
                          bl_system$Vp[, observed, drop = FALSE] %*% U)
    B <- bl_system$B
    if (!is.finite(bl_system$BtB) || bl_system$BtB < .Machine$double.eps) {
      stop("the next-endpoint design is degenerate", call. = FALSE)
    }
    sensitivity <- gls <- NULL
    if (inverse == "gls" || compute_covariance) {
      # The future row is an unknown, not an additional independent measurement.
      observed_system <- bl_system
      observed_system$V <- bl_system$V[, observed, drop = FALSE]
      observed_system$Vp <- bl_system$Vp[, observed, drop = FALSE]
      observed_system$support <- intersect(bl_system$support, observed)
      sensitivity <- build_ic_noise_sensitivity(observed_system, U, tt, p, J_u, noise_sd, dt)
      if (!is.null(state_cov)) {
        columns <- as.vector(outer(observed_system$support, (seq_len(D) - 1L) * M, `+`))
        sensitivity$state_cov <- state_cov[columns, columns, drop = FALSE]
      }
    }
    if (inverse == "gls") {
      gls <- build_ic_gls_weights(sensitivity, ridge)
      projection <- solve(gls$PhitWPhi, gls$PhitW)
    } else {
      projection <- t(build_endpoint_matrix(B, D)) / bl_system$BtB
    }
    project <- function(rhs) as.numeric(projection %*% as.vector(rhs))
    projected_observed <- project(r_observed)

    phi_endpoint <- bl_system$boundary_sign * bl_system$phi_endpoint
    V_next <- bl_system$V[, M + 1L]
    Vp_next <- bl_system$Vp[, M + 1L]
    c2 <- dt^2 / 12
    endpoint_rhs <- function(u, order = em_order) {
      input <- matrix(c(p, u, prediction_time), ncol = 1L)
      f_derivs <- list(as.vector(f_(input)), as.vector(dF_dt_(input)),
                       as.vector(d2F_dt2_(input)), as.vector(d3F_dt3_(input)))
      A <- matrix(0, 5L, D)
      A[1:3, ] <- c2 * g_coeffs(f_derivs, u, 1L)
      if (order == 4L) A <- A - dt^4 / 720 * g_coeffs(f_derivs, u, 3L)
      EM <- phi_endpoint %*% A
      r_next <- -dt * (outer(V_next, f_derivs[[1L]]) + outer(Vp_next, u))
      r_next - EM
    }
    update <- function(u) projected_observed + project(endpoint_rhs(u))

    u_initial <- as.numeric(U[M, ])
    if (initial_guess == "euler") {
      u_initial <- u_initial + dt * as.vector(f_(matrix(c(p, u_initial, tt[M]), ncol = 1L)))
    }
    if (length(u_initial) != D || any(!is.finite(u_initial))) {
      stop("the initial guess must be finite with one value per state", call. = FALSE)
    }
    magnitude <- function(x) sqrt(sum(x * x))
    u <- best_u <- u_initial
    mapped <- update(u)
    defect <- if (all(is.finite(mapped))) magnitude(mapped - u) else Inf
    best_defect <- defect
    history <- list(u)
    defect_history <- defect
    iters <- 0L
    converged <- FALSE
    certificate <- list(converged = FALSE, defect = best_defect, correction = Inf, rcond = NA_real_)
    diverged <- !is.finite(defect)
    for (iter in seq_len(max_iter)) {
      if (diverged) break
      iters <- iter
      u_next <- mapped
      mapped <- update(u_next)
      defect <- if (all(is.finite(mapped))) magnitude(mapped - u_next) else Inf
      history[[iter + 1L]] <- u_next
      defect_history[iter + 1L] <- defect
      if (is.finite(defect) && defect < best_defect) {
        best_u <- u_next
        best_defect <- defect
      }
      u <- u_next
      if (!is.finite(defect)) {
        diverged <- TRUE
        break
      }
      if (defect <= tol * max(1, magnitude(u))) {
        certificate <- certify_prediction_root(u, update, endpoint_rhs, projection, tol)
        if (certificate$converged) {
          converged <- TRUE
          best_u <- u
          best_defect <- certificate$defect
          break
        }
      }
    }
    if (!converged && max_iter > 0L) {
      certificate <- certify_prediction_root(best_u, update, endpoint_rhs, projection, tol)
      converged <- certificate$converged
      best_defect <- certificate$defect
    }

    covariance <- list(prediction_jacobian = NULL, cov_past_next = NULL,
                       cov_u = NULL, cov_state = NULL, mse_proxy = NULL)
    covariance_status <- if (!compute_covariance) "disabled" else "not_converged"
    if ((compute_covariance || select_radius || compute_mse_proxy) && converged) {
      covariance <- build_prediction_covariance(
        U, tt, p, J_uu, noise_sd, dt, observed_system, sensitivity, gls,
        r_observed, endpoint_rhs, best_u, ridge, state_cov = state_cov,
        mse_proxy = select_radius || compute_mse_proxy, hess_cache = hess_cache)
      covariance_status <- "computed"
    }

    c(list(u_hat = best_u, U_hat = rbind(U, matrix(best_u, nrow = 1L)), tt_hat = tt_hat,
         u_initial = u_initial, initial_guess = initial_guess,
         prediction_time = prediction_time, u_history = do.call(rbind, history),
         n_observed = n_observed, forecast_step = forecast_step,
         radius_selection = if (select_radius) radius_selection else "fixed",
         r_bl_min = r_bl_min, radius_guard_limited = radius_guard_limited,
         defect_history = defect_history, iters = iters,
         converged = converged, diverged = diverged, defect = best_defect,
         root_correction = certificate$correction, root_rcond = certificate$rcond,
         r_bl = r_bl, n_bl = n_bl, K_bl = bl_system$n_equations,
         K_int = bl_system$n_interior, em_order = em_order, inverse = inverse,
         covariance_status = covariance_status,
         weak_residual_cov = if (is.null(gls)) NULL else gls$Omega,
         weak_residual_ridge = if (is.null(gls)) NULL else gls$ridge), covariance)
  }

  if (!select_radius) {
    fit <- fit_radius(r_bl, n_bl)
    fit$design <- NULL
    fit$radius_search <- "fixed"
    fit$radius_search_reason <- "fixed_design"
    fit$radius_candidates <- 1L
    fit$radius_score_jump <- FALSE
    if (!compute_covariance) {
      fit[c("prediction_jacobian", "cov_past_next", "cov_u", "cov_state")] <- list(NULL, NULL, NULL, NULL)
      fit$covariance_status <- "disabled"
    }
    return(fit)
  }
  if (is.null(r_bl_grid)) r_bl_grid <- c(4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 80, 100)
  if (!is.numeric(r_bl_grid) || !length(r_bl_grid) || any(!is.finite(r_bl_grid)) ||
      any(r_bl_grid < 2 | r_bl_grid != as.integer(r_bl_grid))) {
    stop("r_bl_grid must contain integer radii of at least 2", call. = FALSE)
  }
  # Retain interior equations in every candidate whenever the grid permits.
  selection_max <- if (isTRUE(include_interior) && M >= 6L) floor((M - 2L) / 2L) else r_bl_max
  radii <- sort(unique(pmin(as.integer(r_bl_grid), selection_max)))
  if (radius_selection == "guarded_mse") {
    # Predictions add no independent observations. Keep boundary support on
    # the original data instead of letting a small radius use only forecasts.
    # Restrict the caller's capped grid; do not silently add candidate radii.
    requested_min <- forecast_step + 4L
    r_bl_min <- min(max(radii), requested_min)
    radius_guard_limited <- r_bl_min < requested_min
    radii <- radii[radii >= r_bl_min]
  }
  rows <- vector("list", length(radii))
  best <- NULL
  best_obj <- Inf
  evaluate <- function(indices) {
    for (i in indices) {
      candidate <- tryCatch(fit_radius(radii[i], radii[i]), error = function(err) err)
      valid <- !inherits(candidate, "error") && isTRUE(candidate$converged) &&
        is.finite(candidate$mse_proxy$obj)
      rows[[i]] <<- data.frame(r_bl = radii[i], n_bl = radii[i],
        obj = if (valid) candidate$mse_proxy$obj else NA_real_,
        var_obj = if (valid) candidate$mse_proxy$var_obj else NA_real_,
        bias_obj = if (valid) candidate$mse_proxy$bias_obj else NA_real_,
        status = if (inherits(candidate, "error")) conditionMessage(candidate) else
          if (!candidate$converged) "not_converged" else if (!valid) "nonfinite_proxy" else "ok")
      if (valid && (candidate$mse_proxy$obj < best_obj ||
          (candidate$mse_proxy$obj == best_obj && candidate$r_bl < best$r_bl))) {
        best <<- candidate
        best_obj <<- candidate$mse_proxy$obj
      }
    }
  }
  search_used <- "full"
  search_reason <- if (is.null(radius_hint)) "no_hint" else if (radius_search == "full") "requested" else "refresh"
  indices <- seq_along(radii)
  local <- radius_search == "adaptive" && !is.null(radius_hint) &&
    (forecast_step - 1L) %% radius_refresh != 0L
  if (local) {
    center <- which.min(abs(radii - radius_hint))
    indices <- seq.int(max(1L, center - 1L), min(length(radii), center + 1L))
    search_used <- "adaptive"
    search_reason <- "local"
  }
  evaluate(indices)
  score_jump <- function() !is.null(radius_score) &&
    best_obj > radius_score_ratio * max(radius_score, .Machine$double.eps)
  if (local && length(indices) < length(radii)) {
    failed <- any(vapply(rows[indices], function(row) row$status != "ok", logical(1L)))
    edge <- !is.null(best) && ((best$r_bl == radii[min(indices)] && min(indices) > 1L) ||
      (best$r_bl == radii[max(indices)] && max(indices) < length(radii)))
    ill_conditioned <- !is.null(best) && best$root_rcond < 1e-6
    reasons <- c(if (failed) "candidate_failure", if (edge) "local_edge",
                 if (score_jump()) "score_jump", if (ill_conditioned) "conditioning")
    if (length(reasons)) {
      evaluate(setdiff(seq_along(radii), indices))
      search_used <- "full"
      search_reason <- paste(reasons, collapse = "+")
    }
  }
  if (is.null(best)) stop("no candidate radius produced a converged prediction with a finite MSE proxy", call. = FALSE)
  best$design <- do.call(rbind, rows)
  best$radius_search <- search_used
  best$radius_search_reason <- search_reason
  best$radius_candidates <- nrow(best$design)
  best$radius_score_jump <- score_jump()
  if (!compute_covariance) {
    best[c("prediction_jacobian", "cov_past_next", "cov_u", "cov_state")] <- list(NULL, NULL, NULL, NULL)
    best$covariance_status <- "disabled"
  }
  best
}

#' Predict multiple steps beyond the supplied trajectory
#'
#' Repeatedly applies \code{predict_next}, carrying the full joint state
#' covariance into every subsequent weak-residual covariance and GLS solve.
#' The supplied model parameters remain fixed throughout the forecast.
#'
#' @details
#' At step \eqn{k}, all supplied states and the preceding \eqn{k-1}
#' predictions enter the next solve. The default is GLS with interior
#' equations. With \code{n_bl = NULL}, the MSE proxy selects a radius anew
#' at every step using the extended trajectory and its current covariance.
#' Adaptive search is the default: it checks the previous radius's nearest
#' candidate and its neighbors, refreshes the full grid every five steps,
#' and expands to the full grid on failed candidates, an outward local
#' minimum, poor root conditioning, or a large proxy increase. Set
#' \code{radius_search = "full"} to evaluate every eligible radius each time.
#' The default \code{radius_selection = "guarded_mse"} keeps radii at least
#' \eqn{k+4} at step \eqn{k} beyond the original observations, limited by the
#' largest capped candidate. Set \code{radius_selection = "mse"} to use the
#' original unrestricted candidate search. The score and full covariance
#' propagation are the same for both rules.
#' Supplying \code{n_bl} fixes the design according to \code{predict_next}.
#' Each solve uses the requested last-value or Euler initialization.
#'
#' If the current trajectory has covariance \eqn{\Sigma_k} and the next
#' prediction Jacobian is \eqn{G_k}, the covariance is augmented by the
#' blocks \eqn{C_k=\Sigma_k G_k^T} and \eqn{P_k=G_k\Sigma_k G_k^T}.
#' These blocks retain correlations between original data and predictions,
#' and between different forecast steps. The next weak-residual covariance
#' is \eqn{L_{k+1}\Sigma_{k+1}L_{k+1}^T}. Appended predictions introduce
#' no independent observation noise. The covariance is therefore generally
#' singular; marginal variances cannot replace the full matrix in recursion.
#'
#' To continue a previous forecast or an observation update, supply its
#' \code{U_hat}, \code{tt_hat}, and \code{cov_state} as \code{U}, \code{tt},
#' and \code{state_cov}, and pass its \code{n_observed}. This preserves the
#' original observation boundary, so the guard does not restart at step one.
#' Within one call this count is carried forward automatically. If omitted,
#' all rows of the initial \code{U} are assumed observed.
#' Also pass the last step's \code{r_bl} as \code{radius_hint} and
#' \code{mse_proxy$obj} as \code{radius_score} to reproduce uninterrupted
#' adaptive selection. Omitting the hint starts with a full-grid search.
#' This function forecasts without future observations;
#' use \code{update_prediction} when a new observation arrives.
#'
#' The covariance approximation is the same as in \code{predict_next}:
#' each derivative holds the selected radius and incoming covariance fixed.
#' Radius-switching effects, data dependence of upstream covariance
#' estimates, parameter uncertainty, and deterministic discretization bias
#' are not included. Consequently, the reported SDs are local propagated
#' uncertainty estimates, not calibrated long-horizon prediction intervals.
#' Full one-step results are retained, so memory grows with the horizon.
#' A failed or unconverged solve raises an error identifying the step;
#' an unconverged point is never used as input to a later prediction.
#'
#' @param U Numeric matrix (M x D) of supplied states, as in
#'   \code{predict_next}.
#' @param n_steps Positive integer number of equally spaced future points.
#' @param radius_search \code{"adaptive"} (default) reuses the previous
#'   radius to guide the candidate search, with full-grid recovery;
#'   \code{"full"} evaluates every eligible candidate at every step.
#' @param radius_refresh Positive integer full-grid refresh interval for
#'   adaptive search, default 5. See \code{predict_next} for details.
#' @param ... Uniquely named arguments to \code{predict_next}, including
#'   model evaluators, \code{tt}, fixed \code{p}, \code{J_u}, and
#'   \code{sigma} or \code{state_cov}. Optional solver settings are applied
#'   at every step. Use \code{radius_selection} to choose guarded (default)
#'   or unguarded MSE selection, and \code{n_observed} when continuing a
#'   forecast. Covariance computation is required for recursion;
#'   \code{compute_covariance = FALSE} is rejected. Thus even OLS requires
#'   \code{J_u} and an input covariance or noise SD.
#' @return A list with:
#' \describe{
#'   \item{U_future, tt_future}{The n_steps by D matrix of forecast states
#'     and their times, excluding the original supplied trajectory.}
#'   \item{u_hat, U_hat, tt_hat}{The final forecast state, complete augmented
#'     trajectory, and complete time grid.}
#'   \item{prediction_sd}{n_steps by D matrix of marginal forecast SDs.}
#'   \item{cov_future}{Joint covariance of all forecast states, ordered as
#'     \code{as.vector(U_future)} (time varies fastest within components).}
#'   \item{cov_past_future}{(M D) by (n_steps D) cross-covariance of the
#'     supplied states and all forecast states in their respective vector
#'     orders. The supplied states may include previous forecasts or updates.}
#'   \item{cov_state}{Full covariance in \code{as.vector(U_hat)} order,
#'     suitable for continuing the forecast.}
#'   \item{steps}{List of complete \code{predict_next} results, one per step,
#'     including candidate radius tables and weak-residual covariances.}
#'   \item{diagnostics}{One row per step, with time, selected radius,
#'     number of boundary-layer and interior tests, iteration count,
#'     fixed-point defect, and MSE proxy (NA when not computed). Also includes
#'     \code{forecast_step} relative to the original observations,
#'     \code{r_bl_min}, \code{radius_guard_limited}, root certificates, actual
#'     radius search and reason, candidate counts, and proxy-jump flag.}
#'   \item{n_observed, radius_selection}{Original observed-prefix length
#'     and selection rule used. Pass \code{n_observed} when continuing the
#'     returned trajectory with either predictor.}
#'   \item{radius_search}{Requested search strategy; actual searches and
#'     recovery reasons are recorded in each step and the diagnostics.}
#'   \item{n_steps, converged}{Requested horizon and TRUE. Results are
#'     returned only when every step converges with a finite covariance.}
#' }
#' @seealso \code{\link{predict_next}}, \code{\link{update_prediction}}
#' @export
predict_steps <- function(U, n_steps, ..., radius_search = c("adaptive", "full"),
                          radius_refresh = 5L) {
  radius_search <- match.arg(radius_search)
  if (!is.numeric(n_steps) || length(n_steps) != 1L || !is.finite(n_steps) ||
      n_steps < 1 || n_steps > .Machine$integer.max || n_steps != floor(n_steps)) {
    stop("n_steps must be a positive integer", call. = FALSE)
  }
  n_steps <- as.integer(n_steps)
  args <- list(...)
  if (length(args) && (is.null(names(args)) || anyNA(names(args)) ||
                      any(!nzchar(names(args))) || anyDuplicated(names(args)))) {
    stop("supply predictor arguments by unique names", call. = FALSE)
  }
  if (!is.null(args$compute_covariance) && !isTRUE(args$compute_covariance)) {
    stop("predict_steps requires compute_covariance = TRUE to propagate the joint state covariance",
         call. = FALSE)
  }
  args$compute_covariance <- TRUE
  args$radius_search <- radius_search
  args$radius_refresh <- radius_refresh
  args$U <- U
  steps <- vector("list", n_steps)
  for (k in seq_len(n_steps)) {
    fit <- tryCatch(do.call(predict_next, args), error = function(err) {
      stop(sprintf("prediction step %d of %d failed: %s", k, n_steps, conditionMessage(err)),
           call. = FALSE)
    })
    if (!isTRUE(fit$converged)) {
      stop(sprintf("prediction step %d of %d at time %.10g did not converge (%d iterations, defect %.6g)",
                   k, n_steps, fit$prediction_time, fit$iters, fit$defect), call. = FALSE)
    }
    if (is.null(fit$cov_state) || any(!is.finite(fit$cov_state))) {
      stop(sprintf("prediction step %d of %d did not produce a finite state covariance", k, n_steps),
           call. = FALSE)
    }
    steps[[k]] <- fit
    args$U <- fit$U_hat
    args$tt <- fit$tt_hat
    args$state_cov <- fit$cov_state
    args$n_observed <- fit$n_observed
    args$radius_hint <- fit$r_bl
    args$radius_score <- if (is.null(fit$mse_proxy)) NULL else fit$mse_proxy$obj
    # Keep the caller's design settings. In particular, do not set n_bl
    # to the selected value: NULL requests a fresh MSE selection next time.
  }

  M <- nrow(U)
  D <- ncol(U)
  N <- M + n_steps
  future_rows <- M + seq_len(n_steps)
  past_columns <- as.vector(outer(seq_len(M), (seq_len(D) - 1L) * N, `+`))
  future_columns <- as.vector(outer(future_rows, (seq_len(D) - 1L) * N, `+`))
  cov_future <- fit$cov_state[future_columns, future_columns, drop = FALSE]
  U_future <- fit$U_hat[future_rows, , drop = FALSE]
  prediction_sd <- matrix(sqrt(pmax(0, diag(cov_future))), n_steps, D,
                          dimnames = dimnames(U_future))
  diagnostics <- do.call(rbind, lapply(seq_len(n_steps), function(k) {
    result <- steps[[k]]
    data.frame(step = k, time = result$prediction_time,
      r_bl = result$r_bl, n_bl = result$n_bl, K_int = result$K_int,
      forecast_step = result$forecast_step, r_bl_min = result$r_bl_min,
      radius_guard_limited = result$radius_guard_limited,
      iters = result$iters, defect = result$defect,
      root_correction = result$root_correction, root_rcond = result$root_rcond,
      radius_search = result$radius_search, radius_search_reason = result$radius_search_reason,
      radius_candidates = result$radius_candidates, radius_score_jump = result$radius_score_jump,
      mse_proxy = if (is.null(result$mse_proxy)) NA_real_ else result$mse_proxy$obj)
  }))
  list(U_future = U_future, tt_future = fit$tt_hat[future_rows],
    u_hat = fit$u_hat, U_hat = fit$U_hat, tt_hat = fit$tt_hat,
    prediction_sd = prediction_sd, cov_future = cov_future,
    cov_past_future = fit$cov_state[past_columns, future_columns, drop = FALSE],
    cov_state = fit$cov_state, steps = steps, diagnostics = diagnostics,
    n_steps = n_steps, n_observed = fit$n_observed,
    radius_selection = fit$radius_selection, radius_search = radius_search, converged = TRUE)
}

#' Combine a fixed-point prediction with a new observation
#'
#' Predicts the next state using \code{predict_next}, then updates that state
#' when its independent noisy observation is available. Propagates the data
#' dependence of the gain into the covariance of the updated trajectory.
#'
#' @details
#' Let \eqn{X=\mathrm{vec}(U)}, \eqn{z=h(X)}, and \eqn{y=u_++\eta}, where
#' \eqn{\eta} has zero mean, covariance \eqn{R}, and is independent of \eqn{X}.
#' For prediction error \eqn{e=z-u_+}, write \eqn{P=\mathrm{Cov}(e)},
#' \eqn{b=\mathrm{E}[e]}, and \eqn{M=\mathrm{E}[ee^T]=P+bb^T}.
#' Within the fixed-gain blending rule \eqn{z^+=z+K(y-z)}, minimizing expected
#' squared error gives \eqn{K=M(M+R)^{-1}}. With the true moments, its error
#' second moment is \eqn{M-M(M+R)^{-1}M}. The implementation substitutes local
#' covariance and bias estimates; it does not guarantee this oracle risk
#' reduction. \code{gain = "mse"} uses \code{cov_u} plus the outer product of
#' the summed truncation and statistical bias estimates from the MSE proxy.
#' \code{gain = "variance"} uses \code{cov_u} alone.
#'
#' The selected radius is held fixed when differentiating. For fixed incoming
#' state covariance and fixed \eqn{R}, the gain derivative is
#' \deqn{\partial_j K=(I-K)(\partial_j\widehat M)(\widehat M+R)^{-1}.}
#' Here \eqn{\widehat P=G\Sigma G^T}, so
#' \deqn{\partial_j\widehat M=(\partial_jG)\Sigma G^T+
#' G\Sigma(\partial_jG)^T+(\partial_j\widehat b)\widehat b^T+
#' \widehat b(\partial_j\widehat b)^T,}
#' with the last two terms omitted for the variance gain. This includes
#' second derivatives of the nonlinear fixed-point predictor through
#' \eqn{\partial_jG}.
#' The derivative of \eqn{\widehat M} is evaluated by centered differences
#' of complete predictor reruns, including its GLS weights, converged fixed
#' point, covariance, and (for the MSE gain) bias proxy. With innovation
#' \eqn{\nu=y-z} and prediction Jacobian \eqn{G=D_Xh},
#' \deqn{(J_X)_j=(I-K)G_j+(\partial_jK)\nu,\qquad J_y=K.}
#' These give the local covariance estimates
#' \deqn{\widehat C_{past,+}=\Sigma J_X^T,\qquad
#'        \widehat P_+=J_X\Sigma J_X^T+KRK^T,\qquad
#'        \widehat C_{y,+}=RK^T.}
#' For a smooth fixed-design map and Gaussian \eqn{X} and \eqn{y}, Stein's
#' identity replaces the Jacobians in the cross-covariances by their
#' expectations to give exact population cross-covariances. The returned
#' variance remains first order. Radius-switching effects, dependence of an
#' incoming covariance estimate on the data, and non-Gaussian recursive inputs
#' are subject to the same limitations as \code{predict_next}.
#'
#' The gain-derivative term changes the covariance, not the updated state.
#' \code{cov_u_frozen_gain} and \code{cov_past_next_frozen_gain} omit only
#' this term for comparison. The supplied observation does not enter radius
#' selection or the prediction. The returned trajectory replaces the last
#' predicted row with the updated state. Pass its \code{cov_state} to the
#' next prediction to propagate the updated correlations through the weak
#' residual. The new observation is not appended again as another state.
#' Carry the returned \code{n_observed} into subsequent predictions too.
#' An update immediately following the observed prefix extends its length
#' by one. If unobserved predictions intervene, the prefix stays unchanged;
#' the radius guard conservatively remains anchored to that original prefix.
#'
#' @param U Matrix of past states supplied to \code{predict_next}.
#' @param observation Numeric length-D observation at the next prediction time.
#' @param ... Named arguments to \code{predict_next}, including the model
#'   evaluators, \code{tt}, supplied parameters \code{p}, \code{J_u}, and
#'   \code{sigma} or \code{state_cov}. Prediction covariance is always computed;
#'   the MSE proxy is requested when needed for the gain.
#' @param observation_sd Positive scalar or length-D SD of the independent
#'   new measurement noise. Defaults to the predictor's \code{sigma} if neither
#'   this nor \code{observation_cov} is supplied.
#' @param observation_cov Optional symmetric positive definite D by D
#'   covariance of the new measurement noise. Supply either this or
#'   \code{observation_sd}, not both. Cross-component correlations are allowed.
#' @param gain \code{"mse"} (default) uses the estimated error second moment;
#'   \code{"variance"} ignores the prediction bias proxy.
#' @param derivative_step Positive relative central-difference step for the
#'   gain moment, scaled by \code{max(1, abs(U[j]))}. Defaults to \code{1e-3}.
#' @param compute_covariance Compute the gain derivative and update covariance.
#'   Defaults to \code{TRUE}. Setting \code{FALSE} returns the identical point
#'   update without the costly gain-derivative reruns.
#' @return A list with:
#' \describe{
#'   \item{u_hat, U_hat, tt_hat}{Updated state, augmented trajectory, and grid.}
#'   \item{n_observed}{Length of the contiguous prefix with observations
#'     available, for the radius guard in subsequent predictions.}
#'   \item{prediction}{The original \code{predict_next} result, including the
#'     selected design and pre-observation covariance.}
#'   \item{observation, observation_cov, innovation}{The supplied new data,
#'     its noise covariance, and \eqn{y-z}.}
#'   \item{gain, gain_method, gain_moment, gain_bias}{The gain matrix, requested
#'     gain rule, estimated error second moment, and bias proxy used.}
#'   \item{gain_jacobian, gain_derivative_term}{D by D by (M D) gain derivative
#'     and the D by (M D) contribution \eqn{(\partial_j K)\nu}.}
#'   \item{update_jacobian}{D by (M D) derivative of the updated state with
#'     respect to \code{as.vector(U)}. Its derivative with respect to the new
#'     observation is \code{gain}.}
#'   \item{cov_u, cov_past_next, cov_observation_update}{Updated state variance,
#'     past-update cross-covariance, and observation-update cross-covariance.}
#'   \item{cov_state}{Augmented covariance in \code{as.vector(U_hat)} order,
#'     suitable for the next \code{predict_next} call.}
#'   \item{cov_u_frozen_gain, cov_past_next_frozen_gain}{Covariance estimates
#'     omitting the gain derivative for comparison.}
#'   \item{derivative_steps, covariance_status}{Perturbation sizes (zero for
#'     unsupported observations), and \code{"computed"} or \code{"disabled"}.
#'     Covariance and Jacobian outputs are \code{NULL} when disabled.}
#' }
#' @seealso \code{\link{predict_next}}
#' @export
update_prediction <- function(U, observation, ..., observation_sd = NULL,
                               observation_cov = NULL, gain = c("mse", "variance"),
                               derivative_step = 1e-3, compute_covariance = TRUE) {
  gain <- match.arg(gain)
  if (!is.logical(compute_covariance) || length(compute_covariance) != 1L ||
      is.na(compute_covariance)) stop("compute_covariance must be TRUE or FALSE", call. = FALSE)
  if (length(derivative_step) != 1L || !is.finite(derivative_step) || derivative_step <= 0) {
    stop("derivative_step must be positive", call. = FALSE)
  }
  args <- list(...)
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))) || anyDuplicated(names(args)))) {
    stop("supply predictor arguments by unique names", call. = FALSE)
  }
  args$U <- U
  args$compute_covariance <- TRUE
  args$compute_mse_proxy <- gain == "mse"
  prediction <- do.call(predict_next, args)
  if (!isTRUE(prediction$converged) || is.null(prediction$cov_u)) {
    stop("the prediction must converge before an observation update", call. = FALSE)
  }
  M <- nrow(U)
  D <- ncol(U)
  if (!is.numeric(observation) || length(observation) != D || any(!is.finite(observation))) {
    stop("observation must contain one finite value per state", call. = FALSE)
  }
  observation <- as.numeric(observation)
  if (!is.null(observation_sd) && !is.null(observation_cov)) {
    stop("supply observation_sd or observation_cov, not both", call. = FALSE)
  }
  if (is.null(observation_cov)) {
    if (is.null(observation_sd)) observation_sd <- args$sigma
    if (!is.numeric(observation_sd) || !length(observation_sd) %in% c(1L, D) ||
        any(!is.finite(observation_sd)) || any(observation_sd <= 0)) {
      stop("supply positive observation_sd or observation_cov", call. = FALSE)
    }
    observation_cov <- diag(rep_len(observation_sd, D)^2, D)
  } else {
    if (!is.matrix(observation_cov) || !is.numeric(observation_cov) ||
        !identical(dim(observation_cov), c(D, D)) || any(!is.finite(observation_cov))) {
      stop("observation_cov must be a finite D by D numeric matrix", call. = FALSE)
    }
    if (max(abs(observation_cov - t(observation_cov))) > 1e-10 * max(abs(observation_cov))) {
      stop("observation_cov must be symmetric", call. = FALSE)
    }
    observation_cov <- (observation_cov + t(observation_cov)) / 2
    if (inherits(tryCatch(chol(observation_cov), error = function(err) err), "error")) {
      stop("observation_cov must be positive definite", call. = FALSE)
    }
  }
  moment <- function(fit) {
    bias <- if (gain == "mse") fit$mse_proxy$truncation_bias + fit$mse_proxy$statistical_bias else numeric(D)
    if (length(bias) != D || any(!is.finite(bias))) stop("the gain bias proxy must be finite", call. = FALSE)
    list(M = fit$cov_u + tcrossprod(bias), bias = bias)
  }
  estimate <- moment(prediction)
  S_chol <- chol(estimate$M + observation_cov)
  S_solve <- function(rhs) backsolve(S_chol, backsolve(S_chol, rhs, transpose = TRUE))
  K <- t(S_solve(t(estimate$M)))
  I_minus_K <- diag(D) - K
  innovation <- observation - prediction$u_hat
  u_hat <- as.vector(prediction$u_hat + K %*% innovation)
  U_hat <- prediction$U_hat
  U_hat[M + 1L, ] <- u_hat
  out <- list(u_hat = u_hat, U_hat = U_hat, tt_hat = prediction$tt_hat,
    n_observed = if (prediction$n_observed == M) M + 1L else prediction$n_observed,
    prediction = prediction, observation = observation, observation_cov = observation_cov,
    innovation = innovation, gain = K, gain_method = gain,
    gain_moment = estimate$M, gain_bias = estimate$bias,
    gain_jacobian = NULL, gain_derivative_term = NULL, update_jacobian = NULL,
    cov_u = NULL, cov_past_next = NULL, cov_observation_update = NULL, cov_state = NULL,
    cov_u_frozen_gain = NULL, cov_past_next_frozen_gain = NULL,
    derivative_steps = NULL, covariance_status = "disabled")
  if (!compute_covariance) return(out)

  # Freeze the selected design, but recompute GLS weights and the complete
  # error-moment estimate for every data perturbation. Incoming covariance,
  # parameters, measurement covariance and noise scales remain fixed.
  args$r_bl <- prediction$r_bl
  args$n_bl <- prediction$n_bl
  args$tol <- min(if (is.null(args$tol)) 1e-9 else args$tol, 1e-12)
  bl <- build_ic_bl_system(prediction$tt_hat, prediction$r_bl, prediction$n_bl,
    include_interior = prediction$K_int > 0L,
    interior_stride = if (is.null(args$interior_stride)) 1L else args$interior_stride,
    endpoint = "final")
  support <- intersect(bl$support, seq_len(M))
  columns <- as.vector(outer(support, (seq_len(D) - 1L) * M, `+`))
  K_jacobian <- array(0, c(D, D, M * D))
  derivative_steps <- numeric(M * D)
  for (j in columns) {
    step <- derivative_step * max(1, abs(U[j]))
    plus <- minus <- NULL
    for (attempt in 1:4) {
      up <- down <- args
      up$U[j] <- U[j] + step
      down$U[j] <- U[j] - step
      plus <- do.call(predict_next, up)
      minus <- do.call(predict_next, down)
      if (plus$converged && minus$converged) break
      step <- step / 2
    }
    if (!isTRUE(plus$converged) || !isTRUE(minus$converged)) {
      stop(sprintf("gain derivative perturbation at input %d did not converge", j), call. = FALSE)
    }
    dM <- (moment(plus)$M - moment(minus)$M) / (2 * step)
    K_jacobian[, , j] <- t(S_solve(t(I_minus_K %*% dM)))
    derivative_steps[j] <- step
  }
  correction <- matrix(0, D, M * D)
  for (j in columns) correction[, j] <- matrix(K_jacobian[, , j], D, D) %*% innovation
  frozen <- I_minus_K %*% prediction$prediction_jacobian
  J_update <- frozen + correction
  old <- as.vector(outer(seq_len(M), (seq_len(D) - 1L) * (M + 1L), `+`))
  next_rows <- seq_len(D) * (M + 1L)
  Sigma <- prediction$cov_state[old, old, drop = FALSE]
  C_update <- Sigma %*% t(J_update)
  C_frozen <- Sigma %*% t(frozen)
  measurement_contribution <- K %*% observation_cov %*% t(K)
  P_update <- J_update %*% C_update + measurement_contribution
  P_frozen <- frozen %*% C_frozen + measurement_contribution
  P_update <- (P_update + t(P_update)) / 2
  cov_state <- prediction$cov_state
  cov_state[old, next_rows] <- C_update
  cov_state[next_rows, old] <- t(C_update)
  cov_state[next_rows, next_rows] <- P_update
  out$gain_jacobian <- K_jacobian
  out$gain_derivative_term <- correction
  out$update_jacobian <- J_update
  out$cov_u <- P_update
  out$cov_past_next <- C_update
  out$cov_observation_update <- observation_cov %*% t(K)
  out$cov_state <- cov_state
  out$cov_u_frozen_gain <- (P_frozen + t(P_frozen)) / 2
  out$cov_past_next_frozen_gain <- C_frozen
  out$derivative_steps <- derivative_steps
  out$covariance_status <- "computed"
  out
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
