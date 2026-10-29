# Endpoint estimation uses boundary-layer weak equations with Euler-Maclaurin
# quadrature corrections. Shared helpers below also support parameter fitting.

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
#' @seealso \code{\link{estimate_IC}} for the initial-condition interface.
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
