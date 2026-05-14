#' Build a joint WENDy problem from a dataset.
#'
#' Like build_wendy_problem() but all residual/covariance functions take
#' (U, p, tt) instead of just (p), and b is returned as build_compute_b(U, tt)
#' rather than a fixed vector, so it updates as U changes.
#'
#' @param wendy_data  List with elements U (mp1 x D), tt (mp1-vector or column
#'   matrix), var (mp1 x D variance matrix).
#' @param f_,J_u,J_uu,J_up,J_p,J_pp,J_upp  Callable symbolic derivatives.
#' @param dF_dt_,d2F_dt2_,d3F_dt3_  Trajectory time-derivatives of f_, used only
#'   when control$include_boundary_layer = TRUE for the EM correction. May be
#'   NULL otherwise.
#' @param J       Number of parameters.
#' @param lip     Logical; TRUE when f is linear in parameters.
#' @param sig     Torch scalar: estimated noise SD.
#' @param device  Torch device.
#' @param control Control list (uses $test_fun_type, $include_boundary_layer,
#'   $use_interp_uncertainty, etc.).
#' @return A list of class "WENDyProblemJoint".
#' @keywords internal
build_wendy_problem_joint <- function(wendy_data, f_, J_u, J_uu, J_up, J_p, J_pp, J_upp,
                                      dF_dt_ = NULL, d2F_dt2_ = NULL, d3F_dt3_ = NULL,
                                      J, lip, sig, device, control) {
  U   <- wendy_data$U
  tt  <- as.vector(wendy_data$tt)
  dt  <- mean(diff(tt))
  var <- wendy_data$var
  D   <- ncol(U)
  mp1 <- nrow(U)

  if (is.null(control$radius_min_time)) control$radius_min_time <- 2 * dt
  if (is.null(control$radius_max_time)) control$radius_max_time <- floor((mp1 - 1) / 2) * dt

  if (sig$numel() != D) {
    sig <- sig * torch::torch_ones(D, dtype = torch::torch_float64(), device = device)
  }

  tf <- if (control$test_fun_type == "SSL") {
    build_full_test_function_matrices_ssl(U, tt, control)
  } else {
    build_full_test_function_matrices_msg(U, tt, control, control$compute_svd)
  }
  K   <- nrow(tf$V)

  V  <- torch::torch_tensor(tf$V,       dtype = torch::torch_float64(), device = device)
  Vp <- torch::torch_tensor(tf$V_prime, dtype = torch::torch_float64(), device = device)

  F_ <- build_F_joint(f_, J, device)

  G  <- build_G_matrix_joint(V, U, tt, F_, J, device)
  g0 <- torch::torch_mm(V, F_(U, rep(0, J), tt))$reshape(c(-1))

  # --- Boundary-layer augmentation ----------------------------------------
  # When control$include_boundary_layer is TRUE, the test-function builder
  # returns extra fields describing the BL rows; assemble the bl spec for
  # build_g_joint. Interior rows of phi_t1/phi_tM are padded with zeros so the
  # boundary outer products only contribute on BL rows.
  bl <- NULL
  if (isTRUE(control$include_boundary_layer) && !is.null(tf$K_bl) && tf$K_bl > 0L) {
    if (is.null(dF_dt_) || is.null(d2F_dt2_) || is.null(d3F_dt3_)) {
      stop("control$include_boundary_layer = TRUE requires dF_dt_, d2F_dt2_, d3F_dt3_.")
    }
    K_int      <- tf$K_interior
    K_bl       <- tf$K_bl
    bdry_scale <- tf$bdry_scale %||% 1.0     # 1 for MSG, 1/dt for SSL
    phi_t1_full <- c(rep(0, K_int), bdry_scale * tf$bl_phi_t1[, "phi0"])
    phi_tM_full <- c(rep(0, K_int), bdry_scale * tf$bl_phi_tM[, "phi0"])
    em_fn <- build_em_correction(tf$bl_phi_t1, tf$bl_phi_tM,
                                 f_, dF_dt_, d2F_dt2_, d3F_dt3_,
                                 dt, scale = bdry_scale, device = device)
    bl <- list(
      K_interior    = K_int,
      K_bl          = K_bl,
      phi_t1        = torch::torch_tensor(phi_t1_full, dtype = torch::torch_float64(), device = device),
      phi_tM        = torch::torch_tensor(phi_tM_full, dtype = torch::torch_float64(), device = device),
      em_correction = em_fn
    )
  }

  # Always use the nonlinear g and b so that both update correctly when U changes.
  # The linear shortcut bakes in G at construction-time U and breaks joint U-optimization.
  g <- build_g_joint(V, F_, bl = bl)

  build_compute_b <- function(U, tt) {
    U_t <- torch::torch_tensor(U, dtype = torch::torch_float64(), device = device)
    -1 * torch::torch_mm(Vp, U_t)$reshape(c(-1))
  }

  # Always use the nonlinear Jp_r so gradients w.r.t. p are correct for any U.
  Jp_r <- build_Jp_r_joint(J_p, K, J, V, device)

  Hp_r <- build_Hp_r_joint(J_pp, K, J, V, device)

  L0 <- build_L0_joint(K, D, mp1, Vp, device)

  L <- if (!lip) build_L_joint(      J_u,  K, V, L0, J, device)
       else      build_L_linear_joint(J_u,  K, V, L0, J, device)

  Jp_L <- if (!lip) build_Jp_L_joint(      J_up, K, J, V, device)
          else      build_Jp_L_linear_joint(J_u,  K, V, J, device)

  Hp_L <- build_Hp_L_joint(J_upp, K, J, V, device)

  W <- if (isTRUE(control$use_interp_uncertainty) && !is.null(var)) {
    var_vec <- c(t(var))
    torch::torch_diag(torch::torch_tensor(var_vec, dtype = torch::torch_float64(), device = device))
  } else {
    NULL
  }

  structure(
    list(
      # Data
      U = U, tt = tt, var = var,
      # Test-function tensors and dimensions
      V = V, Vp = Vp, K = K, D = D, mp1 = mp1, J = J, min_radius = tf$min_radius, rc = tf$radius_c,
      # Minimum radius selection (MSG)
      min_radius_errors = tf$min_radius_errors, min_radius_radii = tf$min_radius_radii,
      # Change point radius selection (SSL)
      rc_errors = tf$rc_errors, rc_radii = tf$rc_radii,
      # Residual components
      F_ = F_, g = g, g0 = g0, build_compute_b = build_compute_b, G = G,
      Jp_r = Jp_r, Hp_r = Hp_r,
      # Covariance factor components
      L0 = L0, L = L, Jp_L = Jp_L, Hp_L = Hp_L,
      # Variance weights
      W = W,
      # Boundary-layer spec (NULL when include_boundary_layer is FALSE)
      bl = bl,
      # Stored for J_u_wnll construction in the system builder
      sig = sig, J_u_fn = J_u, J_uu_fn = J_uu
    ),
    class = "WENDyProblemJoint"
  )
}


#' Assemble a list of joint WENDy problems into a system.
#'
#' Returns a plain list with everything the optimizers need.
#' All function closures take (U, p, tt); wnll/J_wnll/H_wnll take (U, p, tt, b).
#' b is not stored here — call build_compute_b(U, tt) before each evaluation.
#'
#' @param wendy_problems  List of WENDyProblemJoint objects.
#' @param lip             Logical; TRUE when f is linear in parameters.
#' @param diag_reg        Diagonal regularisation added to S.
#' @param use_interp_uncertainty  Logical; if TRUE build block-diagonal W.
#' @param device          Torch device.
#' @keywords internal
build_wendy_system_joint <- function(wendy_problems, lip, diag_reg, use_interp_uncertainty, device) {
  J      <- wendy_problems[[1]]$J
  D      <- wendy_problems[[1]]$D
  mp1    <- wendy_problems[[1]]$mp1
  K_list <- sapply(wendy_problems, `[[`, "K")
  K      <- sum(K_list)
  M      <- length(wendy_problems)

  if (M == 1L) {
    prob      <- wendy_problems[[1]]
    g         <- prob$g
    build_compute_b <- prob$build_compute_b
    G         <- prob$G
    Jp_r      <- prob$Jp_r
    Hp_r      <- prob$Hp_r
    L         <- prob$L
    Jp_L      <- prob$Jp_L
    Hp_L      <- prob$Hp_L
    W         <- prob$W
  } else {
    g_fns <- lapply(wendy_problems, `[[`, "g")
    g     <- local({ gs <- g_fns
      function(U, p, tt) torch::torch_cat(lapply(gs, function(gi) gi(U, p, tt)))
    })

    build_compute_b_fns <- lapply(wendy_problems, `[[`, "build_compute_b")
    build_compute_b <- local({ fns <- build_compute_b_fns
      function(U, tt) torch::torch_cat(lapply(fns, function(fn) fn(U, tt)))
    })

    G <- torch::torch_cat(lapply(wendy_problems, `[[`, "G"), dim = 1L)

    Jp_r_fns <- lapply(wendy_problems, `[[`, "Jp_r")
    Jp_r     <- local({
      fns <- Jp_r_fns
      function(U, p, tt) torch::torch_cat(lapply(fns, function(fn) fn(U, p, tt)), dim = 1L)
    })

    Hp_r_fns <- lapply(wendy_problems, `[[`, "Hp_r")
    Hp_r     <- local({
      fns <- Hp_r_fns
      function(U, p, tt) torch::torch_cat(lapply(fns, function(fn) fn(U, p, tt)), dim = 1L)
    })

    L    <- build_L_block_joint(   lapply(wendy_problems, `[[`, "L"),    K_list, D,    device)
    Jp_L <- build_Jp_L_block_joint(lapply(wendy_problems, `[[`, "Jp_L"), K_list, D, J, device)
    Hp_L <- build_Hp_L_block_joint(lapply(wendy_problems, `[[`, "Hp_L"), K_list, D, J, device)

    W <- if (isTRUE(use_interp_uncertainty) && !is.null(wendy_problems[[1]]$W)) {
      var_vec <- unlist(lapply(wendy_problems, function(pr) c(t(pr$var))))
      torch::torch_diag(torch::torch_tensor(var_vec, dtype = torch::torch_float64(), device = device))
    } else {
      NULL
    }
  }

  S      <- build_S_joint(L, W, diag_reg = diag_reg)
  Jp_S   <- build_J_S_joint(L, Jp_L, J, W, diag_reg = diag_reg)
  wnll   <- build_wnll_joint(S, g, K, D)
  J_wnll <- build_J_wnll_joint(S, Jp_S, Jp_r, g, J)
  H_wnll <- if (!lip) {
    build_H_wnll_joint(       S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, J, W, diag_reg = diag_reg)
  } else {
    build_H_wnll_linear_joint(S, Jp_S, L, Jp_L,       Jp_r,       g, J, W, diag_reg = diag_reg)
  }

  # Analytic gradients w.r.t. U; only supported for single-interpolant systems.
  J_u_wnll <- if (M == 1L) {
    prob <- wendy_problems[[1]]
    build_J_u_wnll_joint(
      J_u_fn  = prob$J_u_fn,
      J_uu_fn = prob$J_uu_fn,
      S = S, g = g, L = L,
      K = K, J_num = J, V = prob$V, Vp = prob$Vp,
      device = device
    )
  } else {
    NULL
  }

  J_u_rss <- if (M == 1L) {
    prob <- wendy_problems[[1]]
    build_J_u_rss_joint(
      J_u_fn = prob$J_u_fn,
      g = g, K = K, J_num = J,
      V = prob$V, Vp = prob$Vp,
      device = device
    )
  } else {
    NULL
  }

  J_sig_wnll <- if (M == 1L) {
    prob <- wendy_problems[[1]]
    build_J_sig_wnll_joint(
      J_u_fn = prob$J_u_fn,
      S = S, g = g, L = L,
      K = K, J_num = J, V = prob$V, Vp = prob$Vp,
      W = W, diag_reg = diag_reg,
      device = device
    )
  } else {
    NULL
  }

  list(
    g = g, build_compute_b = build_compute_b, G = G,
    Jp_r = Jp_r, Hp_r = Hp_r,
    L = L, Jp_L = Jp_L, Hp_L = Hp_L, W = W,
    S = S, Jp_S = Jp_S,
    wnll = wnll, J_wnll = J_wnll, H_wnll = H_wnll,
    J_u_wnll = J_u_wnll, J_u_rss = J_u_rss, J_sig_wnll = J_sig_wnll,
    K = K, D = D, J = J
  )
}
