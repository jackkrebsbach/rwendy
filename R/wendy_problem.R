#' Build a WENDy problem from a dataset.
#'
#' WENDyProblem is constructed from a single dataset and later assembled
#' into a system by build_wendy_system().
#'
#' @param wendy_data  List with elements U (mp1 x D), tt (mp1-vector or column
#'   matrix), var (mp1 x D variance matrix).
#' @param f_,J_u,J_up,J_p,J_pp,J_upp  Callable symbolic derivatives.
#' @param J       Number of parameters.
#' @param lip     Logical; TRUE when f is linear in parameters.
#' @param sig     Numeric vector: estimated noise SD (length 1 or D).
#' @param control Control list (uses $test_fun_type, $use_interp_uncertainty, etc.).
#' @return A list of class "WENDyProblem".
#' @keywords internal
build_wendy_problem <- function(wendy_data, f_, J_u, J_up, J_p, J_pp, J_upp, J, lip, sig, control,
                                dF_dt_ = NULL, d2F_dt2_ = NULL, d3F_dt3_ = NULL,
                                dF_dt_p_ = NULL, d2F_dt2_p_ = NULL, d3F_dt3_p_ = NULL) {
  U   <- wendy_data$U
  tt  <- as.vector(wendy_data$tt)
  dt  <- mean(diff(tt))
  var <- wendy_data$var
  D   <- ncol(U)
  mp1 <- nrow(U)

  if (is.null(control$radius_min_time)) control$radius_min_time <- 2 * dt
  if (is.null(control$radius_max_time)) control$radius_max_time <- floor((mp1 - 1) / 2) * dt

  if (length(sig) != D) sig <- rep(sig[1], D)

  tf <- if (control$test_fun_type == "SSL") {
    build_full_test_function_matrices_ssl(U, tt, control)
  } else {
    build_full_test_function_matrices_msg(U, tt, control, control$compute_svd)
  }
  # The test-function builders return matrices in the raw (sum) convention:
  # V %*% F is the bare sum sum_m phi(t_m) f(t_m), with no quadrature weight.
  # Apply the dt trapezoid weight here, once, so V/Vp represent the weak
  # integrals (int phi f dt, int phi' u dt). Keeping dt out of the builders
  # separates "sample the test functions" from "weight the quadrature".
  # A uniform scaling of V/Vp leaves the GLS estimate and the Fisher covariance
  # invariant (r = g - b and S scale together, the log-det shifts by a
  # p-independent constant), so this only sets the absolute units of g/b/S/wnll.
  K  <- nrow(tf$V)
  V  <- tf$V * dt
  Vp <- tf$V_prime * dt

  # Boundary-layer integration-by-parts term get added to Vp. BL test functions
  # phi do NOT vanish at the endpoints, so the weak form u' = f leaves a boundary
  # term that makes the unbiased BL residual
  #   V F + Vp U + phi(t1) u(t1) - phi(tM) u(tM)   (+ EM quadrature defect in g).
  # The boundary term is linear in U just like Vp U, so we add it as extra entries
  # of Vp at the two endpoint columns:  (Vp_aug U)[k] = (Vp U)[k] + phi_k(t1) u(t1)
  # - phi_k(tM) u(tM). This folds it into BOTH the mean residual (b = -Vp U) AND
  # the covariance (L0 = kron(sig, Vp)). Modelling that channel matters: the
  # endpoint values are noisy observations shared across all BL rows, so when the
  # BL rows carry weight (coarse grids) the unmodelled correlated boundary noise
  # otherwise biases the GLS. phi values are physical units (no dt), matching Vp*dt.
  # bl_rows = K_interior + 1:K_bl. For SSL / MSG-without-SVD the BL rows are the
  # appended block (interior rows have phi(t1)=phi(tM)=0 and are untouched). For
  # MSG-with-SVD the BL rows are mixed into the SVD pool: there the builder reports
  # K_interior=0, K_bl=K so bl_rows=1:K (every mode carries its exact endpoint
  # value, computed by the SVD linear map), and modes with negligible boundary
  # value add ~0.
  if (!is.null(tf$K_bl) && tf$K_bl > 0L) {
    bl_rows <- tf$K_interior + seq_len(tf$K_bl)
    Vp[bl_rows, 1]   <- Vp[bl_rows, 1]   + tf$bl_phi_t1[, "phi0"]
    Vp[bl_rows, mp1] <- Vp[bl_rows, mp1] - tf$bl_phi_tM[, "phi0"]
  }

  F_ <- build_F(U, tt, f_, J)
  G  <- build_G_matrix(V, U, tt, F_, J)
  g0 <- as.vector(V %*% F_(rep(0, J)))

  g <- if (!lip) build_g(V, F_) else build_g_linear(G)

  b_raw <- -as.vector(Vp %*% U)
  b     <- if (!lip) b_raw else b_raw - g0

  Jp_r <- if (!lip) build_Jp_r(J_p, K, D, J, mp1, V, U, tt)
          else      build_Jp_r_linear(G)

  Hp_r <- build_Hp_r(J_pp, K, D, J, mp1, V, U, tt)

  # Euler-Maclaurin defect for the BL rows. Unlike the boundary term it depends
  # on p (through f and its total time derivatives at the boundaries), so it
  # augments the residual via g(p) AND its Jacobian Jp_r. Without it the
  # trapezoidal BL integral carries an O(h^2) bias that grows with dt and
  # corrupts the estimate on coarse grids. The Jacobian dEM/dp is ANALYTIC: the
  # defect is linear in F^(0..3), so its p-derivative is the same combination of
  # the symbolic p-Jacobians dF^(m)/dp (build_em_jacobian) -- no finite
  # differences. The Hessian Hp_r is left as the VF part: the converged estimate
  # is fixed by where the (EM-correct) gradient vanishes, and the reported Fisher
  # covariance uses Jp_r, so only optimizer step quality is affected.
  em_active <- !is.null(tf$K_bl) && tf$K_bl > 0L &&
               !is.null(dF_dt_)   && !is.null(d2F_dt2_)   && !is.null(d3F_dt3_) &&
               !is.null(dF_dt_p_) && !is.null(d2F_dt2_p_) && !is.null(d3F_dt3_p_)
  if (em_active) {
    K_int  <- tf$K_interior
    bl_idx <- K_int + seq_len(tf$K_bl)
    em_fn  <- build_em_correction(tf$bl_phi_t1, tf$bl_phi_tM,
                                  f_, dF_dt_, d2F_dt2_, d3F_dt3_, dt)
    em_jac_fn <- build_em_jacobian(tf$bl_phi_t1, tf$bl_phi_tM,
                                   J_p, dF_dt_p_, d2F_dt2_p_, d3F_dt3_p_, dt, D, J)

    g_base    <- g
    Jp_r_base <- Jp_r
    g <- function(p) g_base(p) + as.vector(rbind(matrix(0, K_int, D), em_fn(U, p, tt)))
    Jp_r <- function(p) {
      Jem <- array(0, c(K, D, J))
      Jem[bl_idx, , ] <- em_jac_fn(U, p, tt)          # K_bl x D x J -> padded K x D x J
      Jp_r_base(p) + matrix(Jem, K * D, J)            # row layout (d-1)*K + k, matches Jp_r
    }
  }

  L0 <- build_L0(K, D, mp1, Vp, sig)

  L <- if (!lip) build_L(       U, tt, J_u, K, V, L0, sig, J)
       else      build_L_linear(U, tt, J_u, K, V, L0, sig, J)

  Jp_L <- if (!lip) build_Jp_L(       U, tt, J_up, K, J, D, V, sig)
          else      build_Jp_L_linear(U, tt, J_u,  K, V, L0, sig, J)

  Hp_L <- build_Hp_L(U, tt, J_upp, K, J, D, V, sig)

  # W is the (mp1*D) diagonal of the interpolation-uncertainty weighting,
  # stored as a vector to avoid materializing a dense diagonal matrix.
  # Column-major flatten of var (mp1 x D) gives W in the (m, b) order, m fast,
  # matching the column index convention for L.
  W <- if (isTRUE(control$use_interp_uncertainty) && !is.null(var)) {
    as.vector(var)
  } else {
    NULL
  }

  structure(
    list(
      U = U, tt = tt, var = var,
      V = V, Vp = Vp, K = K, D = D, mp1 = mp1, J = J,
      min_radius = tf$min_radius, rc = tf$radius_c,
      min_radius_errors = tf$min_radius_errors, min_radius_radii = tf$min_radius_radii,
      rc_errors = tf$rc_errors, rc_radii = tf$rc_radii,
      F_ = F_, g = g, g0 = g0, b = b, G = G,
      Jp_r = Jp_r, Hp_r = Hp_r,
      L0 = L0, L = L, Jp_L = Jp_L, Hp_L = Hp_L,
      W = W,
      # TRUE when the boundary-layer EM correction was folded into g/Jp_r. It
      # makes the residual NONLINEAR in p even for linear-in-p f, so the linear
      # solvers (irls/ols on the constant G) would silently drop it -- downstream
      # routing must use the nonlinear solvers (nirls/nols) when this is set.
      em_active = em_active
    ),
    class = "WENDyProblem"
  )
}


#' Assemble a WENDy problem into system
#'
#' Returns a plain list with everything the optimizers need:
#'   g, b, G, Jp_r, Hp_r, L, Jp_L, Hp_L, W,
#'   S, Jp_S, wnll, J_wnll, H_wnll, K, D, J.
#'
#' @param wendy_problem   A WENDyProblem object.
#' @param lip             Logical; TRUE when f is linear in parameters.
#' @param diag_reg        Diagonal regularisation added to S.
#' @param use_interp_uncertainty  Logical; if TRUE pass through the problem's
#'   W vector (variance weights).
#' @keywords internal
build_wendy_system <- function(wendy_problem, lip, diag_reg, use_interp_uncertainty) {
  J   <- wendy_problem$J
  D   <- wendy_problem$D
  K   <- wendy_problem$K

  g    <- wendy_problem$g
  b    <- wendy_problem$b
  G    <- wendy_problem$G
  Jp_r <- wendy_problem$Jp_r
  Hp_r <- wendy_problem$Hp_r
  L    <- wendy_problem$L
  Jp_L <- wendy_problem$Jp_L
  Hp_L <- wendy_problem$Hp_L
  W    <- if (isTRUE(use_interp_uncertainty)) wendy_problem$W else NULL

  S      <- build_S(L, W, diag_reg = diag_reg)
  Jp_S   <- build_J_S(L, Jp_L, J, K, D, W, diag_reg = diag_reg)
  wnll   <- build_wnll(S, g, b, K, D)
  J_wnll <- build_J_wnll(S, Jp_S, Jp_r, g, b, J)
  H_wnll <- if (!lip) {
    build_H_wnll(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J, W, diag_reg = diag_reg)
  } else {
    build_H_wnll_linear(S, Jp_S, L, Jp_L, Jp_r, g, b, J, W, diag_reg = diag_reg)
  }

  list(
    g = g, b = b, G = G,
    Jp_r = Jp_r, Hp_r = Hp_r,
    L = L, Jp_L = Jp_L, Hp_L = Hp_L, W = W,
    S = S, Jp_S = Jp_S,
    wnll = wnll, J_wnll = J_wnll, H_wnll = H_wnll,
    K = K, D = D, J = J,
    em_active = isTRUE(wendy_problem$em_active)
  )
}
