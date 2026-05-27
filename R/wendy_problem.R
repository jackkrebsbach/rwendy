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
build_wendy_problem <- function(wendy_data, f_, J_u, J_up, J_p, J_pp, J_upp, J, lip, sig, control) {
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
  K  <- nrow(tf$V)
  V  <- tf$V
  Vp <- tf$V_prime

  F_ <- build_F(U, tt, f_, J)
  G  <- build_G_matrix(V, U, tt, F_, J)
  g0 <- as.vector(V %*% F_(rep(0, J)))

  g <- if (!lip) build_g(V, F_) else build_g_linear(G)

  b_raw <- -as.vector(Vp %*% U)
  b     <- if (!lip) b_raw else b_raw - g0

  Jp_r <- if (!lip) build_Jp_r(J_p, K, D, J, mp1, V, U, tt)
          else      build_Jp_r_linear(G)

  Hp_r <- build_Hp_r(J_pp, K, D, J, mp1, V, U, tt)

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
      W = W
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
    K = K, D = D, J = J
  )
}
