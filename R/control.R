#' Control parameters for solveWendy
#'
#' Constructs the list of control parameters used by [solveWendy()].
#' Defaults shown in the usage are used for any argument not supplied.
#'
#' @param optimize Logical. Run IRLS or MLE optimization to estimate parameters.
#' @param estimate_IC Logical. Estimate the initial condition.
#' @param estimate_trajectory Logical. Estimate the state trajectory using RTS.
#' @param noise_sd Numeric. Standard deviation of the noise, if known;
#'   \code{NA} to estimate it from the data.
#' @param compute_svd Logical. Compute the SVD on the test function matrices (MSG).
#' @param diag_reg Numeric. Regularization for the covariance of the weak
#'   residual, for stability.
#' @param max_iterates Integer. Maximum number of iterates for WENDy-IRLS.
#' @param test_fun_type Character. Test function type: Multi-scale Global
#'   (\code{"MSG"}) or Single-scale Local (\code{"SSL"}).
#' @param radius_params Numeric vector. Radius factors to use in the MSG test
#'   functions.
#' @param test_fun Character. \code{"phi"} for the \eqn{C^\infty} bump test
#'   function \eqn{\Phi(t; r, n)}, or \code{"psi"} for the piecewise polynomial
#'   test function \eqn{\Psi(t; r, p)}.
#' @param S Integer. Euler-Maclaurin series order expansion.
#' @param p Integer. Parameters in the piecewise polynomial test function
#'   \eqn{\Psi(t; r, p)}.
#' @param include_boundary_layer Logical. Augment SSL test-function matrices
#'   with boundary-layer rows + EM correction.
#' @param n_bl Integer. Boundary-layer test functions per side for the SSL
#'   augmentation.
#' @param radius_min_time Numeric. User-set lower restriction on the test
#'   function radius.
#' @param radius_max_time Numeric. User-set upper restriction on the test
#'   function radius.
#' @param k_max Integer. The maximum number of test functions.
#' @param max_test_fun_condition_number Numeric. Truncate the SVD of the test
#'   function matrices based on this condition number.
#' @param min_test_fun_info_number Numeric. Cumulative sum of the singular
#'   values retained.
#' @param min_number_points Integer. Target number of data points after
#'   interpolating.
#' @param max_points_interp Integer. Maximum number of points in the data to
#'   interpolate.
#' @param interpolation_method Character. One of \code{"gp"},
#'   \code{"poly_ls_N"}, \code{"spline"}, \code{"linear"}, \code{"cubic"},
#'   \code{"loess"}, or \code{"kernel"}.
#' @param fixed_radius Integer. Fix the base test-function radius, bypassing
#'   auto-selection.
#' @param use_interp_uncertainty Logical. If \code{TRUE}, weight the covariance
#'   by interpolation uncertainty \eqn{W_{ii} = var_{ii} / \sigma^2}.
#' @param smoother Character. State smoother used by
#'   \code{estimate_trajectory}: \code{"gp"} (Matern 5/2) or \code{"erts"}
#'   (EKF/RTS).
#' @param apply_fn Function. Custom apply for multistart, e.g.
#'   \code{parallel::mclapply} or \code{future.apply::future_lapply};
#'   \code{NULL} uses \code{lapply}.
#'
#' @return A named list of control parameters for [solveWendy()].
#' @export
wendy_control <- function(
  optimize = TRUE,
  estimate_IC = FALSE,
  estimate_trajectory = FALSE,
  noise_sd = NA,
  compute_svd = TRUE,
  diag_reg = 10e-10,
  max_iterates = 100,
  test_fun_type = "MSG",
  radius_params = 2^(0:3),
  test_fun = "phi",
  S = 1,
  p = 16,
  include_boundary_layer = FALSE,
  n_bl = NULL,
  radius_min_time = NULL,
  radius_max_time = NULL,
  k_max = 200,
  max_test_fun_condition_number = 1e4,
  min_test_fun_info_number = 0.95,
  min_number_points = 256,
  max_points_interp = 25,
  interpolation_method = NULL,
  fixed_radius = NULL,
  use_interp_uncertainty = TRUE,
  smoother = "erts",
  apply_fn = NULL
) {
  list(
    optimize = optimize,
    estimate_IC = estimate_IC,
    estimate_trajectory = estimate_trajectory,
    noise_sd = noise_sd,
    compute_svd = compute_svd,
    diag_reg = diag_reg,
    max_iterates = max_iterates,
    test_fun_type = test_fun_type,
    radius_params = radius_params,
    test_fun = test_fun,
    S = S,
    p = p,
    include_boundary_layer = include_boundary_layer,
    n_bl = n_bl,
    radius_min_time = radius_min_time,
    radius_max_time = radius_max_time,
    k_max = k_max,
    max_test_fun_condition_number = max_test_fun_condition_number,
    min_test_fun_info_number = min_test_fun_info_number,
    min_number_points = min_number_points,
    max_points_interp = max_points_interp,
    interpolation_method = interpolation_method,
    fixed_radius = fixed_radius,
    use_interp_uncertainty = use_interp_uncertainty,
    smoother = smoother,
    apply_fn = apply_fn
  )
}

# Internal default control parameters; see wendy_control() for documentation.
default_control <- wendy_control()
