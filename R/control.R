#' Default control for solveWendy
#' @export
default_control <- list(
  optimize = TRUE, # boolean: run IRLS or MLE optimzation to estimate parameters
  estimate_IC = FALSE, # boolean: Estimate the initial condition
  estimate_trajectory = FALSE, # Estimate the state trajectory using RTS 
  noise_sd = NA, #  sd of the noise if known
  compute_svd = TRUE, # Compute the SVD on the test function matrices (MSG)
  diag_reg = 10e-10, # Regularization for covariance of the weak residual for stability
  max_iterates = 100, # Maximum number of iterates for WENDy-IRLS
  test_fun_type = "MSG",  # Test function type: Multi-scale Global (MSG) or Single-scale Local (SSL)
  radius_params = 2^(0:3), # Radius factors to use in the MSG test functions
  test_fun = "phi", # phi 𝜱(t; r, n)  𝐂^{∞} Bump test function or psi 𝚿(t; r, p) Piecewise polynomial test function 
  S = 1,  # Euler-Maclaurin series order expansion
  p = 16, # parameters in 𝚿(t; r, p) Piecewise polynomial test function
  include_boundary_layer = FALSE, # augment SSL test-function matrices with boundary-layer rows + EM correction
  n_bl = NULL, # integer: Boundary Layer test functions per side for the SSL augmentation
  radius_min_time = NULL,  #  User set restriction on test function radius
  radius_max_time = NULL,  # User set restriction on test function radius
  k_max = 200, # The max number of test functions
  max_test_fun_condition_number = 1e4, # Truncate SVD of the test function matrices based on the condition number
  min_test_fun_info_number = 0.95, # Cumulative sum of the singular values
  min_number_points = 256, # integer: target number of data points after interpolating
  max_points_interp = 25, # integer: maximum number of points in the data to interpolate 
  interpolation_method = NULL,  # single string: "gp", "poly_ls_N", "spline", "linear", "cubic", "loess", or "kernel"
  fixed_radius = NULL, # integer: fix the base test-function radius, bypassing auto-selection
  use_interp_uncertainty = TRUE, # logical: if TRUE weight covariance by interpolation uncertainty W_ii = var_ii / sigma^2
  smoother = "erts", # state smoother used by estimate_trajectory  "gp" (Matern 5/2) or "erts" (EKF/RTS)
  apply_fn = NULL # function: custom apply for multistart, e.g. parallel::mclapply or future.apply::future_lapply; NULL -> lapply
)