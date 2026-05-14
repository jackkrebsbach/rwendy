#' Default control for solveWendy
#' @export
default_control <- list(
  optimize = TRUE, # boolean: run IRLS or MLE optimzation for parameters 
  estimate_u0 = FALSE, # Estimate the initial condition
  estimate_U_star = FALSE, # Estimate the state using wendy_erts
  noise_sd = NA, # User can supply sd of the noise if known
  compute_svd = TRUE, # Compute the SVD on the test function matrices (MSG)
  diag_reg = 10e-10, # Regularization for covariance of the weak residual for stability
  max_iterates = 100, # Maximum number of iterates for WENDy-IRLS
  test_fun_type = "MSG",  # Multi-scale Global (MSG) or Single-scale Local (SSL)
  radius_params = 2^(0:3), # Radius factors to use in the MSG test functions
  test_fun = "phi", # phi 𝜱(t; r, n)  𝐂^{∞} Bump test function or psi 𝚿(t; r, p) Piecewise polynomial test function 
  S = 1,  # Euler-Maclaurin series order expansion
  p = 16, # parameters in 𝚿(t; r, p) Piecewise polynomial test function
  radius_min_time = NULL,  #  User set restriction on test function radius
  radius_max_time = NULL,  # User set restriction on test function radius
  k_max = 200, # The max number of test functions
  max_test_fun_condition_number = 1e4, # Truncate SVD of the test function matrices based on the condition number
  min_test_fun_info_number = 0.95, # Cumulative sum of the singular values
  min_number_points = 256, # integer: target number of data points after interpolating
  max_points_interp = 25, # integer: only interpolate data when number of data points is less than this
  interpolation_method = NULL,  # "gp", "poly_ls_N", "spline", "linear", "cubic", "loess", or "kernel"
  fixed_radius = NULL, # integer: fix the base test-function radius, bypassing aujko-selection
  use_interp_uncertainty = TRUE, # logical: if TRUE weight covariance by interpolation uncertainty W_ii = var_ii / sigma^2
  smoother = "erts", # JOINT: state smoother — "gp" (Matern 5/2) or "erts" (EKF/RTS)
  gn_alpha = 0,  # JOINT: Tikhonov weight penalising state departure from smoother estimate
  gn_method = "std",  # JOINT (MSG only): "std" -> gn() with adaptive state penalty; "gls" -> gls_gn() weighted by S(p)^{-1}
  two_stage = TRUE, # JOINT: after MSG gn, refine boundary rows with a lazily-built SSL+BL gn_boundary pass
  include_boundary_layer = FALSE, # JOINT: augment test-function matrices with boundary-layer rows + EM correction (MSG only)
  n_bl = 5, # integer: number of BL test functions per side; NULL uses heuristic max(1, floor(radius/5))
  device = "cpu", # If GPUs are available use "cuda"; resolved to a torch device lazily inside solveWendy
  apply_fn = NULL # function: custom apply for multistart, e.g. parallel::mclapply or future.apply::future_lapply; NULL -> lapply
)