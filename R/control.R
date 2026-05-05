default_control <- list(
  optimize = TRUE,
  estimate_u0 = TRUE,
  estimate_U_star = TRUE,
  noise_sd = NA,
  compute_svd = TRUE,
  diag_reg = 10e-10, # Regularization for covariance of the weak residual for stability
  max_iterates = 100, # Maximum number of iterates for WENDy-IRLS
  S = 1,  # Euler-Maclaurin series order expansion
  p = 16, # parameters in 𝚿(t; r, p) Piecewise polynomial test function
  test_fun = "psi", # or phi 𝚿(t; r, p) Piecewise polynomial test function or 𝜱(t; r, n) = exp(η/ [(1 - (t/r)²)]₊) 𝐂^{∞} Bump test function
  test_fun_type = "MSG",  # Multi-scale Global (MSG) or Single-scale Local (SSL)
  radius_params = 2^(0:3), # Radius factors to use in the MSG test functions
  radius_min_time = NULL,  #  User set restriction on test function radius
  radius_max_time = NULL,  # User set restriction on test function radius
  k_max = 200, # The max number of test functions
  max_test_fun_condition_number = 1e4, 
  min_test_fun_info_number = 0.95,
  min_number_points = 256,            
  max_points_interp = 25,           # integer: only interpolate data when number of data points is less than this
  interpolation_method = NULL,  # "poly_ls_N", "spline", "linear", "cubic", "loess", or "kernel"
  fixed_radius = NULL,              # integer: fix the base test-function radius, bypassing auto-selection
  use_interp_uncertainty = TRUE,               # logical: if TRUE weight covariance by interpolation uncertainty W_ii = var_ii / sigma^2
  device = torch::torch_device("cpu"), # If GPUs are available use cuda (this speed up is most appropriate for high dimensional data)
  apply_fn = NULL    # function: custom apply for multistart, e.g. parallel::mclapply or future.apply::future_lapply; NULL -> lapply
)
