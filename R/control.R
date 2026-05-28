#' Default control for solveWendy
#' @export
default_control <- list(
  optimize = TRUE, # boolean: run IRLS or MLE optimzation for parameters
  estimate_IC = FALSE, # Estimate the initial condition
  estimate_trajectory = FALSE, # Estimate the state trajectory using wendy_erts
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
  max_points_interp = 25, # integer: maximum number of points in the data to interpolate 
  interpolation_method = NULL,  # single string: "gp", "poly_ls_N", "spline", "linear", "cubic", "loess", or "kernel"
  fixed_radius = NULL, # integer: fix the base test-function radius, bypassing aujko-selection
  use_interp_uncertainty = TRUE, # logical: if TRUE weight covariance by interpolation uncertainty W_ii = var_ii / sigma^2
  smoother = "erts", # state smoother used by estimate_trajectory -- "gp" (Matern 5/2) or "erts" (EKF/RTS)
  apply_fn = NULL, # function: custom apply for multistart, e.g. parallel::mclapply or future.apply::future_lapply; NULL -> lapply
  include_boundary_layer = FALSE, # augment SSL test-function matrices with boundary-layer rows + EM correction
  n_bl = NULL, # integer: BL test functions per side for the SSL augmentation; NULL uses the heuristic max(3, ceiling(radius/4)). Note: the IC estimate always uses this heuristic and ignores this control.
  joint_basis_rank = "auto", # integer, "auto", or "full": rank of the JOINT state basis (SVD modes of the SSL+BL matrix). "auto" keeps the data modes whose energy exceeds joint_rank_tau * noise floor (adapts to n, noise, and signal complexity -- prevents the large-n overfitting a full basis suffers); "full" = no cap (Kb = mp1); integer = top-r smooth modes
  joint_rank_tau = 0, # numeric: threshold (in units of the per-mode noise floor sum_d sigma_d^2) for keeping a basis mode when joint_basis_rank = "auto". Larger = sparser basis / more smoothing
  joint_basis_p = 16, # integer: shape p of the piecewise-polynomial test functions psi(t;r,p)=(1-(t/r)^2)^p used to build the JOINT state basis (before the SVD). Higher = sharper/narrower; lower = wider/flatter
  joint_state_warmstart = FALSE, # logical: before the joint solve, re-initialise the state by optimising c alone with p held at p0 (ODE-consistent warm start), which can reduce joint iterations
  joint_deriv_pen = 1e-4, # numeric >= 0: smoothness penalty rho*||uhat^(q)||^2 on a time-derivative of the basis expansion (via analytic test-function derivatives). 0 = off
  joint_deriv_order = 2, # 1 or 2: derivative order for joint_deriv_pen. 2 = curvature ||uhat''||^2 (penalises wiggle, not trend -- recommended); 1 = ||uhat'||^2 (biases toward flat)
  joint_deriv_ode = FALSE, # logical: make joint_deriv_pen ODE-aware -- penalise ||uhat^(q) - ODE^(q)(uhat,p)||^2 (q=1: f; q=2: dF/dt) instead of raw ||uhat^(q)||^2, so genuine ODE curvature is not penalised
  joint_lambda = 1000, # numeric, or "auto": weight on the weak-residual term in the JOINT loss L(c,p) = lambda*||weak||^2 + ||uhat - U_obs||^2_{Sigma_eps^-1}. "auto" sets lambda = joint_lambda_scale / mean(diag(S(p_irls))), i.e. a dimensionless weight on the weak term whitened by its mean variance (scales like 1/sigma^2; more regularisation at low noise)
  joint_lambda_scale = 300 # numeric: dimensionless weight B used when joint_lambda = "auto". B ~ 300 puts JOINT in the state-reconstruction plateau across logistic/LV/Lorenz and noise 2-20%
)