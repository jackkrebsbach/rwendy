#' Default control for solveWendy
#' @export
default_control <- list(
  optimize = TRUE, # boolean: run IRLS or MLE optimzation for parameters
  estimate_IC = FALSE, # Estimate the initial condition
  estimate_trajectory = FALSE, # Estimate the state trajectory using ERTS 
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
  # ---- JOINT ----
  joint_basis_rank = "auto", # integer, "auto", or "full": rank of the JOINT state basis (SVD modes of the SSL+BL matrix). "auto" keeps the data modes whose energy exceeds joint_rank_tau * noise floor (adapts to n, noise, and signal complexity -- prevents the large-n overfitting a full basis suffers); "full" = no cap (Kb = mp1); integer = top-r smooth modes
  joint_rank_tau = 0, # numeric: threshold (in units of the per-mode noise floor sum_d sigma_d^2) for keeping a basis mode when joint_basis_rank = "auto". Larger = sparser basis / more smoothing
  joint_basis_p = 16, # integer: shape p of the piecewise-polynomial test functions psi(t;r,p)=(1-(t/r)^2)^p used to build the JOINT state basis (before the SVD). Higher = sharper/narrower; lower = wider/flatter
  joint_state_warmstart = FALSE, # logical: before the joint solve, re-initialise the state by optimising c alone with p held at p0 (ODE-consistent warm start), which can reduce joint iterations
  joint_deriv_pen = 1e-4, # numeric >= 0: smoothness penalty rho*||uhat^(q)||^2 on a time-derivative of the basis expansion (via analytic test-function derivatives). 0 = off
  joint_deriv_order = 2, # 1 or 2: derivative order for joint_deriv_pen. 2 = curvature ||uhat''||^2 (penalises wiggle, not trend -- recommended); 1 = ||uhat'||^2 (biases toward flat)
  joint_deriv_ode = FALSE, # logical: make joint_deriv_pen ODE-aware -- penalise ||uhat^(q) - ODE^(q)(uhat,p)||^2 (q=1: f; q=2: dF/dt) instead of raw ||uhat^(q)||^2, so genuine ODE curvature is not penalised
  joint_lambda = 1000, # numeric, or "auto": weight on the weak-residual term in the JOINT loss L(c,p) = lambda*||weak||^2 + ||uhat - U_obs||^2_{Sigma_eps^-1}. "auto" sets lambda = joint_lambda_scale / mean(diag(S(p_irls))), i.e. a dimensionless weight on the weak term whitened by its mean variance (scales like 1/sigma^2; more regularisation at low noise)
  joint_lambda_scale = 300, # numeric: dimensionless weight B used when joint_lambda = "auto". B ~ 300 puts JOINT in the state-reconstruction plateau across logistic/LV/Lorenz and noise 2-20%
  # ---- JSTATE controls (method = "JSTATE": optimise uhat directly + soft-rank-cap projection penalty) ----
  jstate_lambda = 1000, # numeric or "auto": weight on the weak-residual term in the JSTATE loss. "auto" uses the same heuristic as JOINT (jstate_lambda_scale / mean(diag(S(p_irls))))
  jstate_lambda_scale = 300, # numeric: dimensionless weight when jstate_lambda = "auto"
  jstate_pen = 200, # numeric >= 0: weight rho on the out-of-basis penalty rho*||(I - B^T B) uhat||^2. As rho -> infty, uhat is forced into span(B^T) (soft rank cap). 0 = off
  jstate_data_pen = 1, # numeric >= 0: weight mu on the data-fit term mu*||(uhat - U_obs)/sigma||^2. Default 1 = the natural noise-whitened data fit. Increase to pull uhat tighter to the data (risks tracking noise); decrease to let the weak residual + smoothness penalty dominate
  jstate_basis_rank = "auto", # integer >= 1, or "auto" (default): K, number of smooth modes kept in B. "auto" = min(K_eff, jstate_basis_rank_max), where K_eff is the well-conditioned rank of V_ssl (sv_k > sv_max / jstate_basis_cond_max). Integer overrides auto (still capped at K_eff for safety). At small K + large rho the trivial weak-residual solution (uhat=0, p=0 when F(0,p)=0) can win the global optimum
  jstate_basis_rank_max = NULL, # integer or NULL: optional hard cap on the basis rank for "auto" mode, mirroring k_max for MSG. NULL (default) = no cap (K_eff drives everything). MSG's k_max=200 was tuned for MSG's needs; JSTATE's basis must REPRESENT the trajectory, so capping below K_eff (especially at large mp1) starves the basis and trips the trivial-basin failure. Set explicitly only if you need to limit computation
  jstate_basis_cond_max = 1e6, # numeric: condition-number cap on V_ssl's singular values for "auto" mode. Mirrors max_test_fun_condition_number but at 1e6 (vs 1e4 for MSG) -- JSTATE needs more smooth modes than MSG keeps because the basis must REPRESENT the trajectory (MSG only needs to project the weak residual). Modes with sv_k <= sv_max / cond_max are near-null-space and would inject wiggly directions into B
  jstate_basis_radius = "auto", # integer (samples), or "auto": radius for the SSL psi(t;r,p) test functions whose SVD defines B. "auto" lets SSL pick its own rc via compute_r_c_hat from U_sys
  jstate_basis_p = 16, # integer: shape p of the SSL psi(t;r,p)=(1-(t/r)^2)^p test functions used to build B (with BL augmentation)
  jstate_state_warmstart = TRUE# logical: before the joint solve, optimise uhat alone with p held at p0 (ODE-consistent warm start)
)