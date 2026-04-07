#' @importFrom stats quantile median predict fft lm.fit shapiro.test
#' @importFrom utils modifyList tail
#' @importFrom trust trust
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @importFrom FME modCost modFit
#' @importFrom numbers mGCD bernoulli_numbers
NULL

#' Check for required suggested packages
#' @keywords internal
check_suggested_packages <- function() {

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("Package 'torch' is required but not installed.\n",
         "Please install it manually with:\n",
         "  install.packages('torch')\n",
         "  torch::install_torch()",
         call. = FALSE)
  }

  if (!requireNamespace("symengine", quietly = TRUE)) {
    stop("Package 'symengine' is required but not installed.\n",
         "Please install it manually with:\n",
         "  install.packages('symengine')",
         call. = FALSE)
  }
}

#' Parameter Estimation for ODE Systems
  #'
#' This function estimates parameters of a system of ordinary differential equations (ODEs)
#' The method leverages symbolic derivatives of the ODE right-hand side and a
#' trust-region optimization algorithm.
#'
#' @param f A function of the form \code{f(u, p, t)} defining the ODE right-hand side,
#'   where \code{u} is the state vector, \code{p} is the parameter vector,
#'   and \code{t} is the time variable. 
#' @param p0 Numeric vector or matrix. Initial guess for the parameters for when. Used in MLE or nonlinear least squares solvers.
#' @param U Numeric matrix. Rows represent observed states at time points in \code{tt}.
#'   Columns correspond to state variables.
#' @param tt Numeric vector. Time points corresponding to the rows of \code{U}.
#' @param method String "MLE" | "IRLS" | "OLS" \code{U}.
#'
#' @details
#' The procedure:
#' \itemize{
#'   \item Builds symbolic expressions and derivatives of the ODE system (\eqn{f}, Jacobians, Hessians).
#'   \item Constructs weak negative log-likelihood and functionals based on observed (noisy) data.
#' }
#'
#'
#' @export
solveWendy <- function(f, U, tt, p0 = NULL, noise_dist = c("addgaussian", "lognormal"), method = c("IRLS", "MLE", "OLS", "OE", "HYBRID"), control = NULL){

  check_suggested_packages()

  noise_dist <- match.arg(noise_dist)
  method <- match.arg(method)

  default_control <- list(
    optimize = TRUE,
    noise_sd = NA,
    compute_svd = TRUE,
    diag_reg = 10e-10,
    max_iterates = 200,
    S = 1,  # Euler-Maclaurin series order expansion
    p = 16, # parameters in 𝚿(t; r, p) Piecewise polynomial test function
    test_fun_type = "MSG",  # Multi-scale Global (MSG) or Single-scale Local (SSL)
    radius_params = 2^(0:3),
    radius_min_time = 0.1,
    radius_max_time = 5.0,
    k_max = 200,
    max_test_fun_condition_number = 1e4,
    min_test_fun_info_number = 0.95,
    min_number_points = 256,            
    max_points_interp = 25,           # integer: only interpolate data when number of data points is less than this
    interpolation_method = NULL,  # "poly_ls_N", "spline", "linear", "cubic", "loess", or "kernel"
    fixed_radius = NULL,              # integer: fix the base test-function radius, bypassing auto-selection
    use_interp_uncertainty = TRUE,               # logical: if TRUE weight covariance by interpolation uncertainty W_ii = var_ii / sigma^2
    device = torch::torch_device("cpu") # If GPUs are available use cuda (this speed up is most appropriate for high dimensional data)
  )
  
  if(!is.null(control)) {
    control <- modifyList(default_control, control)
  } else {
    control <- default_control
  }

  if(noise_dist == "lognormal"){
    data <- preprocess_data(U, tt) # remove time points with zeros and take log of the data
    U <- data$U
    tt <- data$tt
  }
  
  U_orig   <- U
  tt_orig  <- as.vector(tt)

  device  <- control$device
  methods <- control$interpolation_method

  # Compute symbolic variables, functions, and gradients of the r.h.s. u̇ = f(p,u,t)
  J      <- detect_n_params(f)
  u_expr <- do.call(c, lapply(1:ncol(U), function(i) symengine::S(paste0("u", i))))
  p_expr <- do.call(c, lapply(seq_len(J), function(i) symengine::S(paste0("p", i))))
  t_expr <- symengine::S("t")

  f_orig_expr <- f(u_expr, p_expr, t_expr)
  f_expr <- switch(noise_dist,
                    lognormal = lognormal_transform(f_orig_expr),
                    f_orig_expr
                   )

  J_u_sym   <- compute_symbolic_jacobian(f_expr, u_expr)
  J_t_sym   <- compute_symbolic_jacobian(f_expr, c(t_expr))
  J_up_sym  <- compute_symbolic_jacobian(J_u_sym, p_expr)
  J_p_sym   <- compute_symbolic_jacobian(f_expr, p_expr)
  J_pp_sym  <- compute_symbolic_jacobian(J_p_sym, p_expr)
  J_upp_sym <- compute_symbolic_jacobian(J_up_sym, p_expr)

  # A function is linear/affine in p iff all second-order mixed partials
  # ∂²f / ∂pᵢ ∂pⱼ are identically zero.
  lip <- all(as.character(symengine::Vector(array(J_pp_sym))) == "0")

  vars <- c(p_expr, u_expr, t_expr)

  # Callable functions of p, u, and t
  f_    <- build_fn(f_expr,    vars)  # f(p,u,t)
  J_t   <- build_fn(J_t_sym, vars)    # ∇ₜf(p,u,t)
  J_u   <- build_fn(J_u_sym,   vars)  # ∇ᵤf(p,u,t)
  J_up  <- build_fn(J_up_sym,  vars)  # ∇ₚ∇ᵤf(p,u,t)
  J_p   <- build_fn(J_p_sym,   vars)  # ∇ₚf(p,u,t)
  J_pp  <- build_fn(J_pp_sym,  vars)  # ∇ₚ∇ₚf(p,u,t)
  J_upp <- build_fn(J_upp_sym, vars)  # ∇ₚ∇ₚ∇ᵤf(p,u,t)


  estimated_sd <- if (!is.na(control$noise_sd)) {
    control$noise_sd
  } else if (nrow(U) >= 20) {
    estimate_std(U, k = 6)
  } else {
    # Standard SD estimation unreliable for sparse data; use poly LS residuals
    degree               <- detect_max_state_order(f_orig_expr, u_expr) + if (noise_dist == "lognormal") 0L else 1L
    degree_interpolation <- paste0("poly_ls_", degree)
    U_fit <- interpolate_to_grid(U_orig, tt_orig, tt_orig, degree_interpolation, substitute_data = FALSE, control = control)$U
    df    <- nrow(U_orig) - (degree + 1L)  # n - p unbiased estimator of the variance
    sqrt(sum((U_orig - U_fit)^2) / df)
  }
  
  sig <- torch::torch_tensor(estimated_sd, dtype = torch::torch_float64(), device = device)

  # When nrow(U) < max_points_interp  we interpolate the sparse data.
  # Degree rationale: Gaussian uses max_order+1 because integrating a degree-d RHS
  # raises the trajectory degree by 1. For lognormal the log transform reduces the
  # effective degree by 1, so max_order suffices. In both cases, if the
  # noise-to-signal ratio is low (<= 0.1) the observations are clean enough that
  # linear interpolation between them suffices.
  if (nrow(U) < control$max_points_interp && is.null(control$interpolation_method)) {
    nsr <- min(estimated_sd / apply(U, 2, sd))
    if (nsr <= 0.1) {
      control$interpolation_method <- "linear"
    } else {
      max_order     <- detect_max_state_order(f_orig_expr, u_expr)
      interp_degree <- if (noise_dist == "lognormal") max_order else max_order + 1L
      control$interpolation_method <- paste0("poly_ls_", interp_degree)
    }
  }

  if (is.null(control$interpolation_method) && nrow(U) >= control$max_points_interp) {
    wendy_data <- list(none = list(U = U, tt = tt, var = matrix(1.0, nrow = nrow(U), ncol = ncol(U))))
  } else {
    methods <- control$interpolation_method
    wendy_data <- setNames(lapply(methods, function(m) {interpolate_data(U, tt, m, control, sigma = estimated_sd)}), methods)
  }

  D <- ncol(wendy_data[[1]]$U)

  wendy_problems <- lapply(wendy_data, function(d) {
    build_wendy_problem(d, f_, J_u, J_up, J_p, J_pp, J_upp, J, lip, sig, device, control)
  })

  system <- build_wendy_system(wendy_problems, lip, control$diag_reg, control$use_interp_uncertainty, device)

  p1 <- wendy_problems[[1]]

  res <- list()

  res$wnll    <- system$wnll
  res$J_wnll  <- system$J_wnll
  res$H_wnll  <- system$H_wnll
  res$g       <- system$g
  res$g0      <- p1$g0
  res$G       <- system$G
  res$b       <- system$b
  res$f       <- f_
  res$J_p     <- J_p
  res$J_u     <- J_u
  res$J_upp   <- J_upp
  res$S       <- system$S
  res$Jp_S    <- system$Jp_S
  res$Jp_r    <- system$Jp_r
  res$F_      <- p1$F_
  res$f_sym   <- f_expr
  res$L       <- system$L
  res$sig     <- sig
  res$V       <- p1$V
  res$V_prime <- p1$Vp
  res$min_radius    <- p1$min_radius
  res$wendy_problems      <- wendy_problems         
  res$wendy_data   <- wendy_data      
  res$U             <- p1$U             
  res$tt            <- wendy_data[[1]]$tt
  res$var           <- p1$var
  res$wendy_methods  <- names(wendy_data)
  res$W <- system$W

  class(res) <- "wendy"
  attr(res, "call") <- match.call()
  attr(res, "method") <- method
  attr(res, "noise_dist") <- noise_dist
  attr(res, "n_obs") <- p1$mp1
  attr(res, "n_params") <- length(p0)
  attr(res, "n_states") <- D

  if(!control$optimize) return(res)

  g      <- system$g
  b      <- system$b
  G      <- system$G
  L      <- system$L
  W      <- system$W
  Jp_r   <- system$Jp_r
  S      <- system$S
  wnll   <- system$wnll
  J_wnll <- system$J_wnll
  H_wnll <- system$H_wnll
  
  # Output error needs an initial guess
  if(is.null(p0) & method == "OE"){
    stop("Output Error Nonlinear optimization selected, an initial guess, p0, is required")
  }
  # If linear in parameters we do a starting guess with weak ordinary least squares
  if(!lip & is.null(p0)){
    stop("Problem is nonlinear in parameters, an initial guess must be supplied p0")
  } 

  data <- switch(method,
                     # Ordinary Least Squares
                     OLS = if(!lip){
                        nols(g, as.array(b$contiguous()), L, Jp_r, p0, reg = 10e-10)
                      } else {
                        ols(as.array(G$contiguous()), as.array(b$contiguous()), L)
                     }, 
                     # Iterative Reweighted Least Squares
                     IRLS = if(!lip){
                          nirls(g, as.array(b$contiguous()), L, Jp_r, p0, W = W, max_its = control$max_iterates)
                        } else{
                          irls(as.array(G$contiguous()), as.array(b$contiguous()), L, W = W, max_its = control$max_iterates)
                      }, 
                     # Maximum Likelihood Estimation
                     MLE = { 
                        if(lip & is.null(p0)){
                          p0 <- ols(as.array(G$contiguous()), as.array(b$contiguous()), L)$p
                        } 
                       mle(p0, wnll, J_wnll, H_wnll, S, Jp_r, control) 
                     },
                     # Output Error
                       OE = output_error(f, U, tt, p0),
                     # Hybrid MLE + OE
                     HYBRID = {
                        if(lip & is.null(p0)){
                          p0 <- ols(as.array(G$contiguous()), as.array(b$contiguous()), L)$p
                        } 
                       data <- mle(p0, wnll, J_wnll, H_wnll, S, Jp_r, control) 
                       output_error(f, U, tt, data$p)
                     }
                  )
  res$data <- data
  res$phat <- data$p 

  return(res)
}