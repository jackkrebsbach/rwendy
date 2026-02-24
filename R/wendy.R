#' @importFrom stats quantile median predict fft lm.fit shapiro.test
#' @importFrom utils modifyList tail
#' @importFrom trust trust
#' @importFrom minpack.lm nls.lm nls.lm.control
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
#' @param p0 Numeric vector. Initial guess for the parameters. Used in MLE or nonlinear least squares solvers.
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
solveWendy <- function(f, p0, U, tt, lip = FALSE, noise_dist = c("addgaussian", "lognormal"),
            method = c("IRLS", "MLE", "OLS"), control = NULL){

  check_suggested_packages()

  noise_dist <- match.arg(noise_dist)
  method <- match.arg(method)

  default_control <- list(
    optimize = TRUE,
    noise_sd = NA,
    compute_svd = TRUE,
    diag_reg = 1e-10,
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
    min_number_points = 25,
    interpolation_method = "linear",  # "spline", "linear", "cubic", "cubic_ls", "loess", or "kernel"
    fixed_radius = NULL,              # integer: fix the base test-function radius, bypassing auto-selection
    device = torch::torch_device("cpu") # If GPUs are available
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
  
  # --- Interpolation: one (U_m, tt) pair per method, all on the same grid ---
  U_orig      <- U
  tt_orig     <- as.vector(tt)
  methods     <- control$interpolation_method
  interp_list <- setNames(lapply(methods, function(m) interpolate_data(U, tt, m, control)), methods)
  tt          <- interp_list[[1]]$tt  # all methods share the same target grid

  device <- control$device

  sig <- if (is.na(control$noise_sd)) {
    if (nrow(U) < 20) {
      U_cubic_fit <- interpolate_to_grid(U_orig, tt_orig, tt_orig, "cubic_ls", substitute_data = FALSE)
      estimated_sd <- sqrt(mean((U_orig - U_cubic_fit)^2))
    } else {
      estimated_sd <- estimate_std(U, k = 6)
    }
    torch::torch_tensor(
      estimated_sd,
      dtype = torch::torch_float64(),
      device = device
    )
  } else {
    torch::torch_tensor(
      control$noise_sd,
      dtype = torch::torch_float64(),
      device = device
    )
  }

  # Build test function matrices per interpolant, then stack V and V' row-wise
  build_tf_matrices <- function(U_m, tt_m) {
    if (control$test_fun_type == "SSL") {
      build_full_test_function_matrices_ssl(U_m, tt_m, control)
    } else {
      build_full_test_function_matrices_msg(U_m, tt_m, control, control$compute_svd)
    }
  }

  tf_list <- lapply(interp_list, function(d) build_tf_matrices(d$U, d$tt))

  V_mat  <- do.call(rbind, lapply(tf_list, `[[`, "V"))
  Vp_mat <- do.call(rbind, lapply(tf_list, `[[`, "V_prime"))

  V  <- torch::torch_tensor(V_mat,  dtype = torch::torch_float64(), device = device) # 𝚽 or 𝚿
  Vp <- torch::torch_tensor(Vp_mat, dtype = torch::torch_float64(), device = device) # 𝚽̇' or 𝚿'

  min_radius <- tf_list[[1]]$min_radius

  J      <- length(p0)
  D      <- ncol(interp_list[[1]]$U)
  mp1    <- nrow(interp_list[[1]]$U)
  K_list <- sapply(tf_list, function(tf) nrow(tf$V))
  K      <- sum(K_list)            # total test functions across all interpolants

  u_expr <- do.call(c, lapply(1:D, function(i) symengine::S(paste0("u", i))))
  p_expr <- do.call(c, lapply(1:length(p0), function(i) symengine::S(paste0("p", i))))
  t_expr <- symengine::S("t")

  f_expr <- switch(noise_dist, # Symbolic representation of the r.h.s.
                    lognormal = lognormal_transform(f(u_expr, p_expr, t_expr)),
                    f(u_expr, p_expr, t_expr)
                   )
  # Compute symbolic gradients of the r.h.s. u̇ = f(p,u,t)
  J_u_sym <- compute_symbolic_jacobian(f_expr, u_expr)
  J_up_sym <- compute_symbolic_jacobian(J_u_sym, p_expr)
  J_p_sym <- compute_symbolic_jacobian(f_expr, p_expr)
  J_pp_sym <- compute_symbolic_jacobian(J_p_sym, p_expr)
  J_upp_sym <- compute_symbolic_jacobian(J_up_sym, p_expr)

  vars <- c(p_expr, u_expr ,t_expr)

  # Callable functions of p, u, and t
  f_ <- build_fn(f_expr, vars)        # f(p,u,t) 
  J_u <- build_fn(J_u_sym, vars)      # ∇ᵤf(p,u,t)      
  J_up <- build_fn(J_up_sym, vars)    # ∇ₚ∇ᵤf(p,u,t)
  J_p <- build_fn(J_p_sym, vars)      # ∇ₚf(p,u,t)
  J_pp <- build_fn(J_pp_sym, vars)    # ∇ₚ∇ₚf(p,u,t)
  J_upp <- build_fn(J_upp_sym, vars)  # ∇ₚ∇ₚ∇ᵤf(p,u,t)

  # Per-interpolant F_m(p): evaluates the ODE rhs on each smoothed dataset
  F_list <- lapply(interp_list, function(d) build_F(d$U, d$tt, f_, J, device))
  F_ <- F_list[[1]]  # kept for backward compat in res$F_

  # V tensors per interpolant (reused below to avoid re-converting)
  V_tensors  <- lapply(tf_list, function(tf)
    torch::torch_tensor(tf$V,       dtype = torch::torch_float64(), device = device))
  Vp_tensors <- lapply(tf_list, function(tf)
    torch::torch_tensor(tf$V_prime, dtype = torch::torch_float64(), device = device))

  # If linear in parameters the function g(p) is an affine transformation Gp + g0 = g(p).
  # In practice we move g0 to the l.h.s. of the linear system  b - g0 = Gp
  # G: rbind(G_1, ..., G_M)
  G <- torch::torch_cat(
    lapply(seq_along(interp_list), function(i)
      build_G_matrix(V_tensors[[i]], interp_list[[i]]$U, interp_list[[i]]$tt, F_list[[i]], J, device)),
    dim = 1L)

  # g0 = cat(V_m F_m(0))
  g0 <- torch::torch_cat(lapply(seq_along(interp_list), function(i)
    torch::torch_mm(V_tensors[[i]], F_list[[i]](rep(0, J)))$reshape(c(-1))))

  # g(p) = cat(V_m F_m(p))  or  Gp  for linear case
  g <- if (!lip) {
    local({
      Vt <- V_tensors; Fl <- F_list
      function(p) torch::torch_cat(lapply(seq_along(Fl), function(i)
        torch::torch_matmul(Vt[[i]], Fl[[i]](p))$reshape(c(-1))))
    })
  } else {
    build_g_linear(G, device)
  }

  # b = cat(-V'_m U_m) # b = -𝚽'U
  b <- torch::torch_cat(lapply(seq_along(interp_list), function(i)
    -1 * torch::torch_mm(Vp_tensors[[i]],
           torch::torch_tensor(interp_list[[i]]$U, dtype = torch::torch_float64(), device = device))$reshape(c(-1))))
  b <- if (!lip) b else b - g0

  # Jp_r: cat(Jp_r_m(p), dim=1)  or  G  for linear case
  Jp_r_fns <- lapply(seq_along(interp_list), function(i)
    build_Jp_r(J_p, K_list[i], D, J, mp1, V_tensors[[i]], interp_list[[i]]$U, interp_list[[i]]$tt, device))
  Jp_r <- if (!lip) {
    local({ fns <- Jp_r_fns
      function(p) torch::torch_cat(lapply(fns, function(fn) fn(p)), dim = 1L) })
  } else {
    build_Jp_r_linear(G)
  }

  # Hp_r: cat(Hp_r_m(p), dim=1)
  Hp_r_fns <- lapply(seq_along(interp_list), function(i)
    build_Hp_r(J_pp, K_list[i], D, J, mp1, V_tensors[[i]], interp_list[[i]]$U, interp_list[[i]]$tt, device))
  Hp_r <- local({ fns <- Hp_r_fns
    function(p) torch::torch_cat(lapply(fns, function(fn) fn(p)), dim = 1L) })

  # L0, L, Jp_L, Hp_L: per-interpolant then block-diagonal
  # S = L L^T is block-diagonal — no cross-covariance between interpolants
  L0_list <- lapply(seq_along(interp_list), function(i)
    build_L0(K_list[i], D, mp1, Vp_tensors[[i]], sig, device))

  L_fns <- lapply(seq_along(interp_list), function(i) {
    if (!lip) build_L(      interp_list[[i]]$U, interp_list[[i]]$tt, J_u,  K_list[i], V_tensors[[i]], L0_list[[i]], sig, J, device)
    else      build_L_linear(interp_list[[i]]$U, interp_list[[i]]$tt, J_u, K_list[i], V_tensors[[i]], L0_list[[i]], sig, J, device)
  })

  Jp_L_fns <- lapply(seq_along(interp_list), function(i) {
    if (!lip) build_Jp_L(      interp_list[[i]]$U, interp_list[[i]]$tt, J_up, K_list[i], J, D, V_tensors[[i]], sig, device)
    else      build_Jp_L_linear(interp_list[[i]]$U, interp_list[[i]]$tt, J_u,  K_list[i], V_tensors[[i]], L0_list[[i]], sig, J, device)
  })

  Hp_L_fns <- lapply(seq_along(interp_list), function(i)
    build_Hp_L(interp_list[[i]]$U, interp_list[[i]]$tt, J_upp, K_list[i], J, D, V_tensors[[i]], sig, device))

  L0   <- build_L0_block(  L0_list,   K_list, D, mp1,    device)
  L    <- build_L_block(   L_fns,     K_list, D, mp1,    device)
  Jp_L <- build_Jp_L_block(Jp_L_fns,  K_list, D, mp1, J, device)
  Hp_L <- build_Hp_L_block(Hp_L_fns,  K_list, D, mp1, J, device)

  S <- build_S(L, diag_reg = control$diag_reg) # Covariance of the weak residual S(p)
  Jp_S <- build_J_S(L, Jp_L, J, K, D, diag_reg = control$diag_reg) # Jacobian of Covariance of the weak residual ∇ₚS(p)

  wnll <- build_wnll(S, g, b, K, D) # Negative log likelihood of the weak form residual (wnll)
  J_wnll <- build_J_wnll(S, Jp_S, Jp_r, g, b, J) # Jacobian of wnll 
  # Hessian of the wnll
  # When linear in parameters there are terms are guaranteed to be zero so we define a new function 
  H_wnll <- if (!lip) build_H_wnll(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J, diag_reg = control$diag_reg) else build_H_wnll_linear(S, Jp_S, L, Jp_L, Jp_r, g, b, J, diag_reg = control$diag_reg)


  res <- list()

  res$wnll <- wnll
  res$J_wnll <- J_wnll
  res$H_wnll <- H_wnll
  res$g <- g
  res$g0 <- g0
  res$G <- G
  res$b <- b
  res$f <- f_
  res$J_p <- J_p
  res$J_u <- J_u
  res$J_upp <- J_upp
  res$S <- S
  res$Jp_S <- Jp_S
  res$Jp_r <- Jp_r
  res$F_ <- F_
  res$f_sym <- f_expr
  res$L <- L
  res$sig <- sig
  res$V <- V
  res$V_prime <- Vp
  res$min_radius <- min_radius
  res$interp_list <- interp_list  # all interpolated datasets
  res$U  <- interp_list[[1]]$U   # first interpolant for backward compat
  res$tt <- tt
  res$interp_methods <- control$interpolation_method

  class(res) <- "wendy"
  attr(res, "call") <- match.call()
  attr(res, "method") <- method
  attr(res, "noise_dist") <- noise_dist
  attr(res, "n_obs") <- mp1
  attr(res, "n_params") <- length(p0)
  attr(res, "n_states") <- D

  if(!control$optimize) return(res)

  data <- switch(method,
                     OLS = if(!lip){ 
                        nols(g, as.array(b$contiguous()), L, Jp_r, p0, reg = 1e-10)
                      } else {
                        ols(as.array(G$contiguous()), as.array(b$contiguous()), L) 
                     },
                     IRLS = if(!lip){
                          nirls(g, as.array(b$contiguous()), L, Jp_r, p0, max_its = control$max_iterates)
                        } else{
                          irls(as.array(G$contiguous()), as.array(b$contiguous()), L, max_its = control$max_iterates)
                      }, # IRLS WENDy / NIRLS WENDy
                     MLE =  mle(p0, wnll, J_wnll, H_wnll, S, Jp_r, control)
                  )
  res$data <- data
  res$phat <- data$p 

  return(res)
}