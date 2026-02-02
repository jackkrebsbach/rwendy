library(symengine)
library(trust)
library(minpack.lm)
library(stats)
library(numbers)

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
  
  cat("Solving WENDy Problem... \n\n")
  noise_dist <- match.arg(noise_dist)
  method <- match.arg(method)

  default_control <- list(
    optimize = TRUE,
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
    min_number_points = 25,
    interpolation_method = "linear",  # "spline" or "linear"
    device = torch::torch_device("cpu") # If GPUs are available
  )
  
  if(!is.null(control)) {
    control <- modifyList(default_control, control)
  } else {
    control <- default_control
  }
  
  # Time spacing must be uniform for WENDy to work
  diff_dt <- diff(as.vector(tt))
  dt <- mean(diff_dt)

  if (max(abs(diff(tt) - dt)) > sqrt(.Machine$double.eps)) {
    cat("Non uniform spacing detected, interpolating data using", control$interpolation_method, "...\n")
    n <- max(floor((max(tt) - min(tt)) / dt), control$min_number_points)
    tt_new <- seq(min(tt), max(tt), length.out = n)

    U <- switch(control$interpolation_method,
      spline = {
        fits <- apply(U, 2, function(col) smooth.spline(tt, col))
        sapply(fits, function(fit) predict(fit, tt_new)$y)
      },
      linear = apply(U, 2, function(col) approx(tt, col, xout = tt_new)$y),
      stop("Unknown interpolation_method: ", control$interpolation_method)
    )

    tt <- matrix(tt_new, ncol = 1)
  }

  if (nrow(U) < control$min_number_points) {
    cat("Warning: Number of time points (", nrow(U), ") is less than min_number_points (",
        control$min_number_points, "). Interpolating to meet minimum requirement...\n", sep = "")

    tt_vec <- as.vector(tt)
    tt_new <- seq(min(tt_vec), max(tt_vec), length.out = control$min_number_points)

    U <- switch(control$interpolation_method,
      spline = {
        fits <- apply(U, 2, function(col) smooth.spline(tt_vec, col))
        sapply(fits, function(fit) predict(fit, tt_new)$y)
      },
      linear = apply(U, 2, function(col) approx(tt_vec, col, xout = tt_new)$y),
      stop("Unknown interpolation_method: ", control$interpolation_method)
    )

    tt <- matrix(tt_new, ncol = 1)
  }

  if(noise_dist == "lognormal"){
    data <- preprocess_data(U, tt) # remove time points with zeros and take log of the data
    U <- data$U
    tt <- data$tt
  }

  torch::torch_set_default_dtype(torch::torch_float64())

  device <- control$device
  sig <- torch::torch_tensor(estimate_std(U, k = 6), dtype = torch::torch_float64(), device = device)

  test_fun_matrices <- if(control$test_fun_type == "SSL"){ 
    build_full_test_function_matrices_ssl(U, tt, control) # Single Scale Local
  } else {
    build_full_test_function_matrices_msg(U, tt, control, control$compute_svd) # Multi Scale Global
  }

  V <- torch::torch_tensor(test_fun_matrices$V, dtype = torch::torch_float64(), device = device) # 𝚽 or 𝚿
  Vp <- torch::torch_tensor(test_fun_matrices$V_prime, dtype = torch::torch_float64(), device = device) # 𝚽̇' or 𝚿'

  min_radius <- test_fun_matrices$min_radius

  J <- length(p0) # Number of parameters
  D <- ncol(U)    # Dimension of system
  mp1 <- nrow(U)  # Number of time points
  K <- nrow(V)    # Number of test functions

  u_expr <- do.call(c, lapply(1:ncol(U), function(i) symengine::S(paste0("u", i))))
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

  F_<- build_F(U, tt, f_, J, device) # F(p,U,t)

  # If linear in parameters the function g(p) is an affine transformation Gp + g0 = g(p).
  # In practice we move g0 to the l.h.s. of the linear system  b - g0 = Gp
  G <- build_G_matrix(V, U, tt, F_, J, device)
  g0 <- torch::torch_mm(V, F_(rep(0,J)))$reshape(c(-1))
  g <- ifelse(!lip, build_g(V, F_),
                   build_g_linear(G, device)
                  )

  b <- -1 * torch::torch_mm(Vp, torch::torch_tensor(U, dtype = torch::torch_float64(), device = device))$reshape(c(-1)) # b = -𝚽'U 
  b <- if (!lip) b else b - g0

  Jp_r <- ifelse(!lip, build_Jp_r(J_p, K, D, J, mp1, V, U, tt, device), # ∇ₚr(p) =  ∇ₚ(g(p) - b) =  ∇ₚg(p) Jacobian of the weak residual
                       build_Jp_r_linear(G)
                      )
  Hp_r <- build_Hp_r(J_pp, K, D, J, mp1, V, U, tt, device) #  ∇ₚ∇ₚr(p) Hessian of the weak residual

  L0 <- build_L0(K, D, mp1, Vp, sig, device) # L0 in factorization of Covariance matrix S = LLᵀ, L = L₁(p) + L₀
  L <- ifelse(!lip, build_L(U, tt, J_u, K, V, L0, sig, J, device),
                    build_L_linear(U, tt, J_u, K, V, L0, sig, J, device)
                  )

  Jp_L <- ifelse(!lip, build_Jp_L(U, tt, J_up, K, J, D, V, sig, device), # Jacobian of covariance factor ∇ₚL
                       build_Jp_L_linear(U, tt, J_u, K, V, L0, sig, J, device)
                      )
  Hp_L <- build_Hp_L(U, tt, J_upp, K, J, D, V, sig, device) # ∇ₚ∇ₚL Hessian of covariance factor

  S <- build_S(L, diag_reg = control$diag_reg) # Covariance of the weak residual S(p)
  Jp_S <- build_J_S(L, Jp_L, J, K ,D) # Jacobian of Covariance of the weak residual ∇ₚS(p)

  wnll <- build_wnll(S, g, b, K, D) # Negative log likelihood of the weak form residual (wnll)
  J_wnll <- build_J_wnll(S, Jp_S, Jp_r, g, b, J) # Jacobian of wnll 
  # Hessian of the wnll
  # When linear in parameters there are terms are guaranteed to be zero so we define a new function 
  H_wnll <- ifelse(!lip, build_H_wnll(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J), 
                         build_H_wnll_linear(S, Jp_S, L, Jp_L, Jp_r, g, b, J)
                        )


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

  if(!control$optimize) return(res)

  data <- switch(method,
                     OLS = if(!lip){ 
                        nols(g, as.array(b$contiguous()), L, Jp_r, p0, reg = 10e-10)
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

  
 cat("Done solving WENDy Problem \n\n")

  return(res)
}
