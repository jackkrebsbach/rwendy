library(deSolve)
library(symengine)
library(trust)
library(uGMAR)
library(trustOptim)

test_params <- list(
  radius_params = 2^(0:3),
  radius_min_time = 0.1,
  radius_max_time = 8.0,
  k_max = 200,
  max_test_fun_condition_number = 1e4,
  min_test_fun_info_number = 0.95
)

#' Parameter Estimation for ODE Systems via Maximum Likelihood Estimation
  #'
#' This function estimates parameters of a system of ordinary differential equations (ODEs)
#' The method leverages symbolic derivatives of the ODE right-hand side and a
#' trust-region optimization algorithm.
#'
#' @param f A function of the form \code{f(u, p, t)} defining the ODE right-hand side,
#'   where \code{u} is the state vector, \code{p} is the parameter vector,
#'   and \code{t} is the time variable. Must return symbolic expressions.
#' @param p0 Numeric vector. Initial guess for the parameters.
#' @param U Numeric matrix. Rows represent observed states at time points in \code{tt}.
#'   Columns correspond to state variables.
#' @param tt Numeric vector. Time points corresponding to the rows of \code{U}.
#' @param method String "MLE" | "IRLS" \code{U}.
#'
#' @details
#' The procedure:
#' \itemize{
#'   \item Builds symbolic expressions and derivatives of the ODE system (\eqn{f}, Jacobians, Hessians).
#'   \item Constructs weak negative log-likelihood and functionals based on observed data.
#' }
#'
#'
#' @export
solveWendy <- function(f, p0, U, tt, noise_dist = "addgaussian", lip = F, method = "MLE", optimize = T, compute_svd = T){

  if(noise_dist == "lognormal"){
    data <- preprocess_data(U, tt)
    U <- data$U
    tt <- data$tt
  }

  sig <- estimate_std(U, k = 6)

  test_fun_matrices <- build_full_test_function_matrices(U, tt, test_params, compute_svd)

  V <- test_fun_matrices$V
  Vp <- test_fun_matrices$V_prime

  min_radius <- test_fun_matrices$min_radius

  J <- length(p0)
  D <- ncol(U)
  mp1 <- nrow(U)
  K <- nrow(V)

  u <- do.call(c, lapply(1:ncol(U), function(i) S(paste0("u", i))))
  p <- do.call(c, lapply(1:length(p0), function(i) S(paste0("p", i))))
  t <- S("t")

  f_expr <- switch(noise_dist, lognormal = lognormal_transform(f(u, p, t)), f(u, p, t))

  J_u_sym <- compute_symbolic_jacobian(f_expr, u)
  J_up_sym <- compute_symbolic_jacobian(J_u_sym, p)
  J_p_sym <- compute_symbolic_jacobian(f_expr, p)
  J_pp_sym <- compute_symbolic_jacobian(J_p_sym, p)
  J_upp_sym <- compute_symbolic_jacobian(J_up_sym, p)

  vars <- c(p, u ,t)

  f_ <- build_fn(f_expr, vars)

  J_u <- build_fn(J_u_sym, vars)
  J_up <- build_fn(J_up_sym, vars)
  J_p <- build_fn(J_p_sym, vars)
  J_pp <- build_fn(J_pp_sym, vars)
  J_upp <- build_fn(J_upp_sym, vars)


  F_ <-build_F(U, tt, f_)
  G <- build_G_matrix(V, U, tt, F_, J)
  g <- build_g(V, F_)

  b <- -1 * as.vector(Vp %*% U)
  g0 <- as.vector(V %*% F_(rep(0,J))) # Lip -> Gp + g0 = g(p)
  b1 <- b - g0

  Jp_r <- build_Jp_r(J_p, K, D, J, mp1, V, U, tt)
  Hp_r <- build_Hp_r(J_pp, K, D, J, mp1, V, U, tt)

  L <- build_L(U, tt, J_u, K, V, Vp, sig)
  Jp_L <- build_Jp_L(U, tt, J_up, K, J, D, V, sig)
  Hp_L <- build_Hp_L(U, tt, J_upp, K, J, D, V, sig)

  S <- build_S(L)
  Jp_S <- build_J_S(L, Jp_L, J, K ,D)

  wnll <- build_wnll(S, g, b, K, D)
  J_wnll <- build_J_wnll(S, Jp_S, Jp_r, g, b, J)
  H_wnll <- build_H_wnll(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J)

  objfun <- function(p) {
      f <- wnll(p)
      g <- J_wnll(p)
      h <- H_wnll(p)
    list(value = f, gradient = g, hessian = h)
  }

  res <- list()

  res$wnll <- wnll
  res$J_wnll <- J_wnll
  res$H_wnll <- H_wnll
  res$g <- g
  res$G <- G
  res$b <- b
  res$b1 <- b1
  res$f <- f_
  res$J_p <- J_p
  res$J_p <- J_p
  res$S <- S
  res$Jp_r <- Jp_r
  res$F_ <- F_
  res$f_sym <- f_expr
  res$L <- L
  res$sig <- sig
  res$V <- V
  res$V_prime <- Vp
  res$min_radius <- min_radius

  if(!optimize) return(res)

  res$phat <- switch(method,
                     IRLS = irls(G, b1, L)$p, # IRLS original WENDy
                     trust::trust(objfun, p0, rinit = 25, rmax = 200, blather = TRUE)$argument # Maximum likelihood estimation
                     )
  return(res)
 }







