#' @importFrom MASS ginv
#' @importFrom stats quantile median predict fft mvfft lm.fit shapiro.test sd rnorm optim setNames
#' @importFrom utils modifyList tail capture.output
#' @importFrom trust trust
#' @importFrom minpack.lm nls.lm
NULL


#' Parameter Estimation for ODE Systems
#'
#' Estimates parameters of a system of ordinary differential equations (ODEs).
#' Leverages symbolic derivatives of the ODE right-hand side and a
#' trust-region optimization algorithm.
#'
#' @param f A function of the form \code{f(u, p, t)} defining the ODE right-hand side,
#'   where \code{u} is the state vector, \code{p} is the parameter vector,
#'   and \code{t} is the time variable. When \code{NULL}, the right-hand side
#'   is discovered from the data with weak-form SINDy (see [solveWSINDy()])
#'   and then refined through the normal WENDy pipeline; the returned object
#'   carries the discovery fit in \code{$wsindy}.
#' @param p0 Numeric vector or matrix. Initial parameter guess. Used by MLE,
#'   nonlinear least squares, and OE solvers; computed automatically when \code{NULL}.
#' @param U Numeric matrix. Rows are observed states at the time points in \code{tt};
#'   columns are state variables.
#' @param tt Numeric vector. Time points corresponding to rows of \code{U}.
#' @param noise_dist One of \code{"addgaussian"} (default) or \code{"lognormal"}.
#' @param method One of \code{"IRLS"}, \code{"OLS"}, \code{"MLE"},
#'   \code{"OE"}, or \code{"HYBRID"}.
#' @param control Named list of control parameters, as produced by
#'   [wendy_control()]; entries supplied here override the defaults.
#'   See [wendy_control()] for the documentation of each option.
#'
#' @details
#' \itemize{
#'   \item Builds symbolic expressions and derivatives of the ODE system (\eqn{f}, Jacobians, Hessians).
#'   \item Constructs weak negative log-likelihood and functionals based on observed (noisy) data.
#' }
#'
#' @return An object of class \code{wendy} with elements including \code{phat}
#'   (the parameter estimate), \code{wnll}/\code{J_wnll}/\code{H_wnll} (callables
#'   for the weak NLL and its derivatives), and the assembled \code{V},
#'   \code{V_prime}, \code{S}, and \code{L} components.
#'
#' @examples
#' # Logistic growth: u' = p1*u - p2*u^2 (linear in p; p_star = (1, 1/10)).
#' f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)
#' p_star  <- c(1, 1/10)
#' u0      <- c(0.1)
#' npoints <- 128
#' t_eval  <- seq(0, 10, length.out = npoints)
#'
#' modelODE <- function(tvec, state, parameters) {
#'   list(as.vector(f(state, parameters, tvec)))
#' }
#' sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)
#'
#' # Additive Gaussian noise at 5% of the RMS signal level.
#' set.seed(8675309 + 1)
#' U_vec    <- as.vector(sol[, -1])
#' noise_sd <- 0.05 * sqrt(mean(U_vec^2))
#' U  <- sol[, 2, drop = FALSE] + rnorm(npoints, sd = noise_sd)
#' tt <- sol[, 1, drop = FALSE]
#'
#' res <- solveWendy(
#'   f, U, tt,
#'   method     = "IRLS",
#'   noise_dist = "addgaussian",
#'   control    = list(estimate_trajectory = FALSE, optimize = TRUE,
#'                     test_fun_type   = "SSL",  estimate_IC = TRUE)
#' )
#' res$phat                       # ~ c(1, 0.1)
#' rel_err(res$phat, p_star)
#'
#' @export
solveWendy <- function(f = NULL, U, tt, p0 = NULL, noise_dist = c("addgaussian", "lognormal"), method = c("IRLS", "MLE", "OLS", "OE", "HYBRID"), control = NULL){
  noise_dist <- match.arg(noise_dist)
  method <- match.arg(method)

  if (!is.null(control)) {
    unknown <- setdiff(names(control), names(default_control))
    if (length(unknown) > 0) {
      stop("Unknown control parameter(s): ", paste(unknown, collapse = ", "),
           ". See ?wendy_control for valid options.")
    }
    control <- modifyList(default_control, control)
  } else {
    control <- default_control
  }

  # No rhs supplied: discover the model structure with WSINDy, then run the
  # discovered (linear-in-p) f through the normal WENDy pipeline below.
  ws <- NULL
  if (is.null(f)) {
    if (noise_dist == "lognormal") {
      stop("WSINDy discovery (f = NULL) assumes additive Gaussian noise; ",
           "noise_dist = \"lognormal\" is not supported in this mode. ",
           "Supply f explicitly, or use noise_dist = \"addgaussian\".")
    }
    ws <- solveWSINDy(U, tt, control = control)
    if (all(ws$W == 0)) {
      stop("WSINDy discovered no terms in any equation. ",
           "Try increasing wsindy_poly_deg, lowering the noise, or adjusting ",
           "wsindy_lambdas. See ?solveWSINDy.")
    }
    f <- ws$f
    if (is.null(p0)) p0 <- ws$p0
  }

  U_orig   <- U
  tt_orig  <- as.vector(tt)

  if(noise_dist == "lognormal"){
    data <- preprocess_data(U, tt) # remove time points with zeros and take log of the data
    U <- data$U
    tt <- data$tt
  }

  # Compute symbolic variables, functions, and gradients of the r.h.s. u̇ = f(p,u,t)
  J <- detect_n_params(f)
  u_expr <- sym_build(lapply(1:ncol(U), function(i) sym_symbol(paste0("u", i))))
  p_expr <- sym_build(lapply(seq_len(J), function(i) sym_symbol(paste0("p", i))))
  t_expr <- sym_symbol("t")

  f_orig_expr <- f(u_expr, p_expr, t_expr)
  f_expr <- switch(noise_dist, lognormal = lognormal_transform(f_orig_expr), f_orig_expr)

  J_p_sym  <- compute_symbolic_jacobian(f_expr, p_expr)
  J_pp_sym <- compute_symbolic_jacobian(J_p_sym, p_expr)

  # A function is linear/affine in p iff all second-order mixed partials
  # ∂²f / ∂pᵢ ∂pⱼ are identically zero.
  lip <- all(sym_strings(J_pp_sym) == "0")

  vars <- c(p_expr, u_expr, t_expr)

  # Callable functions needed by all methods
  f_ <- build_fn(f_expr,   vars)  # f(p,u,t)
  J_p <- build_fn(J_p_sym, vars)  # ∇ₚf(p,u,t)

  if (method != "OE") {
    J_u_sym   <- compute_symbolic_jacobian(f_expr, u_expr)
    J_t_sym   <- compute_symbolic_jacobian(f_expr, c(t_expr))
    J_up_sym  <- compute_symbolic_jacobian(J_u_sym, p_expr)
    J_upp_sym <- compute_symbolic_jacobian(J_up_sym, p_expr)
    J_uu_sym  <- compute_symbolic_jacobian(J_u_sym, u_expr)

    J_t   <- build_fn(J_t_sym,   vars)  # ∇ₜf(p,u,t)
    J_u   <- build_fn(J_u_sym,   vars)  # ∇ᵤf(p,u,t)
    J_uu  <- build_fn(J_uu_sym,  vars)  # ∇ᵤ∇ᵤf(p,u,t)
    J_up  <- build_fn(J_up_sym,  vars)  # ∇ₚ∇ᵤf(p,u,t)
    J_pp  <- build_fn(J_pp_sym,  vars)  # ∇ₚ∇ₚf(p,u,t)
    J_upp <- build_fn(J_upp_sym, vars)  # ∇ₚ∇ₚ∇ᵤf(p,u,t)

    dF_sym   <- compute_symbolic_total_time_deriv(f_expr,  u_expr, f_expr, t_expr)
    d2F_sym  <- compute_symbolic_total_time_deriv(dF_sym,  u_expr, f_expr, t_expr)
    d3F_sym  <- compute_symbolic_total_time_deriv(d2F_sym, u_expr, f_expr, t_expr)

    # Total derivatives with respect to time
    dF_dt_   <- build_fn(dF_sym,  vars)  # d f/dt   (total, along trajectory)
    d2F_dt2_ <- build_fn(d2F_sym, vars)  # d²f/dt²  (total)
    d3F_dt3_ <- build_fn(d3F_sym, vars)  # d³f/dt³  (total)
  }

  estimated_sd <- if (!is.na(control$noise_sd)) {
    control$noise_sd
  } else if (nrow(U) >= 20) {
    estimate_std(U, k = 6)
  } else {
    # Standard SD estimation unreliable for sparse data; use poly LS residuals
    degree <- detect_max_state_order(f_orig_expr, u_expr) + if (noise_dist == "lognormal") 0L else 1L
    degree_interpolation <- paste0("poly_ls_", degree)
    U_fit <- interpolate_to_grid(U_orig, tt_orig, tt_orig, degree_interpolation, substitute_data = FALSE, control = control)$U
    df <- nrow(U_orig) - (degree + 1L)  # n - p unbiased estimator of the variance
    sqrt(sum((U_orig - U_fit)^2) / df)
  }

  D <- ncol(U)
  res <- list()

  if (method != "OE") {
    sig <- estimated_sd

    # When nrow(U) < max_points_interp we interpolate the sparse data.
    # Degree rationale: Gaussian uses max_order+1 because integrating a degree-d RHS
    # raises the trajectory degree by 1. For lognormal the log transform reduces the
    # effective degree by 1, so max_order suffices. In both cases, if the
    # noise-to-signal ratio is low (<= 0.1) the observations are clean enough that
    # linear interpolation between them suffices.
    if (nrow(U) < control$max_points_interp && is.null(control$interpolation_method)) {
      nsr <- min(estimated_sd / apply(U, 2, sd))
      if (nsr <= 0.15) {
        control$interpolation_method <- "linear"
      } else {
        max_order     <- detect_max_state_order(f_orig_expr, u_expr)
        interp_degree <- if (noise_dist == "lognormal") max_order else max_order + 1L
        control$interpolation_method <- paste0("poly_ls_", interp_degree)
      }
    }

    if (is.null(control$interpolation_method) && nrow(U) >= control$max_points_interp) {
      wendy_data <- list(U = U, tt = tt, var = NULL)
    } else {
      if (length(control$interpolation_method) != 1L) {
        stop("control$interpolation_method must be a single method string.")
      }
      wendy_data <- interpolate_data(U, tt, control$interpolation_method, control, sigma = estimated_sd)
    }

    D <- ncol(wendy_data$U)

    wendy_problem <- build_wendy_problem(wendy_data, f_, J_u, J_up, J_p, J_pp, J_upp, J, lip, sig, control)

    system <- build_wendy_system(wendy_problem, lip, control$diag_reg, control$use_interp_uncertainty)

    p1 <- wendy_problem

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
    res$min_radius  <- p1$min_radius
    res$rc          <- p1$rc
    res$rc_errors   <- p1$rc_errors
    res$rc_radii    <- p1$rc_radii
    res$wendy_problem  <- wendy_problem
    res$wendy_data     <- wendy_data
    res$U              <- p1$U
    res$tt             <- wendy_data$tt
    res$var            <- p1$var
    res$W              <- system$W
  }

  opt_ctx <- list(
    system = if (method != "OE") system else NULL,
    method = method,
    f = f,
    U_orig = U_orig,
    tt_orig = tt_orig,
    U_processed  = U, # for lognormal: log-transformed + filtered; else same as U_orig
    tt_processed = as.vector(tt),
    # U/tt that match the V, V_prime grid (post-interpolation if it happened)
    U_sys = if (method != "OE") wendy_data$U else NULL,
    tt_sys = if (method != "OE") as.vector(wendy_data$tt) else NULL,
    lip = lip,
    f_orig_expr = f_orig_expr,
    u_expr = u_expr,
    estimated_sd = estimated_sd,
    f_ = f_,
    J_p = J_p,
    J_u   = if (method != "OE") J_u   else NULL,
    J_uu  = if (method != "OE") J_uu  else NULL,
    J_up  = if (method != "OE") J_up  else NULL,
    J_pp  = if (method != "OE") J_pp  else NULL,
    J_upp = if (method != "OE") J_upp else NULL,
    J_t = if (method != "OE") J_t else NULL,
    dF_dt_   = if (method != "OE") dF_dt_   else NULL,
    d2F_dt2_ = if (method != "OE") d2F_dt2_ else NULL,
    d3F_dt3_ = if (method != "OE") d3F_dt3_ else NULL,
    sig = if (method != "OE") sig else NULL,
    V = if (method != "OE") p1$V          else NULL,
    V_prime = if (method != "OE") p1$Vp         else NULL,
    min_radius = if (method != "OE") p1$min_radius else NULL,
    J = J,
    D = D,
    noise_dist  = noise_dist,
    control = control
  )

  class(res) <- "wendy"
  attr(res, "call") <- match.call()
  attr(res, "method") <- method
  attr(res, "noise_dist") <- noise_dist
  attr(res, "n_obs") <- length(tt)
  attr(res, "n_params") <- J
  attr(res, "n_states") <- D

  res$opt_ctx <- opt_ctx
  res$wsindy  <- ws  # NULL when f was supplied by the user

  if (!control$optimize) return(res)

  if (is.matrix(p0) && nrow(p0) > 1) {
    # Multistart: run optimization from each row of p0
    apply_fn <- control$apply_fn %||% lapply

    all_results <- apply_fn(seq_len(nrow(p0)), function(i) {
      .run_wendy_optimization(opt_ctx, p0[i, ])
    })

    objectives <- sapply(all_results, `[[`, "objective")
    best <- all_results[[which.min(objectives)]]

    res$phat <- best$p
    res$data <- best$data
    res$multistart_results <- all_results
    res$multistart_objectives <- objectives
  } else {
    # Single optimization
    p0_vec <- if (is.matrix(p0)) p0[1, ] else p0
    result  <- .run_wendy_optimization(opt_ctx, p0_vec)
    res$phat <- result$p
    res$data <- result$data
  }
  
  boundary_state <- if (control$estimate_IC && method != "OE") {
     if (identical(control$test_fun_type, "SSL")) {
      r_c_bl <- res$rc
    } else {
      r_c_bl <- compute_r_c_hat(U, tt, control$S, control$p)$rc
    }
    # n_bl deliberately omitted: estimate_IC uses its own max(3, ceiling(r_c/8))
    # heuristic. control$n_bl drives only the SSL boundary-layer augmentation.
    estimate_IC(U, f_, dF_dt_, d2F_dt2_, d3F_dt3_, tt, res$phat, r_c_bl,
                J_u = J_u, sigma = estimated_sd)
  } else NULL

  state <- if (control$estimate_trajectory && method != "OE") {
    if (identical(control$smoother, "erts") && !is.null(res$phat)) {
      # Parameter covariance for ERTS process-noise input: Fisher form
      # (Gᵀ S⁻¹ G)⁻¹.
      C_hat <- tryCatch({
        Sp <- res$S(res$phat)
        Gp <- if (lip) res$G else res$Jp_r(res$phat)
        R  <- chol(Sp)
        solve(crossprod(Gp, backsolve(R, forwardsolve(t(R), Gp))))
      }, error = function(e) NULL)
      # Only seed the filter from u0hat when u0 was actually optimised in
      # estimate_IC (cov_u0 non-NULL); otherwise u0hat is just the noisy obs.
      u0_init <- if (!is.null(boundary_state$cov_u0)) boundary_state$u0hat else NULL
      P0_init <- boundary_state$cov_u0
      wendy_erts(U, f_, J_u, tt, res$phat, control, sigma = estimated_sd,
                 u0_init   = u0_init,
                 P0_init   = P0_init,
                 param_cov = C_hat,
                 J_p       = J_p)
    } else {
      gp_smooth(U, tt, sigma2_n = as.numeric(sig)^2)
    }
  } else NULL

  res$boundary_state <- boundary_state
  res$u0hat <- boundary_state$u0hat
  res$state <- state

  return(res)
}
