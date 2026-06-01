# Null-coalescing operator: return lhs unless it is NULL, then return rhs
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

#' Print method for wendy objects
#'
#' @param x A wendy object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @method print wendy
#' @export
print.wendy <- function(x, ...) {
  cat("WENDy Model Fit\n")
  cat("===============\n\n")
  cat("Method:", attr(x, "method"), "\n")
  cat("Noise distribution:", attr(x, "noise_dist"), "\n")
  cat("Number of observations:", attr(x, "n_obs"), "\n")
  cat("Number of parameters:", attr(x, "n_params"), "\n")
  cat("Number of state variables:", attr(x, "n_states"), "\n\n")
  
  if(!is.null(x$phat)) {
    cat("Estimated parameters:\n")
    print(x$phat)
  }
  
  invisible(x)
}

#' Summary method for wendy objects
#'
#' Computes summary statistics for a fitted WENDy model including
#' the ODE system, parameter estimates, and residual diagnostics.
#'
#' @param object A wendy object returned by \code{\link{solveWendy}}
#' @param ... Additional arguments (ignored)
#' @return An object of class \code{summary.wendy} containing:
#'   \item{call}{The original function call}
#'   \item{method}{Optimization method used}
#'   \item{noise_dist}{Noise distribution assumed}
#'   \item{n_obs}{Number of observations}
#'   \item{n_params}{Number of parameters}
#'   \item{n_states}{Number of state variables}
#'   \item{ode_system}{Formatted ODE system as a string}
#'   \item{parameters}{Data frame of parameter estimates}
#'   \item{residuals}{Summary statistics of weak residuals}
#' @method summary wendy
#' @export
summary.wendy <- function(object, ...) {

  summ <- list()
  summ$call <- attr(object, "call")
  summ$method <- attr(object, "method")
  summ$noise_dist <- attr(object, "noise_dist")
  summ$n_obs <- attr(object, "n_obs")
  summ$n_params <- attr(object, "n_params")
  summ$n_states <- attr(object, "n_states")

  summ$ode_system <- format_ode_system(object$f_sym)

  method <- attr(object, "method")

  if(!is.null(object$phat)) {
    phat <- object$phat

    if (method %in% c("OE", "HYBRID")) {
      # OE / HYBRID: the final estimate comes from output_error, so use the
      # scaled covariance from modFit, restricted to the parameter block
      # (excluding fitted initial conditions).  For HYBRID, WENDy only
      # provides the starting point; the OE estimator's curvature is reported.
      J <- attr(object, "n_params")
      param_cov <- if (!is.null(object$data$cov)) object$data$cov[1:J, 1:J] else NULL
    } else {
      # Fisher form (Gᵀ S(p̂)⁻¹ G)⁻¹: inverse Fisher information of the
      # linear-Gaussian residual model and the asymptotic cov of the IRLS
      # estimator.
      G  <- object$Jp_r(phat)
      Sp <- object$S(phat)
      R  <- chol(Sp)
      param_cov <- solve(crossprod(G, backsolve(R, forwardsolve(t(R), G))))
    }

    std_errors <- if (!is.null(param_cov)) sqrt(diag(param_cov)) else rep(NA, length(phat))

    summ$param_cov <- param_cov
    summ$phat <- phat
    summ$parameters <- data.frame(
      Parameter = paste0("p", seq_along(phat)),
      Estimate  = as.numeric(phat),
      Std.Error = std_errors
    )
  }

  if(!is.null(object$data)) {
    if (method %in% c("OE", "HYBRID")) {
      summ$convergence <- object$data$converged
      summ$iterations  <- object$data$iterations
    } else {
      summ$wnll_value <- as.numeric(object$wnll(object$phat))

      r <- object$g(object$phat) - object$b
      summ$residuals <- list(
        min    = min(as.numeric(r)),
        q1     = quantile(as.numeric(r), 0.25),
        median = median(as.numeric(r)),
        q3     = quantile(as.numeric(r), 0.75),
        max    = max(as.numeric(r))
      )

      if (!is.null(object$data$convergence)) {
        summ$convergence <- object$data$convergence
      }
      if (!is.null(object$data$iterations)) {
        summ$iterations <- object$data$iterations
      }
    }
  }

  if(!is.null(object$min_radius)) {
    summ$min_radius <- object$min_radius
  }

  class(summ) <- "summary.wendy"
  return(summ)
}

#' Print method for summary.wendy objects
#'
#' @param x A summary.wendy object
#' @param digits Number of significant digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @method print summary.wendy
#' @export
print.summary.wendy <- function(x, digits = 4, ...) {
  cat("WENDy Model Summary\n")
  cat("===================\n\n")
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat("Method:", x$method, "\n")
  cat("Noise distribution:", x$noise_dist, "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("State variables:", x$n_states, "\n")
  cat("Parameters:", x$n_params, "\n\n")
  
  # ODE System
  if(!is.null(x$ode_system)) {
    cat("ODE System:\n")
    cat(x$ode_system, "\n\n")
  }
  
  # Parameter estimates
  if(!is.null(x$parameters)) {
    cat("Parameter Estimates:\n")
    print(x$parameters, digits = digits, row.names = FALSE)
    cat("\n")
  }

  # Parameter covariance
  if(!is.null(x$param_cov)) {
    cat("Parameter Covariance:\n")
    J <- nrow(x$param_cov)
    rownames(x$param_cov) <- paste0("p", seq_len(J))
    colnames(x$param_cov) <- paste0("p", seq_len(J))
    print(round(x$param_cov, digits))
    cat("\n")
  }

  # Residuals
  if(!is.null(x$residuals)) {
    cat("Weak Residuals:\n")
    cat(sprintf("   Min: %.*f\n", digits, x$residuals$min))
    cat(sprintf("    Q1: %.*f\n", digits, x$residuals$q1))
    cat(sprintf("Median: %.*f\n", digits, x$residuals$median))
    cat(sprintf("    Q3: %.*f\n", digits, x$residuals$q3))
    cat(sprintf("   Max: %.*f\n", digits, x$residuals$max))
    cat("\n")
  }
  
  # Fit statistics
  if(!is.null(x$wnll_value)) {
    cat(sprintf("Weak Negative Log-Likelihood: %.*f\n", digits, x$wnll_value))
  }
  
  # Convergence
  if(!is.null(x$convergence)) {
    conv_str <- if (is.logical(x$convergence)) {
      ifelse(isTRUE(x$convergence), "successful", "failed")
    } else {
      ifelse(x$convergence == 0, "successful", paste("failed (code:", x$convergence, ")"))
    }
    cat("\nConvergence:", conv_str, "\n")
  }
  if(!is.null(x$iterations)) {
    cat("Iterations:", x$iterations, "\n")
  }
  
  invisible(x)
}

#' Format ODE system as human-readable string
#'
#' @param f_sym Symbolic expression vector from symengine
#' @param format Output format: "latex" (default) or "ascii"
#' @return Formatted string representation of the ODE system
#' @keywords internal
format_ode_system <- function(f_sym, format = c("latex", "ascii")) {
  if(is.null(f_sym)) return(NULL)

  format <- match.arg(format)
  D <- length(f_sym)

  if(format == "latex") {
    eqs <- sapply(1:D, function(i) {
      lhs <- paste0("\\frac{d u_", i, "}{dt}")
      rhs <- symengine::codegen(f_sym[i], "latex")
      paste0("  ", lhs, " &= ", rhs)
    })
    paste0("\\begin{align}\n", paste(eqs, collapse = " \\\\\n"), "\n\\end{align}")
  } else {
    eqs <- sapply(1:D, function(i) {
      lhs <- paste0("du", i, "/dt")
      rhs <- as.character(f_sym[i])
      paste0("  ", lhs, " = ", rhs)
    })
    paste(eqs, collapse = "\n")
  }
}
  
#' Extract coefficients from a wendy object
#'
#' @param object A wendy object
#' @param ... Additional arguments (ignored)
#' @return Numeric vector of estimated parameters
#' @method coef wendy
#' @export
coef.wendy <- function(object, ...) {
  object$phat
}

#' Extract residuals from a wendy object
#'
#' @param object A wendy object
#' @param ... Additional arguments (ignored)
#' @return Numeric vector of weak residuals
#' @method residuals wendy
#' @export
residuals.wendy <- function(object, ...) {
  as.numeric(object$g(object$phat) - object$b)
}

#' Extract corrected residuals from a wendy object
#'
#' @param object A wendy object
#' @param ... Additional arguments (ignored)
#' @return Numeric vector of weak residuals
#' @export
residuals_weighted <- function(object, ...) {
  S_mat <- object$S(object$phat)
  r     <- object$b - object$g(object$phat)
  RT    <- t(chol(S_mat))
  as.numeric(forwardsolve(RT, r))
}

#' Calculate the relative error of two vectors
#'
#' @param x Vector
#' @param y Vector
#' @return Scaler 
#' @export
rel_err <- function(x,y){
  norm(x - y, type = "2") / norm(y, type = "2")
}
