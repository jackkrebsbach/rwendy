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
    param_cov_ols <- NULL   # only the GLS branch has a weak-form G to sandwich

    if (method %in% c("OE", "HYBRID")) {
      # OE / HYBRID: the final estimate comes from output_error, so use the
      # scaled nonlinear-LS covariance, restricted to the parameter block
      # (excluding fitted initial conditions).  For HYBRID, WENDy only
      # provides the starting point; the OE estimator's curvature is reported.
      J <- attr(object, "n_params")
      param_cov <- if (!is.null(object$data$cov)) {
        object$data$cov[1:J, 1:J]
      } else {
        NULL
      }
    } else {
      # GLS Fisher form (Gᵀ S(p̂)⁻¹ G)⁻¹ 
      G <- object$Jp_r(phat)
      Sp <- object$S(phat)
      R <- chol(Sp) # S = RᵀR
      Gt <- backsolve(R, G, transpose = TRUE) # R⁻ᵀG, whitened Jacobian
      param_cov <- chol2inv(chol(crossprod(Gt))) # (G̃ᵀG̃)⁻¹ = (GᵀS⁻¹G)⁻¹ 

      # OLS sandwich estimator G† S(p̂) (G^†)^T: G† = (GᵀG)^-1Gᵀ
      A <- crossprod(G) # GᵀG, J×J
      inner <- crossprod(G, Sp %*% G) # GᵀSG, J×J
      Ainv <- chol2inv(chol(A)) # SPD
      param_cov_ols <- Ainv %*% inner %*% Ainv
    }

    std_errors <- if (!is.null(param_cov)) sqrt(diag(param_cov)) else rep(NA, length(phat))

    summ$param_cov <- param_cov
    summ$param_cov_ols <- param_cov_ols
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
#' @param f_sym Symbolic expression vector from the active symbolic backend
#' @param format Output format: "latex" (default) or "ascii"
#' @return Formatted string representation of the ODE system
#' @keywords internal
format_ode_system <- function(f_sym, format = c("latex", "ascii")) {
  if(is.null(f_sym)) return(NULL)

  format <- match.arg(format)
  D <- sym_length(f_sym)

  if(format == "latex") {
    eqs <- sapply(1:D, function(i) {
      lhs <- paste0("\\frac{d u_", i, "}{dt}")
      rhs <- sym_latex(sym_elt(f_sym, i))
      paste0("  ", lhs, " &= ", rhs)
    })
    paste0("\\begin{align}\n", paste(eqs, collapse = " \\\\\n"), "\n\\end{align}")
  } else {
    rhs_all <- sym_strings(f_sym)
    eqs <- sapply(1:D, function(i) {
      lhs <- paste0("du", i, "/dt")
      paste0("  ", lhs, " = ", rhs_all[i])
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
  .stop_if_no_weak_form(object)
  as.numeric(object$g(object$phat) - object$b)
}

# method = "OE" never builds the weak form, so g/b/S are NULL and the callers
# below would die with "attempt to apply non-function".
.stop_if_no_weak_form <- function(object) {
  if (is.null(object$g) || is.null(object$b) || is.null(object$S)) {
    stop("No weak residual on this wendy object. It is not built for ",
         "method = \"OE\"; use summary() for the output-error fit instead.",
         call. = FALSE)
  }
}

#' Extract corrected residuals from a wendy object
#'
#' @param object A wendy object
#' @param ... Additional arguments (ignored)
#' @return Numeric vector of weak residuals
#' @export
residuals_weighted <- function(object, ...) {
  .stop_if_no_weak_form(object)
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

# Shared drawing for the radius-selection diagnostics: strip a log axis the data
# cannot carry, merge caller graphics args over the defaults, overlay any extra
# series, mark the selected radius and assemble the legend. `extra` is a list of
# list(y, label, col, pch, lty).
.plot_radius_curve <- function(x, y, selected, ix, log, defaults, dots,
                               series_label, selected_label,
                               extra = list(), headroom = 0) {
  ys <- c(y, unlist(lapply(extra, `[[`, "y")))
  if (grepl("y", log, fixed = TRUE) && any(!is.finite(ys) | ys <= 0)) {
    log <- gsub("y", "", log, fixed = TRUE)
  }
  logy <- grepl("y", log, fixed = TRUE)

  if (headroom > 0) {
    yl <- range(ys[is.finite(ys)])
    yl[2] <- if (logy) yl[2] * (yl[2] / yl[1])^headroom else yl[2] + headroom * diff(yl)
    defaults$ylim <- yl
  }

  plot_args <- utils::modifyList(
    c(list(x = x, y = y, type = "b", pch = 16, col = "#1f77b4", log = log), defaults),
    dots)
  do.call(graphics::plot, plot_args)

  for (e in extra) {
    graphics::lines(x, e$y, type = "b", pch = e$pch, lty = e$lty, col = e$col)
  }

  graphics::abline(v = selected, col = "#d62728", lty = 2, lwd = 1.5)
  graphics::points(x[ix], y[ix], col = "#d62728", pch = 19, cex = 1.7)

  graphics::legend(
    "topright",
    legend = c(series_label,
               vapply(extra, `[[`, "", "label"),
               sprintf(selected_label, selected)),
    col    = c(plot_args$col, vapply(extra, `[[`, "", "col"), "#d62728"),
    pch    = c(plot_args$pch, vapply(extra, function(e) e$pch, 0), 19),
    lty    = c(1, vapply(extra, function(e) e$lty, 0), 2),
    bty    = "n"
  )
  invisible(plot_args)
}

#' Plot the test-function radius selection diagnostic
#'
#' Visualises the data-driven choice of test-function support radius made when a
#' WENDy problem is built. WENDy sweeps a range of radii, evaluates an
#' integration-error proxy at each, and selects the change-point ("elbow") of
#' that error curve. This function draws the error curve against the radius and
#' marks the selected radius.
#'
#' Two selection routines are covered, dispatched automatically from what the
#' fit stored:
#' \itemize{
#'   \item \strong{MSG} (\code{test_fun_type = "MSG"}): the minimum radius from
#'     \code{find_min_radius_int_error} (using the \eqn{\Phi} bump test
#'     function) -- stored in \code{object$min_radius} with the error curve in
#'     \code{object$wendy_problem$min_radius_errors} /
#'     \code{min_radius_radii}.
#'   \item \strong{SSL} (\code{test_fun_type = "SSL"}): the critical radius
#'     \eqn{r_c} from \code{compute_r_c_hat} (using the \eqn{\Psi}
#'     piecewise-polynomial test function) -- stored in \code{object$rc} with
#'     the error curve in \code{object$rc_errors} / \code{object$rc_radii}.
#' }
#'
#' Radius selection is skipped (and nothing is plotted) for \code{method = "OE"}
#' and whenever \code{control$fixed_radius} is set.
#'
#' @param object A \code{wendy} object returned by \code{\link{solveWendy}}.
#' @param log Character passed to \code{\link[graphics]{plot}} controlling which
#'   axes use a log scale. Defaults to \code{"y"} because the error decays over
#'   several orders of magnitude; automatically dropped if any error is
#'   non-positive.
#' @param ... Additional graphical parameters forwarded to
#'   \code{\link[graphics]{plot}} (e.g. \code{main}, \code{xlab}, \code{col});
#'   any supplied value overrides the default.
#' @return Invisibly, a list with the swept \code{radii}, the \code{errors}
#'   curve, the selected \code{radius}, its index \code{ix} into \code{radii},
#'   and the selection \code{type} (\code{"MSG"} or \code{"SSL"}).
#' @export
plot_radius_selection <- function(object, log = "y", ...) {
  if (!inherits(object, "wendy")) {
    stop("`object` must be a wendy object returned by solveWendy().")
  }

  wp <- object$wendy_problem

  # SSL r_c change-point (psi test function) vs. MSG minimum radius (phi/psi).
  # The two paths populate disjoint slots, so presence of an error curve tells
  # us which selection ran. Prefer the top-level fields, falling back to the
  # wendy_problem the fit carries.
  ssl_errors <- object$rc_errors %||% wp$rc_errors
  ssl_radii  <- object$rc_radii  %||% wp$rc_radii
  msg_errors <- object$min_radius_errors %||% wp$min_radius_errors
  msg_radii  <- object$min_radius_radii  %||% wp$min_radius_radii

  if (!is.null(ssl_errors) && !is.null(ssl_radii)) {
    type     <- "SSL"
    radii    <- as.numeric(ssl_radii)
    errors   <- as.numeric(ssl_errors)
    selected <- object$rc %||% wp$rc
    main_def <- expression("SSL radius selection (" * Psi * " test functions, " * r[c] * ")")
  } else if (!is.null(msg_errors) && !is.null(msg_radii)) {
    type     <- "MSG"
    radii    <- as.numeric(msg_radii)
    errors   <- as.numeric(msg_errors)
    selected <- object$min_radius %||% wp$min_radius
    main_def <- expression("MSG minimum-radius selection (" * Phi * " test functions)")
  } else {
    stop("No radius-selection diagnostics found on this wendy object. ",
         "Radius selection is skipped for method = \"OE\" and when ",
         "control$fixed_radius is set.")
  }

  if (is.null(selected) || length(selected) == 0L) {
    stop("Selected radius is missing from the wendy object.")
  }

  # The selected radius is one of the swept radii; recover its index.
  ix <- which.min(abs(radii - selected))

  .plot_radius_curve(
    radii, errors, selected, ix, log,
    defaults = list(xlab = "Test-function support radius (grid points)",
                    ylab = "Integration error",
                    main = main_def),
    dots = list(...),
    series_label   = "Integration error",
    selected_label = "Selected radius = %g")

  invisible(list(
    type   = type,
    radii  = radii,
    errors = errors,
    radius = selected,
    ix     = ix
  ))
}

#' Plot the boundary-layer radius selection from estimate_IC
#'
#' Draws the a-priori \eqn{u_0} MSE proxy over the swept \code{r_bl} grid and
#' marks the selected radius. With \code{show_terms}, the variance part
#' \eqn{\sum_d \mathrm{Var}_d/\sigma_d^2} and the bias part
#' \eqn{\sum_d (\Delta^{EM}_d + b_d)^2/\sigma_d^2} are overlaid.
#'
#' @param object A \code{wendy} object from \code{\link{solveWendy}}, or the list
#'   returned by \code{\link{estimate_IC}}.
#' @param log Character passed to \code{\link[graphics]{plot}}; dropped if any
#'   plotted value is non-positive.
#' @param show_terms Logical; overlay the variance and bias parts. Off by
#'   default: the bias part runs orders of magnitude below the objective, so
#'   including it in the y-range flattens the minimum out of view.
#' @param ... Additional graphical parameters forwarded to
#'   \code{\link[graphics]{plot}}.
#' @return Invisibly, a list with the \code{design} table, the selected
#'   \code{r_bl}, and its index \code{ix}.
#' @export
plot_IC_radius_selection <- function(object, log = "y", show_terms = FALSE, ...) {
  bs       <- object$boundary_state %||% object
  design   <- bs$design
  selected <- bs$r_bl

  if (is.null(design) || is.null(selected)) {
    stop("No IC design-selection table found. estimate_IC only sweeps r_bl when ",
         "inverse = \"gls\" and n_bl is NULL.")
  }

  d   <- design[is.finite(design$obj), , drop = FALSE]
  var <- d$var_obj
  bia <- d$obj - var
  ix  <- which.min(abs(d$r_bl - selected))

  extra <- if (show_terms) list(
    list(y = var, label = "variance part", col = "#2ca02c", pch = 1, lty = 3),
    list(y = bia, label = "bias part",     col = "#9467bd", pch = 2, lty = 3)
  ) else list()

  # headroom: the objective is U-shaped in r_bl, so both upper corners are
  # occupied and the legend would otherwise land on the curve.
  .plot_radius_curve(
    d$r_bl, d$obj, selected, ix, log,
    defaults = list(xlab = "Boundary-layer radius",
                    ylab = expression("MSE proxy"),
                    main = expression("IC boundary-layer radius selection")),
    dots = list(...),
    series_label   = "MSE proxy (obj)",
    selected_label = "Selected r_bl = %g",
    extra = extra, headroom = 0.28)

  invisible(list(design = d, r_bl = selected, ix = ix))
}
