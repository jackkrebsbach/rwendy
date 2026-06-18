# ODE integration backend.
#
# WENDy's weak-form estimation never integrates an ODE; a numerical solver is
# only needed for the forward-simulation paths (the Output-Error method and
# simulating example data). We prefer deSolve when it is installed because its
# lsoda auto-switching is the most robust choice (especially on stiff systems),
# but deSolve ships compiled Fortran and so cannot run under webR. To keep the
# package installable everywhere including webR 
.wendy_state <- new.env(parent = emptyenv())

.warn_pracma_fallback <- function() {
  if (!isTRUE(.wendy_state$warned_pracma)) {
    warning(
      "wendy: 'deSolve' is not installed; using the pure-R 'pracma' fallback ",
      "for ODE integration (less robust on stiff systems). ",
      "Install it with install.packages('deSolve') for the recommended solver.",
      call. = FALSE
    )
    .wendy_state$warned_pracma <- TRUE
  }
}

wendy_ode <- function(func, times, y0, parms) {
  times <- as.numeric(times)
  n <- length(times)
  D <- length(y0)

  # deSolve (robust, auto-stiff)
  if (requireNamespace("deSolve", quietly = TRUE)) {
    sol <- tryCatch({
      utils::capture.output(
        res <- suppressWarnings(
          deSolve::ode(y = y0, times = times, func = func, parms = parms)
        ),
        type = "output"
      )
      res
    }, error = function(e) NULL)

    if (is.null(sol) || nrow(sol) != n) return(NULL)
    # deSolve returns cbind(time, states), drop the time column.
    out <- unname(sol[, -1, drop = FALSE])
    if (anyNA(out) || any(!is.finite(out))) return(NULL)
    return(out)
  }

  # pure-R fallback: pracma (works in webR)
  .warn_pracma_fallback()

  # pracma::ode45 takes f(t, y) and integrates over a scalar span on its own
  # adaptive grid, so we strip the parms argument via a closure and interpolate
  # the result back onto the requested observation times.
  deriv <- function(t, y) as.numeric(func(t, y, parms)[[1]])
  sol <- tryCatch(
    pracma::ode45(deriv, t0 = times[1], tfinal = times[n], y0 = y0),
    error = function(e) NULL
  )
  if (is.null(sol) || length(sol$t) < 2) return(NULL)

  ys <- sol$y
  if (is.null(dim(ys))) ys <- matrix(ys, ncol = D)  # D == 1 returns a vector

  out <- matrix(NA_real_, nrow = n, ncol = D)
  for (d in seq_len(D)) {
    out[, d] <- stats::spline(sol$t, ys[, d], xout = times)$y
  }
  if (anyNA(out) || any(!is.finite(out))) return(NULL)
  out
}

# Startup recommendation to install deSolve
# Skipped on webR (R.version$os == "emscripten"), deSolve does not work
# the pracma fallback is the expected path.
.onAttach <- function(libname, pkgname) {
  if (!identical(R.version$os, "emscripten") &&
      !requireNamespace("deSolve", quietly = TRUE)) {
    packageStartupMessage(
      "wendy: 'deSolve' is not installed. Core weak-form estimation works ",
      "without it, but install.packages('deSolve') is recommended for robust ",
      "ODE simulation (used by the Output-Error method and data simulation)."
    )
  }
}