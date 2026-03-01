#' @importFrom stats smooth.spline approx spline loess lm.fit ksmooth predict
#' @importFrom utils head tail
NULL

# Build a polynomial design matrix of the given degree for time vector tt.
# Returns a matrix with columns [1, tt, tt^2, ..., tt^degree].
# @keywords internal
poly_design <- function(tt, degree) {
  do.call(cbind, lapply(0:degree, function(k) tt^k))
}

# Parse the polynomial degree from a method string.
# Supports named forms:  "linear_ls"    -> 1
#                        "quadratic_ls" -> 2
#                        "cubic_ls"     -> 3
#                        "quartic_ls"   -> 4
#                        "quintic_ls"   -> 5
# and a numeric suffix:  "poly_ls_N"    -> N  (any positive integer)
# @keywords internal
poly_ls_degree <- function(method) {
  named <- c(linear = 1L, quadratic = 2L, cubic = 3L,
             quartic = 4L, quintic = 5L, sextic = 6L,
             septic = 7L, octic = 8L, nonic = 9L, decic = 10L)
  m <- regmatches(method, regexpr("^poly_ls_(\\d+)$", method))
  if (length(m) == 1L) return(as.integer(sub("^poly_ls_", "", m)))
  prefix <- sub("_ls$", "", method)
  if (prefix %in% names(named)) return(named[[prefix]])
  stop("Cannot parse polynomial degree from method: '", method, "'")
}

# Apply one interpolation method to map rows of U from tt_vec onto tt_target.
# Returns an interpolated matrix with nrow = length(tt_target), ncol = D.
# @keywords internal
interpolate_to_grid <- function(U, tt_vec, tt_target, method, substitute_data = TRUE) {
  if (grepl("_ls$", method)) {
    degree <- poly_ls_degree(method)
    X      <- poly_design(tt_vec,    degree)
    X_new  <- poly_design(tt_target, degree)
    result <- apply(U, 2, function(col) drop(X_new %*% lm.fit(X, col)$coefficients))
  } else {
    result <- switch(method,
      spline = {
        fits <- apply(U, 2, function(col) smooth.spline(tt_vec, col))
        sapply(fits, function(fit) predict(fit, tt_target)$y)
      },
      linear = apply(U, 2, function(col) approx(tt_vec, col, xout = tt_target)$y),
      cubic  = apply(U, 2, function(col) spline(tt_vec, col, xout = tt_target, method = "natural")$y),
      loess  = apply(U, 2, function(col) predict(loess(col ~ tt_vec), newdata = tt_target)),
      kernel = {
        bw <- diff(range(tt_vec)) / sqrt(length(tt_vec))
        apply(U, 2, function(col) ksmooth(tt_vec, col, kernel = "normal", bandwidth = bw, x.points = tt_target)$y)
      },
      stop("Unknown interpolation_method: ", method)
    )
  }

  # Ensure result is always a matrix (e.g. sapply drops dims for D=1)
  if (!is.matrix(result)) result <- matrix(result, nrow = length(tt_target), ncol = ncol(U))

  # Substitute original observed values at their exact time points
  if(substitute_data){
    tol <- sqrt(.Machine$double.eps)
    for (i in seq_along(tt_vec)) {
      j <- which(abs(tt_target - tt_vec[i]) < tol)
      if (length(j) == 1L) result[j, ] <- U[i, ]
    }
  }

  result
}

# Compute per-point variance for an interpolated grid.
# Returns a numeric vector of length(tt_target):
#   - 1  at positions matching an original observed time (tt_obs)
#   - 1 + leverage  at purely interpolated positions for LS methods
#   - 1 + mean_d( (se_d / residual_scale_d)^2 )  for loess
#   - 1  at all positions for other methods
# U_obs is required for loess (the original observed data matrix).
# @keywords internal
compute_interpolation_variance <- function(tt_obs, tt_target, method, U_obs = NULL) {
  tol        <- sqrt(.Machine$double.eps)
  is_obs     <- vapply(tt_target, function(t) any(abs(tt_obs - t) < tol), logical(1))
  var_vec    <- rep(1.0, length(tt_target))

  if (!any(!is_obs)) return(var_vec)

  tt_interp <- tt_target[!is_obs]

  if (grepl("_ls$", method)) {
    degree <- poly_ls_degree(method)
    X   <- poly_design(tt_obs,    degree)
    X_n <- poly_design(tt_interp, degree)
    XtXinv <- tryCatch(
      solve(crossprod(X)),
      error = function(e) {
        s   <- svd(X)
        thr <- max(s$d) * .Machine$double.eps * max(dim(X))
        d2  <- ifelse(s$d > thr, 1 / s$d^2, 0)
        s$v %*% diag(d2, nrow = length(d2)) %*% t(s$v)
      }
    )
    var_vec[!is_obs] <- 1 + rowSums((X_n %*% XtXinv) * X_n)

  } else if (method == "loess" && !is.null(U_obs)) {
    # predict.loess(se=TRUE) returns $se (SE of fitted value) and $residual.scale.
    # Normalized prediction variance per column: 1 + (se / residual.scale)^2.
    # Average across state dimensions for a single per-time-point scalar.
    # Note: sapply returns a length(tt_interp) x D matrix when length(tt_interp) > 1,
    # but a length-D vector when length(tt_interp) == 1.  Use is.matrix to distinguish.
    se_sq_cols <- sapply(seq_len(ncol(U_obs)), function(d) {
      fit  <- loess(U_obs[, d] ~ tt_obs)
      pred <- predict(fit, newdata = data.frame(tt_obs = tt_interp), se = TRUE)
      (pred$se / pred$residual.scale)^2
    })
    avg_se_sq <- if (is.matrix(se_sq_cols)) rowMeans(se_sq_cols) else mean(se_sq_cols)
    var_vec[!is_obs] <- 1 + avg_se_sq
  }

  var_vec
}

# Apply one named interpolation method to (U, tt), handling non-uniform spacing
# and minimum-number-of-points densification.  Returns list(U, tt, var).
# var is a numeric vector of per-time-point variances relative to the observed
# data noise (var = 1 at original observation times, >= 1 elsewhere).
# @keywords internal
interpolate_data <- function(U, tt, method, control) {
  tt_obs  <- as.vector(tt)          # original observed times (saved before any modification)
  U_obs   <- U                      # original observed data  (saved before any modification)
  diff_dt <- diff(tt_obs)
  dt      <- mean(diff_dt, na.rm = TRUE)

  if (max(abs(diff(tt_obs) - dt)) > sqrt(.Machine$double.eps)) {
    cat("Non uniform spacing detected, interpolating data using", method, "...\n")
    n      <- max(floor((max(tt) - min(tt)) / dt), control$min_number_points)
    tt_new <- seq(min(tt), max(tt), length.out = n)
    U      <- interpolate_to_grid(U, tt_obs, tt_new, method)
    tt     <- matrix(tt_new, ncol = 1)
  }

  if (nrow(U) < control$min_number_points) {
    tt_vec   <- as.vector(tt)
    tt_dense <- tt_vec
    while (2 * length(tt_dense) - 1 <= control$min_number_points * 2) {
      mids     <- (head(tt_dense, -1) + tail(tt_dense, -1)) / 2
      tt_dense <- sort(c(tt_dense, mids))
    }
    U  <- interpolate_to_grid(U, tt_vec, tt_dense, method)
    tt <- matrix(tt_dense, ncol = 1)
  }

  var <- compute_interpolation_variance(tt_obs, as.vector(tt), method, U_obs)
  list(U = U, tt = tt, var = var)
}
