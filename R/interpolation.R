#' @importFrom stats smooth.spline approx spline loess lm predict ksmooth
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
  named <- c(linear = 1L, quadratic = 2L, cubic = 3L, quartic = 4L, quintic = 5L, sextic = 6L, septic = 7L, octic = 8L, nonic = 9L, decic = 10L)
  m <- regmatches(method, regexpr("^poly_ls_(\\d+)$", method))
  if (length(m) == 1L) {
    return(as.integer(sub("^poly_ls_", "", m)))
  }
  prefix <- sub("_ls$", "", method)
  if (prefix %in% names(named)) return(named[[prefix]])
  stop("Cannot parse polynomial degree from method: '", method, "'")
}

# Fit one interpolation model for a single state dimension (column).
# Trains on (tt_obs, y) and predicts at tt_target.
# Returns list(fit, var):
#   fit — numeric vector of length(tt_target), predicted values.
#   var — numeric vector of length(tt_target), prediction-interval variance:
#         - model-based methods (*_ls, loess, spline): se.fit^2 + residual.scale^2
#         - linear: sigma^2 * (1 + (1-t)^2 + t^2) when sigma is supplied, else NA
#         - cubic, kernel: NA (no residual model available)
# sigma: optional scalar noise SD used for exact interpolants that cannot
#        estimate their own residual scale.
# @keywords internal
fit_col <- function(y, tt_obs, tt_target, method, sigma = NULL) {

  if (grepl("_ls$", method)) {
    degree <- poly_ls_degree(method)
    X      <- poly_design(tt_obs,    degree)
    df_obs <- as.data.frame(X[, -1L, drop = FALSE])
    colnames(df_obs) <- paste0("V", seq_len(ncol(df_obs)))
    X_new  <- poly_design(tt_target, degree)
    df_new <- as.data.frame(X_new[, -1L, drop = FALSE])
    colnames(df_new) <- colnames(df_obs)

    model  <- lm(y ~ ., data = cbind(y = y, df_obs))
    pred   <- predict(model, newdata = df_new, se.fit = TRUE, interval = "prediction")
    pred_var <- pred$se.fit^2 #+ pred$residual.scale^2

    list(fit = pred$fit[, "fit"], var = pred_var)

  } else if (method == "spline") {
    fit    <- smooth.spline(tt_obs, y)
    yhat   <- fit$y          
    sigma2 <- sum((y - yhat)^2) / max(length(y) - fit$df, 1)
    list(fit = predict(fit, tt_target)$y, var = rep(sigma2, length(tt_target)))

  } else if (method == "loess") {
    fit  <- loess(y ~ tt_obs)
    pred <- predict(fit, newdata = data.frame(tt_obs = tt_target), se = TRUE)
    list(fit = pred$fit, var = pred$se^2 + pred$residual.scale^2)

  } else if (method == "linear") {
    fit <- approx(tt_obs, y, xout = tt_target)$y
    var <- if (!is.null(sigma)) {
      t_vals <- vapply(tt_target, function(x) {
        i <- max(1L, min(findInterval(x, tt_obs, rightmost.closed = TRUE), length(tt_obs) - 1L))
        (x - tt_obs[i]) / (tt_obs[i + 1L] - tt_obs[i])
      }, numeric(1L))
      sigma^2 * (1 + (1 - t_vals)^2 + t_vals^2)
    } else {
      rep(NA_real_, length(tt_target))
    }
    list(fit = fit, var = var)
  } else if (method == "cubic") {
    list(fit = spline(tt_obs, y, xout = tt_target, method = "natural")$y, var = rep(NA_real_, length(tt_target)))
  } else if (method == "kernel") {
    bw <- diff(range(tt_obs)) / sqrt(length(tt_obs))
    list(fit = ksmooth(tt_obs, y, kernel = "normal", bandwidth = bw, x.points = tt_target)$y,
         var = rep(NA_real_, length(tt_target)))
  } else {
    stop("Unknown interpolation_method: ", method)
  }
}

# Apply one interpolation method across all D columns of U.
# Returns list(U, var):
#   U   — interpolated matrix, nrow = length(tt_target), ncol = D.
#   var — numeric vector of length(tt_target), prediction-interval variance
#         averaged across state dimensions.  Falls back to 1 for exact
#         interpolants that carry no residual model and sigma is not supplied.
# sigma: optional scalar noise SD forwarded to fit_col for exact interpolants.
# @keywords internal
interpolate_to_grid <- function(U, tt_vec, tt_target, method, substitute_data = TRUE, sigma = NULL) {
  cols <- lapply(seq_len(ncol(U)), function(d) {
    fit_col(U[, d], tt_vec, tt_target, method, sigma = sigma)
  })
  U_new   <- do.call(cbind, lapply(cols, `[[`, "fit"))
  var_mat <- do.call(cbind, lapply(cols, `[[`, "var"))

  if (all(is.na(var_mat))) {
    var_vec <- rep(1.0, length(tt_target))
  } else {
    var_vec <- rowMeans(var_mat, na.rm = TRUE)
  }

  if (!is.matrix(U_new)){
    U_new <- matrix(U_new, nrow = length(tt_target), ncol = ncol(U))
  }

  if (substitute_data) {
    tol <- sqrt(.Machine$double.eps)
    for (i in seq_along(tt_vec)) {
      j <- which(abs(tt_target - tt_vec[i]) < tol)
      if (length(j) == 1L) U_new[j, ] <- U[i, ]
    }
  }

  list(U = U_new, var = var_vec)
}

build_scale_vec <- function(tt_target, tt_obs, var_vec) {
  tol       <- sqrt(.Machine$double.eps)
  scale_vec <- 1 + var_vec
  for (i in seq_along(tt_obs)) {
    j <- which(abs(tt_target - tt_obs[i]) < tol)
    if (length(j) == 1L) scale_vec[j] <- 1
  }
  scale_vec
}

interpolate_data <- function(U, tt, method, control, sigma = NULL) {
  tt_obs  <- as.vector(tt)
  U_obs   <- U

  diff_dt <- diff(tt_obs)
  dt      <- mean(diff_dt, na.rm = TRUE)

  if (max(abs(diff(tt_obs) - dt)) > sqrt(.Machine$double.eps)) {
    cat("Non uniform spacing detected, interpolating data using", method, "...\n")
    n      <- max(floor((max(tt) - min(tt)) / dt), control$min_number_points)
    tt_new <- seq(min(tt), max(tt), length.out = n)
    U      <- interpolate_to_grid(U, tt_obs, tt_new, method, sigma = sigma)$U
    tt     <- matrix(tt_new, ncol = 1)
  }

  if (nrow(U) < control$min_number_points) {
    tt_vec   <- as.vector(tt)
    tt_dense <- tt_vec
    while (2 * length(tt_dense) - 1 <= control$min_number_points * 2) {
      mids     <- (head(tt_dense, -1) + tail(tt_dense, -1)) / 2
      tt_dense <- sort(c(tt_dense, mids))
    }
    U  <- interpolate_to_grid(U, tt_vec, tt_dense, method, sigma = sigma)$U
    tt <- matrix(tt_dense, ncol = 1)
  }

  result    <- interpolate_to_grid(U_obs, tt_obs, as.vector(tt), method, substitute_data = FALSE, sigma = sigma)
  scale_vec <- build_scale_vec(as.vector(tt), tt_obs, result$var)
  list(U = U, tt = tt, var = result$var, scale = scale_vec)
}
