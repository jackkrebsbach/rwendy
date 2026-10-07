# Bounds are always specified in physical parameter/state units. A bounded
# state block uses a diagonal coordinate map, so box-constrained optimizers
# enforce the requested bounds rather than bounding correlated GP z values.
resolve_parameter_bound_alias <- function(value, alias, name, alias_name) {
  if (!is.null(value) && !is.null(alias)) {
    stop("Supply only ", name, "; ", alias_name, " is its compatibility alias.")
  }
  if (is.null(value)) alias else value
}

validate_joint_gp_bounds <- function(lower, upper, size, label) {
  expand <- function(x, default, side) {
    if (is.null(x)) {
      return(rep(default, size))
    }
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x) %in% c(1L, size) || anyNA(x)) {
      stop(label, " ", side, " must be a numeric scalar or a vector of length ", size, " without NA.")
    }
    rep_len(x, size)
  }
  lower_bounds <- expand(lower, -Inf, "lower bound")
  upper_bounds <- expand(upper, Inf, "upper bound")
  if (any(lower_bounds >= upper_bounds)) stop("Each ", label, " lower bound must be smaller than its upper bound.")
  list(lower = lower_bounds, upper = upper_bounds, bounded = is.finite(lower_bounds) | is.finite(upper_bounds))
}

apply_state_bounds_to_coordinates <- function(coordinates, bounds, scales = NULL) {
  if (is.null(bounds)) {
    return(coordinates)
  }
  if (is.null(scales)) scales <- rep(1, ncol(coordinates$mu))
  for (component in which(bounds$bounded)) {
    coordinates$mu[, component] <- 0
    coordinates$L[[component]] <- scales[component]
  }
  coordinates
}

build_optimizer_bounds <- function(coordinates, lower, upper, state_bounds = NULL) {
  n_components <- ncol(coordinates$mu)
  n_grid <- nrow(coordinates$mu)
  state_lower <- rep(-Inf, n_grid * n_components)
  state_upper <- rep(Inf, n_grid * n_components)
  if (!is.null(state_bounds)) {
    for (component in which(state_bounds$bounded)) {
      scale <- coordinates$L[[component]]
      if (length(scale) != 1L || scale <= 0) stop("Bounded states require a positive diagonal coordinate map.")
      indices <- (match(component, coordinates$order) - 1L) * n_grid + seq_len(n_grid)
      state_lower[indices] <- (state_bounds$lower[component] - coordinates$mu[, component]) / scale
      state_upper[indices] <- (state_bounds$upper[component] - coordinates$mu[, component]) / scale
    }
  }
  list(lower = c(lower / coordinates$pscale, state_lower), upper = c(upper / coordinates$pscale, state_upper))
}

pack_joint_gp_coordinates <- function(p, U, coordinates) {
  c(p / coordinates$pscale, unlist(lapply(coordinates$order, function(component) {
    coordinate_map <- coordinates$L[[component]]
    centered <- U[, component] - coordinates$mu[, component]
    if (length(coordinate_map) == 1L) centered / coordinate_map else forwardsolve(coordinate_map, centered)
  })))
}

unpack_joint_gp_coordinates <- function(theta, coordinates) {
  n_parameters <- length(coordinates$pscale)
  n_grid <- nrow(coordinates$mu)
  U <- coordinates$mu
  for (k in seq_along(coordinates$order)) {
    component <- coordinates$order[k]
    state_coordinates <- theta[n_parameters + (k - 1L) * n_grid + seq_len(n_grid)]
    coordinate_map <- coordinates$L[[component]]
    increment <- if (length(coordinate_map) == 1L) {
      coordinate_map * state_coordinates
    } else {
      as.vector(coordinate_map %*% state_coordinates)
    }
    U[, component] <- U[, component] + increment
  }
  list(p = theta[seq_len(n_parameters)] * coordinates$pscale, U = U)
}

# Project initial states into the requested bounds before fitting.
clip_state_to_bounds <- function(U, bounds) {
  if (is.null(bounds)) {
    return(U)
  }
  for (component in which(bounds$bounded)) {
    U[, component] <- pmax(bounds$lower[component], pmin(bounds$upper[component], U[, component]))
  }
  U
}
