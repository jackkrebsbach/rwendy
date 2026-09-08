# Latent-component estimation with a fixed weak operator and separate GP fit.
#
# Exactly one latent component. Stage 1 carries the latent grid values
# in PHYSICAL coordinates with a scalar preconditioner and no latent prior
# penalty, which is what keeps the joint MAP over state and amplitude from being
# unbounded below. Stage 2 fits a GP to the stage-1 curve. Stage 3 restarts from
# the stage-1 curve, using that GP only as the PRECONDITIONER -- its penalty is
# off by default. An invertible preconditioner changes neither the
# physical objective nor the set of admissible grid states. See latent_penalty.

.jgl_split <- function(U) {
  full <- colSums(is.na(U)) == 0L
  empty <- colSums(!is.na(U)) == 0L
  if (any(!full & !empty))
    stop("Each component column must be fully observed or entirely NA.")
  if (!any(full)) stop("At least one component must be observed.")
  if (sum(empty) != 1L) stop("Exactly one latent component is supported.")
  if (any(!is.finite(U[,full,drop=FALSE]))) stop("Observed components must be finite.")
  list(observed = which(full), latent = which(empty))
}

# The deviation from the study pipeline: latent start is the row mean of the
# observed columns, carried to the grid.
.jgl_initial <- function(U, tt, grid, obs, lat, latent_initial = NULL) {
  init <- matrix(0, length(grid), ncol(U))
  for (d in obs) init[, d] <- stats::approx(tt, U[, d], xout = grid, rule = 2)$y
  curve <- rowMeans(init[, obs, drop = FALSE])
  observed_rms <- sqrt(colMeans(init[,obs,drop=FALSE]^2))
  if (sqrt(mean(curve^2)) <= sqrt(.Machine$double.eps)*max(observed_rms))
    curve <- init[,obs[which.max(observed_rms)]]
  if (!is.null(latent_initial))
    curve <- stats::approx(tt,latent_initial,xout=grid,rule=2)$y
  for (d in lat) init[, d] <- curve
  init
}

# GP posterior mean and its time derivative on a grid. The derivative is a
# central difference of a SMOOTH deterministic curve -- the data noise is
# already absorbed into alpha -- measured accurate to 1.4e-5 relative, three
# orders below the noise floor it feeds.
.jgl_gp_mean <- function(g, tt)
  as.vector(g$mean + .jgp_gp_cross(g, tt, g$tt) %*% g$alpha)

.jgl_gp_deriv <- function(g, tt) {
  h <- g$span * 1e-4
  (.jgl_gp_mean(g, tt + h) - .jgl_gp_mean(g, tt - h)) / (2 * h)
}

# NOT USED by default -- .jgl_initial (row mean) is the initializer. Retained
# because it is the only p-informed option available: it needs a parameter
# vector, and its accuracy degrades roughly linearly in the error of that
# vector, crossing the row mean at about a factor of two off. With p0 defaulting
# to ones there is no reason to trust it, so it is not called.
# Generic latent initializer: at every grid point solve the OBSERVED rows of
# f(U,p,t) = du/dt for the latent value by least squares, with du/dt from the GP
# posterior mean. Newton on the symbolic Jacobian, so nonlinearity in the latent
# needs no per-system algebra. Where the observed rows do not determine the
# latent -- Lorenz z is multiplied by x, so it is unidentified wherever x ~ 0 --
# those points are filled by interpolating the ones that are determined.
.jgl_invert <- function(prep, p, steps = 5L, tol = .05) {
  obs <- prep$observed; d <- prep$latent[1]; J <- prep$J; M <- length(prep$grid)
  U <- prep$pilot
  for (k in seq_along(obs)) U[, obs[k]] <- .jgl_gp_mean(prep$gps[[k]], prep$grid)
  dmu <- vapply(prep$gps, .jgl_gp_deriv, numeric(M), tt = prep$grid)
  a <- U[, d]; good <- rep(TRUE, M)
  for (it in seq_len(steps)) {
    U[, d] <- a
    input <- rbind(matrix(p, J, M), t(U), prep$grid)
    F <- prep$model$jet[[1]](input)
    jac <- array(prep$model$jet_jac[[1]](input), c(M, prep$D, J + prep$D))
    r <- F[, obs, drop = FALSE] - dmu
    G <- matrix(jac[, obs, J + d], M, length(obs))
    den <- rowSums(G * G)
    good <- sqrt(den) > tol * max(sqrt(den))
    if (!any(good)) return(prep$pilot)
    step <- ifelse(den > 0, rowSums(r * G) / pmax(den, .Machine$double.eps), 0)
    a[good] <- a[good] - step[good]
  }
  if (!all(good))
    a <- stats::approx(prep$grid[good], a[good], xout = prep$grid, rule = 2)$y
  U[, d] <- a
  U
}

.jgl_reference <- function(init, lat) {
  lapply(lat, function(d) {
    scale <- max(stats::sd(init[, d]), .1 * sqrt(mean(init[, d]^2)), 1e-6)
    list(scale = scale, mean = mean(init[, d]) / scale)
  })
}

# Stationary Matern 5/2 at a fixed fraction of the normalized domain. Radius
# bounds are grid/domain based, never chosen using latent truth.
.jgl_hyp <- function(prep, radius = 1/16) {
  span <- diff(range(prep$tt))
  bounds <- c(2 * mean(diff(prep$grid)) / span, 2)
  a0 <- stats::qlogis(pmin(.99, pmax(.01, log(radius / bounds[1]) / diff(log(bounds)))))
  lapply(seq_along(prep$latent), function(k)
    c(prep$reference[[k]], list(origin = min(prep$tt), span = span,
      bounds = bounds, a = a0, m = prep$reference[[k]]$mean, logtau = 0)))
}

.jgl_kernel <- function(hyp, tt, ctl)
  .jgp_correlation((tt - hyp$origin) / hyp$span, hyp$a, hyp$bounds, NULL, ctl)$R

# Prior state: observed components use their fitted GP's PRIOR covariance on the
# grid, the latent component the pilot kernel.
.jgl_state <- function(prep, hyps) {
  M <- length(prep$grid); D <- prep$D
  mu <- matrix(0, M, D); L <- vector("list", D)
  for (k in seq_along(prep$gps)) {
    d <- prep$observed[k]; g <- prep$gps[[k]]
    mu[, d] <- g$mean
    L[[d]] <- t(.jgp_chol(.jgp_gp_cross(g, prep$grid, prep$grid), prep$ctl$gp_jitter)$R)
  }
  for (k in seq_along(hyps)) {
    d <- prep$latent[k]; h <- hyps[[k]]
    mu[, d] <- h$scale * h$m
    L[[d]] <- h$scale * exp(h$logtau) *
      t(chol(.jgl_kernel(h, prep$grid, prep$ctl) + diag(prep$ctl$gp_jitter, M)))
  }
  list(mu = mu, L = L, mean = prep$pilot, include_gp_prior = TRUE)
}

# Fixed posterior coordinates precondition the observed data+GP block. The
# objective still contains the original observations and the original prior;
# only the coordinates the optimizer moves in change.
.jgl_cache <- function(prep) {
  M <- length(prep$grid)
  hyp <- .jgl_hyp(prep)
  prior_state <- .jgl_state(prep, hyp); st <- prior_state
  for (k in seq_along(prep$observed)) {
    o <- prep$observed[k]; L <- st$L[[o]]
    A <- prep$H %*% L / prep$noise[k]
    ch <- chol(diag(M) + crossprod(A)); Q <- backsolve(ch, diag(M))
    st$L[[o]] <- t(chol(tcrossprod(L %*% Q)))
    st$mu[, o] <- prior_state$mu[, o] + L %*% .jgp_solve(ch,
      crossprod(A, (prep$Y[, k] - prep$H %*% prior_state$mu[, o]) / prep$noise[k]))
  }
  list(M = M, d = prep$latent[1], st = st, prior_state = prior_state,
       scale = hyp[[1]]$scale, mean = hyp[[1]]$mean,
       pscale = pmax(abs(prep$p0), .1))
}

.jgl_pack <- function(p, U, st) c(p, unlist(lapply(seq_len(ncol(U)),
  function(d) forwardsolve(st$L[[d]], U[, d] - st$mu[, d]))))

# Full SVD basis of the screened bump pool.
.jgl_modes <- function(prep) {
  ctl <- prep$ctl; design <- prep$weak$design; ni <- nrow(design$interior)
  keep <- which(colSums(abs(prep$weak$maps$interior)) > 0)
  bas <- .jgp_basis_map(.jgp_test_rows(prep$grid, design$interior[keep, , drop = FALSE],
    0L, ctl$bump_eta), ctl$basis_tol, TRUE)
  full <- matrix(0, nrow(bas$map), ni); full[, keep] <- bas$map
  list(map = full, rank = nrow(bas$map), design = design,
       boundary = prep$weak$maps$boundary,
       singular = bas$singular_values[bas$retained])
}

# Component-ordered fits and prior factors for the weak build. In stage 1 the
# latent entry is an INTERPOLATION operator for quadrature only; no prior
# penalty is attached to it.
.jgl_weak_fits <- function(prep) {
  D <- prep$D; hyp <- .jgl_hyp(prep)
  fits <- vector("list", D)
  for (k in seq_along(prep$observed)) fits[[prep$observed[k]]] <- prep$gps[[k]]
  for (k in seq_along(prep$latent)) {
    d <- prep$latent[k]; h <- hyp[[k]]
    fits[[d]] <- list(origin = h$origin, span = h$span,
      tau2 = (h$scale * exp(h$logtau))^2, radius_coef = h$a,
      radius_bounds = h$bounds, control = prep$ctl,
      mean = h$scale * h$m, noise2 = prep$noise[1]^2)
  }
  L <- lapply(fits, function(f) {
    K <- .jgp_gp_cross(f, prep$grid, prep$grid)
    t(.jgp_chol((K + t(K)) / 2, prep$ctl$gp_jitter, "prior")$R)
  })
  list(fits = fits, state = list(L = L))
}

.jgl_build <- function(prep, mo, idx, order = prep$ctl$weak_quad_order) {
  ctl <- prep$ctl; D <- prep$D
  scaling <- .jgp_ode_scaling(matrix(rep(prep$units,each=length(prep$tt)),
    length(prep$tt)),prep$tt,ctl)
  scaling$scale_source <- if (is.null(ctl$ode_component_scale)) prep$scaling$scale_source else "supplied"
  prep$scaling <- scaling; prep$units <- scaling$rms
  fg <- .jgl_weak_fits(prep)
  maps <- list(interior = mo$map[idx, , drop = FALSE], boundary = mo$boundary)
  weak <- .jgp_gauss_pair(prep$grid, mo$design, fg$fits, fg$state, ctl, maps, order)
  weak$parameter_scale <- if (is.null(ctl$weak_design_scale)) pmax(abs(prep$p0),.1) else ctl$weak_design_scale
  if (length(weak$parameter_scale) != prep$J) stop("weak_design_scale must have one value per parameter.")
  weak$state_scale <- prep$units
  gram <- .jgp_ode_metric(weak, D, ctl)
  wt <- .jgp_ode_weights(gram, weak, D, ctl)
  prep$weak <- weak
  prep$weights <- .jgp_scale_weights(wt, prep$scaling$precision, weak$K)
  prep$extension_fits <- fg$fits
  prep$extension_state <- fg$state
  prep
}

# Stage 1: observed z blocks followed by one physical latent w block. There
# are no dummy latent z coordinates and no latent prior term.
.jgl_stage1_objective <- function(prep,cache) {
  coordinates <- cache$st
  coordinates$mu[,cache$d] <- 0
  coordinates$L[[cache$d]] <- cache$scale
  coordinates$order <- c(prep$observed,cache$d)
  coordinates$pscale <- cache$pscale
  .jgl_objective(prep,cache,coordinates)
}

# Both latent stages use exactly the same residual engine as observed fits.
.jgl_objective <- function(prep,cache,coordinates,prior=NULL,penalty=FALSE) {
  coordinates <- .jgp_bounded_coordinates(coordinates,attr(prep$ctl,"state_bounds"),prep$units)
  state <- cache$prior_state; components <- prep$observed
  if (penalty) {
    state$L[[cache$d]] <- prior$L; state$mu[,cache$d] <- prior$mean
    components <- c(components,cache$d)
  }
  evaluate <- .jgp_joint_objective(prep$Y,prep$H,prep$weak,prep$model,prep$weights,
    prep$noise,prep$ctl$lambda,coordinates,state,components,prep$observed)
  fn <- function(par,jacobian=TRUE,scalar=TRUE) {
    v <- evaluate(par,jacobian,scalar)
    cc <- v$contributions
    v$contributions <- c(data=unname(cc["data"]),
      observed_gp=sum(v$prior_contributions[as.character(prep$observed)]),
      ode=sum(cc[c("interior","boundary")]))
    if (penalty) v$contributions <- c(v$contributions,
      latent_gp=unname(v$prior_contributions[as.character(cache$d)]))
    if (!is.null(prior)) {
      ix <- prep$J+(match(cache$d,coordinates$order)-1L)*cache$M+seq_len(cache$M)
      v$latent_z <- if (length(coordinates$L[[cache$d]])==1L)
        sum(forwardsolve(prior$L,v$U[,cache$d]-prior$mean)^2)/2 else sum(par[ix]^2)/2
    }
    v
  }
  attr(fn,"coordinates") <- coordinates
  attr(fn,"bounds") <- .jgp_coordinate_bounds(coordinates,prep$lower,prep$upper,
    attr(prep$ctl,"state_bounds"))
  fn
}

.jgl_stage1 <- function(prep, cache, multiplier, maxit = 5000L,
                        solver = "nlminb", gradient_tol = 1e-4) {
  J <- length(prep$p0); M <- cache$M; nz <- M * length(prep$observed); d <- cache$d
  U <- prep$pilot; p <- prep$p0 * multiplier
  fn <- .jgl_stage1_objective(prep, cache)
  bounds <- attr(fn,"bounds")
  par <- .jgp_pack_coordinates(p,U,attr(fn,"coordinates"))
  par <- pmax(bounds$lower,pmin(bounds$upper,par))
  # Warm phase on the parameters and the latent only, holding the observed
  # blocks at their posterior-preconditioned start.
  warm <- .jgp_scalar_optimize(par, fn, min(75L,maxit), solver,
                       active = c(seq_len(J), J + nz + seq_len(M)),
                       lower=bounds$lower,upper=bounds$upper)
  par <- warm$par
  a <- .jgp_scalar_optimize(par,fn,maxit,solver,lower=bounds$lower,upper=bounds$upper)
  # Preserve the physical weak Jacobian in the stage-1 record; only the final
  # point needs this matrix, never the scalar optimizer's intermediate steps.
  res <- fn(a$par,TRUE,FALSE)
  check <- .jgl_diagnostics(prep,cache,res$p,res$U,gradient_tol)
  list(value = res$value, p = res$p, U = res$U, Ju = res$Ju,
       contributions = res$contributions, multiplier = multiplier,theta=a$par,
       converged = check$stationary, stationarity = check,
       optimizer = a[setdiff(names(a),"par")], warm = warm[setdiff(names(warm),"par")],
       projected_gradient = max(abs(res$gradient)))
}

.jgl_diagnostics <- function(prep,cache,p,U,tolerance,prior=NULL,penalty=FALSE,rank=FALSE) {
  state <- cache$prior_state; components <- prep$observed
  if (penalty) {
    state$L[[cache$d]] <- prior$L; state$mu[,cache$d] <- prior$mean
    components <- c(components,cache$d)
  }
  .jgp_fit_diagnostics(U,p,prep$Y,prep$H,prep$noise,prep$weak,prep$model,
    prep$weights,prep$ctl,state,components,observed=prep$observed,
    lower=prep$lower,upper=prep$upper,tolerance=tolerance,rank=rank)
}

# Stage 2: fit the latent prior to the stage-1 curve. The nugget absorbs the
# stage-1 roughness and is NOT carried into the prior.
.jgl_stage2 <- function(prep, u_latent) {
  fit <- .jgp_gp_fit(prep$grid, u_latent, NULL, prep$ctl)
  K <- .jgp_gp_cross(fit, prep$grid, prep$grid); K <- (K + t(K)) / 2
  ch <- .jgp_chol(K, prep$ctl$gp_jitter, "fitted latent prior")
  list(fit = fit, mean = fit$mean, tau = sqrt(fit$tau2), L = t(ch$R),
       smooth = .jgp_gp_predict(fit, prep$grid)$mean,
       nugget_sd = sqrt(fit$noise2), jitter = ch$added)
}

# Stage 3: joint solve with the fitted latent prior FROZEN.
.jgl_stage3_objective <- function(prep,cache,prior,penalty=FALSE) {
  coordinates <- cache$st
  coordinates$L[[cache$d]] <- prior$L; coordinates$mu[,cache$d] <- prior$mean
  coordinates$order <- seq_len(prep$D); coordinates$pscale <- cache$pscale
  .jgl_objective(prep,cache,coordinates,prior,penalty)
}

.jgl_stage3 <- function(prep, cache, prior, p0, U0, maxit = 5000L,
                        penalty = FALSE, solver = "nlminb", gradient_tol = 1e-4) {
  J <- length(prep$p0); M <- cache$M; D <- prep$D; d <- cache$d
  st <- cache$st; st$L[[d]] <- prior$L; st$mu[, d] <- prior$mean
  fn <- .jgl_stage3_objective(prep, cache, prior, penalty)
  bounds <- attr(fn,"bounds")
  par <- .jgp_pack_coordinates(p0,U0,attr(fn,"coordinates"))
  a <- .jgp_scalar_optimize(par,fn,maxit,solver,lower=bounds$lower,upper=bounds$upper)
  res <- fn(a$par)
  check <- .jgl_diagnostics(prep,cache,res$p,res$U,gradient_tol,prior,penalty)
  list(value = res$value, p = res$p, U = res$U, contributions = res$contributions,
       latent_z = res$latent_z, converged = check$stationary, stationarity = check,
       theta = a$par,
       optimizer = a[setdiff(names(a),"par")],
       projected_gradient = max(abs(res$gradient)))
}

# Screened design pool and frozen component units.
.jgl_prepare <- function(f, U, tt, p0, lower, upper, noise_sd, control,
                         screen = FALSE, latent_initial = NULL) {
  sp <- .jgl_split(U); obs <- sp$observed; lat <- sp$latent
  D <- ncol(U); J <- length(p0)
  Y <- U[, obs, drop = FALSE]
  sd_obs <- if (is.null(noise_sd)) NULL else rep_len(noise_sd, length(obs))
  gps <- lapply(seq_along(obs), function(k) .jgp_gp_fit(tt, Y[, k],
    if (is.null(sd_obs)) NULL else sd_obs[k], control))
  grid <- .jgp_grid(tt, control)
  fine <- seq(min(grid), max(grid), length.out = 2L * length(grid) - 1L)
  model <- .jgp_model(f, D, J, 0L)
  scale_pilot <- .jgl_initial(U, tt, grid, obs, lat)
  init <- .jgl_initial(U, tt, grid, obs, lat, latent_initial)
  initf <- .jgl_initial(U, tt, fine, obs, lat, latent_initial)
  bounds <- attr(control,"state_bounds")
  if (!is.null(bounds)) {
    # Smooth bounded observed starts before interpolation: high-order extension
    # of raw noise can leave a model's domain even when its grid values do not.
    for (k in seq_along(obs)) if (bounds$bounded[obs[k]]) {
      init[,obs[k]] <- .jgl_gp_mean(gps[[k]],grid)
      initf[,obs[k]] <- .jgl_gp_mean(gps[[k]],fine)
    }
    init <- .jgp_clip_state(init,bounds)
    initf <- .jgp_clip_state(initf,bounds)
  }
  design <- if (is.null(control$weak_radii)) .jgp_design_pool(grid, control) else
    .jgp_test_design(tt,gps,control)
  ni <- nrow(design$interior); nb <- nrow(design$boundary)
  raw <- .jgp_weak(grid, design, control, list(interior = diag(ni), boundary = diag(nb)))
  rf <- .jgp_weak(fine, design, control, raw$maps)
  err <- if (screen) .jgp_design_errors(raw, rf, init, initf, p0, model, pmax(abs(p0), 1)) else NULL
  # This screen rejects test functions whose TRAPEZOID quadrature is inaccurate,
  # and it is evaluated at the pilot. With a deliberately crude latent pilot the
  # residual is large everywhere and almost nothing passes, which collapses the
  # mode pool. The public Gauss path keeps the full pool and checks the actual
  # integration rule, residual, Jacobian and Gram independently.
  keep <- if (screen)
    which(err$relative[seq_len(ni)] <= control$weak_design_quad_tol) else seq_len(ni)
  if (!length(keep)) keep <- which(design$interior$radius == max(design$radii))
  bas <- .jgp_basis_map(.jgp_test_rows(grid, design$interior[keep, , drop = FALSE],
    0L, control$bump_eta), control$basis_tol, TRUE)
  size <- .jgp_design_svd_size(bas$singular_values[bas$retained],
    sum(bas$singular_values), control$weak_design_info)
  map <- matrix(0, size$count, ni)
  map[, keep] <- bas$map[seq_len(size$count), , drop = FALSE]
  bm <- .jgp_basis_map(.jgp_test_rows(grid, design$boundary, 0L, control$bump_eta),
    control$basis_tol)
  weak <- .jgp_weak(grid, design, control, list(interior = map, boundary = bm))
  # Freeze the default scale policy independently of a supplied initial curve.
  # This remains a heuristic for unknown latent units; explicit component
  # scales are supported and reported for comparisons across physical units.
  units <- sqrt(colMeans(scale_pilot^2)); units[obs] <- sqrt(colMeans(Y^2))
  if (any(!is.finite(units) | units <= 0)) stop("Unidentified pilot component scale.")
  scaling <- .jgp_ode_scaling(matrix(rep(units,each=length(tt)),length(tt)),tt,control)
  scaling$scale_source <- if (is.null(control$ode_component_scale)) "observed_rms_and_latent_initializer_heuristic" else "supplied"
  units <- scaling$rms
  gram <- .jgp_ode_metric(weak, D, control)
  wt <- .jgp_scale_weights(.jgp_ode_weights(gram, weak, D, control),
    scaling$precision, weak$K)
  prep <- list(Y = Y, tt = tt, grid = grid, H = .jgp_H(grid, tt), gps = gps,
    noise = sqrt(vapply(gps, `[[`, numeric(1), "noise2")), model = model,
    weak = weak, weights = wt, units = units, scaling = scaling, pilot = init, D = D, J = J,
    observed = obs, latent = lat, ctl = control, p0 = p0, lower = lower,
    upper = upper, f = f)
  prep$reference <- .jgl_reference(init, lat)
  prep
}

#' Compatibility wrapper for latent WENDyGP estimation
#'
#' @description Calls [solveWendyGP()] with formulation="latent". New code
#'   should use that unified entry point.
#'
#' @param f Function f(u,p,t), compatible with WENDy's symbolic machinery.
#' @param U Numeric observation matrix with exactly one entirely NA column.
#'   Every observed column must be fully finite, with at least six observations.
#' @param tt Strictly increasing finite observation times.
#' @param p0 Initial parameters, or NULL to infer their count and start at ones.
#' @param noise_sd Positive scalar or per-observed-component noise SD, or NULL.
#' @param control Overrides from [wendygp_control()]. This formulation requires
#'   matern52, gp_gauss, observed GP priors, and test_gram or identity weighting.
#'   It uses the full numerical test span. weak_radii and weak_design_radii
#'   control explicit tests and the multiscale pool respectively. The observed
#'   extension is Lagrange and the latent extension is GP. ode_units and
#'   ode_component_scale control the ODE metric. Supply component scales when
#'   the latent units cannot be inferred from the observed components; the
#'   default latent scale is an initializer-based heuristic, not unit invariant.
#' @param starts Finite parameter multipliers. Stationary completed runs are
#'   preferred, followed by finite feasible stage-3 estimates, then stage-1
#'   fallbacks. Runs in the preferred group are ranked by objective value;
#'   all runs are retained. With latent_penalty=TRUE
#'   only one start is supported, so different fitted priors are not compared.
#' @param modes NULL for the full numerical span, or that span's size. Truncation
#'   is not supported for the unpenalized latent stage.
#' @param screen Must be FALSE. The obsolete trapezoid pilot screen is not used;
#'   the actual Gauss residual, Jacobian and Gram are checked instead.
#' @param latent_penalty Add a frozen latent GP quadratic in stage 3. FALSE uses
#'   its fitted covariance only as an invertible preconditioner and starts from
#'   the stage-1 estimate, preserving its physical stationarity. TRUE starts
#'   from the smoothed stage-1 curve for the new penalized objective.
#' @param maxit Iteration cap per optimization phase, default 5000. An explicit
#'   argument takes precedence over control$maxit; otherwise that override is used.
#' @param gradient_tol Maximum Jacobian-column-scaled projected gradient in
#'   physical coordinates, for stages 1 and 3. Defaults to 1e-4; an explicit
#'   argument takes precedence over control$gtol. No objective normalization is used.
#' @param solver Scalar analytic-gradient optimizer, nlminb or lbfgsb.
#' @param latent_initial Optional finite latent starting curve at observation
#'   times. It does not change the frozen component scales. By default the row
#'   mean of observations is used, falling back to the observed column of
#'   largest RMS if the row mean cancels numerically.
#' @inheritParams solveWendyGP
#' @return A jointgp object (also inheriting wendygp_latent) with parameter/state estimates, all start
#'   records, native optimizer exits, physical stationarity, local Jacobian
#'   rank, and separate grid-resolution and quadrature diagnostics. A finite
#'   feasible stage-3 estimate is returned even without stationarity, with a
#'   warning and converged=FALSE. Stage 1 is returned only if no usable stage-3
#'   estimate exists. If neither stage has a finite feasible estimate, a
#'   wendygp_latent_numerical_error carries all runs.
#'   Local stationarity and rank do not establish global identifiability.
#' @export
solveWendyGPLatent <- function(f,U,tt,p0=NULL,noise_sd=NULL,control=NULL,
                               starts=1,modes=NULL,screen=FALSE,latent_penalty=FALSE,
                               maxit=5000L,gradient_tol=1e-4,
                               solver=c("nlminb","lbfgsb"),latent_initial=NULL,
                               parameter_lower=NULL,parameter_upper=NULL,
                               state_lower=NULL,state_upper=NULL,lower=NULL,upper=NULL) {
  # Compatibility only: validation, preparation and fitting live in solveWendyGP.
  args <- list(f=f,U=U,tt=tt,p0=p0,noise_sd=noise_sd,control=control,
    formulation="latent",starts=starts,modes=modes,screen=screen,
    latent_penalty=latent_penalty,solver=match.arg(solver),latent_initial=latent_initial,
    parameter_lower=parameter_lower,parameter_upper=parameter_upper,
    lower=lower,upper=upper,state_lower=state_lower,state_upper=state_upper)
  if (!missing(maxit)) args$maxit <- maxit
  if (!missing(gradient_tol)) args$gradient_tol <- gradient_tol
  answer <- do.call(solveWendyGP,args)
  answer$call <- match.call()
  answer
}

.jgp_latent_solve <- function(f, U, tt, p0 = NULL, noise_sd = NULL, control = NULL,
                               starts = 1, modes = NULL, screen = FALSE,
                               latent_penalty = FALSE, maxit = 5000L,
                               gradient_tol = 1e-4,
                               solver = c("nlminb", "lbfgsb"), latent_initial = NULL,
                               lower=NULL,upper=NULL,state_lower=NULL,state_upper=NULL) {
  solver <- match.arg(solver)
  supplied <- if (is.null(control)) list() else control
  if (missing(maxit) && "maxit" %in% names(supplied)) maxit <- supplied$maxit
  if (missing(gradient_tol) && "gtol" %in% names(supplied)) gradient_tol <- supplied$gtol
  # Validate unsupported choices before the observed control resolver can
  # coerce them to another integration scheme.
  if (!is.null(supplied$kernel) && !identical(supplied$kernel,"matern52"))
    stop("Latent estimation currently supports kernel='matern52' only.")
  if (isFALSE(supplied$include_gp_prior)) stop("Latent estimation requires observed GP priors.")
  if (identical(supplied$ode_weighting,"gp_delta"))
    stop("Latent estimation supports ode_weighting='test_gram' or 'identity', not 'gp_delta'.")
  if (!is.null(supplied$weak_radius_method) && supplied$weak_radius_method != "svd")
    stop("Latent estimation requires weak_radius_method='svd'; use weak_radii for explicit tests.")
  control <- do.call(wendygp_control, supplied)
  if (control$weak_integration != "gp_gauss") stop("Latent estimation requires weak_integration='gp_gauss'.")
  if (!is.null(control$weak_design_budget)) stop("Latent estimation uses the full test span; weak_design_budget is unsupported.")
  if (!is.matrix(U) || !is.numeric(U) || nrow(U)<6L || ncol(U)<2L)
    stop("U must be a numeric matrix with at least six rows and two components.")
  if (!is.numeric(tt) || length(tt)!=nrow(U) || any(!is.finite(tt)) || any(diff(tt)<=0))
    stop("tt must be strictly increasing, finite, and match nrow(U).")
  sp <- .jgl_split(U); D <- ncol(U)
  if (!is.numeric(starts) || !length(starts) || any(!is.finite(starts)))
    stop("starts must contain finite parameter multipliers.")
  if (!is.numeric(maxit) || length(maxit)!=1L || !is.finite(maxit) || maxit<1 ||
      maxit>.Machine$integer.max || maxit!=as.integer(maxit))
    stop("maxit must be one positive integer.")
  if (!is.numeric(gradient_tol) || length(gradient_tol)!=1L || !is.finite(gradient_tol) || gradient_tol<=0)
    stop("gradient_tol must be one positive finite number.")
  if (!isFALSE(screen)) stop("screen must be FALSE: trapezoid screening is not used by the latent Gauss formulation.")
  if (!is.logical(latent_penalty) || length(latent_penalty)!=1L || is.na(latent_penalty))
    stop("latent_penalty must be TRUE or FALSE.")
  if (latent_penalty && length(starts)>1L)
    stop("latent_penalty=TRUE requires one start; different stage-2 priors cannot be ranked as one fixed objective.")
  if (!is.null(modes) && (!is.numeric(modes) || length(modes)!=1L || !is.finite(modes) || modes<1 ||
      modes>.Machine$integer.max || modes!=as.integer(modes)))
    stop("modes must be NULL or a positive integer equal to the full span.")
  if (!is.null(noise_sd) && (!is.numeric(noise_sd) || !length(noise_sd)%in%c(1L,length(sp$observed)) ||
      any(!is.finite(noise_sd)) || any(noise_sd<=0)))
    stop("noise_sd must be positive, scalar or per observed component.")
  if (!is.null(latent_initial) && (!is.numeric(latent_initial) || length(latent_initial)!=nrow(U) || any(!is.finite(latent_initial))))
    stop("latent_initial must be finite with one value per observation time.")
  J <- if (is.null(p0)) detect_n_params(f) else length(p0)
  if (J<1L) stop("Supply p0 so the number of parameters is known.")
  if (is.null(p0)) p0 <- rep(1,J)
  if (!is.numeric(p0) || any(!is.finite(p0))) stop("p0 must be a finite numeric vector.")
  bounds <- .jgp_bound_pair(lower,upper,J,"parameter")
  lower <- bounds$lower; upper <- bounds$upper
  p0 <- pmax(lower,pmin(upper,p0))
  attr(control,"state_bounds") <- .jgp_bound_pair(state_lower,state_upper,D,"state")
  expected_extension <- ifelse(seq_len(D)%in%sp$observed,"lagrange","gp")
  if ("weak_extension" %in% names(supplied) &&
      !identical(supplied$weak_extension,"lagrange") &&
      !identical(supplied$weak_extension,expected_extension))
    stop("The latent formulation uses Lagrange for observed components and GP for the latent component.")
  control$weak_extension <- expected_extension
  control$em_order <- 0L; control$weak_design_info <- 1
  control$maxit <- maxit; control$gtol <- gradient_tol
  prep0 <- .jgl_prepare(f,U,tt,p0,lower,upper,noise_sd,control,FALSE,latent_initial)
  mo <- .jgl_modes(prep0)
  if (!is.null(modes) && modes!=mo$rank)
    stop("Latent estimation requires the full numerical span (",mo$rank," modes).")
  idx <- seq_len(mo$rank)
  prep <- .jgl_build(prep0,mo,idx); cache <- .jgl_cache(prep)
  accuracy <- function(U,p) .jgp_gauss_accuracy(U,p,prep$extension_state,prep$weak,
    prep$model,prep$extension_fits,prep$weights,prep$H,control,
    observed=prep$observed,observation_times=tt)
  initial_accuracy <- accuracy(prep$pilot,p0)
  runs <- lapply(starts,function(mult) {
    run <- list(multiplier=mult,status="error",stage1=NULL,prior=NULL,stage3=NULL)
    tryCatch({
      run$stage1 <- .jgl_stage1(prep,cache,mult,maxit,solver,gradient_tol)
      run$stage1$accuracy <- accuracy(run$stage1$U,run$stage1$p)
      if (!run$stage1$converged) {
        run$status <- "not_stationary"
        run$message <- "Stage 1 did not reach physical stationarity."
        return(run)
      }
      run$prior <- .jgl_stage2(prep,run$stage1$U[,cache$d])
      if (run$prior$fit$convergence!=0L) {
        run$status <- "gp_not_converged"
        run$message <- "Stage-2 GP fitting did not converge."
        return(run)
      }
      U0 <- run$stage1$U
      if (latent_penalty) U0[,cache$d] <- run$prior$smooth
      run$stage3 <- .jgl_stage3(prep,cache,run$prior,run$stage1$p,U0,maxit,
        penalty=latent_penalty,solver=solver,gradient_tol=gradient_tol)
      if (!run$stage3$converged) {
        run$status <- "not_stationary"
        run$message <- "Stage 3 did not reach physical stationarity."
        return(run)
      }
      run$status <- "stationary"
      run
    },error=function(e) {run$message <- conditionMessage(e); run})
  })
  usable <- function(x) !is.null(x) && is.finite(x$value) &&
    all(is.finite(c(x$p,x$U))) && isTRUE(x$stationarity$feasible)
  valid <- which(vapply(runs,function(r)identical(r$status,"stationary"),logical(1)))
  selected_stage <- vapply(runs,function(r) {
    if (usable(r$stage3)) 3L else
      if (usable(r$stage1)) 1L else 0L
  },integer(1))
  # Return the final optimization's estimate even when it is nonstationary or
  # has a higher objective than stage 1. Keep stage 1 as a numerical fallback.
  stage3_available <- which(selected_stage==3L)
  eligible <- if (length(valid)) valid else if (length(stage3_available))
    stage3_available else which(selected_stage==1L)
  if (!length(eligible)) stop(structure(list(
    message="No finite, feasible latent estimate was produced. Inspect condition$runs for numerical failures.",
    call=NULL,runs=runs),class=c("wendygp_latent_numerical_error","error","condition")))
  values <- vapply(eligible,function(k) runs[[k]][[paste0("stage",selected_stage[k])]]$value,numeric(1))
  chosen <- eligible[which.min(values)]
  best <- runs[[chosen]]; final_stage <- selected_stage[chosen]
  final <- best[[paste0("stage",final_stage)]]
  applied_penalty <- latent_penalty && final_stage==3L
  check <- .jgl_diagnostics(prep,cache,final$p,final$U,
    gradient_tol,best$prior,applied_penalty,rank=TRUE)
  acc <- accuracy(final$U,final$p)
  pipeline_converged <- identical(best$status,"stationary") && check$stationary
  reason <- if (pipeline_converged) final$optimizer$reason else
    paste0(best$message," Returning the stage-",final_stage," estimate.")
  if (!pipeline_converged) warning(reason," See fit$runs; converged=FALSE.",call.=FALSE)
  grid_passed <- initial_accuracy$passed && best$stage1$accuracy$passed && acc$passed
  weak_passed <- acc$quadrature_passed
  if (!grid_passed) {
    msg <- "Latent fit observation-operator error exceeds grid_tol; inspect initial_quadrature and final_grid."
    if (control$grid_action=="error") stop(msg) else warning(msg,call.=FALSE)
  }
  if (!weak_passed) warning(.jgp_quadrature_message(acc,control),call.=FALSE)
  evaluate <- if (final_stage==3L) .jgl_stage3_objective(prep,cache,best$prior,applied_penalty) else
    .jgl_stage1_objective(prep,cache)
  structure(list(phat=final$p,U_hat=final$U,
    U_obs_hat=prep$H%*%final$U[,prep$observed,drop=FALSE],
    U_stage1=best$stage1$U,U_smooth=best$prior$smooth,
    observed=prep$observed,latent=prep$latent,tt=prep$grid,tt_obs=tt,
    Y=prep$Y,noise_sd=prep$noise,gp=prep$gps,latent_prior=best$prior,
    lambda=control$lambda,ode_weighting=control$ode_weighting,ode_units=control$ode_units,
    objective=final$value,converged=pipeline_converged,
    optimizer=final$optimizer,contributions=final$contributions,
    runs=runs,multiplier=best$multiplier,selected_start=chosen,
    diagnostics=list(final_stage=final_stage,pipeline_converged=pipeline_converged,
      stage_status=best$status,modes=length(idx),available=mo$rank,
      weak_rows=prep$weak$K*D,state_variables=length(final$U),
      stage1_variables=J+length(best$stage1$U),
      extension=acc$extension,extension_passed=acc$extension_passed,
      weak_extension=control$weak_extension,ode_scaling=prep$scaling,
      stationarity=check,scaled_gradient=check$scaled_gradient,gradient_tol=gradient_tol,
      grid_passed=grid_passed,weak_grid_passed=weak_passed,
      quadrature_passed=acc$quadrature_passed,
      initial_quadrature=initial_accuracy,stage1_grid=best$stage1$accuracy,final_grid=acc,
      gp_converged=vapply(prep$gps,function(g)g$convergence==0,logical(1)),
      stage1_value=best$stage1$value,stage1_converged=best$stage1$converged,
      latent_nugget_sd=best$prior$nugget_sd,latent_tau=best$prior$tau,
      latent_penalty=applied_penalty,latent_penalty_requested=latent_penalty,latent_z=final$latent_z),
    problem=list(evaluate=function(theta,jacobian=TRUE,scalar=FALSE) evaluate(theta,jacobian,scalar),theta=final$theta,
      coordinates=attr(evaluate,"coordinates"),H=prep$H,weak=prep$weak,model=prep$model,
      weights=prep$weights,control=control,lower=lower,upper=upper,
      state_bounds=attr(control,"state_bounds")),
    include_gp_prior=TRUE,weak_integration="gp_gauss",
    convergence_reason=reason,iterations=final$optimizer$iterations,
    control=control,call=match.call()),class=c("jointgp","wendygp_latent"))
}
