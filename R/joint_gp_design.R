# Multiscale SVD and experimental information-directed weak spaces. No GP radius
# is read by the geometric pool. All candidates share one pilot state/parameter, covariance
# model, observation likelihood, prior, grid, and parameter scaling.

.jgp_design_rows <- function(physical, K, D) {
  unlist(lapply(seq_len(D), function(d) (d - 1L) * K + physical), use.names = FALSE)
}

.jgp_design_pool <- function(grid, control, radii = NULL, shift = 0) {
  span <- diff(range(grid)); a <- min(grid); b <- max(grid)
  if (is.null(radii)) radii <- control$weak_design_radii
  if (is.null(radii)) {
    # Geometric ladder from 0.4*span down to the finest rung still at least one
    # grid spacing wide. A test narrower than h has support containing at most
    # one state node, so it carries no information the state can represent --
    # measured, going below h degraded three of six latent cells and collapsed
    # one by 1000x. The old fixed fractions span*c(1/32,...,0.4) put the finest
    # test at (m-1)/32 grid spacings, so the pool got COARSER relative to the
    # grid as the grid was refined; this is grid adaptive instead.
    h <- span / max(1L, length(grid) - 1L)
    radii <- .4 * span / 2^(0:max(0L, floor(log2(.4 * span / h))))
  }
  radii <- sort(unique(radii))
  if (any(radii >= span / 2)) stop("weak_design_radii must be smaller than half the span.")
  interior <- do.call(rbind, lapply(radii, function(r) {
    if (is.null(control$weak_design_centers)) {
      # No subsampling or added off-grid endpoints. The small tolerance only
      # admits exact support endpoints affected by floating-point roundoff.
      candidates <- grid + shift*r
      tol <- 32 * .Machine$double.eps * max(1, abs(a), abs(b), span)
      centers <- candidates[candidates >= a+r-tol & candidates <= b-r+tol]
    } else {
      n <- max(3L, min(control$weak_design_centers, ceiling((span - 2*r)/(r/4)) + 1L))
      centers <- seq(a+r, b-r, length.out = n)
      if (n > 2L) centers[2:(n-1L)] <- pmax(a+r, pmin(b-r, centers[2:(n-1L)] + shift*r))
    }
    data.frame(center = centers, radius = rep(r, length(centers)))
  }))
  if (!nrow(interior)) stop("No admissible interior grid centers; refine the working grid.")
  br <- if (is.null(control$bl_radii)) max(radii) else control$bl_radii
  if (any(br >= span/2)) stop("bl_radii must be smaller than half the span.")
  boundary <- if (control$include_bl) do.call(rbind, lapply(br, function(r) {
    off <- seq(0, .4*r, length.out = control$bl_count)
    data.frame(center = c(a+off, b-off), radius = r)
  })) else data.frame(center = numeric(), radius = numeric())
  list(interior = interior, boundary = boundary, radii = radii, boundary_radii = br,
       radius_selection = list(method = "geometric", shift = shift,
         center_placement = if (is.null(control$weak_design_centers)) "grid" else "capped",
         grid_points = length(grid)))
}

# Compare the complete weak residual AND parameter sensitivities on nested
# quadratures. Normalize by integral magnitudes, not measurement noise or a
# residual that should be zero at a solution. Screening is local at the pilot.
.jgp_design_errors <- function(weak, fine, U, fine_U, p, model, pscale) {
  v <- .jgp_weak_eval(U, p, weak, model)
  vf <- .jgp_weak_eval(fine_U, p, fine, model)
  F <- model$jet[[1]](rbind(matrix(p, model$J, nrow(U)), t(U), weak$tt))
  magnitude <- abs(weak$V) %*% abs(F) + abs(weak$Vp) %*% abs(U)
  delta <- matrix(v$r - vf$r, weak$K, model$D)
  dG <- sweep(v$Jp - vf$Jp, 2, pscale, "*")
  sensitivity <- matrix(sqrt(rowSums(dG^2)), weak$K, model$D)
  for (d in seq_len(model$D)) magnitude[,d] <- pmax(magnitude[,d],
    .Machine$double.eps * max(magnitude[,d]), .Machine$double.xmin)
  list(base = v, delta = as.vector(delta),
       relative = apply(pmax(abs(delta), sensitivity) / magnitude, 1, max))
}

# Ju C_U|data,GP Ju' without subtracting nearly equal covariance matrices.
# In GP-whitened state coordinates C_z = (I + A'A)^-1, A = H L / sigma.
# Its square root follows from the small observation-space SVD of A.
.jgp_design_nuisance <- function(Ju, state, H, noise) {
  m <- nrow(state$mu); out <- matrix(0, nrow(Ju), nrow(Ju))
  for (d in seq_along(state$L)) {
    L <- state$L[[d]]; A <- H %*% L / noise[d]
    decomp <- svd(A, nu = 0, nv = min(dim(A)))
    B <- Ju[, (d-1L)*m + seq_len(m), drop = FALSE] %*% L
    shrink <- sqrt(1 / (1 + decomp$d^2)) - 1
    B <- B + sweep(B %*% decomp$v, 2, shrink, "*") %*% t(decomp$v)
    out <- out + tcrossprod(B)
  }
  (out + t(out))/2
}

.jgp_design_prepare <- function(tt, grid, fits, state, H, p, model, control, pscale,
                                noise, shift = 0) {
  design <- .jgp_design_pool(grid, control, shift = shift)
  ni <- nrow(design$interior); nb <- nrow(design$boundary)
  raw <- .jgp_weak(grid, design, control,
    list(interior = diag(ni), boundary = diag(nb)))
  ft <- seq(min(grid), max(grid), length.out = 2L*length(grid)-1L)
  fU <- do.call(cbind, lapply(fits, function(g) .jgp_gp_predict(g, ft)$mean))
  fine <- .jgp_weak(ft, design, control, raw$maps)
  check <- .jgp_design_errors(raw, fine, state$mean, fU, p, model, pscale)
  keep_raw <- which(check$relative[seq_len(ni)] <= control$weak_design_quad_tol)
  fallback <- !length(keep_raw)
  if (fallback) keep_raw <- which(design$interior$radius == max(design$radii))
  basis <- .jgp_basis_map(.jgp_test_rows(grid, design$interior[keep_raw,,drop=FALSE],
                                   0L, control$bump_eta), control$basis_tol, spectrum = TRUE)
  q <- basis$map
  map <- matrix(0, nrow(q), ni); map[,keep_raw] <- q
  bmap <- .jgp_basis_map(.jgp_test_rows(grid, design$boundary, 0L, control$bump_eta), control$basis_tol)
  weak <- .jgp_weak(grid, design, control, list(interior = map, boundary = bmap))
  fine <- .jgp_weak(ft, design, control, weak$maps)
  check_modes <- .jgp_design_errors(weak, fine, state$mean, fU, p, model, pscale)
  v <- check_modes$base
  S <- .jgp_ode_metric(weak, model$D, control, v$Ju, state$Sigma)
  sd <- sqrt(pmax(diag(S), .Machine$double.xmin))
  standard_error <- apply(matrix(abs(check_modes$delta) / sd, weak$K, model$D), 1, max)
  keep_modes <- which(check_modes$relative[seq_len(weak$ni)] <= control$weak_design_quad_tol &
    standard_error[seq_len(weak$ni)] <= control$weak_grid_tol)
  if (!length(keep_modes)) { keep_modes <- 1L; fallback <- TRUE }
  physical <- c(keep_modes, weak$ni + seq_len(weak$nb))
  rows <- .jgp_design_rows(physical, weak$K, model$D)
  maps <- list(interior = weak$maps$interior[keep_modes,,drop=FALSE], boundary = weak$maps$boundary)
  selected <- .jgp_weak(grid, design, control, maps)
  Ju <- v$Ju[rows,,drop=FALSE]
  singular_values <- basis$singular_values[basis$retained[keep_modes]]
  singular_total <- sum(basis$singular_values)
  information_available <- sum(singular_values) / singular_total
  information_passed <- information_available + 32*.Machine$double.eps >= control$weak_design_info
  list(weak = selected, G = v$Jp[rows,,drop=FALSE], Ju = Ju,
    r = v$r[rows], Omega = S[rows,rows,drop=FALSE],
    ode_precision = .jgp_ode_scaling(do.call(cbind, lapply(fits, `[[`, "y")), tt, control)$precision,
    singular_values = singular_values, singular_total = singular_total,
    include_gp_prior = !isFALSE(state$include_gp_prior),
    nuisance = if (!isFALSE(state$include_gp_prior)) .jgp_design_nuisance(Ju, state, H, noise) else NULL,
    data_jacobian = if (isFALSE(state$include_gp_prior)) kronecker(diag(1/noise,length(noise)),H) else NULL,
    screening = list(passed = !fallback && (!is.null(control$weak_design_budget) || information_passed),
      center_placement = design$radius_selection$center_placement,
      grid_points = length(grid), raw_count = ni, raw_retained = length(keep_raw),
      mode_count = weak$ni, mode_retained = length(keep_modes),
      raw_relative_error = check$relative[seq_len(ni)],
      mode_relative_error = check_modes$relative[seq_len(weak$ni)],
      mode_standardized_error = standard_error[seq_len(weak$ni)],
      retained_raw = keep_raw, retained_modes = keep_modes,
      singular_values = basis$singular_values, information_available = information_available,
      information_target = control$weak_design_info, information_passed = information_passed,
      tolerance = control$weak_design_quad_tol, fallback = fallback))
}

# MSG information-number convention: cumulative singular values, not their
# squares. Never renormalize away modes lost to numerical/accuracy screening.
.jgp_design_svd_size <- function(values, total, target) {
  if (!length(values) || !is.finite(total) || total <= 0)
    stop("No finite positive singular-value information is available.")
  fraction <- cumsum(values) / total
  reached <- which(fraction + 32*.Machine$double.eps >= target)
  list(count = if (length(reached)) reached[1L] else length(values),
       fractions = fraction, target_met = length(reached) > 0L)
}

# Profile completely free grid-state increments out of the stacked data/ODE
# Jacobian. This is a generalized Schur complement even when the data-only
# state Hessian is singular. It uses no GP covariance or ranking ridge.
.jgp_free_state_information <- function(A, B) {
  scales <- pmax(sqrt(colSums(B^2)), .Machine$double.xmin)
  s <- svd(sweep(B,2,scales,"/"),nu=min(dim(B)),nv=0)
  tol <- max(dim(B)) * .Machine$double.eps
  keep <- which(s$d > max(s$d)*tol)
  if (length(keep)==nrow(B)) residual <- A*0 else {
    Q <- s$u[,keep,drop=FALSE]
    residual <- A-Q%*%crossprod(Q,A)
  }
  list(matrix=crossprod(residual),state_rank=length(keep),
       state_nullity=ncol(B)-length(keep),state_rank_tolerance=tol)
}

# EXACT Schur complement of the local joint Gauss-Newton matrix (up to the
# existing covariance rank projection), computed in residual rather than state
# dimension. Parameter scales are fixed from user inputs, never truth.
.jgp_design_information <- function(space, selected, model, lambda, pscale, control) {
  weak <- space$weak
  physical <- c(sort(selected), weak$ni + seq_len(weak$nb))
  rows <- .jgp_design_rows(physical, weak$K, model$D)
  layout <- list(K = length(physical), ni = length(selected), nb = weak$nb)
  tryCatch({
    wt <- .jgp_ode_weights(space$Omega[rows,rows,drop=FALSE], layout, model$D, control)
    if (!is.null(space$ode_precision)) wt <- .jgp_scale_weights(wt, space$ode_precision, layout$K)
    A <- sqrt(lambda) * sweep(wt$W %*% space$G[rows,,drop=FALSE], 2, pscale, "*")
    free <- NULL
    if (isFALSE(space$include_gp_prior)) {
      B <- rbind(space$data_jacobian,sqrt(lambda)*wt$W%*%space$Ju[rows,,drop=FALSE])
      free <- .jgp_free_state_information(rbind(matrix(0,nrow(space$data_jacobian),ncol(A)),A),B)
      info <- free$matrix
    } else {
      M <- diag(nrow(wt$W)) + lambda * wt$W %*%
        space$nuisance[rows,rows,drop=FALSE] %*% t(wt$W)
      R <- chol((M+t(M))/2)
      info <- crossprod(forwardsolve(t(R), A))
    }
    c(list(ok = TRUE, matrix = info, eigenvalues = pmax(0, eigen(info, symmetric=TRUE, only.values=TRUE)$values),
      covariance_rank = nrow(wt$W),include_gp_prior=!isFALSE(space$include_gp_prior)),
      if (!is.null(free)) free[setdiff(names(free),"matrix")])
  }, error = function(e) list(ok = FALSE, message = conditionMessage(e),
    matrix = matrix(0, model$J, model$J), eigenvalues = rep(0, model$J)))
}

.jgp_design_select <- function(space, method, model, lambda, pscale, control, extra = 0L) {
  count <- space$weak$ni
  spectral <- .jgp_design_svd_size(space$singular_values, space$singular_total,
                                  control$weak_design_info)
  requested <- if (is.null(control$weak_design_budget)) spectral$count else control$weak_design_budget
  budget <- min(count, requested + extra)
  full <- .jgp_design_information(space, seq_len(count), model, lambda, pscale, control)
  # This floor is ONLY a deterministic ranking device for singular candidate
  # information matrices. It never changes Omega, the posterior, or the solver.
  floor <- max(max(full$eigenvalues)*1e-8, .Machine$double.eps)
  metric <- function(x) if (x$ok) -sum(1/(x$eigenvalues + floor)) else -Inf
  history <- list()
  if (method == "svd") selected <- seq_len(budget) else {
    selected <- seq_len(min(budget, control$weak_design_coverage))
    while (length(selected) < budget) {
      remaining <- setdiff(seq_len(count), selected)
      infos <- lapply(remaining, function(j)
        .jgp_design_information(space, c(selected,j), model, lambda, pscale, control))
      scores <- vapply(infos, metric, numeric(1)); at <- which.max(scores)
      selected <- c(selected, remaining[at])
      history[[length(history)+1L]] <- data.frame(step = length(selected),
        mode = remaining[at], score = scores[at], min_eigenvalue = min(infos[[at]]$eigenvalues))
    }
  }
  selected <- sort(selected)
  info <- .jgp_design_information(space, selected, model, lambda, pscale, control)
  maps <- list(interior = space$weak$maps$interior[selected,,drop=FALSE], boundary = space$weak$maps$boundary)
  diagnostics <- list(method = method, selected = selected, budget = budget,
    requested_budget = requested + extra,
    size_rule = if (is.null(control$weak_design_budget)) "singular_value_fraction" else "fixed_budget",
    information_target = control$weak_design_info,
    information_fraction = sum(space$singular_values[selected]) / space$singular_total,
    information_target_met = sum(space$singular_values[selected]) / space$singular_total +
      32*.Machine$double.eps >= control$weak_design_info,
    singular_values = space$singular_values, singular_total = space$singular_total,
    spectral_budget = spectral$count,
    coverage_modes = min(budget, control$weak_design_coverage),
    parameter_scale = pscale, information = info, full_information = full,
    score_floor = floor, score = metric(info), history = if (length(history)) do.call(rbind, history) else NULL,
    screening = space$screening, local_precision_only = TRUE)
  design <- space$weak$design; design$radius_selection <- diagnostics
  weak <- .jgp_weak(space$weak$tt, design, control, maps)
  physical <- c(selected, space$weak$ni + seq_len(space$weak$nb))
  rows <- .jgp_design_rows(physical, space$weak$K, model$D)
  list(weak = weak, Omega = space$Omega[rows,rows,drop=FALSE], diagnostics = diagnostics)
}

.jgp_design_fit <- function(Y, tt, fits, model, control, lower, upper, state, H,
                            pilot, spec, initial_error) {
  weak <- spec$weak; D <- model$D; J <- model$J; m <- nrow(state$mu)
  noise <- sqrt(vapply(fits, `[[`, numeric(1), "noise2"))
  metric <- .jgp_fit_metric(spec$Omega, weak, Y, tt, control)
  weights <- metric$weights
  theta0 <- c(pilot, .jgp_state_start(state))
  objective <- .jgp_objective(Y, H, state, weak, model, weights, noise, control$lambda)
  coordinates <- state$coordinates; coordinates$pscale <- rep(1,J)
  bounds <- .jgp_coordinate_bounds(coordinates,lower,upper,state$bounds)
  result <- .jgp_optimize(theta0,objective,control,bounds$lower,bounds$upper)
  point <- .jgp_unpack(result$par, state, J); final <- objective(result$par, TRUE)
  stationarity <- .jgp_fit_diagnostics(point$U,point$p,Y,H,noise,weak,model,weights,control,
    state,if (state$include_gp_prior) seq_len(D) else integer(),lower=lower,upper=upper)
  accuracy <- .jgp_solution_accuracy(point$U, point$p, state, weak, model, fits, weights, H, control)
  initial <- objective(theta0, FALSE)
  rho <- function(r) sum((weights$Wi %*% r)^2) / weights$interior$rank
  structure(list(phat = point$p, U_hat = point$U, U_obs_hat = H %*% point$U,
    tt = weak$tt, tt_obs = tt, Y = Y, gp = fits, noise_sd = noise,
    initial = list(p = pilot, U = state$mean, theta = theta0), lambda = control$lambda,
    ode_weighting = .jgp_ode_mode(control),
    ode_units = metric$scaling$units,
    include_gp_prior = state$include_gp_prior,
    objective = result$objective, contributions = final$contributions, converged = stationarity$stationary,
    optimizer = list(method=result$method,converged=result$converged,reason=result$reason,
      iterations=result$iterations,evaluations=result$evaluations),
    convergence_reason = result$reason, iterations = result$iterations,
    diagnostics = list(rho_initial = rho(initial$raw), rho_final = rho(final$raw),
      ode_weighting = .jgp_ode_mode(control),
      ode_scaling = metric$scaling,
      include_gp_prior = state$include_gp_prior,
      state_coordinates = if (state$prior_whitened) "gp_whitened" else "physical_or_mixed",
      gp_uncertainty_propagated = .jgp_ode_mode(control) == "gp_delta",
      rho_is_gp_standardized = .jgp_ode_mode(control) == "gp_delta", ode_metric_frozen = TRUE,
      radius_selection = spec$diagnostics, stationarity = stationarity,
      scaled_gradient = stationarity$scaled_gradient,
      design_passed = !identical(spec$diagnostics$screening$passed, FALSE),
      gradient = stationarity$gradient, interior_rank = weights$interior$rank, boundary_rank = weights$boundary$rank,
      gp_jitter = state$jitter, gp_converged = vapply(fits, function(g) g$convergence == 0, logical(1)),
      grid_passed = accuracy$passed,
      grid_fixed = TRUE, grid_refinements = 0L,
      grid_acceptance = "observation_operator",
      weak_grid_passed = max(initial_error) <= control$weak_grid_tol && accuracy$weak_passed,
      initial_error = initial_error, final_grid = accuracy, history = result$history,
      covariance_frozen = TRUE, hyperparameters_frozen = TRUE),
    problem = list(evaluate = objective, theta = result$par, state = state, H = H,
      weak = weak, model = model, Omega = metric$Omega, Omega_raw = spec$Omega, weights = weights,
      control = control, lower = lower, upper = upper,state_bounds=state$bounds)), class = "jointgp")
}

# Controlled groups are used by validation to compare all spaces on a common
# grid. A single requested arm is also available through solveWendyGP.
.jgp_design_group <- function(Y, tt, p0, fits, model, control, lower, upper,
    arms = c("gp", "svd", "sensitivity")) {
  pscale <- control$weak_design_scale
  if (is.null(pscale)) pscale <- pmax(abs(p0), 1)
  if (length(pscale) != model$J) stop("weak_design_scale must have one entry per parameter.")
  grid <- .jgp_grid(tt, control)
  noise <- sqrt(vapply(fits, `[[`, numeric(1), "noise2"))
  state <- .jgp_state(fits, grid, control); H <- .jgp_H(grid, tt)
  pd <- .jgp_design_pool(grid, control, diff(range(tt))*c(.25,.4))
  pw <- .jgp_weak(grid, pd, control)
  pilot <- .jgp_initialize(state$mean, p0, pw, model, state$Sigma, control, lower, upper)
  space <- .jgp_design_prepare(tt, grid, fits, state, H, pilot, model, control, pscale, noise)
  shifted <- if ("sensitivity_shifted" %in% arms)
    .jgp_design_prepare(tt, grid, fits, state, H, pilot, model, control, pscale, noise, shift=.125) else NULL
  specs <- initial_errors <- list()
  interpolation <- vapply(fits, .jgp_interpolation_error, numeric(1), tt=grid, obs=tt, H=H)
  for (arm in arms) {
    if (arm == "gp") {
      cc <- control; cc$weak_radius_method <- "gp"
      w <- .jgp_weak(grid, .jgp_test_design(tt, fits, cc), cc)
      v <- .jgp_weak_eval(state$mean, pilot, w, model)
      specs[[arm]] <- list(weak=w, Omega=.jgp_ode_metric(w,model$D,control,v$Ju,state$Sigma),
                          diagnostics=list(method="gp_geometry_common_pilot"))
    } else specs[[arm]] <- .jgp_design_select(if (arm=="sensitivity_shifted") shifted else space,
      if (arm=="svd") "svd" else "sensitivity", model, control$lambda, pscale, control,
      extra=if (arm=="sensitivity_enriched") 2L else 0L)
    s <- specs[[arm]]
    wt <- .jgp_fit_metric(s$Omega, s$weak, Y, tt, control)$weights
    initial_errors[[arm]] <- if (control$lambda > 0)
      .jgp_refinement(state, pilot, s$weak, model, fits, wt, control) else c(interior=0,boundary=0)
  }
  initial_weak_ok <- vapply(initial_errors, function(x) max(x)<=control$weak_grid_tol, logical(1))
  initial_ok <- setNames(rep(max(interpolation)<=control$grid_tol,length(arms)),arms)
  answers <- lapply(arms, function(arm) tryCatch({
    fit <- .jgp_design_fit(Y,tt,fits,model,control,lower,upper,state,H,pilot,specs[[arm]],initial_errors[[arm]])
    fit$diagnostics$grid_passed <- fit$diagnostics$grid_passed && max(interpolation)<=control$grid_tol
    fit$diagnostics$initial_interpolation <- interpolation
    fit
  }, error=identity)); names(answers) <- arms
  passed <- vapply(answers, function(x) !inherits(x,"error") && x$diagnostics$grid_passed, logical(1))
  history <- list(list(grid=length(grid),stage="final",initial_passed=initial_ok,
                       initial_weak_passed=initial_weak_ok,passed=passed))
  # Shifted, unselected weak tests are a coverage diagnostic, not independent
  # validation observations or a chi-square calibration test.
  hd <- .jgp_design_pool(grid, control, diff(range(tt))*c(.125,.25,.4), shift=-.125)
  hw <- .jgp_weak(grid, hd, control)
  hv <- .jgp_weak_eval(state$mean, pilot, hw, model)
  holdout <- tryCatch(.jgp_ode_weights(.jgp_ode_metric(hw,model$D,control,hv$Ju,state$Sigma),
                                      hw,model$D,control), error=identity)
  if (!inherits(holdout,"error")) holdout <- .jgp_scale_weights(holdout,
    .jgp_ode_scaling(Y, tt, control)$precision, hw$K)
  for (arm in arms) if (!inherits(answers[[arm]],"error")) {
    fit <- answers[[arm]]; fit$diagnostics$grid_history <- history
    if (!inherits(holdout,"error")) {
      fit$diagnostics$heldout_tests <- tryCatch({
        rr <- .jgp_weak_eval(fit$U_hat,fit$phat,hw,model,FALSE)$r
        acc <- .jgp_solution_accuracy(fit$U_hat,fit$phat,state,hw,model,fits,holdout,H,control)
        list(rms=sqrt(mean((holdout$W%*%rr)^2)), grid_passed=acc$passed, weak_grid_passed=acc$weak_passed,
             accuracy=acc, independent_data=FALSE)
      }, error=function(e) list(error=conditionMessage(e), grid_passed=FALSE, independent_data=FALSE))
    } else fit$diagnostics$heldout_tests <- list(error=conditionMessage(holdout),
      grid_passed=FALSE, independent_data=FALSE)
    answers[[arm]] <- fit
  }
  list(fits=answers, pilot=pilot, grid=grid, history=history,
       screening=space$screening, parameter_scale=pscale)
}
