# Opt-in staged GP projection. This module does not change solveWendyGP.
# Fit actual observations once; the quadrature nodes are never observations.
# The nonlinear solves minimize a posterior Mahalanobis distance subject to
# retained weak equations, NOT a tempered data/GP/ODE penalty.

wendygp_staged_control <- function(...) {
  defaults <- list(grid_min = 64L, grid = NULL, pilot = "gp_mean",
    weak_radii = NULL, weak_info = .95, weak_budget = NULL, bump_eta = 9,
    include_bl = TRUE, bl_radius = NULL, bl_count = 1L, boundary_points = NULL,
    quad_order = 12L, quad_tol = 1e-5, basis_tol = 1e-10,
    posterior_tol = 1e-12, rank_tol = 1e-11, constraint_tol = 1e-7,
    stationarity_tol = 1e-6, maxit = 100L, warn = TRUE)
  supplied <- list(...)
  if (is.null(names(supplied)) && length(supplied)) stop("Controls must be named.")
  if (any(!names(supplied) %in% names(defaults)) || anyDuplicated(names(supplied)))
    stop("Unknown or duplicate staged-GP control.")
  z <- modifyList(defaults, supplied, keep.null = TRUE)
  for (nm in c("grid_min", "bl_count", "quad_order", "maxit"))
    if (!is.numeric(z[[nm]]) || length(z[[nm]]) != 1L || !is.finite(z[[nm]]) ||
        z[[nm]] < 1 || z[[nm]] != as.integer(z[[nm]])) stop(nm, " must be a positive integer.")
  if (z$grid_min < 9L || z$quad_order < 2L) stop("Require grid_min >= 9 and quad_order >= 2.")
  for (nm in c("bump_eta", "quad_tol", "basis_tol", "posterior_tol", "rank_tol",
               "constraint_tol", "stationarity_tol"))
    if (!is.numeric(z[[nm]]) || length(z[[nm]]) != 1L || !is.finite(z[[nm]]) ||
        z[[nm]] <= 0) stop(nm, " must be positive and finite.")
  for (nm in c("basis_tol", "posterior_tol", "rank_tol"))
    if (z[[nm]] >= 1) stop(nm, " must be smaller than one.")
  for (nm in c("weak_budget", "boundary_points"))
    if (!is.null(z[[nm]]) && (!is.numeric(z[[nm]]) || length(z[[nm]]) != 1L ||
        !is.finite(z[[nm]]) || z[[nm]] < 1 || z[[nm]] != as.integer(z[[nm]])))
      stop(nm, " must be NULL or a positive integer.")
  if (!is.numeric(z$weak_info) || length(z$weak_info) != 1L ||
      !is.finite(z$weak_info) || z$weak_info <= 0 || z$weak_info > 1)
    stop("weak_info must be in (0,1].")
  if (length(z$pilot) != 1L || is.na(z$pilot) || !z$pilot %in% c("gp_mean", "weak_mle"))
    stop("pilot must be gp_mean or weak_mle; phat supplies an external pilot.")
  for (nm in c("include_bl", "warn"))
    if (!is.logical(z[[nm]]) || length(z[[nm]]) != 1L || is.na(z[[nm]]))
      stop(nm, " must be logical.")
  z
}

.jgstage_blockdiag <- function(blocks) {
  out <- matrix(0, sum(vapply(blocks, nrow, integer(1))), sum(vapply(blocks, ncol, integer(1))))
  i <- j <- 0L
  for (b in blocks) {
    out[i + seq_len(nrow(b)), j + seq_len(ncol(b))] <- b
    i <- i + nrow(b); j <- j + ncol(b)
  }
  out
}

# Full latent covariance conditional on the ORIGINAL observations. No nugget
# is added at prediction/optimization points, and no diagonal-only surrogate.
.jgstage_cross <- function(fit, a, b) {
  Ka <- matern52_matrix(a, fit$tt, fit$sigma2, fit$ell)
  Kb <- matern52_matrix(b, fit$tt, fit$sigma2, fit$ell)
  matern52_matrix(a, b, fit$sigma2, fit$ell) -
    crossprod(forwardsolve(t(fit$L), t(Ka)), forwardsolve(t(fit$L), t(Kb)))
}

.jgstage_posterior <- function(fits, grid, tol) {
  components <- lapply(fits, function(fit) {
    S <- .jgstage_cross(fit, grid, grid); S <- (S + t(S))/2
    e <- eigen(S, symmetric = TRUE); scale <- max(e$values)
    if (!is.finite(scale) || scale <= 0 || min(e$values) < -1e-7 * scale)
      stop("GP posterior covariance is not numerically positive semidefinite.")
    keep <- which(e$values > tol * scale)
    list(Sigma = S, basis = sweep(e$vectors[,keep,drop=FALSE], 2, sqrt(e$values[keep]), "*"),
      inverse_factor = sweep(e$vectors[,keep,drop=FALSE], 2, sqrt(e$values[keep]), "/"),
      rank = length(keep), discarded = length(grid) - length(keep),
      discarded_variance = sum(pmax(e$values[-keep], 0)))
  })
  sizes <- vapply(components, `[[`, integer(1), "rank")
  ends <- cumsum(sizes); starts <- c(0L, head(ends, -1L))
  list(fits = fits, grid = grid, components = components,
    indices = Map(function(a,b) seq.int(a+1L,b), starts, ends), rank = sum(sizes),
    B = .jgstage_blockdiag(lapply(components, `[[`, "basis")),
    mean = do.call(cbind, lapply(fits, function(f) gp_predict(f, grid)$mean)))
}

# Conditional-mean representation u(t)=mu(t)+C(t,I) Sigma_I^dagger B z.
# Residual GP innovations outside the retained inducing-state span are not
# optimized; the grid and discarded ranks are therefore explicit diagnostics.
.jgstage_extension <- function(post, tt) {
  list(mean = do.call(cbind, lapply(post$fits, function(f)
      as.vector(matern52_matrix(tt, f$tt, f$sigma2, f$ell) %*% f$alpha))),
    B = lapply(seq_along(post$fits), function(d)
      .jgstage_cross(post$fits[[d]], tt, post$grid) %*% post$components[[d]]$inverse_factor))
}

.jgstage_values <- function(z, ext, post) {
  U <- ext$mean
  for (d in seq_along(post$fits)) U[,d] <- U[,d] + ext$B[[d]] %*% z[post$indices[[d]]]
  U
}

.jgstage_weak <- function(post, control, maps = NULL, order = control$quad_order) {
  span <- diff(range(post$grid))
  radii <- control$weak_radii
  if (is.null(radii)) radii <- span * c(1/32, 1/16, 1/8, 1/4, .4)
  br <- if (is.null(control$bl_radius)) max(radii) else control$bl_radius
  if (!is.numeric(radii) || !length(radii) || any(!is.finite(radii)) ||
      any(radii <= 0 | radii >= span/2) || length(br) != 1L ||
      !is.finite(br) || br <= 0 || br >= span/2) stop("Invalid physical weak radii.")
  design <- .jgp_design_pool(post$grid, list(weak_design_radii = radii,
    weak_design_centers = NULL, bl_radii = br, include_bl = control$include_bl,
    bl_count = control$bl_count))
  # Original observation locations are integration breaks, not extra states.
  breaks <- sort(unique(c(post$grid, post$fits[[1]]$tt)))
  quad <- .jgp_gauss_rule(breaks, order)
  rows <- rbind(design$interior, design$boundary)
  raw <- .jgp_test_rows(quad$tt, rows, 0L, control$bump_eta)
  ni <- nrow(design$interior); nb <- nrow(design$boundary)
  selection <- NULL
  if (is.null(maps)) {
    X <- sweep(raw, 2, sqrt(quad$weights), "*")
    spec <- .jgp_basis_map(X[seq_len(ni),,drop=FALSE], control$basis_tol, TRUE)
    size <- .jgp_design_svd_size(spec$singular_values[spec$retained],
      sum(spec$singular_values), control$weak_info)
    count <- if (is.null(control$weak_budget)) size$count else min(control$weak_budget, nrow(spec$map))
    maps <- list(interior = spec$map[seq_len(count),,drop=FALSE],
      boundary = .jgp_basis_map(X[ni+seq_len(nb),,drop=FALSE], control$basis_tol))
    selection <- list(raw_rows = ni, modes = count,
      information_fraction = sum(spec$singular_values[spec$retained[seq_len(count)]])/sum(spec$singular_values),
      requested_fraction = control$weak_info, explicit_budget = control$weak_budget)
  }
  transform <- function(x) rbind(maps$interior %*% x[seq_len(ni),,drop=FALSE],
    maps$boundary %*% x[ni+seq_len(nb),,drop=FALSE])
  phi <- transform(raw)
  prime <- transform(.jgp_test_rows(quad$tt, rows, 1L, control$bump_eta))
  endpoint <- transform(.jgp_test_rows(range(post$grid), rows, 0L, control$bump_eta))
  list(tt = quad$tt, V = sweep(phi, 2, quad$weights, "*"),
    Vp = sweep(prime, 2, quad$weights, "*"), endpoint = cbind(-endpoint[,1], endpoint[,2]),
    gram = tcrossprod(sweep(phi, 2, sqrt(quad$weights), "*")),
    ni = nrow(maps$interior), nb = nrow(maps$boundary), K = nrow(phi),
    extension = .jgstage_extension(post, quad$tt),
    edge_extension = .jgstage_extension(post, range(post$grid)),
    maps = maps, design = design, selection = selection, order = order)
}

.jgstage_eval <- function(z, p, weak, post, model, jacobian = TRUE) {
  U <- .jgstage_values(z, weak$extension, post)
  edge <- .jgstage_values(z, weak$edge_extension, post)
  q <- nrow(U); K <- weak$K; D <- model$D; J <- model$J
  input <- rbind(matrix(p, J, q), t(U), weak$tt)
  F <- model$jet[[1]](input)
  r <- as.vector(weak$V %*% F + weak$Vp %*% U - weak$endpoint %*% edge)
  Jz <- Jp <- NULL
  if (jacobian) {
    df <- array(model$jet_jac[[1]](input), c(q, D, J+D))
    Jp <- matrix(weak$V %*% matrix(df[,,seq_len(J),drop=FALSE], q, D*J), K*D, J)
    Jz <- matrix(0, K*D, post$rank)
    for (a in seq_len(D)) for (b in seq_len(D)) {
      block <- weak$V %*% sweep(weak$extension$B[[b]], 1, df[,a,J+b], "*")
      if (a == b) block <- block + weak$Vp %*% weak$extension$B[[b]] -
        weak$endpoint %*% weak$edge_extension$B[[b]]
      Jz[(a-1L)*K+seq_len(K), post$indices[[b]]] <- block
    }
  }
  list(r = r, J = Jz, Jp = Jp)
}

.jgstage_svd <- function(A, tol, full = FALSE) {
  if (!nrow(A) || !ncol(A)) return(list(u = matrix(0,nrow(A),0),
    v = if (full) diag(ncol(A)) else matrix(0,ncol(A),0), d = numeric(), rank = 0L))
  s <- svd(A, nu = min(dim(A)), nv = if (full) ncol(A) else min(dim(A)))
  keep <- which(s$d > max(s$d) * tol)
  list(u = s$u[,keep,drop=FALSE], v = if (full) s$v else s$v[,keep,drop=FALSE],
    d = s$d[keep], rank = length(keep))
}

# SQP with identity objective Hessian: min ||z||^2/2 subject to h(z)=0.
# The adaptive L1 merit multiplier below ONLY globalizes the numerical steps;
# it is not an ODE coefficient in the target. Success requires both feasibility
# and projected stationarity. Rank-deficient constraints are never silently
# dropped from the feasibility test, nor is a least-squares failure called MAP.
.jgstage_project <- function(evaluate, start, control) {
  z <- as.numeric(start); history <- list(); reason <- "iteration_limit"
  converged <- FALSE; merit_weight <- 1
  for (iteration in seq_len(control$maxit)) {
    cur <- evaluate(z, TRUE)
    if (any(!is.finite(c(z,cur$r,cur$J)))) { reason <- "nonfinite"; break }
    s <- .jgstage_svd(cur$J, control$rank_tol)
    row_projection <- if (s$rank) as.vector(s$v %*% crossprod(s$v,z)) else z*0
    stationarity <- sqrt(sum((z-row_projection)^2))/(1+sqrt(sum(z^2)))
    feasibility <- if (length(cur$r)) max(abs(cur$r)) else 0
    history[[iteration]] <- data.frame(iteration = iteration, objective = sum(z^2)/2,
      feasibility = feasibility, stationarity = stationarity, rank = s$rank)
    if (feasibility <= control$constraint_tol && stationarity <= control$stationarity_tol) {
      converged <- TRUE; reason <- "feasible_stationary"; break
    }
    rhs <- as.vector(cur$J %*% z) - cur$r
    target <- if (s$rank) as.vector(s$v %*% (crossprod(s$u,rhs)/s$d)) else z*0
    direction <- target-z
    if (sqrt(sum(direction^2)) <= .Machine$double.eps*(1+sqrt(sum(z^2)))) {
      reason <- "no_resolved_step"; break
    }
    multiplier <- if (s$rank) as.vector(-s$u %*% (crossprod(s$v,target)/s$d)) else 0
    merit_weight <- max(merit_weight, 1+max(abs(multiplier)))
    dr <- as.vector(cur$J %*% direction)
    dc <- sum(ifelse(cur$r == 0, abs(dr), sign(cur$r)*dr))
    slope <- sum(z*direction)+merit_weight*dc
    if (slope >= 0 && dc < 0) {
      merit_weight <- max(merit_weight, 1+2*abs(sum(z*direction))/(-dc))
      slope <- sum(z*direction)+merit_weight*dc
    }
    if (slope >= 0) { reason <- "unresolved_constraints"; break }
    merit <- sum(z^2)/2 + merit_weight*sum(abs(cur$r))
    accepted <- FALSE
    for (backtrack in 0:24) {
      step <- 2^(-backtrack); trial <- z+step*direction
      rr <- tryCatch(evaluate(trial, FALSE)$r, error = function(e) NA_real_)
      if (all(is.finite(rr)) && sum(trial^2)/2 + merit_weight*sum(abs(rr)) <=
          merit + 1e-4*step*slope) { z <- trial; accepted <- TRUE; break }
    }
    if (!accepted) { reason <- "line_search_failed"; break }
  }
  final <- evaluate(z, TRUE); s <- .jgstage_svd(final$J, control$rank_tol)
  projected <- if (s$rank) as.vector(s$v %*% crossprod(s$v,z)) else z*0
  feasibility <- if (length(final$r)) max(abs(final$r)) else 0
  stationarity <- sqrt(sum((z-projected)^2))/(1+sqrt(sum(z^2)))
  converged <- all(is.finite(c(final$r,z))) && feasibility <= control$constraint_tol &&
    stationarity <= control$stationarity_tol
  if (converged) reason <- "feasible_stationary"
  list(z = z, converged = converged, reason = reason, iterations = length(history),
    objective = sum(z^2)/2, feasibility = feasibility, stationarity = stationarity,
    constraint_rank = s$rank, constraints = length(final$r),
    history = if (length(history)) do.call(rbind,history) else data.frame())
}

# In white posterior coordinates, conditioning A z=b gives a minimum-norm
# center zc and an orthonormal null-space factor Z. This is exactly the Gaussian
# Schur-complement conditional, without cancellation in Sigma_BB-Sigma_BI... .
.jgstage_condition <- function(A, b, tol) {
  s <- .jgstage_svd(A, tol, full = TRUE)
  active <- seq_len(s$rank)
  center <- if (s$rank) as.vector(s$v[,active,drop=FALSE] %*% (crossprod(s$u,b)/s$d)) else numeric(ncol(A))
  null <- s$v[,setdiff(seq_len(ncol(A)),active),drop=FALSE]
  list(center = center, basis = null, rank = s$rank,
    compatibility = if (length(b)) max(abs(A %*% center-b)) else 0)
}

.jgstage_noise <- function(Y, tt) {
  if (nrow(Y) >= 20L) return(estimate_std(Y, k = 6L))
  if (nrow(Y) <= 4L) stop("Supply noise_sd when there are at most four observations.")
  # Prespecified sparse-data trend diagnostic, not an ODE-derived trajectory.
  X <- outer((tt-min(tt))/diff(range(tt)), 0:3, "^")
  vapply(seq_len(ncol(Y)), function(d) sqrt(sum(lm.fit(X,Y[,d])$residuals^2)/(nrow(Y)-4L)), numeric(1))
}

# Experimental, deliberately not exported/installed over the production API.
# phat: an externally fitted parameter vector, frozen in BOTH correction stages.
# Otherwise pilot="gp_mean" fits interior moments at the continuous GP mean;
# pilot="weak_mle" uses original weak MLE on actual data, WITHOUT interpolation.
# The latter is currently scalar-component only; supply phat for other baselines.
solveWendyGPStaged <- function(f, Y, tt, p0, phat = NULL, noise_sd = NULL, control = list()) {
  control <- do.call(wendygp_staged_control, control)
  Y <- as.matrix(Y); tt <- as.numeric(tt); p0 <- as.numeric(p0)
  if (!is.numeric(Y) || nrow(Y) != length(tt) || ncol(Y) < 1L || length(tt) < 5L ||
      any(!is.finite(Y)) || any(!is.finite(tt)) || any(diff(tt) <= 0) ||
      !length(p0) || any(!is.finite(p0))) stop("Require complete observed components, increasing times and finite p0.")
  if (!is.null(phat) && (!is.numeric(phat) || length(phat) != length(p0) || any(!is.finite(phat))))
    stop("phat must be a finite parameter vector of the same length as p0.")
  noise_source <- if (is.null(noise_sd)) "independent_data_estimate" else "supplied"
  if (is.null(noise_sd)) noise_sd <- .jgstage_noise(Y,tt)
  if (!is.numeric(noise_sd) || !length(noise_sd) %in% c(1L,ncol(Y)) ||
      any(!is.finite(noise_sd)) || any(noise_sd <= 0)) stop("noise_sd must be positive, scalar or per component.")
  noise_sd <- rep_len(noise_sd,ncol(Y))
  grid <- control$grid
  if (is.null(grid)) grid <- seq(min(tt),max(tt),length.out=max(control$grid_min,length(tt)))
  if (!is.numeric(grid) || length(grid) < 9L || any(!is.finite(grid)) || any(diff(grid) <= 0) ||
      !isTRUE(all.equal(range(grid),range(tt)))) stop("grid must be increasing and span the observation interval.")
  # Exactly one marginal-likelihood GP fit per observed component.
  fits <- lapply(seq_len(ncol(Y)), function(d) gp_fit_1d(tt,Y[,d],sigma2_n=noise_sd[d]^2))
  post <- .jgstage_posterior(fits,grid,control$posterior_tol)
  weak <- .jgstage_weak(post,control)
  fine <- .jgstage_weak(post,control,weak$maps,2L*control$quad_order)
  model <- .jgp_model(f,ncol(Y),length(p0),0L)
  ii <- .jgp_design_rows(seq_len(weak$ni),weak$K,model$D)
  z0 <- numeric(post$rank); pilot_fit <- NULL; pilot_converged <- NA
  method <- if (is.null(phat)) control$pilot else "supplied"
  if (method == "weak_mle") {
    if (ncol(Y) != 1L) stop("For a multi-component original-WENDy pilot, supply phat explicitly.")
    pilot_fit <- solveWendy(f,Y,tt,p0=p0,method="MLE",control=list(noise_sd=noise_sd[1],
      max_points_interp=0L,interpolation_method=NULL,include_boundary_layer=FALSE,
      test_fun_type="MSG",estimate_IC=FALSE,estimate_trajectory=FALSE))
    phat <- as.numeric(pilot_fit$phat)
    pilot_converged <- isTRUE(pilot_fit$data$converged)
  } else if (method == "gp_mean") {
    ev <- function(p) .jgstage_eval(z0,p,weak,post,model,TRUE)
    if (model$affine) {
      v <- ev(p0); G <- v$Jp[ii,,drop=FALSE]
      if (qr(G)$rank < length(p0)) stop("GP-mean interior parameter pilot is rank deficient.")
      phat <- as.vector(qr.solve(G,as.vector(G %*% p0)-v$r[ii]))
      pilot_converged <- TRUE
    } else {
      pilot_fit <- minpack.lm::nls.lm(p0, fn=function(p)ev(p)$r[ii],
        jac=function(p)ev(p)$Jp[ii,,drop=FALSE],
        control=minpack.lm::nls.lm.control(maxiter=control$maxit,maxfev=10000L))
      phat <- as.numeric(pilot_fit$par)
      pilot_converged <- pilot_fit$info %in% 1:4
    }
  }
  if (length(phat) != length(p0) || any(!is.finite(phat))) stop("Nonfinite parameter pilot.")
  # Fixed, data-derived unit conversion ONLY for numerical feasibility tests.
  component_scale <- pmax(sqrt(colMeans(Y^2)),noise_sd)/sqrt(diff(range(tt)))
  row_scale <- rep(component_scale, each=weak$K)
  evaluate <- function(z, jacobian=TRUE, rows=seq_len(weak$K*model$D)) {
    v <- .jgstage_eval(z,phat,weak,post,model,jacobian)
    list(r=v$r[rows]/row_scale[rows],
      J=if(jacobian)sweep(v$J[rows,,drop=FALSE],1,row_scale[rows],"/") else NULL)
  }
  interior <- .jgstage_project(function(z,jacobian)evaluate(z,jacobian,ii),z0,control)
  z_int <- interior$z; z <- z_int
  boundary <- list(converged=NA,reason="disabled",applied=FALSE)
  fixed_rows <- integer(); rc <- 0L
  if (control$include_bl && interior$converged) {
    rc <- control$boundary_points
    if (is.null(rc)) rc <- sum(grid < min(grid)+weak$design$boundary_radii)
    rc <- min(rc, floor((length(grid)-1L)/2L))
    local_fixed <- seq.int(rc+1L,length(grid)-rc)
    fixed_rows <- unlist(lapply(seq_len(model$D),function(d)(d-1L)*length(grid)+local_fixed))
    A <- post$B[fixed_rows,,drop=FALSE]
    conditioned <- .jgstage_condition(A,as.vector(A %*% z_int),
      max(dim(A))*.Machine$double.eps)
    Z <- conditioned$basis; center <- conditioned$center
    boundary <- .jgstage_project(function(w,jacobian) {
      v <- evaluate(as.vector(center+Z %*% w),jacobian)
      if(jacobian)v$J <- v$J %*% Z
      v
    },as.vector(crossprod(Z,z_int-center)),control)
    boundary$conditioning <- conditioned
    boundary$candidate <- as.vector(center+Z %*% boundary$z)
    boundary$fixed_interior_error <- max(abs(A %*% (boundary$candidate-z_int)))
    boundary$applied <- boundary$converged && boundary$fixed_interior_error <=
      control$constraint_tol*max(noise_sd)
    if (boundary$applied) z <- boundary$candidate
  } else if (control$include_bl) boundary$reason <- "interior_not_converged"
  check <- function(zz) {
    a <- .jgstage_eval(zz,phat,weak,post,model,TRUE)
    b <- .jgstage_eval(zz,phat,fine,post,model,TRUE)
    residual_error <- max(abs(a$r-b$r)/row_scale)
    jacobian_error <- sqrt(sum((a$J-b$J)^2))/max(1,sqrt(sum(b$J^2)))
    W <- .jgp_whitener(weak$gram,control$basis_tol)$W
    gram_error <- max(abs(W %*% (fine$gram-weak$gram) %*% t(W)))
    list(residual_error=residual_error,jacobian_error=jacobian_error,gram_error=gram_error,
      fine_feasibility=max(abs(b$r)/row_scale),
      fine_interior_feasibility=max(abs(b$r[ii])/row_scale[ii]),
      passed=max(residual_error,jacobian_error,gram_error)<=control$quad_tol)
  }
  accuracy <- list(initial=check(z0),interior=check(z_int),final=check(z))
  grid_values <- function(zz)matrix(as.vector(post$mean)+post$B %*% zz,nrow=length(grid))
  continuous_grid <- .jgstage_values(z,.jgstage_extension(post,grid),post)
  representation_error <- max(abs(sweep(continuous_grid-grid_values(z),2,noise_sd,"/")))
  frozen_error <- if (length(fixed_rows)) max(abs(as.vector(continuous_grid-
    .jgstage_values(z_int,.jgstage_extension(post,grid),post))[fixed_rows])) else 0
  converged <- interior$converged && (!control$include_bl || isTRUE(boundary$applied)) &&
    all(vapply(accuracy,`[[`,logical(1),"passed")) &&
    representation_error <= control$constraint_tol &&
    frozen_error <= control$constraint_tol*max(noise_sd) && !isFALSE(pilot_converged)
  if (control$warn && !converged)
    warning("Staged GP projection is not validated numerically; inspect stages and diagnostics. No grid refinement or soft-penalty fallback was used.",call.=FALSE)
  structure(list(phat=as.numeric(phat),tt=grid,tt_obs=tt,Y=Y,noise_sd=noise_sd,
    Uhat=grid_values(z),U_obs=.jgstage_values(z,.jgstage_extension(post,tt),post),
    initial=list(U=post$mean,z=z0,phat=as.numeric(phat)),z=z,
    stages=list(interior=interior,boundary=boundary),
    posterior=post,weak=weak,weak_check=fine,model=model,control=control,
    pilot=list(method=method,fit=pilot_fit,converged=pilot_converged),converged=converged,
    diagnostics=list(gp_fits=length(fits),gp_training_points=vapply(fits,function(g)length(g$tt),integer(1)),
      noise_source=noise_source,state_points=length(grid),state_rank=post$rank,
      quadrature_points=length(weak$tt),check_points=length(fine$tt),accuracy=accuracy,
      fixed_interior_rows=fixed_rows,boundary_points=rc,parameter_change=0,
      representation_error_noise_units=representation_error,continuous_frozen_error=frozen_error,
      covariance_scope="original-data posterior; retained inducing-state conditional-mean span",
      uncertainty_scope="GP uncertainty is pre-ODE and conditional on fitted hyperparameters; no corrected-state intervals")),
    class="wendy_gp_staged")
}

predictWendyGPStaged <- function(object, tt, stage=c("final","interior","gp")) {
  stage <- match.arg(stage); tt <- as.numeric(tt)
  if (!inherits(object,"wendy_gp_staged") || any(!is.finite(tt)) ||
      any(tt < min(object$tt) | tt > max(object$tt))) stop("Prediction times must lie in the fitted interval.")
  z <- switch(stage,final=object$z,interior=object$stages$interior$z,gp=object$initial$z)
  .jgstage_values(z,.jgstage_extension(object$posterior,tt),object$posterior)
}
