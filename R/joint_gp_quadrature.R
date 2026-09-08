# Experimental continuous-GP weak integration. Extra quadrature nodes never
# become state variables. No ODE-derived endpoint jets or covariance ridge.
.jgp_gauss_rule <- function(grid, order) {
  j <- seq_len(order - 1L); b <- j / sqrt(4*j^2 - 1)
  A <- matrix(0, order, order); A[cbind(j,j+1L)] <- b; A <- A + t(A)
  e <- eigen(A, symmetric=TRUE); ix <- order(e$values)
  nodes <- e$values[ix]; w <- 2*e$vectors[1,ix]^2 / sum(e$vectors[1,ix]^2)
  h <- diff(grid); mid <- (head(grid,-1L)+tail(grid,-1L))/2
  list(tt=as.vector(t(outer(mid,rep(1,order))+outer(h/2,nodes))),
       weights=as.vector(t(outer(h/2,w))), order=order)
}

.jgp_gauss_raw <- function(tt, design, derivative, eta) {
  rbind(.jgp_test_rows(tt,design$interior,derivative,eta),
        .jgp_test_rows(tt,design$boundary,derivative,eta))
}

# Stencil WIDTH is tied to the fitted GP length scale, not fixed: a stencil
# spanning more than about one length scale is exactly what makes high-order
# interpolation diverge, and the same rule then adapts to grid and problem with
# no user knob. The working grid is uniform, so h is exact.
.jgp_lagrange_order <- function(fit, grid, tt) {
  h <- mean(diff(grid)); m <- length(grid)
  # Each integration cell must carry one polynomial. Choosing the width at
  # each quadrature node introduces stencil switches inside a cell, hence
  # discontinuities that an unsplit Gauss rule resolves poorly.
  mid <- (head(grid,-1L)+tail(grid,-1L))/2
  r <- .jgp_radius((mid-fit$origin)/fit$span, fit$radius_coef, fit$radius_bounds)$r*fit$span
  width <- pmin(2L*pmax(2L, pmin(8L, as.integer(round(r/h/2)))), m)
  width[findInterval(tt,grid,rightmost.closed=TRUE,all.inside=TRUE)]
}

# Centered-stencil Lagrange rows. Weights sum to one, so the constant mean
# offset the GP extension needs is a no-op here and the code shape is shared.
.jgp_lagrange_rows <- function(grid,tt,k) {
  m <- length(grid); E <- matrix(0,length(tt),m)
  interval <- findInterval(tt,grid,rightmost.closed=TRUE)
  # Batch rows with the same stencil width. Each row uses the original product
  # formula, including exact interpolation at grid nodes (no barycentric 0/0).
  for (width in unique(k)) {
    rows <- which(k==width)
    lo <- pmin(pmax(interval[rows]-width%/%2L+1L,1L),m-width+1L)
    ix <- outer(lo,seq_len(width)-1L,"+")
    x <- matrix(grid[ix],nrow=length(rows))
    for (a in seq_len(width)) {
      value <- rep(1,length(rows))
      for (b in setdiff(seq_len(width),a))
        value <- value*((tt[rows]-x[,b])/(x[,a]-x[,b]))
      E[cbind(rows,ix[,a])] <- value
    }
  }
  E
}

.jgp_lagrange_extension <- function(fits, grid, tt) {
  lapply(fits, function(f) .jgp_lagrange_rows(grid, tt, .jgp_lagrange_order(f, grid, tt)))
}

# Per-component choice. Both are fixed representations in the weak operator;
# neither an invertible coordinate change nor an interpolant guarantees that
# an unobserved grid state is identifiable.
.jgp_extension <- function(fits, state, grid, tt, control) {
  if (!length(control$weak_extension) %in% c(1L, length(fits)))
    stop("weak_extension must be scalar or have one entry per component.")
  kind <- rep_len(control$weak_extension, length(fits))
  lapply(seq_along(fits), function(d)
    if (identical(kind[d], "lagrange"))
      .jgp_lagrange_rows(grid, tt, .jgp_lagrange_order(fits[[d]], grid, tt))
    else t(.jgp_solve(t(state$L[[d]]), t(.jgp_gp_cross(fits[[d]], tt, grid)))))
}

# Conditional sd of the state at the quadrature nodes given the grid, under the
# fitted GP. This measures whether the GRID resolves the state between nodes --
# a property of grid and kernel, not of whichever interpolant carries it there --
# so it is useful context for both extensions. It is NOT an error bound for
# the optimized curve or the Lagrange interpolant. k(t,t)=tau2 for both kernels.
.jgp_gauss_resolution <- function(fits, state, grid, tt) {
  vapply(seq_along(fits), function(d) {
    C <- .jgp_gp_cross(fits[[d]], tt, grid)
    S <- .jgp_solve(t(state$L[[d]]), t(C))
    sqrt(max(pmax(fits[[d]]$tau2 - colSums(t(C)*S), 0)))
  }, numeric(1))
}

.jgp_gauss_extension <- function(fits, state, grid, tt) {
  lapply(seq_along(fits),function(d) {
    cross <- .jgp_gp_cross(fits[[d]],tt,grid)
    t(.jgp_solve(t(state$L[[d]]),t(cross)))
  })
}

.jgp_gauss_weak <- function(grid, design, fits, state, control, maps=NULL,
                             order=control$weak_quad_order,whitened=TRUE,resolution=TRUE) {
  quad <- .jgp_gauss_rule(grid,order)
  raw <- .jgp_gauss_raw(quad$tt,design,0L,control$bump_eta)
  nint <- nrow(design$interior); nbraw <- nrow(design$boundary)
  spectrum <- NULL
  if (is.null(maps)) {
    X <- sweep(raw,2,sqrt(quad$weights),"*")
    spectrum <- .jgp_basis_map(X[seq_len(nint),,drop=FALSE],control$basis_tol,TRUE)
    size <- .jgp_design_svd_size(spectrum$singular_values[spectrum$retained],
      sum(spectrum$singular_values),control$weak_design_info)
    count <- if (!is.null(control$weak_radii)) nrow(spectrum$map) else
      if (is.null(control$weak_design_budget)) size$count else
      min(control$weak_design_budget,nrow(spectrum$map))
    if (count < 1L) stop("Continuous weak test pool has no retained interior mode.")
    bmap <- if (nbraw) .jgp_basis_map(X[nint+seq_len(nbraw),,drop=FALSE],control$basis_tol) else matrix(0,0,0)
    maps <- list(interior=spectrum$map[seq_len(count),,drop=FALSE],boundary=bmap)
  }
  ni <- nrow(maps$interior); nb <- nrow(maps$boundary); K <- ni+nb
  transform <- function(x)rbind(maps$interior%*%x[seq_len(nint),,drop=FALSE],
    maps$boundary%*%x[nint+seq_len(nbraw),,drop=FALSE])
  phi <- transform(raw)
  prime <- transform(.jgp_gauss_raw(quad$tt,design,1L,control$bump_eta))
  endpoints <- transform(.jgp_gauss_raw(range(grid),design,0L,control$bump_eta))
  E <- .jgp_extension(fits,state,grid,quad$tt,control)
  Eb <- .jgp_extension(fits,state,grid,range(grid),control)
  w <- list(integration="gp_gauss",tt=grid,quad_tt=quad$tt,quad_weights=quad$weights,
    quad_order=order,ni=ni,nb=nb,K=K,design=design,maps=maps,
    V=sweep(phi,2,quad$weights,"*"),Vp=sweep(prime,2,quad$weights,"*"),
    B=cbind(-endpoints[,1],endpoints[,2]),
    gram=tcrossprod(sweep(phi,2,sqrt(quad$weights),"*")),
    E=E,Eb=Eb,mean=vapply(fits,`[[`,numeric(1),"mean"),em_order=0L,
    extension_kind=control$weak_extension,
    resolved=if (resolution) .jgp_gauss_resolution(fits,state,grid,quad$tt) else NULL)
  w$linear <- lapply(seq_along(fits),function(d)w$Vp%*%E[[d]]-w$B%*%Eb[[d]])
  # Fixed for the whole solve; consumed by the constant df/du blocks. The
  # L-whitened variants are the same factors right-multiplied by the prior
  # Cholesky, so a whitened state Jacobian assembles at identical cost and the
  # optimizer never pays a separate whitening pass. EL = C K^-1 L = C L^-T.
  w$VE <- lapply(E,function(x)w$V%*%x)
  if (whitened) {
  w$EL <- lapply(seq_along(fits),function(d)E[[d]]%*%state$L[[d]])
  w$VEL <- lapply(seq_along(fits),function(d)w$VE[[d]]%*%state$L[[d]])
  w$linearL <- lapply(seq_along(fits),function(d)w$linear[[d]]%*%state$L[[d]])
  }
  if (!is.null(spectrum)) {
    s <- spectrum$singular_values[spectrum$retained]; total <- sum(spectrum$singular_values)
    fraction <- sum(head(s,ni))/total
    w$selection <- list(method=if(is.null(control$weak_radii))"svd" else "explicit",selected=seq_len(ni),budget=ni,
      size_rule=if(!is.null(control$weak_radii))"full_numerical_span" else
        if(is.null(control$weak_design_budget))"singular_value_fraction" else "fixed_budget",
      singular_values=s,singular_total=total,information_fraction=fraction,
      information_target=control$weak_design_info,
      information_target_met=fraction+32*.Machine$double.eps>=control$weak_design_info,
      integration="gp_gauss",screening=list(raw_count=nint,mode_retained=ni,
        center_placement=design$radius_selection$center_placement,
        passed=fraction+32*.Machine$double.eps>=control$weak_design_info))
  }
  w
}

# Hold the selected tests and optimized grid fixed. The companion changes only
# integration nodes, and is reused at the initial and final physical states.
.jgp_gauss_pair <- function(grid,design,fits,state,control,maps=NULL,
                            order=control$weak_quad_order) {
  weak <- .jgp_gauss_weak(grid,design,fits,state,control,maps,order)
  weak$gram_whitener <- .jgp_whitener(weak$gram,control$covariance_tol)$W
  weak$fine <- .jgp_gauss_weak(grid,design,fits,state,control,weak$maps,2L*order,
    whitened=FALSE,resolution=FALSE)
  weak$gram_error <- norm(weak$gram_whitener%*%(weak$fine$gram-weak$gram)%*%t(weak$gram_whitener),"2")
  weak$grid_extension <- .jgp_extension(fits,state,grid,grid,control)
  weak
}

.jgp_gauss_eval <- function(U,p,weak,model,jacobian=TRUE,whiten=FALSE,adjoint=FALSE) {
  m <- nrow(U); D <- model$D; J <- model$J; K <- weak$K; q <- length(weak$quad_tt)
  # whiten selects the coordinates of the STATE Jacobian only. U is physical
  # either way, so the residual always uses the physical extension.
  gE <- if (whiten) weak$EL else weak$E
  gVE <- if (whiten) weak$VEL else weak$VE
  gLin <- if (whiten) weak$linearL else weak$linear
  uq <- matrix(0,q,D); ub <- matrix(0,2L,D)
  for (d in seq_len(D)) {
    centered <- U[,d]-weak$mean[d]
    uq[,d] <- weak$mean[d]+weak$E[[d]]%*%centered
    ub[,d] <- weak$mean[d]+weak$Eb[[d]]%*%centered
  }
  input <- rbind(matrix(p,J,q),t(uq),weak$quad_tt)
  F <- model$jet[[1]](input)
  r <- weak$V%*%F+weak$Vp%*%uq-weak$B%*%ub
  Jp <- Ju <- pullback <- NULL
  if (jacobian) {
    df <- array(model$jet_jac[[1]](input),c(q,D,J+D))
    if (adjoint) {
      # Apply the transposed weak Jacobian directly. For each state b this is
      # E_b' sum_a(df_a/du_b * V' score_a) + linear_b' score_b.
      # No K-by-m derivative blocks or V diag(df) E products are needed.
      pullback <- function(score) {
        score <- matrix(score,K,D)
        force <- crossprod(weak$V,score)
        gp <- as.vector(crossprod(matrix(df[,,seq_len(J),drop=FALSE],q*D,J),
          as.vector(force)))
        gu <- matrix(0,m,D)
        for (b in seq_len(D)) {
          force_b <- rowSums(matrix(df[,,J+b],q,D)*force)
          gu[,b] <- crossprod(gE[[b]],force_b)+crossprod(gLin[[b]],score[,b])
        }
        list(p=gp,U=as.vector(gu))
      }
    } else {
      Vt <- t(weak$V)
      Jp <- matrix(weak$V%*%matrix(df[,,seq_len(J),drop=FALSE],q,D*J),K*D,J)
      Ju <- matrix(0,K*D,m*D)
      # The dense case is V diag(w) E. Scaling rows of V' uses R's column-major
      # recycling directly, avoiding sweep's array permutations. Most blocks
      # never need the dense product. model$ju_zero blocks contribute nothing, and
      # model$ju_const blocks carry one scalar, so they reduce to a multiple of
      # the precomputed weak$VE. Only genuinely state- or time-varying blocks pay
      # the full product; scalar fits use the adjoint above instead.
      for (a in seq_len(D)) for (b in seq_len(D)) {
        if (model$ju_zero[a,b] && a!=b) next
        ir <- (a-1L)*K+seq_len(K); ic <- (b-1L)*m+seq_len(m)
        block <- if (model$ju_zero[a,b]) 0 else
          if (model$ju_const[a,b]) df[1L,a,J+b]*gVE[[b]] else
          crossprod(Vt*df[,a,J+b],gE[[b]])
        Ju[ir,ic] <- if (a==b) block+gLin[[b]] else block
      }
    }
  }
  list(r=as.vector(r),Jp=Jp,Ju=Ju,pullback=pullback)
}

.jgp_quadrature_message <- function(check,control) {
  values <- c(residual=max(check$weak),Jacobian=check$jacobian_relative,Gram=check$gram_error)
  tolerances <- c(control$weak_grid_tol,control$weak_design_quad_tol,control$weak_quad_gram_tol)
  failed <- !is.finite(values) | values>tolerances
  details <- sprintf("%s discrepancy %.3g (tolerance %.3g)",names(values)[failed],
    values[failed],tolerances[failed])
  paste0("Final weak quadrature check failed: ",paste(details,collapse="; "),
    ". Inspect fit$diagnostics$final_grid.",
    if (control$weak_quad_order<32L)
      " A higher control$weak_quad_order refines integration without adding state variables." else "")
}

.jgp_gauss_accuracy <- function(U,p,state,weak,model,fits,weights,H,control,
                                observed=seq_len(model$D), observation_times=NULL) {
  interp <- setNames(numeric(length(observed)), observed)
  stabilization <- extension <- numeric(model$D)
  exact <- all(rowSums(H!=0)==1L)
  # Use the representation stored with the objective, not a caller's default.
  control$weak_extension <- weak$extension_kind
  Eg <- weak$grid_extension
  if (is.null(Eg)) Eg <- .jgp_extension(fits,state,weak$tt,weak$tt,control)
  for(d in seq_len(model$D)) {
    centered <- U[,d]-weak$mean[d]
    stabilization[d] <- max(abs(U[,d]-(weak$mean[d]+Eg[[d]]%*%centered)))
    amplitude <- sqrt(mean(U[,d]^2))
    extension[d] <- if (amplitude > 0) weak$resolved[d]/amplitude else Inf
  }
  if (!exact) for (k in seq_along(observed)) {
    d <- observed[k]
    times <- if (is.null(observation_times)) fits[[d]]$tt else observation_times
    if (length(times) != nrow(H)) stop("Observation times must match H's rows.")
    dc <- control; dc$weak_extension <- rep_len(weak$extension_kind,model$D)[d]
    Eo <- .jgp_extension(fits[d],list(L=state$L[d]),weak$tt,times,dc)[[1]]
    interp[k] <- max(abs(weak$mean[d]+Eo%*%(U[,d]-weak$mean[d])-H%*%U[,d])) /
      sqrt(fits[[d]]$noise2)
  }
  fine <- weak$fine
  if (is.null(fine)) fine <- .jgp_gauss_weak(weak$tt,weak$design,fits,state,
    control,weak$maps,2L*weak$quad_order)
  coarse_eval <- .jgp_gauss_eval(U,p,weak,model)
  fine_eval <- .jgp_gauss_eval(U,p,fine,model)
  delta <- coarse_eval$r-fine_eval$r
  rms <- function(W) if (nrow(W)) sqrt(mean((W%*%delta)^2)) else 0
  unweighted <- c(interior=rms(weights$Wi),boundary=rms(weights$Wb))
  weighted <- sqrt(control$lambda)*unweighted
  # Compare derivatives in declared, fixed physical-coordinate scales. This
  # does not depend on either solver's optimization preconditioner.
  ps <- if (is.null(weak$parameter_scale)) rep(1,length(p)) else weak$parameter_scale
  us <- if (is.null(weak$state_scale)) sqrt(vapply(fits,`[[`,numeric(1),"tau2")) else weak$state_scale
  scales <- c(ps,rep(us,each=nrow(U)))
  jc <- sweep(weights$W%*%cbind(coarse_eval$Jp,coarse_eval$Ju),2,scales,"*")
  jf <- sweep(weights$W%*%cbind(fine_eval$Jp,fine_eval$Ju),2,scales,"*")
  jacobian_relative <- norm(jc-jf,"F")/max(norm(jf,"F"),.Machine$double.eps)
  GW <- weak$gram_whitener
  if (is.null(GW)) GW <- .jgp_whitener(weak$gram,control$covariance_tol)$W
  gram_error <- weak$gram_error
  if (is.null(gram_error)) gram_error <- norm(GW%*%(fine$gram-weak$gram)%*%t(GW),"2")
  quadrature_passed <- all(is.finite(c(weighted,jacobian_relative,gram_error))) &&
    max(weighted)<=control$weak_grid_tol &&
    jacobian_relative<=control$weak_design_quad_tol && gram_error<=control$weak_quad_gram_tol
  extension_passed <- all(is.finite(extension)) && max(extension)<=control$weak_extension_tol
  list(interpolation=interp,stabilization_state=stabilization,
    observed=observed,weak=weighted,weak_unweighted=unweighted,
    jacobian_relative=jacobian_relative,gram_error=gram_error,
    quadrature_passed=quadrature_passed,
    extension=extension,extension_tol=control$weak_extension_tol,
    extension_passed=extension_passed,
    extension_interpretation="GP conditional uncertainty / fitted-state RMS; advisory, not an integration error",
    weak_passed=quadrature_passed,
    passed=all(is.finite(interp)) && max(interp)<=control$grid_tol,
    quad_order=weak$quad_order,check_order=fine$quad_order,
    check_points=length(fine$quad_tt),
    quadrature_points=length(weak$quad_tt),state_variables=length(U))
}

.jgp_gauss_solve <- function(f,Y,tt,p0,fits,control,lower,upper) {
  grid <- .jgp_grid(tt,control);state <- .jgp_state(fits,grid,control)
  H <- .jgp_H(grid,tt);D <- ncol(Y);J <- length(p0)
  model <- .jgp_model(f,D,J,0L)
  design <- if(is.null(control$weak_radii)) .jgp_design_pool(grid,control) else
    .jgp_test_design(tt,fits,control)
  weak <- .jgp_gauss_pair(grid,design,fits,state,control)
  weak$parameter_scale <- if(is.null(control$weak_design_scale))pmax(abs(p0),1) else control$weak_design_scale
  if(length(weak$parameter_scale)!=J)stop("weak_design_scale must have one entry per parameter.")
  pilot <- .jgp_initialize(state$mean,p0,weak,model,state$Sigma,control,lower,upper)
  base <- .jgp_weak_eval(state$mean,pilot,weak,model)
  Omega <- .jgp_ode_metric(weak,D,control,base$Ju,state$Sigma)
  weights <- .jgp_fit_metric(Omega,weak,Y,tt,control)$weights
  initial <- .jgp_gauss_accuracy(state$mean,pilot,state,weak,model,fits,weights,H,control)
  spec <- list(weak=weak,Omega=Omega,diagnostics=weak$selection)
  spec$diagnostics$screening$quadrature_passed <- initial$quadrature_passed
  spec$diagnostics$screening$passed <- spec$diagnostics$screening$passed && initial$quadrature_passed
  fit <- .jgp_design_fit(Y,tt,fits,model,control,lower,upper,state,H,pilot,spec,
                         initial$weak)
  interp <- vapply(fits,.jgp_interpolation_error,numeric(1),tt=grid,obs=tt,H=H)
  fit$weak_integration <- "gp_gauss"
  fit$diagnostics$weak_integration <- "gp_gauss"
  fit$diagnostics$weak_grid_passed <- fit$diagnostics$final_grid$quadrature_passed
  fit$diagnostics$quadrature_passed <- fit$diagnostics$final_grid$quadrature_passed
  fit$diagnostics$extension_passed <- fit$diagnostics$final_grid$extension_passed
  fit$diagnostics$grid_passed <- max(interp)<=control$grid_tol && fit$diagnostics$grid_passed
  fit$diagnostics$initial_quadrature <- initial
  fit$diagnostics$initial_interpolation <- interp
  fit$diagnostics$quadrature_points <- length(weak$quad_tt)
  fit$diagnostics$state_variables <- length(state$mu)
  fit$diagnostics$em_order_applied <- 0L
  fit$diagnostics$gp_delta_scope <- if(control$ode_weighting=="gp_delta")"frozen_grid_covariance" else "not_propagated"
  if(!fit$diagnostics$grid_passed) {
    msg <- "Observation-operator error exceeds grid_tol on the fixed working grid."
    if(control$grid_action=="error")stop(msg) else warning(msg,call.=FALSE)
  }
  if(!fit$diagnostics$quadrature_passed)
    warning(.jgp_quadrature_message(fit$diagnostics$final_grid,control),call.=FALSE)
  if(!isTRUE(weak$selection$information_target_met))
    warning("Requested weak singular-value fraction was not reached; inspect radius_selection.",call.=FALSE)
  fit
}
