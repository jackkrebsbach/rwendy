# Final diagnostics use the shared residual engine in physical (p,U)
# coordinates, independently of each optimization preconditioner.
.jgp_fit_diagnostics <- function(U, p, Y, H, noise, weak, model, weights, control,
                                 prior_state = NULL, prior_components = integer(),
                                 observed = seq_len(model$D),
                                 lower = rep(-Inf, length(p)),
                                 upper = rep(Inf, length(p)),
                                 tolerance = control$gtol, rank = TRUE,
                                 state_bounds = attr(control,"state_bounds")) {
  m <- nrow(U); D <- ncol(U); J <- length(p); n <- nrow(Y)
  evaluate <- .jgp_joint_objective(Y,H,weak,model,weights,noise,control$lambda,
    coordinates=list(mu=matrix(0,m,D),L=rep(list(1),D),order=seq_len(D),pscale=rep(1,J)),
    prior_state=prior_state,prior_components=prior_components,observed=observed)
  final <- evaluate(c(p,as.vector(U)))
  residual <- final$r; jacobian <- final$J
  gradient <- as.vector(crossprod(jacobian,residual))
  theta <- c(p,as.vector(U))
  if (is.null(state_bounds)) state_bounds <- .jgp_bound_pair(NULL,NULL,D,"state")
  lo <- c(lower,rep(state_bounds$lower,each=m)); hi <- c(upper,rep(state_bounds$upper,each=m))
  active_tolerance <- 32*.Machine$double.eps*pmax(1,abs(theta),
    ifelse(is.finite(lo),abs(lo),0),ifelse(is.finite(hi),abs(hi),0))
  projected <- gradient
  projected[theta<=lo+active_tolerance & gradient>0 | theta>=hi-active_tolerance & gradient<0] <- 0
  norms <- sqrt(colSums(jacobian^2))
  scaled <- max(abs(projected)/pmax(norms,1e-8))
  finite <- all(is.finite(c(residual,gradient,norms))) && is.finite(sum(residual^2))
  local <- NULL
  if (rank && finite) {
    singular <- svd(sweep(jacobian,2,pmax(norms,1e-8),"/"),nu=0,nv=0)$d
    cutoff <- sqrt(.Machine$double.eps)
    retained <- if (length(singular) && max(singular)>0) sum(singular>max(singular)*cutoff) else 0L
    local <- list(rank=retained,variables=ncol(jacobian),nullity=ncol(jacobian)-retained,
      singular_values=singular,relative_tolerance=cutoff,
      interpretation="local joint residual Jacobian rank; not global identifiability")
  }
  feasible <- all(theta>=lo-active_tolerance & theta<=hi+active_tolerance)
  list(stationary=finite && feasible && scaled<=tolerance,feasible=feasible,scaled_gradient=scaled,
    gradient=gradient,projected_gradient=projected,column_norms=norms,
    tolerance=tolerance,coordinates="physical",local_rank=local,
    criterion="Jacobian-column-scaled projected gradient")
}
