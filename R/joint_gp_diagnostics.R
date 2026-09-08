# Final diagnostics use the shared residual engine in physical (p,U)
# coordinates, independently of each optimization preconditioner.
.jgp_fit_diagnostics <- function(U, p, Y, H, noise, weak, model, weights, control,
                                 prior_state = NULL, prior_components = integer(),
                                 observed = seq_len(model$D),
                                 lower = rep(-Inf, length(p)),
                                 upper = rep(Inf, length(p)),
                                 tolerance = control$gtol, rank = TRUE) {
  m <- nrow(U); D <- ncol(U); J <- length(p); n <- nrow(Y)
  evaluate <- .jgp_joint_objective(Y,H,weak,model,weights,noise,control$lambda,
    coordinates=list(mu=matrix(0,m,D),L=rep(list(1),D),order=seq_len(D),pscale=rep(1,J)),
    prior_state=prior_state,prior_components=prior_components,observed=observed)
  final <- evaluate(c(p,as.vector(U)))
  residual <- final$r; jacobian <- final$J
  gradient <- as.vector(crossprod(jacobian,residual))
  theta <- c(p,as.vector(U))
  lo <- c(lower,rep(-Inf,m*D)); hi <- c(upper,rep(Inf,m*D))
  projected <- gradient
  projected[theta<=lo & gradient>0 | theta>=hi & gradient<0] <- 0
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
  list(stationary=finite && scaled<=tolerance,scaled_gradient=scaled,
    gradient=gradient,projected_gradient=projected,column_norms=norms,
    tolerance=tolerance,coordinates="physical",local_rank=local,
    criterion="Jacobian-column-scaled projected gradient")
}
