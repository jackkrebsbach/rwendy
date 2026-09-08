# One fixed-covariance least-squares objective, with explicit coordinate maps.
# The scalar adapter uses derivative products: it need not assemble a stacked
# residual Jacobian or transform a dense weak Jacobian at every scalar step.
.jgp_joint_objective <- function(Y,H,weak,model,weights,noise,lambda,coordinates,
                                 prior_state=NULL,prior_components=integer(),
                                 observed=seq_len(model$D),prior_whitened=FALSE) {
  m <- nrow(coordinates$mu); D <- model$D; J <- model$J; n <- nrow(Y)
  order <- coordinates$order; pscale <- coordinates$pscale
  index <- lapply(seq_len(D),function(d) J+(match(d,order)-1L)*m+seq_len(m))
  physical <- function(theta) {
    U <- coordinates$mu
    for (d in seq_len(D)) {
      z <- theta[index[[d]]]; L <- coordinates$L[[d]]
      U[,d] <- U[,d]+if (length(L)==1L) L*z else as.vector(L%*%z)
    }
    list(p=theta[seq_len(J)]*pscale,U=U)
  }
  transform <- function(B,d) {
    L <- coordinates$L[[d]]
    if (length(L)==1L) B*L else B%*%L
  }
  Jdata <- matrix(0,n*length(observed),J+m*D)
  for (k in seq_along(observed))
    Jdata[(k-1L)*n+seq_len(n),index[[observed[k]]]] <- -transform(H,observed[k])/noise[k]
  Jprior <- matrix(0,m*length(prior_components),J+m*D)
  for (k in seq_along(prior_components)) {
    d <- prior_components[k]
    L <- coordinates$L[[d]]
    Jprior[(k-1L)*m+seq_len(m),index[[d]]] <- if (prior_whitened) diag(m) else
      forwardsolve(prior_state$L[[d]],if (length(L)==1L) diag(L,m) else L)
  }
  whiten <- prior_whitened && identical(weak$integration,"gp_gauss")
  ode_scale <- sqrt(lambda); last <- cached <- NULL
  function(theta,jacobian=TRUE,scalar=FALSE) {
    if (identical(theta,last) && (!jacobian || !is.null(cached$gradient)) &&
        (scalar || !jacobian || !is.null(cached$J))) return(cached)
    point <- physical(theta); U <- point$U; p <- point$p
    rd <- sweep(Y-H%*%U[,observed,drop=FALSE],2,noise,"/")
    rp <- matrix(0,m,length(prior_components))
    for (k in seq_along(prior_components)) {
      d <- prior_components[k]
      rp[,k] <- if (prior_whitened) theta[index[[d]]] else
        forwardsolve(prior_state$L[[d]],U[,d]-prior_state$mu[,d])
    }
    v <- .jgp_weak_eval(U,p,weak,model,jacobian,whiten)
    wr <- as.vector(weights$W%*%v$r); ode <- ode_scale*wr
    qi <- weights$interior$rank
    contributions <- c(data=sum(rd^2)/2,gp=sum(rp^2)/2,
      interior=sum(ode[seq_len(qi)]^2)/2,boundary=sum(ode[-seq_len(qi)]^2)/2)
    gradient <- jac <- NULL
    if (jacobian) {
      if (scalar) {
        score <- lambda*crossprod(weights$W,wr)
        gu <- as.vector(crossprod(v$Ju,score))
        gradient <- as.vector(crossprod(Jdata,as.vector(rd))+crossprod(Jprior,as.vector(rp)))
        gradient[seq_len(J)] <- as.vector(crossprod(v$Jp,score))*pscale
        for (d in seq_len(D)) {
          g <- gu[(d-1L)*m+seq_len(m)]; L <- coordinates$L[[d]]
          gradient[index[[d]]] <- gradient[index[[d]]]+if (whiten) g else
            if (length(L)==1L) L*g else as.vector(crossprod(L,g))
        }
      } else {
        JU <- matrix(0,nrow(v$Ju),m*D)
        for (d in seq_len(D)) {
          block <- v$Ju[,(d-1L)*m+seq_len(m),drop=FALSE]
          JU[,index[[d]]-J] <- if (whiten) block else transform(block,d)
        }
        jac <- rbind(Jdata,Jprior,ode_scale*weights$W%*%
          cbind(sweep(v$Jp,2,pscale,"*"),JU))
        gradient <- as.vector(crossprod(jac,c(rd,rp,ode)))
      }
    }
    last <<- theta
    cached <<- list(r=c(rd,rp,ode),J=jac,value=sum(contributions),gradient=gradient,
      p=p,U=U,Ju=if (!whiten) v$Ju else NULL,
      prior=if (prior_whitened) list(rows=length(rd)+seq_len(m*D),parameters=J,weight=1) else NULL,
      contributions=contributions,prior_contributions=setNames(colSums(rp^2)/2,prior_components),raw=v$r)
    cached
  }
}

# The optimizer choice does not change the objective or its derivative source.
.jgp_optimize <- function(par,evaluate,control,lower=rep(-Inf,length(par)),
                          upper=rep(Inf,length(par))) {
  solver <- attr(control,"solver")
  if (is.null(solver)) solver <- "lm"
  if (solver=="lm") {
    out <- .jgp_lm(par,evaluate,control,lower,upper); out$method <- solver
    return(out)
  }
  fn <- function(x) evaluate(x,TRUE,scalar=TRUE)
  out <- .jgp_scalar_optimize(par,fn,control$maxit,solver,lower=lower,upper=upper)
  out$objective <- fn(out$par)$value; out$history <- NULL
  out
}

.jgp_scalar_optimize <- function(par,fn,maxit,solver,active=NULL,
                                 lower=rep(-Inf,length(par)),upper=rep(Inf,length(par))) {
  if (!is.null(active)) {
    full <- par
    sub <- function(a,gradient=FALSE) {
      x <- full; x[active] <- a; y <- fn(x)
      if (gradient) y$gradient[active] else y$value
    }
    r <- .jgp_scalar_optimize(par[active],list(value=sub,
      gradient=function(a)sub(a,TRUE)),maxit,solver,lower=lower[active],upper=upper[active])
    full[active] <- r$par; r$par <- full
    return(r)
  }
  fv <- if (is.function(fn)) function(x)fn(x)$value else fn$value
  gv <- if (is.function(fn)) function(x)fn(x)$gradient else fn$gradient
  par <- pmax(lower,pmin(upper,par))
  if (solver=="nlminb") {
    a <- stats::nlminb(par,fv,gv,lower=lower,upper=upper,
      control=list(iter.max=maxit,eval.max=4*maxit,rel.tol=1e-12))
    list(par=a$par,converged=a$convergence==0,code=a$convergence,reason=a$message,
      iterations=a$iterations,evaluations=a$evaluations,method=solver)
  } else {
    a <- stats::optim(par,fv,gv,method="L-BFGS-B",lower=lower,upper=upper,
      control=list(maxit=maxit,factr=1e7))
    list(par=a$par,converged=a$convergence==0,code=a$convergence,reason=a$message,
      iterations=NA_integer_,evaluations=a$counts,method=solver)
  }
}
