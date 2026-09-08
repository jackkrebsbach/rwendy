# Bounds are always specified in physical parameter/state units. A bounded
# state block uses a diagonal coordinate map, so box-constrained optimizers
# enforce the requested bounds rather than bounding correlated GP z values.
.jgp_parameter_bound_alias <- function(value,alias,name,alias_name) {
  if (!is.null(value) && !is.null(alias))
    stop("Supply only ",name,"; ",alias_name," is its compatibility alias.")
  if (is.null(value)) alias else value
}

.jgp_bound_pair <- function(lower,upper,size,label) {
  expand <- function(x,default,side) {
    if (is.null(x)) return(rep(default,size))
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x)%in%c(1L,size) || anyNA(x))
      stop(label," ",side," must be a numeric scalar or a vector of length ",size," without NA.")
    rep_len(x,size)
  }
  lo <- expand(lower,-Inf,"lower bound"); hi <- expand(upper,Inf,"upper bound")
  if (any(lo>=hi)) stop("Each ",label," lower bound must be smaller than its upper bound.")
  list(lower=lo,upper=hi,bounded=is.finite(lo)|is.finite(hi))
}

.jgp_bounded_coordinates <- function(coordinates,bounds,scales=NULL) {
  if (is.null(bounds)) return(coordinates)
  if (is.null(scales)) scales <- rep(1,ncol(coordinates$mu))
  for (d in which(bounds$bounded)) {
    coordinates$mu[,d] <- 0
    coordinates$L[[d]] <- scales[d]
  }
  coordinates
}

.jgp_coordinate_bounds <- function(coordinates,lower,upper,state_bounds=NULL) {
  D <- ncol(coordinates$mu); m <- nrow(coordinates$mu)
  lo <- rep(-Inf,m*D); hi <- rep(Inf,m*D)
  if (!is.null(state_bounds)) for (d in which(state_bounds$bounded)) {
    L <- coordinates$L[[d]]
    if (length(L)!=1L || L<=0) stop("Bounded states require a positive diagonal coordinate map.")
    ix <- (match(d,coordinates$order)-1L)*m+seq_len(m)
    lo[ix] <- (state_bounds$lower[d]-coordinates$mu[,d])/L
    hi[ix] <- (state_bounds$upper[d]-coordinates$mu[,d])/L
  }
  list(lower=c(lower/coordinates$pscale,lo),upper=c(upper/coordinates$pscale,hi))
}

.jgp_pack_coordinates <- function(p,U,coordinates) {
  c(p/coordinates$pscale,unlist(lapply(coordinates$order,function(d) {
    L <- coordinates$L[[d]]; centered <- U[,d]-coordinates$mu[,d]
    if (length(L)==1L) centered/L else forwardsolve(L,centered)
  })))
}

# Projection is used only to initialize an optimization, never to repair its
# reported solution after fitting.
.jgp_clip_state <- function(U,bounds) {
  if (is.null(bounds)) return(U)
  for (d in which(bounds$bounded)) U[,d] <- pmax(bounds$lower[d],pmin(bounds$upper[d],U[,d]))
  U
}
