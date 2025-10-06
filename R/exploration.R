{
  library(deSolve)
  library(symengine)
  library(trust)
  library(uGMAR)
  library(wendy)
  library(trustOptim)

  source("./R/noise.R")
  source("./R/symbolics.R")
  source("./R/test_functions.R")
  source("./R/weak_residual.R")
  source("./R/wendy.R")

}

{
  f <- function(u, p, t) { c(p[1] * u[1] - p[2] * u[1]^2) }

  logistic <- function(u, p, t) {
    list(p[[1]] * u[[1]] - p[[2]] * u[[1]]^2)
  }

  noise_sd <- 0.05;
  p_star <- c(1, 1);
  u0 <- c(0.01);
  p0 <- c(0.5, 0.5);
  npoints <- 50;
  t_span <- c(0, 10);
  t_eval <- seq(t_span[1], t_span[2], length.out = npoints);
}

{
  modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
  sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)
  U <- matrix(c(sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)), ncol = 1)
  tt <- matrix(sol[, 1], ncol = 1)
}

res <- solveWendy(f, p0, U, tt, compute_svd = TRUE, optimize = TRUE)
test_funs <- WendySolver(logistic, U, p0, tt, log_level = "none", compute_svd_ = FALSE, optimize_ = FALSE)

V <- test_funs$V
V_r <- res$V








