# Logistic Example
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

  modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
  sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)
  U <- matrix(c(sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)), ncol = 1)
  tt <- matrix(sol[, 1], ncol = 1)
}

# Lorenz Example
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

  f <- function(u, p, t) {
    du1 <- p[1] * (u[2] - u[1])
    du2 <- u[1] * (p[2] - u[3]) - u[2]
    du3 <- u[1] * u[2] - p[3] * u[3]
    c(du1, du2, du3)
  }

  lorenz <- function(u, p, t) {
    du1 <- p[[1]] * (u[[2]] - u[[1]])
    du2 <- u[[1]] * (p[[2]] - u[[3]]) - u[[2]]
    du3 <- u[[1]] * u[[2]] - p[[3]] * u[[3]]
    list(du1, du2, du3)
  }

  noise_sd <- 0.05
  p_star <- c(10.0, 28.0, 8.0 / 3.0)
  p0 <- c(12.0, 21, 4.0)
  u0 <- c(2, 1, 1)
  npoints <- 116
  t_span <- c(0, 10)
  t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

  modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

  sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

  noise <- matrix(
    rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
    nrow = nrow(sol)
  )

  U <- sol[, -1] + noise

  tt <- matrix(sol[, 1], ncol = 1)
}

res <- solveWendy(f, p0, U, tt, optimize = T)
test_funs <- WendySolver(lorenz, U, p0, tt, log_level = "info", compute_svd_ = TRUE, optimize_ = FALSE)

V_r <- res$V_prime
V <- test_funs$V_prime






