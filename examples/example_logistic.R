library(deSolve)
library(symengine)
library(trust)
library(uGMAR)
library(trustOptim)

source("./R/noise.R")
source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/weak_residual.R")
source("./R/wendy.R")

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

noise_sd <- 0.05;
p_star <- c(1, 1);
u0 <- c(0.01);
p0 <- c(0.5, 0.5);
npoints <- 200;
t_span <- c(0.005, 10);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Guassian Noise
U <- matrix(c(sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)), ncol = 1)
tt <- matrix(sol[, 1], ncol = 1)

# Log Normal Noise
# noisy <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))
# U <- matrix(noisy, ncol = 1)
# tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, p0, U, tt, method = "IRLS")

sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)

plot(U, cex = 0.5)
points(sol[, 2], cex = 0.5, col = "blue")
points(sol_hat[, 2], cex = 0.5, col = "red")

print(res$phat)

