library(wendy)
library(deSolve)
library(symengine)
library(trust)

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}


noise_ratio <- 0.05;
p_star <- c(1, 1);
u0 <- c(0.01);
p0 <- c(0.5, 0.5);
npoints <- 103;
t_span <- c(0.001, 10);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
noise_sd <- noise_ratio * (norm(sol[,2, drop = F], "F")^2 / npoints)
U <- matrix(c(sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)), ncol = 1)
tt <- matrix(sol[, 1], ncol = 1)

# Log Normal Noise
# noisy <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))
# U <- matrix(noisy, ncol = 1)
# tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, p0, U, tt, method = "MLE")

sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)















