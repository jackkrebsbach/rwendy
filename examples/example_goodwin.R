library(deSolve)
library(symengine)
library(trustOptim)
#library(wendy)

source("./R/noise.R")
source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/weak_residual.R")
source("./R/wendy.R")


f <- function(u, p, t) {
  du1 <- p[1] / (36 + p[2] * u[2]) - p[3]
  du2 <- p[4] * u[1] - p[5]
  c(du1, du2)
}

noise_sd <- 0.05
npoints <- 200
p_star <- c(72, 1, 2, 1, 1)
p0 <- c(70, 1.56, 2.5, 1.75, 0.6)
u0 <- c(7, -10)
t_span <- c(0, 60)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)


modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise

tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, p0, U, tt)

wnll <- res$wnll
J_wnll <- res$J_wnll

p_hat <- res$solution

sol_hat <- deSolve::ode(u0, t_eval, modelODE, p_hat)

plot(U[, c(1, 2)], cex = 0.5)
points(sol_hat[, c(1, 2)], cex = 0.5, col = "red")

print(res$solution)
