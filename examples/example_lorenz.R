library(deSolve)
library(plotly)
library(trust)
library(symengine)

source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/noise.R")
source("./R/optim.R")
source("./R/weak_residual.R")
source("./R/wendy.R")

f <- function(u, p, t) {
  du1 <- p[1] * (u[2] - u[1])
  du2 <- u[1] * (p[2] - u[3]) - u[2]
  du3 <- u[1] * u[2] - p[3] * u[3]
  c(du1, du2, du3)
}

noise_ratio <- 0.05
p_star <- c(10.0, 28.0, 8.0 / 3.0)
p0 <- c(12.0, 21, 4.0)
u0 <- c(2, 1, 1)
npoints <- 256
t_span <- c(0.1, 20)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

noise_sd <- sqrt(noise_ratio * norm(as.array(sol[,-1]), type = "2") / npoints)
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, p0, U, tt, method = "MLE", optimize = T)

sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)[, -1]

plot_ly(
 x = sol[, 2],
 y = sol[, 3],
 z = sol[, 4],
 type = 'scatter3d',
 mode = 'marker',
 marker = list(color = 'blue', size = 3),
 name = "data"
) |>
 add_trace(
   x = sol_hat[, 1],
   y = sol_hat[, 2],
   z = sol_hat[, 3],
   mode = 'lines',
   marker = list(color = 'red', size = 3),
   name = "fit"
 )


