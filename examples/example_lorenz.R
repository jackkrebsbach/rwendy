library(deSolve)
library(plotly)

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

noise_sd <- 0.15
p_star <- c(10.0, 28.0, 4.0)
p0 <- c(12.10, 21, 3.0)
u0 <- c(2, 1, 1)

npoints <- 256
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

res <- solveWendy(f, p0, U, tt)

#phat <- res$solution
phat <- res$argument

sol_hat <- deSolve::ode(u0, t_eval, modelODE, phat)[, -1]

plot_ly(
  x = U[, 1],
  y = U[, 2],
  type = 'scatter',
  mode = 'markers',
  marker = list(color = 'blue', size = 3),
  name = "data"
) |>
  add_trace(
    x = sol_hat[, 1],
    y = sol_hat[, 2],
    type = 'scatter',
    mode = 'markers',
    marker = list(color = 'red', size = 3),
    name = "fit"
  )

print(phat)


