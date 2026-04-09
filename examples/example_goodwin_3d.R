# %%
# library(wendy)
library(deSolve)
library(plotly)


invisible({devtools::load_all()})

f <- function(u, p, t) {
  du1 <- p[1] / (2.15 + p[3] * u[3]^p[4]) - p[2] * u[1]
  du2 <- p[5] * u[1] - p[6] * u[2]
  du3 <- p[7] * u[2] - p[8] * u[3]
  c(du1, du2, du3)
}

npoints <- 512
p_star <- c(3.4884, 0.0969, 1, 10, 0.0969, 0.0581, 0.0969, 0.0775)
p0 <- c(3, 0.1, 4, 12, 0.1, 0.1, 0.1, 0.1)
u0 <- c(0.3617, 0.9137, 1.393)
t_span <- c(0, 80)

modelODE <- function(tvec, state, parameters) {
  list(as.vector(f(state, parameters, tvec)))
}

t_eval <- seq(t_span[1], t_span[2], length.out = npoints)
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

nr <- 0.05
U_vec <- as.vector(sol[,-1])
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)
U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

control <- list(radius_max_time = 20)
res <- solveWendy(f, U, tt, method = "IRLS", control = control)

p_hat <- res$phat

sol_hat <- deSolve::ode(u0, t_eval, modelODE, p_hat)[,-1]

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

