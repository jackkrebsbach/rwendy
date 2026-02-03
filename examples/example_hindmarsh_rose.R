
# %%
library(wendy)
library(deSolve)
library(plotly)

f <- function(u, p, t) {
  du1 <- p[1] * u[2] + p[2] * u[1]^3 + p[3] * u[1]^2 + p[4] * u[3]
  du2 <- p[5] + p[6] * u[1]^2 + p[7] * u[2]
  du3 <- p[8] * u[1] + p[9] + p[10] * u[3]
  c(du1, du2, du3)
}

npoints <- 512
p_star <- c(10, -10, 30, -10, 10, -50, -10, 0.04, 0.0319, -0.01)
p0 <-     c(15, -12, 35, -12, 5 , -60, -12, 0.08, 0.06, -0.04)
u0 <- c(-1.31, -7.6, -0.2)
t_span <- c(0, 10)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

t_eval <- seq(t_span[1], t_span[2], length.out = npoints)
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

nr <- 0.02
U_vec <- as.vector(sol[,-1])
noise_sd <- nr * sqrt(mean(U_vec^2))

noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, p0, U, tt, method = "MLE")

p_hat <- res$phat

U_hat <- deSolve::ode(
  y = u0,
  times = t_eval,
  func = modelODE,
  parms = p_hat,
  method = "lsodes",
  atol = 1e-10,
  rtol = 1e-10
)[, -1]

plot_ly() |>
  add_trace(
    x = sol[, -1][, 1],
    y = sol[, -1][, 2],
    z = sol[, -1][, 3],
    type = 'scatter3d',
    mode = 'lines',
    name = 'True Trajectory'
  ) |>
  add_trace(
    x = U_hat[, 1],
    y = U_hat[, 2],
    z = U_hat[, 3],
    type = 'scatter3d',
    mode = 'lines',
    name = 'Estimated Trajectory',
    line = list(color = 'red', opacity = 0.95)
  ) |>
  add_trace(
    x = U[, 1],
    y = U[, 2],
    z = U[, 3],
    type = 'scatter3d',
    mode = 'markers',
    name = 'Noisy Data',
    marker = list(color = 'black', size = 2, opacity = 0.35)
  ) |>
  layout(
    paper_bgcolor = 'rgba(0,0,0,0)',
    plot_bgcolor = 'rgba(0,0,0,0)',
    title = "Hindmarsh Rose",
    legend = list(
      x = 0.02,
      y = 0.98,
      xanchor = "left",
      yanchor = "top"
    )
  )

cat("pstar:", paste(p_star, collapse = " "), "\n")
cat("p_hat:", paste(format(p_hat, digits = 3, scientific = FALSE), collapse = " "),"\n")
