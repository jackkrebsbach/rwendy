
# %%
# Latent Lorenz. The z component is never observed; it is inferred jointly with
# the parameters from x and y alone. Observed components use a Lagrange
# extension, the latent one a GP extension, because a component with no data
# anchor needs the smoothing the GP conditional mean supplies.
library(deSolve)
library(plotly)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  du1 <- p[1] * (u[2] - u[1])
  du2 <- u[1] * (p[2] - u[3]) - u[2]
  du3 <- u[1] * u[2] - p[3] * u[3]

  c(du1, du2, du3)
}

p_star <- c(10.0, 28.0, 8.0 / 3.0)
u0 <- c(-8, 8, 27)
npoints <- 64
t_span <- c(0, 5)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star,
                    rtol = 1e-12, atol = 1e-12)
truth <- sol[, -1]
tt <- sol[, 1]
  
set.seed(8675309)

nr <- 0.1

U_vec <- as.array(sol[-1])
noise_sd <- nr * sqrt(mean(U_vec^2))
# noise_sd <- 2

noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U_full <- sol[, -1] + noise

names3 <- c("x", "y", "z")

# Hide z entirely. A latent component is an all-NA column; the observed columns
# must stay fully finite.
U <- U_full
U[, 3] <- NA_real_

# stopifnot(all(is.na(U[, 3])), all(is.finite(U[, 1:2])))

# Joint inference. Stage 1 carries the latent grid values with no prior penalty,
# stage 2 fits its prior to the stage-1 curve, stage 3 refits with that prior
# frozen and restarted from the smoothed curve.
time <- system.time({
  fit <- solveWendyGP(f, U, tt, control = list(grid_min = 129L))
})
print(time)

cat(sprintf("\nobserved components: %s   latent: %s\n",
  paste(names3[fit$observed], collapse = ", "), names3[fit$latent]))
cat(sprintf("extension per component: %s\n",
  paste(fit$diagnostics$weak_extension, collapse = " / ")))
cat(sprintf("test modes %d of %d available, %d weak rows, %d state variables\n",
  fit$diagnostics$modes, fit$diagnostics$available,
  fit$diagnostics$weak_rows, fit$diagnostics$state_variables))

# This gate is the GP conditional sd at the quadrature nodes divided by the
# component scale: how much the working grid leaves the state undetermined
# BETWEEN grid points. It is advisory here -- it fails at grid_min = 129 while
# the fit is still good -- and the remedy if you want it clear is a larger
# grid_min. Quadrature errors are checked separately at twice weak_quad_order.

cat(sprintf("grid resolution gate: %.2e (tol %.0e) passed=%s\n",
  max(fit$diagnostics$extension), 5e-3, fit$diagnostics$extension_passed))

cat(sprintf("\np*   = [%s]\n", paste(sprintf("%8.4f", p_star), collapse = ", ")))
cat(sprintf("phat = [%s]\n", paste(sprintf("%8.4f", fit$phat), collapse = ", ")))
cat(sprintf("phat RCE = %.3e\n",
  rel_err(fit$phat, p_star)))

# Score the recovered latent against the truth it never saw.
ref <- deSolve::ode(y = u0, times = fit$tt, func = modelODE, parms = p_star,
                    rtol = 1e-12, atol = 1e-12)[, -1]
z_hat <- fit$U_hat[, fit$latent]
z_true <- ref[, fit$latent]
cat(sprintf("latent z NMSE = %.3e\n", sum((z_hat - z_true)^2) / sum(z_true^2)))

# One figure: component colors stay consistent across truth, data and estimates.
# Only x and y have observation markers; z truth is a reference, never fit data.
component_colors <- c(x = "#0072B2", y = "#D55E00", z = "#009E73")
fig <- plot_ly()
for (d in seq_along(names3)) {
  component <- names3[d]
  color <- component_colors[[component]]
  is_latent <- d %in% fit$latent
  fig <- fig |>
    add_trace(x = fit$tt, y = ref[, d], type = "scatter", mode = "lines",
              name = paste(component, if (is_latent) "truth (never observed)" else "truth"),
              legendgroup = component, line = list(color = color, width = 1.5),
              opacity = 0.65) |>
    add_trace(x = fit$tt, y = fit$U_hat[, d], type = "scatter", mode = "lines",
              name = paste(component, if (is_latent) "inferred" else "fitted"),
              legendgroup = component, line = list(color = color, width = 2.5, dash = "dash"))
  if (d %in% fit$observed) {
    fig <- fig |>
      add_trace(x = tt, y = U[, d], type = "scatter", mode = "markers",
                name = paste(component, "observations"), legendgroup = component,
                marker = list(color = color, size = 4), opacity = 0.55)
  }
}
fig <- fig |> layout(
  title = list(text = paste0("Latent Lorenz: x and y observed; z inferred",
    "<br><sup>Estimated parameters: [", paste(sprintf("%.3f", fit$phat), collapse = ", "),
    "] · True: [", paste(sprintf("%.3f", p_star), collapse = ", "), "]</sup>")),
  xaxis = list(title = "Time"), yaxis = list(title = "State"),
  hovermode = "x unified",
  legend = list(groupclick = "toggleitem"))
print(fig)
