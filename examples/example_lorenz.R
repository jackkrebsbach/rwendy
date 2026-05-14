
# %%
# library(wendy)
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
p0 <- c(12.0, 21, 4.0)
u0 <- c(-8, 10, 27)
npoints <- 256
t_span <- c(0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-12)

nr <- 0.15
U_vec <- as.array(sol[-1])
noise_sd <- nr * sqrt(mean(U_vec^2))

set.seed(8675309)

noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, U, tt, method = "IRLS", control = list(estimate_u0 = FALSE))

# sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)[, -1]

# plot_ly(
#  x = sol[, 2],
#  y = sol[, 3],
#  z = sol[, 4],
#  type = 'scatter3d',
#  mode = 'marker',
#  marker = list(color = 'blue', size = 3),
#  name = "data"
# ) |>
#  add_trace(
#    x = sol_hat[, 1],
#    y = sol_hat[, 2],
#    z = sol_hat[, 3],
#    mode = 'marker',
#    marker = list(color = 'red', size = 3),
#    name = "fit"
#  )

res2 <- solveWendy(f, U, tt, method = "JOINT",
  control = list(
      test_fun_type = "SSL",
      include_boundary_layer = TRUE,
      alpha_gn = 0
    )
  )

res3 <- solveWendy(f, U, tt, method = "JOINT",
  control = list(
      test_fun_type = "SSL",
      include_boundary_layer = FALSE,
      alpha_gn = 0
    )
  )

res4 <- solveWendy(f, U, tt, method = "JOINT",
  control = list(
      test_fun_type = "MSG",
      alpha_gn = 0
    )
  )

cat(sprintf("\np̂₁ = [%s]", paste(sprintf("%.3f", res$phat), collapse = ", ")))
cat(sprintf("\n     rel error = %.4f", wendy::rel_err(res$phat, p_star)))

cat(sprintf("\np̂₂ = [%s]", paste(sprintf("%.3f", res2$phat), collapse = ", ")))
cat(sprintf("\n     rel error = %.4f", wendy::rel_err(res2$phat, p_star)))

cat(sprintf("\np̂₃ = [%s]", paste(sprintf("%.3f", res3$phat), collapse = ", ")))
cat(sprintf("\n     rel error = %.4f", wendy::rel_err(res3$phat, p_star)))

cat(sprintf("\np̂₄ = [%s]", paste(sprintf("%.3f", res4$phat), collapse = ", ")))
cat(sprintf("\n     rel error = %.4f", wendy::rel_err(res4$phat, p_star)))

# cat(sprintf("True u₀ [%s] ", paste(sprintf("%.3f", u0), collapse = ", ")))

# cat(sprintf("\nEstimated  û₀ = [%s]", paste(sprintf("%.3f", res$u0hat), collapse = ", ")))
# cat(sprintf("\n     rel error = %.4f", wendy::rel_err(res$u0hat, u0)))

# cat(sprintf("\nNoisy u₀ = [%s]", paste(sprintf("%.3f", U[1,]), collapse = ", ")))
# cat(sprintf("\n     rel error = %.4f", wendy::rel_err(U[1,], u0)))