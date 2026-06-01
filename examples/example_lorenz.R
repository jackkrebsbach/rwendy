
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
npoints <- 512
t_span <- c(0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-12)

nr <- 0.05
U_vec <- as.array(sol[-1])
noise_sd <- nr * sqrt(mean(U_vec^2))

set.seed(8675309)

noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

time <- system.time({
  res <- solveWendy(f, U, tt, method = "IRLS", control = list(test_fun_type = "MSG", estimate_IC = TRUE, estimate_trajectory = TRUE))
})

print(time)

# resOE <- solveWendy(f, U, tt, p0 = p0, method = "OE")

# cat(sprintf("\nIRLS p̂ = [%s]", paste(sprintf("%.3f", res$phat), collapse = ", ")))
# cat(sprintf("    rel error = %.4f", wendy::rel_err(res$phat, p_star)))

# cat(sprintf("\nOE p̂ = [%s]  rel_err = %.4f",
#             paste(sprintf("%.4f", resOE$phat), collapse = ", "),
#             rel_err(resOE$phat, p_star)))

# cat(sprintf("\nu₀ = [%s]",
#               paste(sprintf("%.4f",u0), collapse = ", ")))

# cat(sprintf("\nû₀ = [%s]",
#             paste(sprintf("%.4f", res$u0hat), collapse = ", ")))

# cat(sprintf("\nERTS û₀ = [%s]  rel_err = %.4f\n",
#             paste(sprintf("%.4f", res$state$U_star[1,]), collapse = ", "),
#             rel_err(res$state$U_star[1,], u0)))


# cat(sprintf("\nOE û₀ = [%s]",
#               paste(sprintf("%.4f", resOE$data$u0), collapse = ", ")))