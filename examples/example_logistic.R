
# %%
# library(wendy)
library(deSolve)
library(uGMAR)
library(devtools)
library(ggplot2)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1/10)
u0 <- c(0.1)
p0 <- c(1.1, 0.2)
npoints <- 10
t_span <- c(0.0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-14)

set.seed(8675309 + 1)

nr <- 0.05
U_vec <- as.vector(sol[,-1])

# Additive Gaussian Noise
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- rnorm(npoints, mean = 0, sd = noise_sd)
U <- sol[, 2, drop = FALSE] + noise
tt <- sol[, 1, drop = FALSE]

# Multiplicative Lognormal Noise
# noise_sd <- nr
# noise <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))

cat(sprintf("σ = %.2f", noise_sd))

p0 = t(matrix(c(10, 10, 1.5, 1.5, 0.5, 0.5, 2, 2, 5, 
                5, 10,10 , 0.25, 10, 0.1, 8, 1.5, 1.5,
                0.5, 0.5, 2, 2, 5, 5, 10,10 , 0.25, 10,
                0.1, 8), nrow = 2))

time <- system.time({
  res <- solveWendy(f, U, tt, method = "MLE", noise_dist = "addgaussian",
                control = list(estimate_U_star = FALSE, optimize = TRUE, test_fun_type = "SSL"))
})

t_eval2 <- seq(t_span[1], t_span[2], length.out = npoints);
sol_hat <- deSolve::ode(u0, t_eval2, modelODE, res$phat)

t_eval_dense <- seq(t_span[1], t_span[2], length.out = npoints);
sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

interp_colors <- c("purple", "darkorange", "forestgreen", "deeppink", "cyan4")
problem_names <- names(res$wendy_data)

plot(tt, U, cex = 1, xlab = "Time", ylab=  "u₁", col="black")
for(i in seq_along(res$wendy_problems)){
  prob <- res$wendy_problems[[i]]
  points(prob$tt, prob$U, cex = 0.75, col = interp_colors[i])
}

lines(t_eval2, sol_hat[,2], cex = 0.25, col = "#1f77b4")
lines(t_eval_dense, sol_true[,2], cex = 0.5, col = "red")

title(paste0("nr: ", nr, "\n n: ", npoints, "\n p̂: ", round(res$phat[1],3), " ", round(res$phat[2], 3)))

legend(
  "bottomright",
  legend = c("data", "inferred trajectory", "true trajectory", problem_names),
  col    = c("black", "#1f77b4", "red", interp_colors[seq_along(problem_names)]),
  pch    = c(1, NA, NA, rep(1, length(problem_names))),
  lty    = c(NA, 1, 1, rep(NA, length(problem_names))),
  xpd    = TRUE,
  bty    = "n",
  cex = 0.8
)

cat(sprintf("\np̂ = [%s]  rel_err = %.4f\n",
            paste(sprintf("%.4f", res$phat), collapse = ", "),
            rel_err(res$phat, p_star)))