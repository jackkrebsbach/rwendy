
# %%
# library(wendy)
library(deSolve)
library(devtools)
library(ggplot2)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1/10)
u0 <- c(0.1)
p0 <- c(0.75, 0.75)
npoints <- 256
t_span <- c(0.0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

set.seed(8675309)

nr <- 0.25
U_vec <- as.vector(sol[,-1])

# Additive Gaussian Noise
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)

# Multiplicative Lognormal Noise
# noise_sd <- nr
# noise <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))

cat(sprintf("σ = %.2f", noise_sd))

U <- matrix(c(noise), ncol = 1)
tt <- sol[, 1, drop = FALSE]

p0 = t(matrix(c(10, 10, 1.5, 1.5, 0.5, 0.5, 2, 2, 5, 
                5, 10,10 , 0.25, 10, 0.1, 8, 1.5, 1.5,
                0.5, 0.5, 2, 2, 5, 5, 10,10 , 0.25, 10,
                0.1, 8), nrow = 2))

time <- system.time({
  res <- solveWendy(f, U, tt, method = "IRLS", noise_dist = "addgaussian", control = list(estimate_U_star = TRUE, optimize = TRUE))
})

res2 <- solveWendy(f, res$state$U_star, tt, method = "OLS", noise_dist = "addgaussian")
res3 <- solveWendy(f, sol[,-1, drop = FALSE], tt, method = "OLS", noise_dist = "addgaussian")

cat(sprintf("\np̂₁ = [%s]", paste(sprintf("%.3f", res$phat), collapse = ", ")))
cat(sprintf("\n     rel error = [%.4f]", wendy::rel_err(res$phat, p_star)))
cat(sprintf("\np̂₂ = [%s]", paste(sprintf("%.3f", res2$phat), collapse = ", ")))
cat(sprintf("\n     rel error = [%.4f]", wendy::rel_err(res2$phat, p_star)))
cat(sprintf("\np* = [%s]", paste(sprintf("%.3f", p_star), collapse = ", ")))
cat(sprintf("\nRelative error: %.1e \nWalltime: %.4f\n ", rel_err(res$phat, p_star), time[1]))

t_eval2 <- seq(t_span[1], t_span[2], length.out = 256);
sol_hat <- deSolve::ode(u0, t_eval2, modelODE, res2$phat)

t_eval_dense <- seq(t_span[1], t_span[2], length.out = 256);
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
# points(tt, U_vec, cex = 0.75, col = "purple")
# points(tt, res$state$U_star, cex = 0.75, col = "green")

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



# u0 <- c(0.1);
# t_span <- c(0.0, 10);
# t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

# modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
# sol_3 <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = res3$phat)

# # cat(sprintf("\n| True state OLS p̂ = [%s]", 
# #     paste(sprintf("%.2f", res2$phat), collapse = ", ")))
# # cat(message(sprintf("\n| Kalman filtered state p̂ = [%s]", 
# #     paste(sprintf("%.2f", res3$phat), collapse = ", "))))

# plot(tt, U, col = "black")
# points(tt, U_vec, cex = 1, xlab = "Time", ylab=  "u₁", col="green")
# points(tt, res$state$U_star, col = "purple")
# lines(tt, sol_3[,2], col = "orange")
# legend(
#   "bottomright",
#   legend = c(
#     "Noisy data",
#     sprintf("True state p̂ = [%s]", paste(sprintf("%.2f", res2$phat), collapse = ", ")),
#     sprintf("Kalman filter state p̂ = [%s]", paste(sprintf("%.2f", res3$phat), collapse = ", ")),
#     "Inferred Trajectory"
#   ),
#   col = c("black","green", "purple", "orange"),
#   pch = c(1, 1, 1),
#   bty = "n",
#   cex = 0.8
# )

# plot(res$wendy_problems[[1]]$min_radius_radii, res$wendy_problems[[1]]$min_radius_errors)
# print(res$min_radius)

# r1 <- residuals(res)
# r2 <- residuals_weighted(res)

# d1 <- density(r1)
# d2 <- density(r2)

# plot(d1, main = "Residuals", xlab = "", xlim=c(min(r2), max(r2)) ,ylim = range(c(d1$y, d2$y)))
# lines(d2, col = "red")
# legend("topright", legend = c("unweighted", "weighted"), col = c("black", "red"), lty = 1, bty = "n")

# rt <- sol[,-1] - U
# r2 <- sol[,-1] - res$state$U_star
# plot(density(rt), xlim =c(min(rt), max(rt)))
# lines(density(rt))