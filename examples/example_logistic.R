
# %%
# library(wendy)
library(deSolve)
library(devtools)
library(ggplot2)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1);
u0 <- c(0.01);
p0 <- c(0.75, 0.75);
npoints <- 1024
t_span <- c(0.0, 10);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# set.seed(86753)

nr <- 0.2
U_vec <- as.vector(sol[,-1])

# Additive Gaussian Noise
# noise_sd <- nr * sqrt(mean(U_vec^2))
# noise <- sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)

# Multiplicative Lognormal Noise
noise_sd <- nr
noise <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))

U <- matrix(c(noise), ncol = 1)
tt <- sol[, 1, drop = FALSE]

# radii <- seq(2, 20)
# errs <- rep(0, length(radii))

# for(i in seq(length(radii))){
#    r <- radii[i]
#   print(res$phat)
#   rel_err <- norm(res$phat - p_star, type = "2")/norm(p_star, type = "2")
#   errs[i] <- rel_err
# }

# plot(radii, errs)

# library(ggplot2)

# df <- data.frame(radius = radii, rel_error = errs)

# ggplot(df, aes(x = radius, y = rel_error)) +
#   geom_line(color = "#378ADD", linewidth = 0.8) +
#   geom_point(color = "#378ADD", fill = "white", shape = 21, size = 3, stroke = 1.2) +
#   scale_x_continuous(breaks = seq(min(radii), max(radii), by = 2)) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
#   labs(
#     title = "Relative coefficient error vs. Radius of test functions",
#     subtitle = "SSL Test Functions",
#     x = "Radius",
#     y = "Relative coefficient error"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title    = element_text(face = "bold", size = 14),
#     plot.subtitle = element_text(color = "gray50", size = 11),
#     panel.grid.minor = element_blank(),
#     axis.line     = element_line(color = "gray80"),
#     plot.margin   = margin(12, 16, 12, 12)
#   )

# t_eval2 <- seq(t_span[1], t_span[2], length.out = 256);
# sol_hat <- deSolve::ode(u0, t_eval2, modelODE, res$phat)

# t_eval_dense <- seq(t_span[1], t_span[2], length.out = 256);
# sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

# interp_colors <- c("purple", "darkorange", "forestgreen", "deeppink", "cyan4")
# problem_names <- names(res$wendy_data)

# plot(tt, U, cex = 1, xlab = "Time", ylab=  "u₁", col="black")
# for(i in seq_along(res$wendy_problems)){
#   prob <- res$wendy_problems[[i]]
#   points(prob$tt, prob$U, cex = 0.75, col = interp_colors[i])
# }
# lines(t_eval2, sol_hat[,2], cex = 0.25, col = "#1f77b4")
# lines(t_eval_dense, sol_true[,2], cex = 0.5, col = "red")
# title(paste0("nr: ", nr,
#             "\n n: ", npoints,
#             "\n p̂: ", round(res$phat[1],3), " ", round(res$phat[2], 3)
#             ))
# legend(
#   "bottomright",
#   legend = c("data", "inferred trajectory", "true trajectory", problem_names),
#   col    = c("black", "#1f77b4", "red", interp_colors[seq_along(problem_names)]),
#   pch    = c(1, NA, NA, rep(1, length(problem_names))),
#   lty    = c(NA, 1, 1, rep(NA, length(problem_names))),
#   xpd    = TRUE,
#   bty    = "n",
#   cex = 0.8
# )

# plot(res$wendy_problems[[1]]$min_radius_radii, res$wendy_problems[[1]]$min_radius_errors)
# print(time)
# print(res$phat)

# p0_multi <- matrix(data = c(1.5, 0.5, 0.5, 0.5, 0, 0.2), ncol = 2)

res <- solveWendy(f, U, tt, p0 = c(0.5, 0.5), method = "OE", noise_dist = "lognormal")
sum <- summary(res)
cov <- solve(res$H_wnll(res$phat))

print(res$phat)

plot(res$wendy_problems[[1]]$min_radius_radii, res$wendy_problems[[1]]$min_radius_errors)
# print("Standard deviation from hessian")
# print(sqrt(diag(cov)))
# print("Standard deviation from wendy estimator")
# print(sqrt(diag(sum$param_cov)))