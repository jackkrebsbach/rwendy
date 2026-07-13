
# %%
# library(wendy)
library(deSolve)
library(MASS)
library(uGMAR)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  u1 <- p[1] * u[1] + p[2] * u[1] * u[2]
  u2 <- p[3] * u[2] + p[4] * u[1] * u[2]
  c(u1, u2)
}

p_star <- c(3, -1, -6, 1)
u0 <- c(1, 1)

p_star <- c(1, -0.1, -1.5, 0.075)
u0 <- c(10,5)
p0 <- c(2, -0.1, -1, 0.25)
npoints <- 128
t_span <- c(0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
nr <- 0.025

# set.seed(8675309 + 1)

U_vec <- as.vector(sol[,-1])
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)
U <- sol[, -1] + noise
# U[,1] <- mean(U[,2])
tt <- matrix(sol[, 1], ncol = 1)

time <- system.time({
res  <- solveWendy(f, U, p0 = p0, tt, method = "IRLS",
 control = list(test_fun_type = "MSG")
)
})

# plot(tt, U[,1], col = adjustcolor("brown", alpha.f = 0.3), cex = 0.5,
#    ylab = "State u₁ & u₂",
#   main = "Lokta Volterra")
# points(tt, U[,2], col = adjustcolor("blue", alpha.f = 0.3), cex = 0.5)

# lines(tt, sol[,-1][,1], cex = 0.5, col = "black")
# lines(tt, sol[,-1][,2], cex = 0.5, col = "black")

# # lines(res$tt, res$U[,1], cex = 0.5, col = "black")
# # lines(res$tt, res$U[,2], cex = 0.5, col = "black")

# cat(sprintf("\np̂_IRLS  = [%s]  rel_err = %.4f",
#             paste(sprintf("%.4f", res$phat), collapse = ", "),
#             rel_err(res$phat, p_star)))

# # cat(sprintf("\np̂_MLE   = [%s]  rel_err = %.4f\n",
# #             paste(sprintf("%.4f", res2$phat), collapse = ", "),
# #             rel_err(res2$phat, p_star)))
plot_radius_selection(res)