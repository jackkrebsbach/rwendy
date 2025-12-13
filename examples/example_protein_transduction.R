
# %%
library(deSolve)

invisible(sapply(list.files("./R", pattern = "\\.R$", full.names = TRUE), source))

f <- function(u, p, t) {
  u1 <- p[1] * u[1] + p[2] * u[1] * u[3] + p[3] * u[4]
  u2 <- p[4] * u[1]
  u3 <- p[5] * u[1] * u[3] + p[6] * u[4] + p[7] * u[5] / (0.3 + u[5])
  u4 <- p[8] * u[1] * u[3] + p[9] * u[4]
  u5 <- p[10] * u[4] + p[11] * u[5] / (0.3 + u[5])
   c(u1, u2, u3, u4, u5)
}

p_star <- c(-0.07, -0.6, 0.35, 0.07, -0.6, 0.05, 0.17, 0.6, -0.35, 0.3, -0.017)
u0 <- c(1, 0, 1, 0, 1)
p0 <- c(-0.01, -0.1, 0.1, 0.01, -0.1, 0.01, 0.1, 0.1, -0.1, 0.1, -0.01)
npoints <- 256
t_span <- c(0.001, 25)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
nr <- 0.005
U_vec <- as.array(sol[-1])
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)
U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

control <- list(radius_max_time = 15)
res <- solveWendy(f, p0, U, tt, control = control, method = "MLE")
sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)

plot(tt, U[,1], cex = 0.5)
points(tt, U[,2], cex = 0.5)
points(tt, U[,3], cex = 0.5)
points(tt, U[,4], cex = 0.5)
points(tt, U[,5], cex = 0.5)
points(sol_hat[,1], sol_hat[,2], cex = 0.5, col = "red")
points(sol_hat[,1], sol_hat[,3], cex = 0.5, col = "red")
points(sol_hat[,1], sol_hat[,4], cex = 0.5, col = "red")
points(sol_hat[,1], sol_hat[,5], cex = 0.5, col = "red")
points(sol_hat[,1], sol_hat[,6], cex = 0.5, col = "red")