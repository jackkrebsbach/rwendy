
# %%
library(wendy)
library(deSolve)
library(devtools)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  u1 <- p[1] * u[2]
  u2 <- p[2] * u[2] + p[3] * u[1] + p[4] * u[1]^3
  c(u1, u2)
}

p_star <- c(1, -0.2, -0.05, -1);
p0 <- c(2, -0.1, -1, -0.25);
u0 <- c(0,2);
npoints <- 128
t_span <- c(0.005, 20);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);
# t_eval <- sort(runif(npoints, min = t_span[1], max = t_span[2])) non-uniform spacing

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
nr <- 0.25
U_vec <- as.vector(sol[,-1])
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, U, tt, method = "MLE", control = list(optimize = TRUE))
sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)

plot(U[,1],U[,2], cex = 0.5)
points(sol_hat[,2], sol_hat[,3], cex = 0.5, col = "red")

plot(res$wendy_problems[[1]]$min_radius_radii, res$wendy_problems[[1]]$min_radius_errors)
abline(v = res$wendy_problems[[1]]$min_radius, col = "red")
