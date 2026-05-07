
# %%
library(deSolve)
library(devtools)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  u1 <- p[1] * u[1] + p[2] * u[1]^3  + p[3] * u[2]
  u2 <- p[4] * u[1] + p[5] + p[6] * u[2]
  c(u1, u2)
}

p_star <- c(3, -3, 3, -1/3, 17/150, 1/15);
u0 <- c(0,0.1);
u0 <- c(-1,1)
p0 <- c(1,-1, 1, -0.2, 0.01, 0.01);
npoints <- 1000
t_span <- c(0.001, 25);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
nr <- 0.05
noise_sd <- sqrt(nr * norm(as.array(sol[,-1]), type = "2") / npoints)
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, U, tt, method = "IRLS", control = list(estimate_U_star = FALSE))
# res <- solveWendy(f, U, tt, p0=p0, method = "OE")
sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)

plot(tt, U[,1], cex = 0.5)
lines(tt, sol[,2], col = "blue")
points(tt, sol_hat[,2], cex = 0.5, col = "red")

print(norm(res$phat - p_star, type = "2") / norm(p_star, type = "2"))