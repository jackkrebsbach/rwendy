
# %%
library(wendy)
library(deSolve)

f <- function(u, p, t) {
  du1 <- -p[1] * u[1] + p[3] * u[2] + u[3] * (p[1] * exp(-p[1] * p[2])) / (1 - exp(-p[1] * p[2]))
  du2 <- p[1] * u[1] - p[3] * u[2] - p[4] * (1 - exp(-p[5] * t^2)) * u[2]
  du3 <- p[4] * (1 - exp(-p[5] * t^2)) * u[2] - (p[1] * exp(-p[1] * p[2]) / (1 - exp(-p[1] * p[2]))) * u[3]
  c(du1, du2, du3)
}

npoints <- 256
p_star = c(0.2, 1.5, 0.074, 0.113, 0.0024)
p0 <- c(0.1, 2.5, 0.1, 0.1, 0.001)
u0 <- c(1, 0, 0)
t_span <- c(0, 50)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

nr <- 0.005
noise_sd <- sqrt(nr)
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)
# Log Normal Noise
U <- sol[, -1] * exp(noise)
tt <- matrix(sol[, 1], ncol = 1)

tt <- sol[, 1]
u1 <- U[,1]
u2 <- U[, 2]
u3 <- U[, 3]

plot(tt, u1, type = "p", pch = 16, cex = 1, col = adjustcolor("black", alpha.f = 0.5), xlab = "Time", ylab = "Proportion", main = "SIR Compartments", ylim=c(0, 1))
points(tt, u2, pch = 16, cex = 1, col = adjustcolor("red", alpha.f = 0.5))
points(tt, u3, pch = 16, cex = 1, col = adjustcolor("blue", alpha.f = 0.5))

legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("black", "red", "blue"), lwd = 1)

res <- solveWendy(f, p0, U, tt, method = "MLE", noise_dist = "lognormal")

sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = res$phat)

lines(sol[, "time"], sol[, 2], col = "black", lwd = 2)
lines(sol[, "time"], sol[, 3], col = "red", lwd = 2)
lines(sol[, "time"], sol[, 4], col = "blue", lwd = 2)