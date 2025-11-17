library(deSolve)
library(symengine)

source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/noise.R")
source("./R/optim.R")
source("./R/weak_residual.R")
source("./R/wendy.R")

f <- function(u, p, t) {
  du1 <- -p[1] * u[1] +
          p[3] * u[2] +
          u[3] * (p[1] * exp(-p[1] * p[2])) / (1 - exp(-p[1] * p[2]))
  du2 <- p[1] * u[1] -
         p[3] * u[2] -
         p[4] * (1 - exp(-p[5] * t^2)) * u[2]
  du3 <- p[4] * (1 - exp(-p[5] * t^2)) * u[2] -
        (p[1] * exp(-p[1] * p[2]) / (1 - exp(-p[1] * p[2]))) * u[3]
  c(du1, du2, du3)
}

noise_sd <- 0.05
npoints <- 256
p_star <- c(1.99, 1.5, 0.074, 0.113, 0.0024)
p0 <- c(0.5, 0.5, 0.1, 0.1, 0.001)
u0 <- c(1, 0, 0)
t_span <- c(0, 50)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)
# U <- sol[, -1] * exp(noise)
U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

tt <- sol[, 1]
u1 <- U[,1]
u2 <- U[, 2]
u3 <- U[, 3]

plot(tt, u1, type = "p", lwd = 0.25, col = "black",
     xlab = "Time", ylab = "State",
     main = "SIR Components")

points(tt, u2, lwd = 0.25, col = "red")
points(tt, u3, lwd = 0.25, col = "blue")

legend("topright", legend = c("u1", "u2", "u3"), col = c("black", "red", "blue"), lwd = 1)

# Log Normal Noise
#noisy <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))
#U <- matrix(noisy, ncol = 1)
#tt <- matrix(sol[, 1], ncol = 1)

#constraints = [(1e-4, 1.0), (1e-4, 2.0), (1e-4, 1.0), (1e-4, 1.0), (1e-4, 1.0)]
# constraints

res <- solveWendy(f, p0, U, tt, method = "MLE", optimize = T)

# upper <-  c(1,2,1,1,1)
# lower <- c(1e-4, 1e-4, 1e-4, 1e-4, 1e-4)
# constraints <- list(lower, upper)
#
# obj <- res$wnll
# J_obj <- res$J_wnll
# H_obj <- res$H_wnll
#
# result <- optimr(
#   par = p0,
#   fn = obj,
#   gr = J_obj,
#   method = "L-BFGS-B",
#   lower = lower,
#   upper = upper,
#   control = list(trace = 1)
# )

# sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)
# plot(U, cex = 0.5)
# points(sol_hat[,2], cex = 0.5, col = "red")

#cat(calc_gradient(p_star, obj))
#cat(J_obj(p_star))
