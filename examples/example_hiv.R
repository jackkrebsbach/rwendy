
# %%
library(wendy)
library(deSolve)

f <- function(u, p, t) {
   eta <- 9e-5 * (1 - 0.9 * cos(pi * t / 1000))

   u1 <- p[1] - p[2] * u[1] - eta * u[1] * u[3] * 100
   u2 <- eta * u[1] * u[3] * 100 - p[3] * u[2] 
   u3 <- p[4] * u[2] - p[5] * u[3] * 100

   c(u1, u2, u3) 
}

npoints <- 100
t_span <- c(0, 20)
p_star <- c(36, 0.108, 0.5, 5, 3)
u0 <- param.true$x0
D <- length(u0)
p0 <- c(38, 0.05, 0.25, 3, 2)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

t_eval <- seq(t_span[1], t_span[2], length.out = npoints)
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

U <- matrix(0, nrow = nrow(sol), ncol = length(u0))
tt <- matrix(sol[, 1], ncol = 1)
U_star <- sol[,-1]

nr <- 0.05

for(d in  seq(D)){
  U_d <- U_star[,d]
  noise_sd <- nr * sqrt(mean(U_d^2))
  noise <- rnorm(npoints, mean = 0, sd = noise_sd)
  U[,d] <- U_d + noise
}

# noise <- matrix(
#   rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = nr),
#   nrow = nrow(sol)
# )
# U <- sol[, -1] + noise

tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, p0, U, tt, method = "MLE", noise_dist = "addgaussian",
        control = list(test_fun_type = "MSG", radius_min_time = 0.01, radius_max_time = 5))

cat("pstar:", p_star, "\nphat :", res$phat)

mpl_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c")
sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)
U_hat <- sol_hat[,-1]
U_star <- sol[,-1]

par(mfrow = c(3, 1), mar = c(4, 4, 2, 1)) 
for(d in seq(D)){
  legend(
    "topright",
    legend = c("Truth", "Noisy observation", "Estimated"),
    col    = c("#000000", mpl_colors[1], mpl_colors[2]),
    pch    = c(16, 16, 16),
  )
  plot(tt, U_star[,d], ylab = paste0("u",d))
  points(tt, U[,d], col = mpl_colors[1])
  points(tt, U_hat[,d], col = mpl_colors[2])
}
