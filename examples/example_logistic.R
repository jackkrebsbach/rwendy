
# %%
library(wendy)
library(deSolve)

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1);
u0 <- c(0.01);
p0 <- c(0.5, 0.5);
npoints <- 10
t_span <- c(0.005, 10);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
nr <- 0.1
U_vec <- as.vector(sol[,-1])
noise_sd <- nr * sqrt(mean(U_vec^2))
U <- matrix(c(sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)), ncol = 1)
tt <- sol[, 1, drop = FALSE]

res <- solveWendy(f, p0, U, tt, lip = FALSE, method = "IRLS",
  control = list(test_fun_type = "MSG", min_number_points = 18))

t_eval2 <- seq(t_span[1], t_span[2], length.out = 256);
sol_hat <- deSolve::ode(u0, t_eval2, modelODE, res$phat)

plot(res$tt, res$U, cex = 0.5, xlab = "Time", ylab=  "u₁", col="purple")
points(tt, U, cex = 0.5, col = "black")
lines(t_eval2, sol_hat[,2], cex = 0.5, col = "#1f77b4")
lines(t_eval, sol[,2], cex = 0.5, col = "red")
title(paste0("nr ", nr, "\n p̂: ", round(res$phat[1],3), " ", round(res$phat[2], 3)))
legend(
  "bottomright",
  legend = c("data", "inferred trajectory", "true trajectory", "interpolated data"),
  col    = c("black", "#1f77b4", "red", "purple"),
  pch    = c(1, NA),          
  lty    = c(NA, 1),          
  xpd    = TRUE,             
  bty    = "n"                
)
print(res$phat)