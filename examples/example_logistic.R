
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
u0 <- c(0.1);
p0 <- c(0.75, 0.75);
npoints <- 32
t_span <- c(0.0, 10);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# set.seed(86753)

nr <- 0.25
U_vec <- as.vector(sol[,-1])

# Additive Gaussian Noise
# noise_sd <- nr * sqrt(mean(U_vec^2))
# noise <- sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)

# Multiplicative Lognormal Noise
noise_sd <- nr
noise <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))

U <- matrix(c(noise), ncol = 1)
tt <- sol[, 1, drop = FALSE]

res <- solveWendy(f, U, tt, method = "IRLS", noise_dist = "lognormal")

t_eval2 <- seq(t_span[1], t_span[2], length.out = 256);
sol_hat <- deSolve::ode(u0, t_eval2, modelODE, res$phat)

t_eval_dense <- seq(t_span[1], t_span[2], length.out = 256);
sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

interp_colors <- c("purple", "darkorange", "forestgreen", "deeppink", "cyan4")
problem_names <- names(res$wendy_data)

plot(tt, U, cex = 1, xlab = "Time", ylab=  "u₁", col="black")
for(i in seq_along(res$wendy_problems)){
  prob <- res$wendy_problems[[i]]
  points(prob$tt, exp(prob$U), cex = 0.75, col = interp_colors[i])
}

lines(t_eval2, sol_hat[,2], cex = 0.25, col = "#1f77b4")
lines(t_eval_dense, sol_true[,2], cex = 0.5, col = "red")
points(tt, exp(res$state$U_star), cex = 1, xlab = "Time", ylab=  "u₁", col="green")

title(paste0("nr: ", nr,
            "\n n: ", npoints,
            "\n p̂: ", round(res$phat[1],3), " ", round(res$phat[2], 3)
            ))
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

cat(rel_err(res$phat, p_star))
cat(rel_err(res$u0hat, u0), "\n")
cat(rel_err(res$state$U_star[1, ], u0), "\n")
cat(rel_err(U[1, ], u0))