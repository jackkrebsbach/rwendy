
# %%
# library(wendy)
library(deSolve)
library(devtools)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1);
u0 <- c(0.01);
p0 <- c(0.5, 0.5);
npoints <-10 
t_span <- c(0.005, 10);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

nr <- 0.5
U_vec <- as.vector(sol[,-1])
noise_sd <- nr * sqrt(mean(U_vec^2))

# Additive Gaussian Noise
noise <- sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)

# Multiplicative Lognormal Noise
# noise <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = nr))

U <- matrix(c(noise), ncol = 1)
tt <- sol[, 1, drop = FALSE]

res <- solveWendy(f, p0, U, tt, lip = TRUE, method = "MLE", noise_dist = "addgaussian",
  control = list(
    min_number_points = 25,
    interpolation_method = "cubic_ls",
    test_fun_type = "MSG"
    )
  )

t_eval2 <- seq(t_span[1], t_span[2], length.out = 256);
sol_hat <- deSolve::ode(u0, t_eval2, modelODE, res$phat)

t_eval_dense <- seq(t_span[1], t_span[2], length.out = 256);
sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

interp_colors <- c("purple", "darkorange", "forestgreen", "deeppink", "cyan4")
interp_names  <- names(res$interp_list)

plot(tt, U, cex = 1, xlab = "Time", ylab=  "u₁", col="black")
for(i in seq_along(res$interp_list)){
  interps <- res$interp_list[[i]]
  points(interps$tt, interps$U, cex = 0.75, col = interp_colors[i])
}
lines(t_eval2, sol_hat[,2], cex = 0.25, col = "#1f77b4")
lines(t_eval_dense, sol_true[,2], cex = 0.5, col = "red")
title(paste0("nr: ", nr,
            "\n n: ", npoints,
            "\n p̂: ", round(res$phat[1],3), " ", round(res$phat[2], 3)
            ))
legend(
  "bottomright",
  legend = c("data", "inferred trajectory", "true trajectory", interp_names),
  col    = c("black", "#1f77b4", "red", interp_colors[seq_along(interp_names)]),
  pch    = c(1, NA, NA, rep(1, length(interp_names))),
  lty    = c(NA, 1, 1, rep(NA, length(interp_names))),
  xpd    = TRUE,
  bty    = "n",
  cex = 0.8
)
