
# %%
# library(wendy)
library(deSolve)
library(uGMAR)
library(devtools)
library(ggplot2)

options(wendy.symbolic_backend = "native")

invisible({devtools::load_all()})

# options(wendy.symbolic_backend = "symengine")

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1/10)
u0 <- c(0.1)
p0 <- c(1.25, 0.25)
npoints <- 128
t_span <- c(0.0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-14)

set.seed(8675309 + 3)

nr <- 0.15
U_vec <- as.vector(sol[,-1])

# Additive Gaussian Noise
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- rnorm(npoints, mean = 0, sd = noise_sd)
U <- sol[, 2, drop = FALSE] + noise
tt <- sol[, 1, drop = FALSE]

# Multiplicative Lognormal Noise
# noise_sd <- nr
# noise <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))

cat(sprintf("σ = %.2f", noise_sd))

time <- system.time({
res <- solveWendy(f, U, tt, method = "IRLS",
  control = list(
    estimate_IC         = TRUE,
    estimate_trajectory = TRUE
))
})

t_eval_dense <- seq(t_span[1], t_span[2], length.out = npoints);
sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

U_star    <- res$state$U_star
u0hat     <- res$u0hat
u0_smooth <- U_star[1, 1]
tt_vec    <- as.vector(tt)

# Total posterior SE from ERTS: conditional-on-p̂ posterior + parameter
# uncertainty folded in via the trajectory sensitivity (law of total variance).
se_post   <- sqrt(res$state$P_smooth[, 1, 1])

band_post_lo  <- U_star[, 1] - 2 * se_post
band_post_hi  <- U_star[, 1] + 2 * se_post

ylim <- range(c(U, band_post_lo, band_post_hi), finite = TRUE)
plot(t_eval_dense, sol_true[, 2], col = "red", type = "l",
     xlab = "Time", ylab = "u₁", ylim = ylim)
polygon(c(tt_vec, rev(tt_vec)), c(band_post_lo, rev(band_post_hi)),
        col = adjustcolor("#1f77b4", alpha.f = 0.25), border = NA)
lines(tt_vec, U_star[, 1], col = "#1f77b4", lwd = 2)
points(tt, U, col = "black", cex = 0.5)
points(tt[1], u0hat,     pch = 2, col = "#1f77b4", cex = 1)
points(tt[1], u0_smooth, pch = 2, col = "#ff7f0e", cex = 1)
points(tt[1], U[1,], pch = 2, col = "green", cex = 1)

title(paste0("nr: ", nr, "\n n: ", npoints, "\n p̂: ", round(res$phat[1],3), " ", round(res$phat[2], 3)))

legend(
  "bottomright",
  legend = c("true trajectory", "smoothed state (ERTS)",
             "±2 SE (ERTS posterior)",
             "û₀ BLDC", "û₀ ERTS ",
              "u₀ Noisy"),
  col    = c("red", "#1f77b4",
             adjustcolor("#1f77b4", alpha.f = 0.6),
             "#1f77b4", "#ff7f0e", "green"),
  pch    = c(NA, NA, 15, 17, 17, 17),
  lty    = c(1,  1,  NA, NA, NA,NA),
  xpd    = TRUE,
  bty    = "n",
  cex = 0.7
)

cat(sprintf("\np̂_IRLS  = [%s]  rel_err = %.4f",
            paste(sprintf("%.4f", res$phat), collapse = ", "),
            rel_err(res$phat, p_star)))