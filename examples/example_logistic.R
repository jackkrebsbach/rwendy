
# %%
# library(wendy)
library(base)
library(MASS)
library(deSolve)
library(numDeriv)
library(uGMAR)
library(devtools)
library(ggplot2)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1)
u0 <- c(0.01)
p0 <- c(1.25, 0.25)
npoints <- 128
t_span <- c(0.0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-14)

# set.seed(8675309 + 16)

nr <- 0.3
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

t_eval_dense <- seq(t_span[1], t_span[2], length.out = 501L)
sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

time <- system.time({
  res <- solveWendy(f = f, U, tt, method = "IRLS",  control = list(estimate_IC = TRUE, estimate_trajectory = TRUE))
})

U_star    <- res$state$U_star
u0hat     <- res$u0hat
u0_smooth <- U_star[1, 1]
tt_vec    <- as.vector(tt)

se_post <- sqrt(res$state$P_smooth[, 1, 1])
band_lo <- U_star[, 1] - 2 * se_post
band_hi <- U_star[, 1] + 2 * se_post

u0_estimates <- c(u0_smooth, u0hat)
u0_se        <- c(se_post[1], sqrt(res$boundary_state$cov_u0[1, 1]))
u0_colors    <- c("#ff7f0e", "#9467bd")

plot_trajectory <- function(xlim, ylim) {
  plot(NA, xlim = xlim, ylim = ylim, xlab = "Time", ylab = "u₁")
  polygon(c(tt_vec, rev(tt_vec)), c(band_lo, rev(band_hi)),
          col = adjustcolor("#1f77b4", alpha.f = 0.25), border = NA)
  lines(t_eval_dense, sol_true[, 2], col = "red")
  lines(tt_vec, U_star[, 1], col = "#1f77b4", lwd = 2)
  points(tt_vec, U, col = "black", cex = 0.75)
}

old_par <- par(mfrow = c(1, 2))

plot_trajectory(range(tt_vec), range(c(U, band_lo, band_hi), finite = TRUE))
points(rep(tt_vec[1], 2), u0_estimates, pch = 17, col = u0_colors, cex = 1.2)
title(paste0("nr: ", nr, "\n n: ", npoints, "\n p̂: ", round(res$phat[1], 3)))
legend(
  "bottomright",
  legend = c("true trajectory", "data", "smoothed state (ERTS)",
             "±2 SE (ERTS posterior)", "û₀ ERTS", "û₀ estimate_IC"),
  col    = c("red", "black", "#1f77b4",
             adjustcolor("#1f77b4", alpha.f = 0.6), u0_colors),
  pch    = c(NA, 1, NA, 15, 17, 17),
  lty    = c(1, NA, 1, NA, NA, NA),
  xpd    = TRUE,
  bty    = "n",
  cex    = 0.7
)

zoom_width <- 0.1 * diff(t_span)
dodge      <- 0.03 * zoom_width
u0_t       <- tt_vec[1] + c(-dodge, dodge)
u0_lo      <- u0_estimates - 2 * u0_se
u0_hi      <- u0_estimates + 2 * u0_se
in_zoom    <- tt_vec <= tt_vec[1] + zoom_width

plot_trajectory(c(tt_vec[1] - 2 * dodge, tt_vec[1] + zoom_width),
                range(c(band_lo[in_zoom], band_hi[in_zoom], u0_lo, u0_hi), finite = TRUE))
arrows(u0_t, u0_lo, u0_t, u0_hi, angle = 90, code = 3, length = 0.04, col = u0_colors, lwd = 1.5)
points(u0_t, u0_estimates, pch = 17, col = u0_colors, cex = 1.2)
title("Initial condition, ±2 SE")

par(old_par)

cat(sprintf("\n\np*                : %s\n", paste(sprintf("%.4f", p_star), collapse = "  ")))
cat(sprintf("p_hat WENDy: %s\n", paste(sprintf("%.4f", res$phat), collapse = "  ")))
cat(sprintf("u0 true          : %.4f\n", u0))
cat(sprintf("u0hat estimate_IC: %.4f (SE %.4f)\n", u0hat, u0_se[2]))
cat(sprintf("u0hat ERTS       : %.4f (SE %.4f)\n", u0_smooth, u0_se[1]))
cat("\nIRLS Wall Time:", time[["elapsed"]],"seconds")