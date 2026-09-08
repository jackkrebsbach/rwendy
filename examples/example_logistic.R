
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
npoints <- 30
t_span <- c(0.0, 10)
# Nonuniform observation times, denser near the start, retaining both endpoints.
t_fraction <- seq(0, 1, length.out = npoints)^2
t_eval <- t_span[1] + diff(t_span) * t_fraction

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-14)

# set.seed(8675309 + 16)

nr <- 0.05
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

time_j <- system.time({ 
  resj <- solveWendyGP(f, U, tt) 
})

U_star    <- res$state$U_star
u0hat     <- res$u0hat
u0_smooth <- U_star[1, 1]
tt_vec    <- as.vector(tt)

t_eval_dense <- seq(t_span[1], t_span[2], length.out = 501L)
sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

U_star    <- res$state$U_star
u0hat     <- res$u0hat
u0_smooth <- U_star[1, 1]
tt_vec    <- as.vector(tt)

# Total posterior SE from ERTS: conditional-on-p̂ posterior + parameter
se_post   <- sqrt(res$state$P_smooth[, 1, 1])

band_post_lo  <- U_star[, 1] - 2 * se_post
band_post_hi  <- U_star[, 1] + 2 * se_post

ylim <- range(c(U, band_post_lo, band_post_hi), finite = TRUE)
plot(t_eval_dense, sol_true[, 2], col = "red", type = "l",
     xlab = "Time", ylab = "u₁", ylim = ylim)
# polygon(c(tt_vec, rev(tt_vec)), c(band_post_lo, rev(band_post_hi)),
        # col = adjustcolor("#1f77b4", alpha.f = 0.25), border = NA)
# lines(tt_vec, U_star[, 1], col = "#1f77b4", lwd = 2)
points(tt, U, col = "black", cex = 0.5)
lines(resj$tt, resj$U_hat[, 1], col = "#2ca02c", lwd = 2, lty = 2)
# points(tt[1], u0hat,     pch = 2, col = "#1f77b4", cex = 1)
# points(tt[1], u0_smooth, pch = 2, col = "#ff7f0e", cex = 1)
# points(resj$tt[1], resj$U_hat[1, 1], pch = 2, col = "#2ca02c", cex = 1.2)
# points(tt[1], U[1,], pch = 2, col = "green", cex = 1.5)

title(paste0("nr: ", nr, "\n n: ", npoints, "\n p̂: ", round(res$phat[1],3), "\n p̂_GP: ", round(resj$phat[1],3), " "))

legend(
  "bottomright",
  legend = c("true trajectory", "smoothed state (ERTS)",
             "±2 SE (ERTS posterior)", "joint GP state",
             "û₀ BLDC", "û₀ ERTS ", "û₀ joint GP",
              "u₀ Noisy"),
  col    = c("red", "#1f77b4",
             adjustcolor("#1f77b4", alpha.f = 0.6), "#2ca02c",
             "#1f77b4", "#ff7f0e", "#2ca02c", "green"),
  pch    = c(NA, NA, 15, NA, 17, 17, 17, 17),
  lty    = c(1,  1,  NA, 2,  NA, NA, NA, NA),
  xpd    = TRUE,
  bty    = "n",
  cex = 0.7
)

res1 <- solveWendy(f = f, U, tt, method = "IRLS", control = list(estimate_IC = FALSE, estimate_trajectory = TRUE))

#three-way comparison
u0_true  <- sol[1, 2]
U_true_g <- deSolve::ode(y = u0, times = resj$tt, func = modelODE,
                         parms = p_star, rtol = 1e-12, atol = 1e-14)[, 2]
u0_gp    <- resj$U_hat[1, 1]
rmse <- function(a, b) sqrt(mean((a - b)^2))

cat(sprintf("\n\np*                : %s\n", paste(sprintf("%.4f", p_star), collapse = "  ")))
cat(sprintf("p_hat WENDy IRLS  : %s\n", paste(sprintf("%.4f", res$phat), collapse = "  ")))
cat(sprintf("p_hat joint GP    : %s\n", paste(sprintf("%.4f", resj$phat), collapse = "  ")))
cat(sprintf("per-parameter rel err (WENDy | GP): %s\n",
            paste(sprintf("%.4f | %.4f", abs(res$phat - p_star)/abs(p_star),
                          abs(resj$phat - p_star)/abs(p_star)), collapse = "   ")))

cat(sprintf("\n%-22s %10s %10s %10s %10s %8s\n",
            "arm", "p rel err", "u0", "u0 rel err", "state RMSE", "time_s"))
cat(sprintf("%-22s %10s %10.5f %10.4f %10.4f %8s\n", "raw data", "-", U[1, 1],
            abs(U[1, 1] - u0_true)/u0_true, rmse(U[, 1], sol[, 2]), "-"))
cat(sprintf("%-22s %10.4f %10.5f %10.4f %10.4f %8.2f\n", "WENDy IRLS + ERTS",
            rel_err(res$phat, p_star), u0_smooth, abs(u0_smooth - u0_true)/u0_true,
            rmse(res$state$U_star[, 1], sol[, 2]), time[["elapsed"]]))
cat(sprintf("%-22s %10s %10.5f %10.4f %10s %8s\n", "  + estimate_IC", "\"",
            res$u0hat, abs(res$u0hat - u0_true)/u0_true, "\"", "\""))
cat(sprintf("%-22s %10.4f %10.5f %10.4f %10.4f %8.2f\n", "joint GP",
            rel_err(resj$phat, p_star), u0_gp, abs(u0_gp - u0_true)/u0_true,
            rmse(resj$U_obs_hat[, 1], sol[, 2]), time_j[["elapsed"]]))

# state RMSE above is on the observation times for every arm; the joint GP also
# carries a denser working grid, scored here against the true trajectory on it
cat(sprintf("\njoint GP state RMSE on its own %d-point grid: %.4f\n",
            length(resj$tt), rmse(resj$U_hat[, 1], U_true_g)))
cat(sprintf("joint GP converged: %s (%s), grid accuracy passed: %s\n",
            resj$converged, resj$convergence_reason, resj$diagnostics$grid_passed))
cat(sprintf("\n(true u0 = %.5f, sigma = %.3f, nr = %.2f, n = %d)\n",
            u0_true, noise_sd, nr, npoints))

# The observation times are fixed, but noise is a fresh random draw each run.
# Repeated seeded runs are needed to compare estimation errors and timings.


cat("\nIRLS Wall Time:", time[["elapsed"]],"seconds")
cat("\nJoint GP Wall Time:", time_j[["elapsed"]], "seconds")
