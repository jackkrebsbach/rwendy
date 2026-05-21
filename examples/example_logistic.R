
# %%
# library(wendy)
library(deSolve)
library(uGMAR)
library(devtools)
library(ggplot2)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1/10)
u0 <- c(0.1)
p0 <- c(1.25, 0.25)
npoints <- 256
t_span <- c(0.0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-14)

set.seed(8675309 + 4)

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

time <- system.time({
  res <- solveWendy(f, U, tt, method = "IRLS",
  control = list(
    estimate_u0     = TRUE,
    estimate_U_star = TRUE,
    smoother        = "erts",
    test_fun_type   = "MSG"
  ))
})

t_eval_dense <- seq(t_span[1], t_span[2], length.out = npoints);
sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

U_star    <- res$state$U_star
u0hat     <- res$u0hat
u0_smooth <- U_star[1, 1]
tt_vec    <- as.vector(tt)

# Posterior SE from ERTS (conditional on p̂) and the wider band including
# parameter uncertainty propagated via S(t) Ĉ S(t)^T.
se_post   <- sqrt(res$state$P_smooth[, 1, 1])
P_total   <- res$state$P_total %||% res$state$P_smooth
se_total  <- sqrt(P_total[, 1, 1])

band_post_lo  <- U_star[, 1] - 2 * se_post
band_post_hi  <- U_star[, 1] + 2 * se_post
band_total_lo <- U_star[, 1] - 2 * se_total
band_total_hi <- U_star[, 1] + 2 * se_total

interp_colors <- c("purple", "darkorange", "forestgreen", "deeppink", "cyan4")
problem_names <- names(res$wendy_data)

ylim <- range(c(U, band_total_lo, band_total_hi), finite = TRUE)
plot(t_eval_dense, sol_true[, 2], col = "red", type = "l",
     xlab = "Time", ylab = "u₁", ylim = ylim)
# Outer band (total: posterior + parameter sensitivity), then inner (ERTS only).
polygon(c(tt_vec, rev(tt_vec)), c(band_total_lo, rev(band_total_hi)),
        col = adjustcolor("#ff7f0e", alpha.f = 0.18), border = NA)
polygon(c(tt_vec, rev(tt_vec)), c(band_post_lo, rev(band_post_hi)),
        col = adjustcolor("#1f77b4", alpha.f = 0.25), border = NA)
lines(tt_vec, U_star[, 1], col = "#1f77b4", lwd = 2)
points(tt[1], u0hat,     pch = 17, col = "#1f77b4", cex = 1.5)
points(tt[1], u0_smooth, pch = 17, col = "#ff7f0e", cex = 1.5)

title(paste0("nr: ", nr, "\n n: ", npoints, "\n p̂: ", round(res$phat[1],3), " ", round(res$phat[2], 3)))

legend(
  "bottomright",
  legend = c("true trajectory", "smoothed state (ERTS)",
             "±2 SE (ERTS posterior)", "±2 SE (incl. parameter unc.)",
             "û₀ (BL fit)", "û₀ (ERTS smoothed)"),
  col    = c("red", "#1f77b4",
             adjustcolor("#1f77b4", alpha.f = 0.6),
             adjustcolor("#ff7f0e", alpha.f = 0.6),
             "#1f77b4", "#ff7f0e"),
  pch    = c(NA, NA, 15, 15, 17, 17),
  lty    = c(1,  1,  NA, NA, NA, NA),
  xpd    = TRUE,
  bty    = "n",
  cex = 0.8
)

cat(sprintf("\np̂ = [%s]  rel_err = %.4f",
            paste(sprintf("%.4f", res$phat), collapse = ", "),
            rel_err(res$phat, p_star)))

cat(sprintf("\nû₀ = [%s]  abs_err = %.4f\n",
            paste(sprintf("%.4f", res$u0hat), collapse = ", "),
            abs(res$u0hat - u0)))

cat(sprintf("\nERTS û₀ = [%s]  rel_err = %.4f\n",
            paste(sprintf("%.4f", res$state$U_star[1,]), collapse = ", "),
            rel_err(res$state$U_star[1,], u0)))

# resOE <- solveWendy(f, U, tt, p0 = p0, method = "OE")

# cat(sprintf("\nOE p̂ = [%s]  rel_err = %.4f",
#             paste(sprintf("%.4f", resOE$phat), collapse = ", "),
#             rel_err(resOE$phat, p_star)))

# cat(sprintf("\nOE û₀ = [%s]  abs_err = %.4f\n",
#             paste(sprintf("%.4f", resOE$data$u0), collapse = ", "),
#             abs(resOE$data$u0 - u0)))