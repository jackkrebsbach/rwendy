
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
npoints <- 120
t_span <- c(0.0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-14)

# set.seed(8675309 + 4)

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
  res1 <- solveWendy(f, U, tt, method = "IRLS",
  control = list(
    estimate_IC         = TRUE,
    estimate_trajectory = TRUE,
    test_fun_type   = "MSG"
  ))
})

t_eval_dense <- seq(t_span[1], t_span[2], length.out = npoints);
sol_true <- deSolve::ode(y = u0, times = t_eval_dense, func = modelODE, parms = p_star)

U_star    <- res$state$U_star
u0hat     <- res$u0hat
u0_smooth <- U_star[1, 1]
tt_vec    <- as.vector(tt)

# Posterior SE from ERTS (conditional on p̂).
# se_post   <- sqrt(res$state$P_smooth[, 1, 1])

# band_post_lo  <- U_star[, 1] - 2 * se_post
# band_post_hi  <- U_star[, 1] + 2 * se_post

# ylim <- range(c(U, band_post_lo, band_post_hi), finite = TRUE)
# plot(t_eval_dense, sol_true[, 2], col = "red", type = "l",
#      xlab = "Time", ylab = "u₁", ylim = ylim)
# polygon(c(tt_vec, rev(tt_vec)), c(band_post_lo, rev(band_post_hi)),
#         col = adjustcolor("#1f77b4", alpha.f = 0.25), border = NA)
# lines(tt_vec, U_star[, 1], col = "#1f77b4", lwd = 2)
# points(tt, U, col = "black", cex = 0.5)
# points(tt[1], u0hat,     pch = 2, col = "#1f77b4", cex = 1.5)
# points(tt[1], u0_smooth, pch = 2, col = "#ff7f0e", cex = 1.5)

# title(paste0("nr: ", nr, "\n n: ", npoints, "\n p̂: ", round(res$phat[1],3), " ", round(res$phat[2], 3)))

# legend(
#   "bottomright",
#   legend = c("true trajectory", "smoothed state (ERTS)",
#              "±2 SE (ERTS posterior)",
#              "û₀ (BL fit)", "û₀ (ERTS smoothed)"),
#   col    = c("red", "#1f77b4",
#              adjustcolor("#1f77b4", alpha.f = 0.6),
#              "#1f77b4", "#ff7f0e"),
#   pch    = c(NA, NA, 15, 17, 17),
#   lty    = c(1,  1,  NA, NA, NA),
#   xpd    = TRUE,
#   bty    = "n",
#   cex = 0.8
# )


# cat(sprintf("\nû₀ = [%s]  abs_err = %.4f",
#             paste(sprintf("%.4f", res$u0hat), collapse = ", "),
#             abs(res$u0hat - u0)))

# cat(sprintf("\nERTS û₀ = [%s]  abs_err = %.4f\n",
#             paste(sprintf("%.4f", res$state$U_star[1,]), collapse = ", "),
#             abs(res$state$U_star[1,] - u0)))

# resOE <- solveWendy(f, U, tt, p0 = p0, method = "OE")

# cat(sprintf("\nOE p̂ = [%s]  rel_err = %.4f",
#             paste(sprintf("%.4f", resOE$phat), collapse = ", "),
#             rel_err(resOE$phat, p_star)))

# cat(sprintf("\nOE û₀ = [%s]  abs_err = %.4f\n",
#             paste(sprintf("%.4f", resOE$data$u0), collapse = ", "),
#             abs(resOE$data$u0 - u0)))

# print(summary(res)$param_cov)

# dt <- mean(diff(tt))
# V <- res$V / dt
# VT <- t(V) 

# I <- V %*% VT

# res <- solveWendy(f, U, tt, method = "IRLS",
#   control = list(
#     optimize = FALSE,
#     test_fun_type   = "SSL",
#     include_boundary_layer = TRUE,
#     n_bl = 10
#   ))

# V <- t(svd(res$V)$v)
# c <- V %*% sol[,-1]

# plot(t(V) %*% c)

# tt_vec <- as.vector(tt)
# mp1    <- nrow(U)

# bl_left_idx  <- which(V[, 1]   != 0)
# bl_right_idx <- which(V[, mp1] != 0)
# bl_idx       <- sort(union(bl_left_idx, bl_right_idx))

# ylim <- range(V[bl_idx, , drop = FALSE], finite = TRUE)
# plot(NA, xlim = range(tt_vec), ylim = ylim,
#      xlab = "Time", ylab = expression(phi[k]),
#      main = sprintf("Boundary-layer test functions (r_c = %d, K_bl = %d)",
#                     res$rc, length(bl_idx)))
# abline(h = 0, col = "grey80")
# for (k in bl_left_idx)  lines(tt_vec, V[k, ], col = "#1f77b4", lwd = 1.5)
# for (k in bl_right_idx) lines(tt_vec, V[k, ], col = "#d62728", lwd = 1.5)
# legend("top", legend = c("left BL", "right BL"), horiz = TRUE,
#        col = c("#1f77b4", "#d62728"), lwd = 1.5, bty = "n")

res <- solveWendy(f, U, tt, method = "JOINT",
                    control = list(
                      joint_rank_tau = 0,
                      joint_deriv_pen = 1e-4,
                      joint_lambda = 1000,
                      joint_deriv_ode = FALSE,
                      joint_basis_rank = "full"
                    )
                )

cat(sprintf("\np̂_JOINT = [%s]  rel_err = %.4f",
            paste(sprintf("%.4f", res$phat), collapse = ", "),
            rel_err(res$phat, p_star)))

cat(sprintf("\np̂_IRLS  = [%s]  rel_err = %.4f\n",
            paste(sprintf("%.4f", res1$phat), collapse = ", "),
            rel_err(res1$phat, p_star)))

state <- res$data$uhat

plot(tt, U[,1], col = adjustcolor("blue", alpha.f = 0.3), cex = 0.5)
lines(tt, state[,1], cex = 0.5, col = "blue")
lines(tt, sol[,-1], cex = 0.5, col = "black")