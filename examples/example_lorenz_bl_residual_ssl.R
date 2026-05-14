# %%
# Validate that the production SSL + boundary-layer system delivers O(h^6)
# convergence of the augmented residual at (p*, U*).
#
# In SSL convention V has no h-factor, so V%*%F is in "sum" form (~ integral/h).
# The boundary terms and EM correction are scaled by bdry_scale = 1/h to match.
# Concretely, sys$g(U,p,tt) - sys$build_compute_b(U,tt) returns the augmented
# residual in SSL sum form: r_SSL[BL] = O(h^5).
# Multiplying by h recovers the integral-form residual: h * r_SSL[BL] = O(h^6).
#
# The script samples on grids aligned to r_time (so the right edge of phi's
# support lands on a grid point) and reports both slopes.

library(deSolve)
library(ggplot2)
library(devtools)
invisible({ devtools::load_all("/Users/jackkrebsbach/Documents/ml/rwendy") })

f <- function(u, p, t) c(
  p[1] * (u[2] - u[1]),
  u[1] * (p[2] - u[3]) - u[2],
  u[1] * u[2] - p[3] * u[3]
)

p_star  <- c(10.0, 28.0, 8.0 / 3.0)
u0_true <- c(-8, 10, 27)
t_span  <- c(0, 10)
D       <- 3L

modelODE <- function(tvec, state, parms) list(as.vector(f(state, parms, tvec)))

# Fix BL support in TIME and choose dt so r_time / dt is an integer.
r_time      <- 1.0
n_per_rtime <- c(10L, 20L, 40L, 80L, 160L, 320L)
dt_vec      <- r_time / n_per_rtime
npoints_vec <- as.integer((t_span[2] - t_span[1]) / dt_vec + 1)

err_sum <- numeric(length(dt_vec))      # ||r_SSL[BL]||_F  -- sum form
err_int <- numeric(length(dt_vec))      # ||h * r_SSL[BL]||_F -- integral form
dts     <- dt_vec

for (i in seq_along(dt_vec)) {
  dt      <- dt_vec[i]
  npoints <- npoints_vec[i]
  t_eval  <- seq(t_span[1], t_span[2], length.out = npoints)
  radius_samples <- as.integer(n_per_rtime[i])      # samples * dt = r_time

  sol <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                      parms = p_star, rtol = 1e-15, atol = 1e-15)
  U  <- sol[, -1]
  tt <- matrix(sol[, 1], ncol = 1)

  res <- solveWendy(f, U, tt, method = "IRLS",
                    control = list(optimize = FALSE, estimate_u0 = FALSE,
                                   test_fun_type = "SSL",
                                   include_boundary_layer = TRUE,
                                   fixed_radius = radius_samples))

  sys    <- res$system_joint
  prob   <- res$wendy_problems_joint[[1]]
  K_int  <- prob$bl$K_interior
  K_bl   <- prob$bl$K_bl
  K_tot  <- K_int + K_bl
  tt_vec <- as.vector(res$tt)

  b_vec <- sys$build_compute_b(U, tt_vec)
  g_vec <- sys$g(U, p_star, tt_vec)
  r_vec <- as.numeric((g_vec - b_vec)$cpu())          # length K_tot * D, row-major

  r_mat <- matrix(r_vec, nrow = K_tot, ncol = D, byrow = TRUE)
  r_bl  <- r_mat[(K_int + 1):(K_int + K_bl), , drop = FALSE]

  err_sum[i] <- norm(r_bl,        "F")
  err_int[i] <- norm(dt * r_bl,   "F")

  cat(sprintf("dt = %.6f  radius = %d  K_int = %d  K_bl = %d  ||r_BL|| = %.3e  ||h*r_BL|| = %.3e\n",
              dt, radius_samples, K_int, K_bl, err_sum[i], err_int[i]))
}

slope_sum <- unname(coef(lm(log(err_sum) ~ log(dts)))[2])
slope_int <- unname(coef(lm(log(err_int) ~ log(dts)))[2])
cat(sprintf("\nSlopes:  sum form (||r_BL||) = %.3f  (expected 5)\n", slope_sum))
cat(sprintf("         integral form (||h*r_BL||) = %.3f  (expected 6)\n", slope_int))

# Plot integral-form residual against O(h^4), O(h^6) reference
h2_ref <- err_int[1] * (dts / dts[1])^2
h4_ref <- err_int[1] * (dts / dts[1])^4
h6_ref <- err_int[1] * (dts / dts[1])^6

df     <- data.frame(dt = dts, error = err_int)
df_ref <- data.frame(
  dt    = rep(dts, 3),
  value = c(h2_ref, h4_ref, h6_ref),
  rate  = factor(rep(c("O(h^2)", "O(h^4)", "O(h^6)"), each = length(dts)),
                 levels = c("O(h^2)", "O(h^4)", "O(h^6)"))
)

ggplot() +
  geom_line(data = df_ref, aes(x = dt, y = value, linetype = rate),
            color = "gray40", alpha = 0.8, linewidth = 0.6) +
  geom_line(data = df, aes(x = dt, y = error), color = "steelblue", linewidth = 0.9) +
  geom_point(data = df, aes(x = dt, y = error),
             color = "steelblue", fill = "white", shape = 21, size = 3, stroke = 1.2) +
  scale_x_log10(breaks = dts, labels = scales::label_scientific()(dts)) +
  scale_y_log10(labels = scales::label_scientific()) +
  scale_linetype_manual(values = c("O(h^2)" = "dashed",
                                   "O(h^4)" = "dotted",
                                   "O(h^6)" = "longdash")) +
  annotation_logticks(sides = "bl") +
  labs(
    title    = sprintf("SSL + BL augmented residual (integral form): slope %.2f", slope_int),
    subtitle = sprintf("Lorenz, p = p*, U = U_true (clean), r_time = %.2f", r_time),
    x        = "Step size  h",
    y        = expression("||" * h %.% r[BL] * "||"[F]),
    linetype = "reference"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())
