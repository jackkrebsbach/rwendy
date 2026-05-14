# %%
devtools::load_all()
library(deSolve)
library(ggplot2)

f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)

p_star  <- c(1.0, 0.1)
u0      <- c(0.1)
n_pts   <- 256
t_span  <- c(0.0, 10.0)
tt      <- seq(t_span[1], t_span[2], length.out = n_pts)

ode_rhs <- function(tvec, state, parms) list(as.vector(f(state, parms, tvec)))
sol     <- deSolve::ode(y = u0, times = tt, func = ode_rhs, parms = p_star)
U_true  <- sol[, -1, drop = FALSE]

noise_ratio <- 0.4
sigma       <- noise_ratio * sqrt(mean(U_true^2))

cat(sprintf("σ = %.4f  (noise ratio = %.2f)\n", sigma, noise_ratio))

n_rep   <- 500
p_joint <- matrix(NA_real_, n_rep, 2)
p_mle   <- matrix(NA_real_, n_rep, 2)

p1_seq    <- seq(0.25, 1.75, length.out = 40)
p2_seq    <- seq(0.0, 0.3, length.out = 40)
grid_base <- expand.grid(p1 = p1_seq, p2 = p2_seq)
wnll_sum  <- numeric(nrow(grid_base))   # accumulates WNLL across realizations
wnll_n    <- 0L                         # number of successfully evaluated surfaces

for (i in seq_len(n_rep)) {
  set.seed(8675309 + i)
  U_noisy <- U_true + matrix(rnorm(length(U_true), 0, sigma), nrow = nrow(U_true))

  res_j <- tryCatch(
    solveWendy(f, U_noisy, tt, method = "JOINT", control = list(gn_method = "gls")),
    error = function(e) NULL
  )
  if (!is.null(res_j) && !anyNA(res_j$phat))
    p_joint[i, ] <- res_j$phat

  res_m <- tryCatch(
    solveWendy(f, U_noisy, tt, method = "MLE"),
    error = function(e) NULL
  )
  if (!is.null(res_m) && !anyNA(res_m$phat)) {
    p_mle[i, ] <- res_m$phat

    # Accumulate WNLL surface from this realization's system
    vals <- mapply(function(p1, p2) {
      tryCatch(as.numeric(res_m$wnll(c(p1, p2))), error = function(e) NA_real_)
    }, grid_base$p1, grid_base$p2)

    if (!anyNA(vals)) {
      wnll_sum <- wnll_sum + vals
      wnll_n   <- wnll_n + 1L
    }
  }

  if (i %% 50 == 0)
    cat(sprintf("  realization %d / %d\n", i, n_rep))
}

p_joint_mean <- colMeans(p_joint, na.rm = TRUE)
p_mle_mean   <- colMeans(p_mle,   na.rm = TRUE)

rce_joint <- apply(p_joint, 1, function(phat) {
  if (anyNA(phat)) return(NA_real_)
  norm(phat - p_star, type = "2") / norm(p_star, type = "2")
})
rce_mle <- apply(p_mle, 1, function(phat) {
  if (anyNA(phat)) return(NA_real_)
  norm(phat - p_star, type = "2") / norm(p_star, type = "2")
})

rce_joint_mean <- mean(rce_joint, na.rm = TRUE)
rce_mle_mean   <- mean(rce_mle,   na.rm = TRUE)

cat(sprintf("\np*          = [%.4f, %.4f]\n", p_star[1], p_star[2]))
cat(sprintf("JOINT mean  = [%.4f, %.4f]  bias = [%.4f, %.4f]  RCE = %.4f\n",
            p_joint_mean[1], p_joint_mean[2],
            p_joint_mean[1] - p_star[1], p_joint_mean[2] - p_star[2],
            rce_joint_mean))
cat(sprintf("MLE   mean  = [%.4f, %.4f]  bias = [%.4f, %.4f]  RCE = %.4f\n",
            p_mle_mean[1], p_mle_mean[2],
            p_mle_mean[1] - p_star[1], p_mle_mean[2] - p_star[2],
            rce_mle_mean))

grid          <- grid_base
grid$wnll     <- wnll_sum / wnll_n   # E[WNLL(p)] over all realizations

wnll_bounds    <- quantile(grid$wnll, c(0.01, 0.99), na.rm = TRUE)
grid$wnll_clip <- pmin(pmax(grid$wnll, wnll_bounds[1]), wnll_bounds[2])

pts <- data.frame(
  p1     = c(p_star[1],     p_joint_mean[1], p_mle_mean[1]),
  p2     = c(p_star[2],     p_joint_mean[2], p_mle_mean[2]),
  method = factor(c("True p*", "JOINT (mean)", "MLE (mean)"),
                  levels = c("True p*", "JOINT (mean)", "MLE (mean)")),
  label  = c(
    sprintf("(%.3f, %.3f)", p_star[1],       p_star[2]),
    sprintf("(%.3f, %.3f)  RCE=%.3f", p_joint_mean[1], p_joint_mean[2], rce_joint_mean),
    sprintf("(%.3f, %.3f)  RCE=%.3f", p_mle_mean[1],   p_mle_mean[2],   rce_mle_mean)
  )
)

ggplot(grid, aes(x = p1, y = p2)) +
  geom_tile(aes(fill = wnll_clip)) +
  scale_fill_viridis_c(option = "plasma", name = "WNLL") +
  geom_point(
    data  = pts,
    aes(x = p1, y = p2, colour = method, shape = method),
    size  = 5, stroke = 1.5
  ) +
  geom_text(
    data  = pts,
    aes(x = p1, y = p2, label = label, colour = method),
    vjust = -1.2, size = 3.5, fontface = "bold", show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("True p*" = "white", "JOINT (mean)" = "#00eeff", "MLE (mean)" = "#ffdd00")
  ) +
  scale_shape_manual(
    values = c("True p*" = 8, "JOINT (mean)" = 17, "MLE (mean)" = 16)
  ) +
  labs(
    title    = sprintf(
      "E[WNLL] — logistic (n = %d, noise ratio = %.2f, %d realizations)",
      n_pts, noise_ratio, n_rep),
    x        = expression(p[1]),
    y        = expression(p[2]),
    colour   = NULL,
    shape    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("examples/logistic_joint_vs_mle_wnll.png", width = 7, height = 6, dpi = 150)