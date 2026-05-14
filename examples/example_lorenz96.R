
# %%
# Lorenz 96 system
# du_i/dt = (u_{i+1} - u_{i-2}) * u_{i-1} - u_i + F,  i = 1..N (cyclic)
# Single unknown parameter: forcing constant F (p_star = 8)

# library(wendy)
library(deSolve)
library(ggplot2)
library(devtools)

invisible(devtools::load_all())

N <- 40  # number of dimensions

f <- function(u, p, t) {
  F_forcing <- p[1]
  du <- numeric(N)
  for (i in seq_len(N)) {
    ip1 <- ((i    ) %% N) + 1   # i+1 mod N
    im1 <- ((i - 2) %% N) + 1   # i-1 mod N
    im2 <- ((i - 3) %% N) + 1   # i-2 mod N
    du[i] <- (u[ip1] - u[im2]) * u[im1] - u[i] + F_forcing
  }
  du
}

p_star  <- c(8.0)
u0      <- rep(p_star, N); u0[1] <- u0[1] + 0.01   # small perturbation from steady state
npoints <- 256
t_span  <- c(0, 4)
tt      <- matrix(seq(t_span[1], t_span[2], length.out = npoints), ncol = 1)

modelODE <- function(tvec, state, parms) list(as.vector(f(state, parms, tvec)))
sol      <- deSolve::ode(y = u0, times = as.vector(tt), func = modelODE, parms = p_star,
                         rtol = 1e-10, atol = 1e-10)
U_true   <- sol[, -1]

set.seed(123)
nr       <- 0.05
noise_sd <- nr * sqrt(mean(U_true^2))
U        <- U_true + matrix(rnorm(npoints * N, sd = noise_sd), nrow = npoints)

cat(sprintf("N = %d   n = %d   nr = %.2f   σ = %.4f\n", N, npoints, nr, noise_sd))

# --- Solve ---
res_irls <- solveWendy(f, U, tt, method = "IRLS")
res_root <- solveWendy(f, U, tt, method = "JOINT")

cat(sprintf("p*    = %.4f\n",             p_star))
cat(sprintf("IRLS  p̂ = %.4f   rel_err = %.4f\n", res_irls$phat, rel_err(res_irls$phat, p_star)))
cat(sprintf("JOINT  p̂ = %.4f   rel_err = %.4f\n", res_root$phat, rel_err(res_root$phat, p_star)))

# --- Reconstructed trajectory from p̂ ---
sol_irls <- deSolve::ode(y = u0, times = as.vector(tt), func = modelODE, parms = res_irls$phat)[, -1]
sol_root <- deSolve::ode(y = u0, times = as.vector(tt), func = modelODE, parms = res_root$phat)[, -1]

# --- Plot first two dimensions ---
df <- data.frame(
  t      = rep(as.vector(tt), 4),
  u      = c(U_true[, 1], U[, 1], sol_irls[, 1], sol_root[, 1]),
  label  = rep(c("True", "Noisy", "IRLS", "JOINT"), each = npoints)
)
df$label <- factor(df$label, levels = c("Noisy", "True", "IRLS", "JOINT"))

ggplot(df, aes(x = t, y = u, colour = label)) +
  geom_point(data = subset(df, label == "Noisy"),
             size = 0.6, alpha = 0.5) +
  geom_line(data = subset(df, label != "Noisy"),
            linewidth = 0.8) +
  scale_colour_manual(values = c(Noisy = "grey70", True = "steelblue",
                                 IRLS = "forestgreen", JOINT = "firebrick")) +
  labs(title = sprintf("Lorenz 96 (N=%d) — u₁ trajectory", N),
       subtitle = sprintf("p* = %.1f  |  IRLS p̂ = %.4f  |  JOINT p̂ = %.4f",
                          p_star, res_irls$phat, res_root$phat),
       x = "Time", y = "u₁", colour = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
