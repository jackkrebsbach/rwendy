## Validate analytical J_wnll and H_wnll against numerical approximations
## using uGMAR::calc_gradient / calc_hessian

library(devtools)
library(deSolve)
library(uGMAR)

invisible(devtools::load_all())

# ── 1. Simulate logistic ODE ─────────────────────────────────────────────────
f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)

p_star  <- c(1, 1)
u0      <- c(0.01)
npoints <- 512
t_span  <- c(0, 10)
t_eval  <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) list(as.vector(f(state, parameters, tvec)))
sol      <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

set.seed(42)
noise_sd <- 0.1 * sqrt(mean(sol[, 2]^2))
U  <- matrix(sol[, 2] + rnorm(npoints, sd = noise_sd), ncol = 1)
tt <- sol[, 1, drop = FALSE]

# ── 2. Build WENDy system (skip optimization, evaluate at p_star) ─────────────
res <- solveWendy(f, p_star, U, tt,
                  lip    = TRUE,
                  method = "MLE",
                  control = list(optimize = FALSE))

wnll   <- res$wnll
J_wnll <- res$J_wnll
H_wnll <- res$H_wnll

# ── 3. Analytical values at p_star ───────────────────────────────────────────
cat("=== Analytical ===\n")
cat("wnll(p_star):  ", wnll(p_star), "\n")

grad_analytical <- J_wnll(p_star)
cat("J_wnll(p_star):", grad_analytical, "\n")

hess_analytical <- H_wnll(p_star)
cat("H_wnll(p_star):\n")
print(hess_analytical)

# ── 4. Numerical approximations via uGMAR ────────────────────────────────────
cat("\n=== Numerical (uGMAR) ===\n")

grad_numerical <- calc_gradient(p_star, wnll)
cat("calc_gradient: ", grad_numerical, "\n")

hess_numerical <- calc_hessian(p_star, wnll)
cat("calc_hessian:\n")
print(hess_numerical)

# ── 5. Differences ───────────────────────────────────────────────────────────
cat("\n=== Differences ===\n")
cat("Max |gradient error|:", max(abs(grad_analytical - grad_numerical)), "\n")
cat("Max |hessian error|: ", max(abs(hess_analytical - hess_numerical)), "\n")

cat("\nRelative gradient error:",
    max(abs((grad_analytical - grad_numerical) / (abs(grad_numerical) + 1e-12))), "\n")
cat("Relative hessian error: ",
    max(abs((hess_analytical - hess_numerical) / (abs(hess_numerical) + 1e-12))), "\n")
