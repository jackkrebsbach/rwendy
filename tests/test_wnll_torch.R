library(deSolve)
library(torch)
library(numbers)
# Source package internals so we can access build_wnll_torch
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)

torch_set_default_dtype(torch_float64())

cat("=== Test: build_wnll vs build_wnll_torch ===\n\n")

# Set up a simple logistic ODE problem
f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

p_star <- c(1, 1)
u0 <- c(0.01)
p0 <- c(0.5, 0.5)
npoints <- 100
t_span <- c(0.005, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) {
  list(as.vector(f(state, parameters, tvec)))
}
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

set.seed(42)
nr <- 0.15
U_vec <- as.vector(sol[, -1])
noise_sd <- nr * sqrt(mean(U_vec^2))
U <- matrix(c(sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)), ncol = 1)
tt <- sol[, 1, drop = FALSE]

# Run solveWendy without optimization to get all the components
res <- solveWendy(f, p0, U, tt, lip = TRUE, method = "MLE",
  noise_dist = "addgaussian",
  control = list(test_fun_type = "SSL", optimize = FALSE))

# Build both versions
wnll_R <- res$wnll  # original (R linear algebra)
wnll_torch <- build_wnll_torch(res$S, res$g, res$b, nrow(res$V), ncol(U))

# Test at multiple parameter values
test_params <- list(
  p0,
  c(1.0, 1.0),
  c(0.1, 0.1),
  c(2.0, 0.5),
  c(0.8, 1.2)
)

cat("Comparing wnll values at different parameter vectors:\n\n")
all_pass <- TRUE

for (i in seq_along(test_params)) {
  p <- test_params[[i]]
  val_R <- wnll_R(p)
  val_torch <- wnll_torch(p)
  diff <- abs(val_R - val_torch)
  rel_err <- diff / max(abs(val_R), 1e-15)
  pass <- rel_err < 1e-6

  cat(sprintf("  p = [%s]: R = %.10f, Torch = %.10f, rel_err = %.2e %s\n",
    paste(p, collapse = ", "), val_R, val_torch, rel_err,
    ifelse(pass, "PASS", "FAIL")))

  if (!pass) all_pass <- FALSE
}

# Benchmark
cat("\n--- Timing comparison (100 evaluations) ---\n")
p_test <- c(0.8, 1.2)

t_R <- system.time({
  for (i in 1:100) wnll_R(p_test)
})

t_torch <- system.time({
  for (i in 1:100) wnll_torch(p_test)
})

cat(sprintf("  R version:     %.4f s\n", t_R["elapsed"]))
cat(sprintf("  Torch version: %.4f s\n", t_torch["elapsed"]))
cat(sprintf("  Speedup:       %.2fx\n", t_R["elapsed"] / t_torch["elapsed"]))

cat(sprintf("\n=== %s ===\n", ifelse(all_pass, "ALL TESTS PASSED", "SOME TESTS FAILED")))
