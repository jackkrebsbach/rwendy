library(deSolve)
library(torch)
library(numbers)
for (ff in list.files("R", full.names = TRUE, pattern = "[.]R$")) source(ff)

torch_set_default_dtype(torch_float64())

cat("=== Test: build_wnll vs build_wnll_torch (Lorenz, D=3, J=3) ===\n\n")

f <- function(u, p, t) {
  du1 <- p[1] * (u[2] - u[1])
  du2 <- u[1] * (p[2] - u[3]) - u[2]
  du3 <- u[1] * u[2] - p[3] * u[3]
  c(du1, du2, du3)
}

p_star <- c(10.0, 28.0, 8.0 / 3.0)
p0 <- c(12.0, 21, 4.0)
u0 <- c(-8, 10, 27)
npoints <- 256
t_eval <- seq(0, 10, length.out = npoints)

modelODE <- function(tvec, state, parameters) {
  list(as.vector(f(state, parameters, tvec)))
}
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

set.seed(42)
nr <- 0.05
U_vec <- as.array(sol[, -1])
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)
U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

# --- Test nonlinear (lip = FALSE) ---
cat("--- Nonlinear case (lip = FALSE) ---\n")
res_nl <- solveWendy(f, p0, U, tt, lip = FALSE, method = "MLE",
  noise_dist = "addgaussian",
  control = list(test_fun_type = "SSL", optimize = FALSE))

K <- nrow(res_nl$V)
D <- ncol(U)
J <- length(p0)
cat(sprintf("  K=%d, D=%d, J=%d, K*D=%d\n", K, D, J, K * D))

wnll_R_nl <- res_nl$wnll
wnll_torch_nl <- build_wnll_torch(res_nl$S, res_nl$g, res_nl$b, K, D)

# Check shapes at a test point
cat("\n  Shape checks at p0:\n")
Sp <- res_nl$S(p0)
cat(sprintf("    S(p) shape: [%s], dtype: %s\n", paste(Sp$shape, collapse=", "), Sp$dtype$.type()))
r <- (res_nl$g(p0) - res_nl$b)$reshape(c(-1))
cat(sprintf("    r shape: [%s], dtype: %s\n", paste(r$shape, collapse=", "), r$dtype$.type()))
L <- linalg_cholesky(Sp)
cat(sprintf("    L shape: [%s], dtype: %s\n", paste(L$shape, collapse=", "), L$dtype$.type()))
r_col <- r$unsqueeze(2L)
cat(sprintf("    r_col shape: [%s]\n", paste(r_col$shape, collapse=", ")))
S_inv_r <- torch_cholesky_solve(r_col, L)
cat(sprintf("    S_inv_r shape: [%s]\n", paste(S_inv_r$shape, collapse=", ")))
S_inv_r_1d <- S_inv_r$squeeze(2L)
cat(sprintf("    S_inv_r squeezed shape: [%s]\n", paste(S_inv_r_1d$shape, collapse=", ")))

test_params <- list(
  p0,
  p_star,
  c(5.0, 14.0, 1.5),
  c(15.0, 35.0, 4.0),
  c(8.0, 25.0, 2.5)
)

cat("\n  Comparing wnll values:\n")
all_pass <- TRUE
for (i in seq_along(test_params)) {
  p <- test_params[[i]]
  val_R <- wnll_R_nl(p)
  val_torch <- wnll_torch_nl(p)
  diff <- abs(val_R - val_torch)
  rel_err <- diff / max(abs(val_R), 1e-15)
  pass <- rel_err < 5e-6
  cat(sprintf("    p=[%s]: R=%.8f, Torch=%.8f, rel_err=%.2e %s\n",
    paste(sprintf("%.1f", p), collapse=", "), val_R, val_torch, rel_err,
    ifelse(pass, "PASS", "FAIL")))
  if (!pass) all_pass <- FALSE
}

# --- Test linear (lip = TRUE) ---
cat("\n--- Linear case (lip = TRUE) ---\n")
res_li <- solveWendy(f, p0, U, tt, lip = TRUE, method = "MLE",
  noise_dist = "addgaussian",
  control = list(test_fun_type = "SSL", optimize = FALSE))

K_li <- nrow(res_li$V)
cat(sprintf("  K=%d, D=%d, J=%d, K*D=%d\n", K_li, D, J, K_li * D))

wnll_R_li <- res_li$wnll
wnll_torch_li <- build_wnll_torch(res_li$S, res_li$g, res_li$b, K_li, D)

cat("\n  Comparing wnll values:\n")
for (i in seq_along(test_params)) {
  p <- test_params[[i]]
  val_R <- wnll_R_li(p)
  val_torch <- wnll_torch_li(p)
  diff <- abs(val_R - val_torch)
  rel_err <- diff / max(abs(val_R), 1e-15)
  pass <- rel_err < 5e-6
  cat(sprintf("    p=[%s]: R=%.8f, Torch=%.8f, rel_err=%.2e %s\n",
    paste(sprintf("%.1f", p), collapse=", "), val_R, val_torch, rel_err,
    ifelse(pass, "PASS", "FAIL")))
  if (!pass) all_pass <- FALSE
}

# --- Timing ---
cat("\n--- Timing (100 evaluations, nonlinear) ---\n")
p_test <- c(8.0, 25.0, 2.5)

t_R <- system.time({ for (i in 1:100) wnll_R_nl(p_test) })
t_torch <- system.time({ for (i in 1:100) wnll_torch_nl(p_test) })

cat(sprintf("  R version:     %.4f s\n", t_R["elapsed"]))
cat(sprintf("  Torch version: %.4f s\n", t_torch["elapsed"]))
cat(sprintf("  Speedup:       %.2fx\n", t_R["elapsed"] / t_torch["elapsed"]))

cat(sprintf("\n=== wnll: %s ===\n", ifelse(all_pass, "ALL TESTS PASSED", "SOME TESTS FAILED")))

# =====================================================
# J_wnll tests
# =====================================================
cat("\n\n=== Test: build_J_wnll vs build_J_wnll_torch (Lorenz) ===\n\n")
all_pass_J <- TRUE

# --- Nonlinear J_wnll ---
cat("--- Nonlinear case (lip = FALSE) ---\n")
J_wnll_R_nl <- res_nl$J_wnll
J_wnll_torch_nl <- build_J_wnll_torch(res_nl$S, res_nl$Jp_S, res_nl$Jp_r, res_nl$g, res_nl$b, J)

# Shape checks
cat("  Shape checks at p0:\n")
J_Sp_check <- res_nl$Jp_S(p0)$contiguous()
J_rp_check <- res_nl$Jp_r(p0)$contiguous()
cat(sprintf("    Jp_S(p) shape: [%s], dtype: %s\n", paste(J_Sp_check$shape, collapse=", "), J_Sp_check$dtype$.type()))
cat(sprintf("    Jp_r(p) shape: [%s], dtype: %s\n", paste(J_rp_check$shape, collapse=", "), J_rp_check$dtype$.type()))

cat("\n  Comparing J_wnll values:\n")
for (i in seq_along(test_params)) {
  p <- test_params[[i]]
  grad_R <- J_wnll_R_nl(p)
  grad_torch <- J_wnll_torch_nl(p)
  max_abs_err <- max(abs(grad_R - grad_torch))
  max_rel_err <- max(abs(grad_R - grad_torch) / pmax(abs(grad_R), 1e-15))
  pass <- max_rel_err < 1e-4
  cat(sprintf("    p=[%s]:\n", paste(sprintf("%.1f", p), collapse=", ")))
  cat(sprintf("      R:     [%s]\n", paste(sprintf("%.8f", grad_R), collapse=", ")))
  cat(sprintf("      Torch: [%s]\n", paste(sprintf("%.8f", grad_torch), collapse=", ")))
  cat(sprintf("      max_rel_err=%.2e %s\n", max_rel_err, ifelse(pass, "PASS", "FAIL")))
  if (!pass) all_pass_J <- FALSE
}

# --- Linear J_wnll ---
cat("\n--- Linear case (lip = TRUE) ---\n")
J_wnll_R_li <- res_li$J_wnll
J_wnll_torch_li <- build_J_wnll_torch(res_li$S, res_li$Jp_S, res_li$Jp_r, res_li$g, res_li$b, J)

cat("  Comparing J_wnll values:\n")
for (i in seq_along(test_params)) {
  p <- test_params[[i]]
  grad_R <- J_wnll_R_li(p)
  grad_torch <- J_wnll_torch_li(p)
  max_rel_err <- max(abs(grad_R - grad_torch) / pmax(abs(grad_R), 1e-15))
  pass <- max_rel_err < 1e-4
  cat(sprintf("    p=[%s]: max_rel_err=%.2e %s\n",
    paste(sprintf("%.1f", p), collapse=", "), max_rel_err, ifelse(pass, "PASS", "FAIL")))
  if (!pass) all_pass_J <- FALSE
}

# --- Timing J_wnll ---
cat("\n--- Timing J_wnll (100 evaluations, nonlinear) ---\n")
t_R_J <- system.time({ for (i in 1:100) J_wnll_R_nl(p_test) })
t_torch_J <- system.time({ for (i in 1:100) J_wnll_torch_nl(p_test) })

cat(sprintf("  R version:     %.4f s\n", t_R_J["elapsed"]))
cat(sprintf("  Torch version: %.4f s\n", t_torch_J["elapsed"]))
cat(sprintf("  Speedup:       %.2fx\n", t_R_J["elapsed"] / t_torch_J["elapsed"]))

cat(sprintf("\n=== J_wnll: %s ===\n", ifelse(all_pass_J, "ALL TESTS PASSED", "SOME TESTS FAILED")))

# =====================================================
# H_wnll tests
# =====================================================
cat("\n\n=== Test: build_H_wnll vs build_H_wnll_torch (Lorenz) ===\n\n")
all_pass_H <- TRUE

# --- Nonlinear H_wnll ---
cat("--- Nonlinear case (lip = FALSE) ---\n")
H_wnll_R_nl <- res_nl$H_wnll
H_wnll_torch_nl <- build_H_wnll_torch(res_nl$S, res_nl$Jp_S, res_nl$L, res_nl$Jp_L, res_nl$Hp_L,
                                        res_nl$Jp_r, res_nl$Hp_r, res_nl$g, res_nl$b, J)

cat("  Comparing H_wnll values:\n")
for (idx in seq_along(test_params)) {
  p <- test_params[[idx]]
  H_R <- H_wnll_R_nl(p)
  H_torch <- H_wnll_torch_nl(p)
  max_abs_err <- max(abs(H_R - H_torch))
  max_rel_err <- max(abs(H_R - H_torch) / pmax(abs(H_R), 1e-15))
  pass <- max_rel_err < 1e-1
  cat(sprintf("    p=[%s]:\n", paste(sprintf("%.1f", p), collapse=", ")))
  cat(sprintf("      R Hessian:\n"))
  for (row in seq_len(J)) cat(sprintf("        [%s]\n", paste(sprintf("%12.4f", H_R[row,]), collapse=", ")))
  cat(sprintf("      Torch Hessian:\n"))
  for (row in seq_len(J)) cat(sprintf("        [%s]\n", paste(sprintf("%12.4f", H_torch[row,]), collapse=", ")))
  cat(sprintf("      max_rel_err=%.2e %s\n", max_rel_err, ifelse(pass, "PASS", "FAIL")))
  if (!pass) all_pass_H <- FALSE
}

# --- Linear H_wnll ---
cat("\n--- Linear case (lip = TRUE) using build_H_wnll_linear_torch ---\n")
H_wnll_R_li <- res_li$H_wnll
H_wnll_torch_li <- build_H_wnll_linear_torch(res_li$S, res_li$Jp_S, res_li$L, res_li$Jp_L,
                                               res_li$Jp_r, res_li$g, res_li$b, J)

cat("  Comparing H_wnll_linear values:\n")
for (idx in seq_along(test_params)) {
  p <- test_params[[idx]]
  H_R <- H_wnll_R_li(p)
  H_torch <- H_wnll_torch_li(p)
  max_rel_err <- max(abs(H_R - H_torch) / pmax(abs(H_R), 1e-15))
  pass <- max_rel_err < 1e-1
  cat(sprintf("    p=[%s]:\n", paste(sprintf("%.1f", p), collapse=", ")))
  cat(sprintf("      R Hessian:\n"))
  for (row in seq_len(J)) cat(sprintf("        [%s]\n", paste(sprintf("%12.4f", H_R[row,]), collapse=", ")))
  cat(sprintf("      Torch Hessian:\n"))
  for (row in seq_len(J)) cat(sprintf("        [%s]\n", paste(sprintf("%12.4f", H_torch[row,]), collapse=", ")))
  cat(sprintf("      max_rel_err=%.2e %s\n", max_rel_err, ifelse(pass, "PASS", "FAIL")))
  if (!pass) all_pass_H <- FALSE
}

# --- Timing H_wnll ---
cat("\n--- Timing H_wnll (20 evaluations, nonlinear) ---\n")
t_R_H <- system.time({ for (i in 1:20) H_wnll_R_nl(p_test) })
t_torch_H <- system.time({ for (i in 1:20) H_wnll_torch_nl(p_test) })

cat(sprintf("  R version:     %.4f s\n", t_R_H["elapsed"]))
cat(sprintf("  Torch version: %.4f s\n", t_torch_H["elapsed"]))
cat(sprintf("  Speedup:       %.2fx\n", t_R_H["elapsed"] / t_torch_H["elapsed"]))

cat("\n--- Timing H_wnll_linear (20 evaluations, linear) ---\n")
t_R_HL <- system.time({ for (i in 1:20) H_wnll_R_li(p_test) })
t_torch_HL <- system.time({ for (i in 1:20) H_wnll_torch_li(p_test) })

cat(sprintf("  R version:     %.4f s\n", t_R_HL["elapsed"]))
cat(sprintf("  Torch version: %.4f s\n", t_torch_HL["elapsed"]))
cat(sprintf("  Speedup:       %.2fx\n", t_R_HL["elapsed"] / t_torch_HL["elapsed"]))

cat(sprintf("\n=== H_wnll: %s ===\n", ifelse(all_pass_H, "ALL TESTS PASSED", "SOME TESTS FAILED")))
